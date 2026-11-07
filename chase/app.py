"""CHASE - Shiny for Python app.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import csv
import dataclasses
import getpass
import hashlib
import io as _io
import math
import os
import re
import time
from datetime import date
from pathlib import Path

import warnings

import numpy as np
import plotly
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shiny.types import SilentException

from . import __version__, help as helppage, plots
from . import render as rh
from .core import alerts, chem, classes, demo, io, normalize, plates, session as sessmod
from .core import cache as chase_cache
from .core import stats as cstats
from .core import triage as tri
from .core.triage import STATUS, STATUSES, TriageBoard

# Newer Shiny versions warn that this way of writing download handlers is old. It still works fine,
# so I hide the warning here.
warnings.filterwarnings("ignore", message="render.download is deprecated")
# The browser loads the files from CHASE. Nothing has to be downloaded from the internet.
WWW = Path(__file__).parent / "www"
PLOTLY_JS = Path(plotly.__file__).parent / "package_data"


def _col_choices(header, allow_none=False):
    """Return select-box choices ``{"0": header0, ...}`` for a header row, optionally with '(none)'."""
    ch = {str(i): (h or f"Column {i + 1}") for i, h in enumerate(header)}
    return ({"-1": "(none)"} | ch) if allow_none else ch


# The colour pickers shown under 'Colors': (name in DEFAULT_COLORS, label next to the picker).
# chase.js finds these inputs by their class, chase-color.
COLOR_KEYS = [("active", "Active"), ("moderate", "Moderate (k = 4)"), ("weak", "Weak"), ("inactive", "Inactive"),
              ("pos", "Positive control"), ("neg", "Negative control"), ("pains", "PAINS mark"), ("covalent", "Covalent mark")]
COLOR_ROW = ('<div class="ctools colorrow" id="colorrow" hidden>'
             + "".join(f'<label class="cpick"><input type="color" class="chase-color" data-key="{k}" value="{plots.DEFAULT_COLORS[k]}" aria-label="{n} colour">{n}</label>'
                       for k, n in COLOR_KEYS)
             + '<button type="button" class="btn btn-sm btn-outline-secondary" data-chase="color_reset">Reset colours</button>'
             + '<span class="hint">Colours apply to the plots, chips and chemotype strips, and are remembered in this browser.</span></div>')
# Colours coming from the browser must look like #rrggbb. Anything else is ignored.
HEX = re.compile(r"^#[0-9a-fA-F]{6}$")


def safe_name(s: str) -> str:
    """Make a string safe for use in a file name."""
    return re.sub(r"[^\w.-]+", "_", s).strip("_") or "campaign"


def step(n, title, *children):
    """Return a numbered step box for the sidebar."""
    return ui.div(ui.h2(ui.span(str(n), class_="n"), title), *children, class_="step")


def card(title, *children, header_extra=(), id=None):
    """Return a card with a title bar (plus optional controls) and a body."""
    return ui.div(ui.div(ui.h3(title), *header_extra, class_="ccard-h"), *children, class_="ccard", id=id)


# ================================================================ user interface
# The sidebar on the left with the four setup steps.
sidebar = ui.sidebar(
    step(1, "Load files",
         ui.input_file("screen_file", "Screening data", accept=[".xlsx", ".xls", ".csv", ".tsv", ".txt"], placeholder="Raw Data, Plate, Well, Batch ID"),
         ui.output_ui("screen_hint"),
         ui.accordion(ui.accordion_panel("Columns",
                                         ui.div(ui.input_select("col_raw", "Raw signal", {}), ui.input_select("col_plate", "Plate", {}), class_="row2"),
                                         ui.div(ui.input_select("col_well", "Well", {}), ui.input_select("col_bid", "Batch ID", {}), class_="row2"),
                                         ui.input_select("sheet", "Sheet", {})), open=False, id="acc_screen"),
         ui.input_file("struct_file", "Structures", accept=[".csv", ".tsv", ".txt", ".smi", ".xlsx"], placeholder="Batch ID + SMILES"),
         ui.output_ui("struct_hint"),
         ui.accordion(ui.accordion_panel("Columns",
                                         ui.input_select("col_sid", "Join on (ID)", {}),
                                         ui.div(ui.input_select("col_smi", "SMILES", {}), ui.input_select("col_name", "Display name", {}), class_="row2")),
                      open=False, id="acc_struct"),
         ui.input_file("session_file", "Session", accept=[".json"], button_label="Open session...", placeholder="CHASE .json file"),
         ui.download_button("save_session", "Save session", class_="btn-sm btn-outline-secondary w-100"),
         ui.input_action_button("load_example", "Load example data", class_="btn-sm btn-outline-secondary w-100 mt-2")),
    step(2, "Controls & normalization",
         ui.div(ui.input_selectize("pos", "Positive (100%)", [], options={"create": True}),
                ui.input_selectize("neg", "Negative (0%)", [], options={"create": True}), class_="row2"),
         ui.output_ui("ctrl_hint"),
         ui.input_radio_buttons("scope", "Normalize", {"plate": "Per plate", "global": "All plates"}, inline=True),
         ui.input_radio_buttons("stat", "Control statistic", {"mean": "Mean", "median": "Median"}, inline=True),
         ui.p("% = 100 * (raw - neg) / (pos - neg)", class_="hint mono"),
         ui.input_action_button("run", "Normalize & plot", class_="btn-primary w-100 mt-2")),
    step(3, "Activity classes",
         ui.input_select("cls_method", "Method", {"kmeans": "k-means on % response", "threshold": "Fixed thresholds", "robust": "Robust z (median + n*MAD)"}),
         ui.panel_conditional("input.cls_method === 'kmeans'", ui.input_radio_buttons("k", "Number of classes (k)", {"2": "2", "3": "3", "4": "4"}, selected="3", inline=True)),
         ui.panel_conditional("input.cls_method === 'threshold'", ui.div(ui.input_numeric("thr_a", "Active >= %", 50), ui.input_numeric("thr_w", "Weak >= %", 20), class_="row2")),
         ui.panel_conditional("input.cls_method === 'robust'", ui.input_numeric("rob_n", "n * robust SD", 3, step=0.5, min=0.5)),
         ui.output_ui("cls_hint")),
    step(4, "Chemotypes",
         ui.input_radio_buttons("chemo_method", "Group by", {"butina": "Similarity", "murcko": "Scaffold"}, inline=True),
         ui.panel_conditional("input.chemo_method === 'butina'", ui.input_slider("sim", "Tanimoto >=", 0.3, 0.9, 0.5, step=0.05)),
         ui.output_ui("chemo_hint"),
         ui.download_button("export_csv", "Export results (.csv)", class_="btn-sm btn-outline-secondary w-100")),
    width=340, open="always", bg="#f3f5f4",
)

# The main area ('bench'): summary tiles, the % response plot, chemotypes next to the detail panel,
# the triage board, activity cliffs, plates (maps, replicates, QC) and data checks.
# The plots start as empty, and chase.js draws into them.
bench = ui.div(
    ui.output_ui("banner"),
    ui.output_ui("errbox"),
    ui.output_ui("stats"),
    card("% response by compound",
         ui.HTML(COLOR_ROW),
         ui.div(id="plot"),
         ui.div("Hover a point to see its structure | drag to zoom | double-click to reset | click a point to pin it | bars above the plot show each plate's Z' (dotted line at 0.5)", class_="plot-foot"),
         header_extra=(ui.output_ui("plot_meta", inline=True), ui.span(class_="spacer"),
                       ui.input_select("flag_mode", None, {"hits": "Flags: hits only", "all": "Flags: all compounds", "none": "Flags: off"}, selected="hits", width="170px"),
                       ui.input_checkbox("show_z", "Z' per plate", True),
                       ui.input_checkbox("show_ctrl", "Control wells", True),
                       ui.tags.button("Colors", class_="btn btn-sm btn-outline-secondary", type="button", onclick="chaseToggle('colorrow', this)", **{"aria-expanded": "false", "aria-controls": "colorrow"}),
                       ui.tags.button("Save JPEG", class_="btn btn-sm btn-outline-secondary", onclick="chaseJPEG('plot', 'CHASE_scatter')"),
                       ui.download_button("plot_csv", "Save data (.csv)", class_="btn-sm btn-outline-secondary"))),
    ui.div(
        card("Chemotypes",
             ui.div(ui.input_select("sort", "Sort", {}), ui.input_numeric("min_size", "Min. size", 2, min=1, step=1),
                    ui.input_checkbox("incl_single", "Include singleton hits", True), class_="ctools"),
             ui.output_ui("chemo_grid"),
             ui.div(ui.input_action_button("more", "Show more", class_="btn-sm btn-outline-secondary"), class_="text-center pb-3"),
             header_extra=(ui.input_radio_buttons("chemo_tab", None, {"act": "Active", "ina": "Inactive"}, inline=True), ui.output_ui("chemo_counts", inline=True))),
        ui.div(ui.output_ui("detail"), class_="ccard", id="detailcard"),
        class_="split"),
    card("Triage board",
         ui.output_ui("tri_notice"),
         ui.output_ui("board"),
         header_extra=(ui.output_ui("tri_meta", inline=True), ui.span(class_="spacer"),
                       ui.input_text("campaign", "Campaign", ""), ui.input_text("tri_user", "Your name", ""),
                       ui.download_button("export_triage", "Export triage (.csv)", class_="btn-sm btn-outline-secondary"),
                       ui.div(ui.input_file("triage_file", None, accept=[".csv"], button_label="Import triage (.csv)...", placeholder=""), class_="trimport"),
                       ui.tags.button("Clear board", class_="btn btn-sm btn-outline-secondary", type="button", **{"data-chase": "tri_clear"}))),
    card("Activity cliffs",
         ui.div(ui.input_slider("cliff_sim", "Tanimoto >=", 0.6, 0.95, 0.7, step=0.05), ui.input_numeric("cliff_delta", "|Delta % response| >=", 25, min=0, step=5),
                ui.input_checkbox("cliff_chemo", "Selected chemotype only", False), class_="ctools"),
         ui.output_ui("cliffs"),
         ui.div(ui.input_action_button("cliff_more", "Show more", class_="btn-sm btn-outline-secondary"), class_="text-center pb-3"),
         header_extra=(ui.output_ui("cliff_meta", inline=True),)),
    ui.navset_card_underline(
        ui.nav_panel("Plate maps",
                     ui.div(ui.span("% response"), ui.span(class_="hmbar"), ui.span("-50  0  100+", class_="mono"),
                            ui.span("Outlined cells are wells without data. Click a plate to inspect it.", class_="hint"), class_="ctools hmlegend"),
                     ui.output_ui("minimaps"),
                     ui.div(ui.output_ui("bighead", inline=True), ui.span(class_="sp"),
                            ui.tags.button("Save JPEG", class_="btn btn-sm btn-outline-secondary", onclick="chaseHeatJPEG('bigmap', 'CHASE_plate')"), class_="bighead"),
                     ui.div(id="bigmap")),
        ui.nav_panel("Replicates",
                     ui.div(ui.input_numeric("rep_t", "Flag replicates differing by more than (points)", 20, min=1, step=5), class_="ctools"),
                     ui.output_ui("rep_stats"), ui.div(id="repplot"), ui.output_ui("rep_list")),
        ui.nav_panel("QC table", ui.output_ui("qc")),
        id="plates_tab", title=ui.h3("Plates", class_="m-0 me-3"),
    ),
    ui.accordion(ui.accordion_panel("Data checks", ui.output_ui("notes")), open=False),
    class_="bench",
)

def _stamp(name: str) -> str:
    """Return a short fingerprint of a file in www/, added to its URL as ?v=...

    Browsers keep a copy of chase.css and chase.js. When the file changes, its fingerprint
    and therefore its URL change too, so the browser loads the new version instead of the old copy.
    """
    return hashlib.sha1((WWW / name).read_bytes()).hexdigest()[:10]


# Two tabs: 'Analysis' and 'How to use'. The page head loads the stylesheet, plotly.js and chase.js.
app_ui = ui.page_navbar(
    ui.nav_panel("Analysis", ui.layout_sidebar(sidebar, bench, fillable=False)),
    ui.nav_panel("How to use", ui.div(ui.HTML(helppage.BODY), ui.tags.section(ui.h3("Software versions"), ui.tags.pre(helppage.versions_text(), class_="formula")),
                                      ui.HTML(helppage.ACK), ui.HTML(helppage.REFS), ui.HTML(helppage.FOOT), class_="help")),
    title=ui.span(ui.span("CHASE", class_="brand"), ui.HTML(helppage.ACRO.replace('class="acro"', 'class="acro d-none d-md-inline"'))),
    id="nav", fillable=False, window_title="CHASE",
    header=ui.head_content(ui.tags.link(rel="stylesheet", href=f"www/chase.css?v={_stamp('chase.css')}"), ui.tags.script(src="plotly/plotly.min.js"), ui.tags.script(src=f"www/chase.js?v={_stamp('chase.js')}")),
)


def consume(f: dict, reader):
    """Read an uploaded file, then delete Shiny's temporary copy straight away (nothing stays on the server)."""
    p = Path(f["datapath"])
    try:
        return reader(p)
    finally:
        try:
            p.unlink(missing_ok=True)
            p.parent.rmdir()          # Shiny puts every upload in its own folder, so remove that too
        except OSError:
            pass


# The name suggested for triage decisions: the user account that is running CHASE.
DEFAULT_USER = getpass.getuser()


def server(input: Inputs, output: Outputs, session: Session):
    """Shiny server function: one call per browser session. All data live in this session's reactive values."""
    # ================================================================ per-session state
    # When the browser tab is closed, throw away this session's caches (structure drawings, alerts).
    session.on_ended(lambda: chase_cache.drop(session.id))
    rv_screen = reactive.value(None)
    rv_struct = reactive.value(None)
    rv_demo = reactive.value(False)
    rv_an = reactive.value(None)
    # rv_ver goes up by one every time the data change (a new analysis or new structures).
    # The calculations below read it, so they know when to start again.
    rv_ver = reactive.value(0)
    rv_err = reactive.value("")
    rv_cpd = reactive.value(None)          # number of the compound pinned in the detail panel
    rv_anchor = reactive.value(None)       # the selected chemotype is 'the one that contains this compound'
    rv_pair = reactive.value(None)
    rv_plate = reactive.value(0)
    rv_shown = reactive.value(24)
    rv_cliff_shown = reactive.value(20)
    # rv_bump goes up by one every time the triage board changes.
    rv_bump = reactive.value(0)
    rv_pending = reactive.value(None)
    rv_note = reactive.value("")           # a short message shown once under the triage board
    tboard = TriageBoard()                 # this tab's decisions: in memory only
    # The last list of similar pairs computed: (data version, threshold used, the pairs).
    pair_cache: dict = {}
    rv_colors = reactive.value(dict(plots.DEFAULT_COLORS))

    async def push_colors(c: dict):
        """Send colours to the browser (updates the colour pickers and the browser's remembered choice)."""
        await session.send_custom_message("chase-colors", {"colors": c})

    @reactive.effect
    @reactive.event(input.chase_colors)
    def _colors_in():
        """Take colour choices from the browser, keeping only valid ``#rrggbb`` values."""
        got = input.chase_colors() or {}
        c = dict(plots.DEFAULT_COLORS)
        c.update({k: v for k, v in got.items() if k in c and isinstance(v, str) and HEX.match(v)})
        if c != rv_colors.get():
            rv_colors.set(c)

    @reactive.effect
    @reactive.event(input.chase_color_reset)
    async def _colors_reset():
        """Restore the default colours."""
        rv_colors.set(dict(plots.DEFAULT_COLORS))
        await push_colors(dict(plots.DEFAULT_COLORS))
    # The campaign name is filled in from the file name until the user types a name of their own.
    camp = {"name": "", "auto": True, "last_auto": ""}

    def bump(rv):
        """Increment a counter value without creating a reactive dependency (used to signal changes)."""
        with reactive.isolate():
            rv.set(rv.get() + 1)

    # ------------------------------------------------------------ files
    def set_screen(sf: io.ScreenFile, label: str):
        """Make ``sf`` the current screening file: fill the column selectors and suggest the control IDs.

        Returns:
            The extracted well table using the guessed columns.
        """
        rv_screen.set(sf)
        g = sf.guess()
        ch = _col_choices(sf.header)
        for k in ("raw", "plate", "well", "bid"):
            ui.update_select(f"col_{k}", choices=ch, selected=str(g[k]))
        ui.update_select("sheet", choices={n: n for n in sf.sheets}, selected=sf.sheet)
        rows = sf.extract(g)
        cands, pos, neg = io.control_candidates(rows)
        # Offer the most frequent IDs as choices for the controls and select the best guess.
        # The user can still type any other ID.
        opts = [c for c, _ in cands]
        ui.update_selectize("pos", choices=opts, selected=pos or None)
        ui.update_selectize("neg", choices=opts, selected=neg or None)
        sf.label = f"{label} | {len(rows):,} wells | {rows['plate'].nunique()} plates"
        sf.ctrl_hint = (f"Suggested from IDs repeated on most plates: {', '.join(f'{c} ({n})' for c, n in cands[:4])}. Check before normalizing."
                        if cands else "No repeated Batch IDs found. Type the control IDs.")
        return rows

    def set_struct(stf: io.StructFile, label: str):
        """Make ``stf`` the current structures file and fill its column selectors."""
        rv_struct.set(stf)
        ui.update_select("col_sid", choices=_col_choices(stf.header), selected=str(stf.id))
        ui.update_select("col_smi", choices=_col_choices(stf.header), selected=str(stf.smi))
        ui.update_select("col_name", choices=_col_choices(stf.header, True), selected=str(stf.name))
        stf.label = f"{label} | {len(stf.body):,} structures"

    def screen_cols(sf):
        """Column mapping chosen in the selectors, falling back to the guess when a selector is not yet valid."""
        g = sf.guess()
        out = {}
        for k in ("raw", "plate", "well", "bid"):
            try:
                # The column selectors hold numbers as text. Before they are filled in, they can be empty or out of date.
                v = int(input[f"col_{k}"]())
            except (TypeError, ValueError, SilentException):
                v = g[k]
            out[k] = v if 0 <= v < len(sf.header) else g[k]
        return out

    @reactive.calc
    def rows():
        """Well table of the current screening file with the chosen columns."""
        sf = rv_screen.get()
        req(sf)
        return sf.extract(screen_cols(sf))

    @reactive.calc
    def struct_map():
        """Batch-ID -> structure mapping of the current structures file with the chosen columns (None if none)."""
        stf = rv_struct.get()
        if stf is None:
            return None
        for attr, iid in (("id", "col_sid"), ("smi", "col_smi"), ("name", "col_name")):
            try:
                v = int(input[iid]())
                if -1 <= v < len(stf.header):
                    setattr(stf, attr, v)
            except (TypeError, ValueError, SilentException):
                pass
        return stf.mapping()

    def clear_selection():
        """Unpin the compound and deselect the chemotype and cliff pair."""
        rv_cpd.set(None)
        rv_anchor.set(None)
        rv_pair.set(None)

    def do_analyze(rows_df, pos, neg, scope, stat, smap, label=None):
        """Normalize, store the analysis and reset selections; report control problems in the error box.

        Returns:
            True on success.
        """
        try:
            an = normalize.analyze(rows_df, pos or "", neg or "", scope, stat, smap)
        except normalize.NormalizationError as e:
            rv_err.set(str(e))
            return False
        rv_err.set("")
        # Remember the file name. It is used for the campaign name and for the names of exported files.
        an.file_name = label or (rv_screen.get().file_name if rv_screen.get() else "screen")
        rv_an.set(an)
        clear_selection()
        rv_plate.set(0)
        bump(rv_ver)
        if camp["auto"]:
            name = "Example screen" if rv_demo.get() else re.sub(r"\.[^.]+$", "", an.file_name.split(" (")[0])
            camp["last_auto"] = name
            ui.update_text("campaign", value=name)
        return True

    def load_demo():
        """Load and analyse the generated example screen and structures."""
        screen, struct = demo.make_demo()
        rv_demo.set(True)
        camp["auto"] = True
        sf = io.ScreenFile.from_sheets({"Example": screen}, "example_screen.xlsx")
        rows_df = set_screen(sf, "example_screen.xlsx (generated)")
        stf = io.StructFile.from_sheets({"s": struct}, "example_structures.csv", set(rows_df["key"]))
        set_struct(stf, "example_structures.csv (generated)")
        ui.update_selectize("pos", selected="CTRL")
        ui.update_selectize("neg", selected="DMSO")
        do_analyze(rows_df, "CTRL", "DMSO", "plate", "mean", stf.mapping(), "example_screen.xlsx")

    # Load the example data once when the app opens, so the user doesn't start with an empty screen.
    started = {"done": False}

    @reactive.effect
    def _start():
        """Show the example data when a session starts (runs once)."""
        if started["done"]:
            return
        started["done"] = True
        with reactive.isolate():
            load_demo()

    @reactive.effect
    @reactive.event(input.load_example)
    def _example():
        """Load example data on request."""
        load_demo()

    def leave_demo():
        """Drop the example data when the user loads their own file."""
        if rv_demo.get():
            rv_demo.set(False)
            rv_screen.set(None)
            rv_struct.set(None)
            rv_an.set(None)
            clear_selection()
            camp["auto"] = True

    @reactive.effect
    @reactive.event(input.screen_file)
    def _screen_upload():
        """Read an uploaded screening file."""
        f = input.screen_file()[0]
        try:
            sheets = consume(f, lambda p: io.read_table(p, f["name"]))
            sf = io.ScreenFile.from_sheets(sheets, f["name"])
        except Exception as e:  # noqa: BLE001
            rv_err.set(f"Could not read the screening file: {e}")
            return
        leave_demo()
        camp["auto"] = True
        rv_an.set(None)
        clear_selection()
        rv_err.set("")
        set_screen(sf, f["name"])

    @reactive.effect
    @reactive.event(input.sheet)
    def _sheet():
        """Switch to another sheet of the screening workbook."""
        sf = rv_screen.get()
        if sf is None or input.sheet() not in sf.sheets or input.sheet() == sf.sheet:
            return
        nsf = io.ScreenFile(file_name=sf.file_name, sheets=sf.sheets)
        nsf.set_sheet(input.sheet())
        set_screen(nsf, sf.file_name)

    @reactive.effect
    @reactive.event(input.struct_file)
    def _struct_upload():
        """Read an uploaded structures file."""
        f = input.struct_file()[0]
        try:
            sheets = consume(f, lambda p: io.read_table(p, f["name"]))
            # Hand over the Batch IDs of the screen. Then the ID column of the structures file can be found by
            # checking which column matches them best.
            keys = set(rv_screen.get().extract(rv_screen.get().guess())["key"]) if (rv_screen.get() and not rv_demo.get()) else None
            stf = io.StructFile.from_sheets(sheets, f["name"], keys)
        except Exception as e:  # noqa: BLE001
            rv_err.set(f"Could not read the structures file: {e}")
            return
        leave_demo()
        rv_err.set("")
        set_struct(stf, f["name"])

    @reactive.effect
    def _relink():
        """Re-attach structures to compounds when the structures file or its columns change."""
        m = struct_map()
        with reactive.isolate():
            an = rv_an.get()
            if an is None:
                return
            # Selections point at compound numbers. Clear them if the structures really changed.
            before = an.cpds["smiles"].tolist()
            normalize.relink(an, m)
            if an.cpds["smiles"].tolist() != before:
                clear_selection()
                bump(rv_ver)

    @reactive.effect
    @reactive.event(input.run)
    def _run():
        """'Normalize & plot' button."""
        sf = rv_screen.get()
        if sf is None:
            rv_err.set("Load a screening file first.")
            return
        do_analyze(rows(), input.pos(), input.neg(), input.scope(), input.stat(), struct_map())

    # ------------------------------------------------------------ analysis chain
    # The analysis runs as a chain of steps, each feeding the next:
    # rv_an (normalized data) -> cl_base (activity classes) -> prep (fingerprints) -> pairs (similar pairs)
    # -> groups (chemotypes) -> chemo (chemotype statistics).
    # Shiny only reruns a step when something it depends on has changed.
    @reactive.calc
    def cl_base():
        """Activity classes from the current method and parameters."""
        an = rv_an.get()
        req(an)
        rv_ver.get()
        k = int(input.k() or 3)
        return classes.assign(an.cpds["mean"].to_numpy(), input.cls_method(), k, input.thr_a(), input.thr_w(), input.rob_n() or 3)

    @reactive.calc
    def cl():
        """Classes plus this user's colours (changing colours does not recompute classes or chemotype statistics)."""
        return dataclasses.replace(cl_base(), colors=dict(rv_colors.get()))

    @reactive.calc
    def prep():
        """Fingerprints and scaffolds of all compounds with structures (None if there are none)."""
        an = rv_an.get()
        req(an)
        rv_ver.get()
        smiles = an.cpds["smiles"].tolist()
        if not any(smiles):
            return None
        with ui.Progress(min=0, max=1) as p:
            p.set(0, message="Reading structures")
            return chem.prepare(smiles, lambda f, m: p.set(f, message=m))

    @reactive.calc
    def flags():
        """PAINS / covalent flags, stored on an.cpds as '; '-joined alert names. Returns the number of flagged compounds."""
        an = rv_an.get()
        req(an)
        rv_ver.get()
        smiles = an.cpds["smiles"].tolist()
        if any(smiles):
            with ui.Progress(min=0, max=1) as p:
                p.set(0, message="Checking PAINS and covalent alerts")
                res = alerts.flag_many(smiles, lambda f, m: p.set(f, message=m))
        else:
            res = [((), ())] * len(smiles)
        an.cpds["pains"] = ["; ".join(a) for a, _ in res]
        an.cpds["covalent"] = ["; ".join(b) for _, b in res]
        return int(sum(1 for a, b in res if a or b))

    @reactive.calc
    def pairs():
        """All similarity pairs needed for clustering and cliffs, cached and reused when only the cut-off rises."""
        pr = prep()
        if pr is None or len(pr.idx) > cstats.PAIR_LIMIT:
            return None
        # The same list of similar pairs is used for the clustering (Tanimoto >= sim) and for the activity
        # cliffs (>= CLIFF_MIN), so we compute it for the lower of the two thresholds. If we already have a
        # list for the same data that was made with the same or a lower threshold, we simply reuse it.
        need = min(input.sim(), cstats.CLIFF_MIN) if input.chemo_method() == "butina" else cstats.CLIFF_MIN
        v = rv_ver.get()
        c = pair_cache.get("p")
        if c and c[0] == v and c[1] <= need + 1e-9:
            return c[2]
        with ui.Progress(min=0, max=1) as p:
            p.set(0, message="Comparing fingerprints")
            P = chem.sim_pairs(pr.fps, need, lambda f, m: p.set(f, message=m))
        pair_cache["p"] = (v, need, P)
        return P

    @reactive.calc
    def groups():
        """Chemotypes from the chosen method (independent of activity classes)."""
        pr = prep()
        if pr is None:
            return None
        an = rv_an.get()
        return cstats.build_chemotypes(pr, an.cpds["smiles"].tolist(), input.chemo_method(), input.sim(), pairs())

    @reactive.calc
    def chemo():
        """Chemotypes with statistics for the current activity classes."""
        g = groups()
        if g is None:
            return None
        c = cl_base()
        cstats.chemo_stats(g, rv_an.get().cpds["mean"].to_numpy(), c.cls, len(c.labels))
        return g

    @reactive.calc
    def sel_chemo_id():
        """ID of the selected chemotype: the one containing the anchor compound (survives re-clustering)."""
        a = rv_anchor.get()
        if a is None:
            return None
        ch = chemo()
        return ch.cpd_chemo.get(a) if ch else None

    @reactive.calc
    def key_to_idx():
        """Compound key -> compound index for the current analysis."""
        an = rv_an.get()
        rv_ver.get()
        return an.key_to_idx if an is not None else {}

    # ------------------------------------------------------------ triage
    @reactive.effect
    def _campaign():
        """Track the campaign name; once the user types one, it is no longer replaced automatically."""
        v = (input.campaign() or "").strip()
        # If the name differs from the one filled in, the user typed it.
        if v and v != camp["last_auto"]:
            camp["auto"] = False
        camp["name"] = v or camp["last_auto"] or "Campaign"

    @reactive.effect
    def _user_default():
        """Fill 'Your name' with the name remembered by the browser, or the computer account name."""
        try:
            # chase.js sends the name this browser remembered from last time, once the page is connected.
            saved = input.chase_saved_user()
        except SilentException:
            saved = None
        with reactive.isolate():
            cur = (input.tri_user() or "").strip()
            if saved and cur in ("", DEFAULT_USER):
                ui.update_text("tri_user", value=saved)
            elif not cur:
                if DEFAULT_USER:
                    ui.update_text("tri_user", value=DEFAULT_USER)

    @reactive.effect
    @reactive.event(input.tri_user, ignore_init=True)
    async def _user_store():
        """Ask the browser to remember 'Your name'."""
        if (input.tri_user() or "").strip():
            await session.send_custom_message("chase-store", {"key": "chase.user", "value": input.tri_user().strip()})

    def campaign_name():
        """Current campaign name (typed, automatic, or 'Campaign')."""
        v = (input.campaign() or "").strip()
        return v or camp["last_auto"] or "Campaign"

    @reactive.calc
    def decisions():
        """Current triage decisions of this tab."""
        rv_bump.get()
        camp["name"] = campaign_name()
        return tboard.decisions()

    @reactive.effect
    async def _dirty():
        """Tell the browser how many changes are unsaved (it warns before the tab is closed)."""
        decisions()
        await session.send_custom_message("chase-dirty", {"n": tboard.unsaved})

    @reactive.effect
    @reactive.event(input.chase_saved)
    def _saved():
        """The board was just downloaded (session or triage CSV): nothing is unsaved any more."""
        tboard.mark_saved()
        rv_note.set("")
        bump(rv_bump)

    @reactive.effect
    @reactive.event(input.chase_tri_clear)
    def _ask_clear():
        """Ask for confirmation before clearing the board."""
        if not len(tboard):
            return
        ui.modal_show(ui.modal(
            ui.p(f"Remove all {len(tboard)} decisions from the board?" + (" Some of them have not been saved yet; export them first if you want to keep them." if tboard.unsaved else "")),
            title="Clear the triage board",
            footer=ui.div(ui.modal_button("Cancel", class_="btn-outline-secondary"), ui.input_action_button("tri_clear_ok", "Clear board", class_="btn-danger ms-2")),
            easy_close=True))

    @reactive.effect
    @reactive.event(input.tri_clear_ok)
    def _clear():
        """Clear the board after confirmation."""
        tboard.clear()
        ui.modal_remove()
        rv_note.set("")
        bump(rv_bump)

    @reactive.effect
    @reactive.event(input.triage_file)
    def _import_triage():
        """Import a triage CSV and merge it into the board (newer decision wins per compound)."""
        f = input.triage_file()[0]
        try:
            ds, bad, cname = tri.from_csv(consume(f, lambda p: p.read_text(encoding="utf-8-sig")))
        except (ValueError, UnicodeDecodeError) as e:
            rv_note.set(f"Could not import {f['name']}: {e}")
            bump(rv_bump)
            return
        was_empty = not len(tboard)
        # If the board was empty, the import just brings back decisions that were saved before, so none of
        # them counts as unsaved. In that case we also take over the campaign name from the file.
        n = tboard.merge(ds)
        if was_empty:
            tboard.mark_saved()
        if was_empty and cname and camp["auto"]:
            camp.update(auto=False, last_auto="")
            ui.update_text("campaign", value=cname)
        rv_note.set(f"Imported {n} decision{'s' if n != 1 else ''} from {f['name']}" + (f"; skipped {len(bad)} unreadable row{'s' if len(bad) != 1 else ''}." if bad else "."))
        bump(rv_bump)

    def user():
        """Name recorded with decisions."""
        return (input.tri_user() or "").strip() or DEFAULT_USER or "anonymous"

    def write(key: str, status: str | None, note: str | None = None):
        """Set (``status``) or remove (``status=None``) the decision for compound ``key``, keeping the existing note unless ``note`` is given."""
        an = rv_an.get()
        if an is None:
            return
        idx = key_to_idx().get(key)
        with reactive.isolate():
            prev = decisions().get(key)
            if status is None:
                tboard.delete(key)
                bump(rv_bump)
                return
            c = an.cpds.loc[idx] if idx is not None else None
            ch = chemo()
            k = cl()
            # Save a copy of the compound's details with the decision (ID, name, %, class, chemotype), so the board
            # and the CSV can be read on their own. If the compound is not in the data loaded now, keep the details
            # from its previous decision.
            doc = dict(key=key, bid=c["bid"] if c is not None else (prev or {}).get("bid", key), name=c["name"] if c is not None else (prev or {}).get("name", ""),
                       status=status, note=note if note is not None else (prev or {}).get("note", ""),
                       pct=round(float(c["mean"]), 2) if c is not None else (prev or {}).get("pct"),
                       cls=k.labels[k.cls[idx]][0] if idx is not None else (prev or {}).get("cls", ""),
                       chemo=(ch.cpd_chemo.get(idx, "") if ch and idx is not None else ""), by=user(), at=time.time() * 1000)
        tboard.set(doc)
        bump(rv_bump)

    @reactive.effect
    @reactive.event(input.chase_tri)
    def _tri():
        """Decision button: set it, or clear it when the same decision is clicked again."""
        d = input.chase_tri()
        key, status = d.get("key"), d.get("status")
        with reactive.isolate():
            cur = decisions().get(key)
        if status == "clear" or (cur and cur["status"] == status):
            write(key, None)
        elif status in STATUS:
            write(key, status)

    @reactive.effect
    @reactive.event(input.chase_triall)
    def _triall():
        """Apply a decision to the selected chemotype: its undecided members, or all if every member is decided."""
        status = input.chase_triall().get("status")
        with reactive.isolate():
            ch = chemo()
            cid = sel_chemo_id()
            if not ch or not cid or status not in STATUS:
                return
            an = rv_an.get()
            keys = [an.cpds.at[i, "key"] for i in ch.by_id[cid].members]
            dec = decisions()
            # Only fill in members that have no decision yet. If all of them already have one, overwrite them all.
            todo = [k for k in keys if k not in dec] or keys
        for k in todo:
            write(k, status)
        ui.notification_show(f"{len(todo)} compounds set to {STATUS[status][1]}", duration=3)

    @reactive.effect
    @reactive.event(input.chase_note)
    def _note():
        """Save an edited note (sent when the note box loses focus)."""
        d = input.chase_note()
        with reactive.isolate():
            cur = decisions().get(d.get("key"))
        if cur and (cur.get("note") or "") != d.get("note", ""):
            write(d["key"], cur["status"], d.get("note", ""))

    # ------------------------------------------------------------ selection events
    def pin_cpd(ci: int, follow_chemo: bool):
        """Pin compound ``ci``; with ``follow_chemo`` also select its chemotype and switch to the right tab."""
        # The selected chemotype is tied to a compound rather than to a chemotype number, so the selection
        # survives when the chemotypes are recalculated with other settings.
        rv_cpd.set(ci)
        if follow_chemo:
            ch = chemo()
            if ch and ci in ch.cpd_chemo:
                rv_anchor.set(ci)
                tab = "act" if ch.by_id[ch.cpd_chemo[ci]].n_hit > 0 else "ina"
                if tab != input.chemo_tab():
                    ui.update_radio_buttons("chemo_tab", selected=tab)

    @reactive.effect
    @reactive.event(input.chase_sel_chemo)
    def _sel_chemo():
        """Chemotype card clicked: select it, or deselect it if already selected."""
        cid = input.chase_sel_chemo().get("id")
        ch = chemo()
        if not ch or cid not in ch.by_id:
            return
        # Clicking the selected card a second time deselects it. Otherwise the chemotype is tied to its most
        # active member, and a pinned compound from another chemotype is unpinned.
        if sel_chemo_id() == cid:
            rv_anchor.set(None)
        else:
            rv_anchor.set(ch.by_id[cid].members[0])
            if rv_cpd.get() is not None and ch.cpd_chemo.get(rv_cpd.get()) != cid:
                rv_cpd.set(None)

    @reactive.effect
    @reactive.event(input.chase_sel_cpd)
    def _sel_cpd():
        """Compound row clicked: pin it."""
        pin_cpd(int(input.chase_sel_cpd()["cpd"]), False)

    @reactive.effect
    @reactive.event(input.chase_click)
    def _click():
        """Point clicked on a plot: pin the compound behind it."""
        d = input.chase_click()
        an = rv_an.get()
        if an is None:
            return
        # Points on the replicate plot carry a compound number. Points on the main plot carry a well number.
        if d["plot"] == "repplot":
            pin_cpd(int(d["idx"]), True)
            return
        w = an.wells.loc[int(d["idx"])]
        if w["type"] == "sample":
            pin_cpd(int(w["cpd"]), True)

    @reactive.effect
    @reactive.event(input.chase_sel_pair)
    async def _sel_pair():
        """Cliff pair clicked: select or deselect it and scroll to the plot."""
        d = input.chase_sel_pair()
        p = (int(d["hi"]), int(d["lo"]))
        rv_pair.set(None if rv_pair.get() == p else p)
        if rv_pair.get():
            await session.send_custom_message("chase-scroll", {"id": "plot", "block": "center"})

    @reactive.effect
    @reactive.event(input.chase_act)
    def _act():
        """Small actions from the detail panel (clear selections, open a chemotype)."""
        a = input.chase_act().get("act")
        if a == "clear_cpd":
            rv_cpd.set(None)
        elif a == "clear_chemo":
            rv_anchor.set(None)
        elif a == "clear_pair":
            rv_pair.set(None)
        elif a == "open_chemo":
            ch = chemo()
            cid = input.chase_act().get("id")
            if ch and cid in ch.by_id:
                rv_anchor.set(ch.by_id[cid].members[0])
                ui.update_radio_buttons("chemo_tab", selected="act" if ch.by_id[cid].n_hit > 0 else "ina")

    @reactive.effect
    @reactive.event(input.chase_board)
    async def _board():
        """Board item clicked: pin that compound if it is in the loaded data."""
        idx = key_to_idx().get(input.chase_board().get("key"))
        if idx is None:
            ui.notification_show("This compound is not in the loaded data.", duration=3)
            return
        pin_cpd(idx, True)
        await session.send_custom_message("chase-scroll", {"id": "detailcard"})

    @reactive.effect
    @reactive.event(input.chase_plate)
    def _plate():
        """Plate thumbnail clicked: show it large."""
        rv_plate.set(int(input.chase_plate()["pi"]))

    @reactive.effect
    @reactive.event(input.chemo_tab, input.sort, input.min_size, input.incl_single, input.chemo_method, input.sim)
    def _reset_shown():
        """Show the first 24 cards again when tab, sorting or filters change."""
        rv_shown.set(24)

    @reactive.effect
    @reactive.event(input.more)
    def _more():
        """Show 24 more chemotype cards."""
        rv_shown.set(rv_shown.get() + 24)

    @reactive.effect
    @reactive.event(input.cliff_sim, input.cliff_delta, input.cliff_chemo)
    def _reset_cliffs():
        """Show the first 20 cliff pairs again when the cliff filters change."""
        rv_cliff_shown.set(20)

    @reactive.effect
    @reactive.event(input.cliff_more)
    def _cliff_more():
        """Show 20 more cliff pairs."""
        rv_cliff_shown.set(rv_cliff_shown.get() + 20)

    @reactive.effect
    @reactive.event(input.chemo_tab)
    def _sort_choices():
        """Offer the sort options of the active or inactive tab."""
        ch = {"hits": "Hits", "enrich": "Enrichment (p)", "rate": "Hit rate", "max": "Max % response", "mean": "Mean % response", "size": "Size"} \
            if input.chemo_tab() == "act" else {"size": "Size", "mean": "Mean % response"}
        ui.update_select("sort", choices=ch, selected=next(iter(ch)))

    # ------------------------------------------------------------ tooltip + plots
    @reactive.effect
    @reactive.event(input.chase_hover)
    async def _hover():
        """Send the structure tooltip for the hovered point (the browser shows it only if it is still current)."""
        h = input.chase_hover()
        an = rv_an.get()
        if an is None:
            return
        with reactive.isolate():
            # For a point on the replicate plot, show the compound's first well.
            # Calling flags() makes sure the PAINS and covalent chips can be shown.
            wi = int(h["idx"])
            if h["plot"] == "repplot":
                wi = an.cpd_wells[wi][0]
            flags()
            html = rh.tooltip(an, cl(), chemo(), wi, decisions())
        await session.send_custom_message("chase-tip", {"nonce": h["nonce"], "html": html})

    @reactive.effect
    async def _plot():
        """Draw or update the % response plot."""
        an = rv_an.get()
        if an is None:
            await session.send_custom_message("chase-plot", {"id": "plot", "fig": None, "empty": "Set the controls in step 2, then press <b>Normalize &amp; plot</b>."})
            return
        cid = sel_chemo_id()
        ch = chemo().by_id.get(cid) if cid else None
        flags()
        # The figure is built here in Python and sent to the browser as JSON. chase.js draws it with
        # Plotly.react, which keeps the user's zoom.
        fig = plots.scatter(an, cl(), input.show_ctrl(), ch, rv_cpd.get(), rv_pair.get(), show_z=input.show_z(), flag_mode=input.flag_mode())
        await session.send_custom_message("chase-plot", {"id": "plot", "fig": fig, "config": plots.CONFIG})

    @reactive.calc
    def dims():
        """Plate layout of the current data."""
        an = rv_an.get()
        req(an)
        rv_ver.get()
        return plates.plate_dims(an.wells)

    @reactive.effect
    async def _bigmap():
        """Draw the large map of the selected plate."""
        an = rv_an.get()
        if an is None:
            await session.send_custom_message("chase-plot", {"id": "bigmap", "fig": None})
            return
        pi = min(rv_plate.get(), len(an.plates) - 1)
        await session.send_custom_message("chase-plot", {"id": "bigmap", "fig": plots.heatmap(an, pi, dims()), "config": plots.CONFIG})

    @reactive.calc
    def rep():
        """Replicate comparison for the current threshold."""
        an = rv_an.get()
        req(an)
        rv_ver.get()
        return plates.replicates(an, float(input.rep_t() or 20))

    @reactive.effect
    async def _repplot():
        """Draw the replicate plot, or explain why there is none."""
        an = rv_an.get()
        r = rep() if an is not None else None
        if r is None or len(r["pts"]) < 3:
            n = len(r["pts"]) if r else 0
            msg = (f"Only {n} compound{' was' if n == 1 else 's were'} tested in more than one well." if n else "No compound was tested in more than one well.") \
                + " The replicate plot needs a duplicate or replicate screen."
            await session.send_custom_message("chase-plot", {"id": "repplot", "fig": None, "empty": msg if an is not None else ""})
            return
        await session.send_custom_message("chase-plot", {"id": "repplot", "fig": plots.replicate_plot(an, cl(), r), "config": plots.CONFIG})

    # ------------------------------------------------------------ outputs
    # Some outputs call flags() first, so the PAINS and covalent columns exist before the HTML is built.
    @render.ui
    def banner():
        """Banner shown while the example data are loaded."""
        if not rv_demo.get():
            return None
        return ui.div(ui.HTML("<b>Example data</b> &nbsp;The plot shows a generated two-plate screen so you can try the tools. Load your screening file and structures in step 1 to analyze your own data."), class_="banner")

    @render.ui
    def errbox():
        """Error message box."""
        e = rv_err.get()
        return ui.div(e, class_="alert-chase") if e else None

    @render.ui
    def screen_hint():
        """Summary line under the screening-file input."""
        sf = rv_screen.get()
        return ui.p(getattr(sf, "label", sf.file_name), class_="filehint") if sf else None

    @render.ui
    def struct_hint():
        """Summary line under the structures-file input."""
        stf = rv_struct.get()
        return ui.p(getattr(stf, "label", stf.file_name), class_="filehint") if stf else None

    @render.ui
    def ctrl_hint():
        """Suggested control IDs."""
        sf = rv_screen.get()
        return ui.p(getattr(sf, "ctrl_hint", "Pick Batch IDs that appear on every plate."), class_="hint") if sf else ui.p("Pick Batch IDs that appear on every plate.", class_="hint")

    @render.ui
    def cls_hint():
        """Explanation of the current class borders."""
        return ui.p(cl().hint, class_="hint") if rv_an.get() is not None else None

    @render.ui
    def chemo_hint():
        """Explanation of the chemotype method."""
        t = "Butina clustering on RDKit Morgan fingerprints (radius 3, 2048 bits)." if input.chemo_method() == "butina" else "RDKit Bemis-Murcko scaffolds: ring systems and linkers, side chains removed."
        ch = chemo() if rv_an.get() is not None else None
        extra = f" {ch.note}" if ch and ch.note else ""
        return ui.p(t + extra, class_="hint")

    @render.ui
    def stats():
        """Summary tiles."""
        an = rv_an.get()
        if an is None:
            return None
        flags()
        return ui.HTML(rh.stats_strip(an, cl(), chemo()))

    @render.ui
    def plot_meta():
        """Number of wells and plates next to the plot title."""
        an = rv_an.get()
        if an is None:
            return None
        return ui.span(f"{int((an.wells['type'] == 'sample').sum()):,} sample wells | {len(an.plates)} plates", class_="meta")

    def chemo_list():
        """Chemotypes for the current tab, filtered by size and sorted."""
        ch = chemo()
        if not ch:
            return []
        tab = input.chemo_tab()
        ms = max(1, int(input.min_size() or 1))
        incl = input.incl_single()
        # The 'Active' tab shows chemotypes with at least one hit (single compounds only if asked for).
        # The 'Inactive' tab shows chemotypes without any hit.
        if tab == "act":
            lst = [c for c in ch.chemos if c.n_hit > 0 and (c.n >= ms or (incl and c.n == 1))]
        else:
            lst = [c for c in ch.chemos if c.n_hit == 0 and c.n >= ms]
        key = input.sort() or ("hits" if tab == "act" else "size")
        # Ways to sort the cards. When two cards are equal, more hits and then a higher maximum come first.
        sorters = {
            "hits": lambda c: (-c.n_top, -c.n_hit, -c.max), "enrich": lambda c: (c.p, -c.n_hit), "rate": lambda c: (-c.rate, -c.n_hit, -c.max),
            "max": lambda c: -c.max, "mean": (lambda c: -c.mean) if tab == "act" else (lambda c: c.mean), "size": lambda c: (-c.n, -c.n_hit)}
        return sorted(lst, key=sorters.get(key, sorters["size"]))

    @render.ui
    def chemo_counts():
        """Number of active and inactive chemotypes next to the tabs."""
        ch = chemo() if rv_an.get() is not None else None
        if not ch:
            return None
        ms = max(1, int(input.min_size() or 1))
        a = sum(1 for c in ch.chemos if c.n_hit > 0 and (c.n >= ms or (input.incl_single() and c.n == 1)))
        b = sum(1 for c in ch.chemos if c.n_hit == 0 and c.n >= ms)
        return ui.span(f"{a:,} active | {b:,} inactive", class_="meta")

    @render.ui
    def chemo_grid():
        """Chemotype cards (or an explanation when there are none)."""
        an = rv_an.get()
        if an is None:
            return None
        if not (an.cpds["smiles"] != "").any():
            msg = "None of the Batch IDs in the screening file match the structures file. Check the 'Join on' column in step 1." if rv_struct.get() else \
                "Load a structures file (Batch ID + SMILES) in step 1 to group compounds into chemotypes."
            return ui.div(msg, class_="empty")
        lst = chemo_list()
        if not lst:
            return ui.div("No chemotype contains an active compound with these settings." if input.chemo_tab() == "act" else "No fully inactive chemotype reaches the minimum size.", class_="empty")
        # Only the first rv_shown cards are drawn. 'Show more' adds another 24.
        n = rv_shown.get()
        flags()
        html = rh.chemo_cards(lst[:n], input.chemo_tab(), sel_chemo_id(), an.cpds["mean"].to_numpy(), cl(), decisions(), an.cpds["key"].tolist(), an.cpds)
        more = f'<p class="hint text-center">{len(lst) - n:,} more below the button</p>' if len(lst) > n else ""
        return ui.HTML(f'<div class="cgrid">{html}</div>{more}')

    @render.ui
    def detail():
        """Detail panel."""
        an = rv_an.get()
        if an is None:
            return ui.div("Results appear here after normalization.", class_="empty")
        flags()
        return ui.HTML(rh.detail(an, cl(), chemo() if (an.cpds["smiles"] != "").any() else None, rv_cpd.get(), sel_chemo_id(), rv_pair.get(), decisions()))

    @render.ui
    def tri_meta():
        """Number of decisions next to the board title."""
        n = len(decisions())
        return ui.span(f"{n} decision{'s' if n != 1 else ''}" if n else "", class_="meta")

    @render.ui
    def tri_notice():
        """Saved / unsaved status and how to keep decisions."""
        decisions()
        n = tboard.unsaved
        msg = ("Decisions live only in this browser tab; CHASE writes nothing to disk. Keep them with <b>Save session</b> or <b>Export triage (.csv)</b>, "
               "and bring them back with <b>Open session</b> or <b>Import triage (.csv)</b>.")
        if n:
            state = f'<span class="chip warn">{n} unsaved change{"s" if n != 1 else ""}</span>'
        else:
            state = '<span class="chip good">saved</span>' if len(tboard) else ""
        extra = rv_note.get()
        return ui.div(ui.HTML(f"{state} {msg}" + (f" <b>{rh.esc(extra)}</b>" if extra else "")), class_="trinotice")

    @render.ui
    def board():
        """Triage board."""
        an = rv_an.get()
        if an is not None:
            flags()
        return ui.HTML(rh.board(decisions(), an, key_to_idx()))

    def cliff_items():
        """Activity-cliff pairs for the current filters."""
        ch = chemo()
        an = rv_an.get()
        only = None
        # 'Selected chemotype only': keep pairs where at least one of the two compounds is in the selected chemotype.
        if input.cliff_chemo() and sel_chemo_id():
            only = set(ch.by_id[sel_chemo_id()].members)
        return cstats.cliff_list(ch, an.cpds["mean"].to_numpy(), float(input.cliff_sim()), float(input.cliff_delta() or 0), only)

    @render.ui
    def cliff_meta():
        """Number of cliff pairs."""
        an = rv_an.get()
        ch = chemo() if an is not None else None
        if not ch or ch.cliffs is None:
            return None
        return ui.span(f"{len(cliff_items()):,} pairs", class_="meta")

    @render.ui
    def cliffs():
        """Activity-cliff list."""
        an = rv_an.get()
        if an is None:
            return None
        if not (an.cpds["smiles"] != "").any():
            return ui.div("Load a structures file to find activity cliffs.", class_="empty")
        ch = chemo()
        if ch is None or ch.cliffs is None:
            return ui.div(ch.cliff_note if ch else "Activity cliffs are not available for this data set.", class_="empty")
        flags()
        lst = cliff_items()
        if not lst:
            return ui.div("No pairs reach these similarity and Delta thresholds. Lower either one to see more.", class_="empty")
        return ui.HTML('<div class="cliffs">' + rh.cliffs(lst, an, cl(), ch, rv_pair.get(), rv_cliff_shown.get()) + "</div>")

    @render.ui
    def minimaps():
        """Plate thumbnails."""
        an = rv_an.get()
        if an is None:
            return None
        return ui.HTML(rh.mini_maps(an, dims(), min(rv_plate.get(), len(an.plates) - 1)))

    @render.ui
    def bighead():
        """Header of the large plate map: Z', S/B and edge versus interior."""
        an = rv_an.get()
        if an is None:
            return None
        pi = min(rv_plate.get(), len(an.plates) - 1)
        p = an.plates[pi]
        es = plates.edge_stats(an.wells, pi, dims())
        edge = f" &nbsp;|&nbsp; Edge wells {es['edge']:.1f}% vs interior {es['inner']:.1f}% (Delta {es['edge'] - es['inner']:+.1f} pts)" if es["nE"] and es["nI"] else ""
        return ui.HTML(f"<b>Plate {rh.esc(p['label'])}</b> &nbsp;{p['nSample']} samples | Z' {rh.fmt(p['z'], 2)} | S/B {rh.fmt(p['sb'], 1)}{edge}")

    @render.ui
    def rep_stats():
        """Replicate statistics tiles."""
        an = rv_an.get()
        if an is None:
            return None
        r = rep()
        if len(r["pts"]) < 3:
            return None
        items = [("Compounds", f"{len(r['pts']):,}"), ("Pearson r", f"{r['r']:.3f}"), ("Spearman rho", f"{r['rho']:.3f}"), ("Median |Delta|", f"{r['med']:.1f} pts"), (f"|Delta| > {r['T']:.0f} pts", f"{len(r['disc']):,}")]
        return ui.HTML('<div class="stats repstats">' + "".join(f'<div class="stat"><span class="k">{k}</span><span class="v">{v}</span></div>' for k, v in items) + "</div>")

    @render.ui
    def rep_list():
        """Compounds whose replicates disagree most."""
        an = rv_an.get()
        if an is None:
            return None
        r = rep()
        if len(r["pts"]) < 3 or not r["disc"]:
            return None
        rows_html = []
        for p in r["disc"][:30]:
            c = an.cpds.loc[p["ci"]]
            ws = an.wells.loc[an.cpd_wells[p["ci"]][:2]]
            loc = " / ".join(f"P{rh.esc(w.plate)}:{rh.esc(w.well)}" for w in ws.itertuples())
            rows_html.append(f'<button class="mrow" type="button" data-chase="sel_cpd" data-cpd="{p["ci"]}">{rh.tile(c["smiles"])}<span class="who"><b>{rh.esc(c["bid"])}</b><span>{loc}</span></span>'
                             f'<span class="pct mono">{p["a"]:.1f} / {p["b"]:.1f}%<span class="sm">Delta {p["d"]:.1f}</span></span></button>')
        return ui.HTML('<div class="replist"><h4>Largest disagreements</h4><div class="rlist">' + "".join(rows_html) + "</div></div>")

    @render.ui
    def qc():
        """Plate QC table."""
        an = rv_an.get()
        if an is None:
            return None
        return ui.HTML(rh.qc_table(an, cl()))

    @render.ui
    def notes():
        """Data checks: empty wells, missing structures, Excel-mangled IDs, unreadable SMILES, flags, replicates."""
        an = rv_an.get()
        if an is None:
            return None
        out = []
        if an.n_empty:
            out.append(f"{an.n_empty:,} wells have a signal but no Batch ID. They are left out of the analysis.")
        nos = an.cpds[an.cpds["smiles"] == ""]
        if rv_struct.get() is not None and len(nos):
            ids = " ".join(f"<code>{rh.esc(b)}</code>" for b in nos["bid"].head(12))
            out.append(f"{len(nos):,} compounds have no SMILES in the structures file: {ids}{' ...' if len(nos) > 12 else ''}")
            # An ID like 274000000 is usually '274E06' that Excel has read as a number in scientific notation.
            sci = [b for b in nos["bid"] if re.fullmatch(r"\d{3,}0{6,}", b)]
            if sci:
                out.append(" ".join(f"<code>{rh.esc(b)}</code>" for b in sci) + " look like IDs that Excel turned into numbers (an ID such as <code>274E06</code> is read as 274 x 10^6). "
                           "Format the Batch ID column as text in the source workbook to keep them.")
        pr = prep()
        if pr is not None and pr.errors:
            out.append(f"{len(pr.errors):,} SMILES could not be read and are not clustered: " + " ".join(f"<code>{rh.esc(an.cpds.at[i, 'bid'])}</code>" for i in list(pr.errors)[:6]))
        r = rep()
        if len(r["pts"]) >= 3 and r["disc"]:
            out.append(f"{len(r['disc']):,} of {len(r['pts']):,} replicated compounds differ by more than {r['T']:.0f} points between replicates. See Plates -> Replicates.")
        if (an.cpds["smiles"] != "").any():
            flags()
            hit = cl().cls > 0
            npn, ncv = int((an.cpds["pains"] != "").sum()), int((an.cpds["covalent"] != "").sum())
            hpn, hcv = int(((an.cpds["pains"] != "") & hit).sum()), int(((an.cpds["covalent"] != "") & hit).sum())
            out.append(f"{npn:,} compounds match a PAINS alert ({hpn} of them hits) and {ncv:,} carry a covalent warhead ({hcv} hits). "
                       "Flagged hits are marked on the plot and in every list; the compound panel names the alerts.")
            bad = alerts.covalent_alerts()[1]
            if bad:
                out.append("Some lines of the covalent alert file could not be read and were skipped: " + ", ".join(f"<code>{rh.esc(b)}</code>" for b in bad[:5]))
        reps = int((an.cpds["n"] > 1).sum())
        if reps:
            out.append(f"{reps:,} compounds were tested in more than one well. Classes and chemotypes use the mean of their wells; the plot shows every well.")
        bor = [p for p in an.plates if p["borrowed"]]
        if bor:
            out.append("Plates " + ", ".join(f"P{rh.esc(p['label'])}" for p in bor) + " lack one of the controls and were normalized with all-plate control values.")
        ch = chemo()
        if ch and ch.note:
            out.append(rh.esc(ch.note))
        if not out:
            out.append("No problems found in the loaded files.")
        return ui.HTML('<div class="notes">' + "".join(f"<div>{n}</div>" for n in out) + "</div>")

    # ------------------------------------------------------------ exports
    def base_name():
        """Base for download file names, from the screening file name."""
        an = rv_an.get()
        n = an.file_name if an is not None else "screen"
        return re.sub(r"[^\w.-]+", "_", re.sub(r"\.[^.]+$", "", n.split(" (")[0]))

    @render.download(filename=lambda: f"{base_name()}_results.csv")
    def export_csv():
        """Download one row per compound with class, chemotype statistics and flags."""
        an = rv_an.get()
        req(an)
        c, ch = cl(), chemo()
        buf = _io.StringIO()
        w = csv.writer(buf)
        flags()
        w.writerow(["batch_id", "compound_name", "smiles", "n_wells", "wells", "pct_response_mean", "pct_response_sd", "activity_class", "chemotype",
                    "chemotype_size", "chemotype_hits", "chemotype_hit_rate", "chemotype_EF", "chemotype_p", "chemotype_q", "pains_alerts", "covalent_alerts"])
        # One row per compound, starting with the most active.
        for i in np.argsort(-an.cpds["mean"].to_numpy()):
            r = an.cpds.loc[i]
            cid = ch.cpd_chemo.get(int(i)) if ch else None
            t = ch.by_id[cid] if cid else None
            wl = " ".join(f"P{x.plate}:{x.well}" for x in an.wells.loc[an.cpd_wells[i]].itertuples())
            w.writerow([r["bid"], r["name"], r["smiles"], r["n"], wl, f"{r['mean']:.2f}", "" if not math.isfinite(r["sd"]) else f"{r['sd']:.2f}", c.labels[c.cls[i]][0],
                        cid or "", t.n if t else "", t.n_hit if t else "", f"{t.rate * 100:.1f}" if t else "", f"{t.ef:.2f}" if t and math.isfinite(t.ef) else "",
                        f"{t.p:.3g}" if t else "", f"{t.q:.3g}" if t else "", r["pains"], r["covalent"]])
        yield buf.getvalue()

    @render.download(filename=lambda: f"{safe_name(camp['name'] or 'campaign')}_triage.csv")
    def export_triage():
        """Download the triage board as CSV."""
        yield tri.to_csv(decisions(), camp["name"])

    @render.download(filename=lambda: f"{base_name()}_CHASE-session_{date.today().isoformat()}.json")
    def save_session():
        """Download the whole session (data, settings, selection, decisions) as JSON."""
        sf = rv_screen.get()
        req(sf)
        an = rv_an.get()
        rws = rows()
        stf = rv_struct.get()
        # The setting names are written in camelCase, like in the HTML version, so session files work in both.
        settings = dict(pos=input.pos(), neg=input.neg(), scope=input.scope(), stat=input.stat(), clsMethod=input.cls_method(), k=int(input.k() or 3),
                        thrA=input.thr_a(), thrW=input.thr_w(), robN=input.rob_n(), chemoMethod=input.chemo_method(), sim=input.sim(), minSize=input.min_size(),
                        inclSingle=input.incl_single(), tab=input.chemo_tab(), showCtrl=input.show_ctrl(), cliffSim=input.cliff_sim(), cliffDelta=input.cliff_delta(),
                        repT=input.rep_t(), campaign=camp["name"], plateSel=rv_plate.get(), showZ=input.show_z(), colors=rv_colors.get(),
                        flagMode=input.flag_mode())
        # The selection is saved by compound ID, because compound numbers can change when the session is opened again.
        sel = {"cpd": an.cpds.at[rv_cpd.get(), "key"] if (an is not None and rv_cpd.get() is not None) else None,
               "chemoRepOf": an.cpds.at[rv_anchor.get(), "key"] if (an is not None and rv_anchor.get() is not None) else None}
        yield sessmod.build(__version__, an.file_name if an is not None else sf.file_name, rws[["raw", "plate", "well", "bid"]].values.tolist(),
                            stf.file_name if stf else None, stf.rows() if stf else None, settings, an is not None, sel, camp["name"], list(decisions().values()))

    @reactive.effect
    @reactive.event(input.session_file)
    def _open_session():
        """Restore a session file: data, settings, selection and triage decisions."""
        f = input.session_file()[0]
        try:
            d = sessmod.parse(consume(f, lambda p: p.read_text(encoding="utf-8")))
        except sessmod.SessionError as e:
            rv_err.set(str(e))
            return
        st = d.get("settings") or {}
        rv_demo.set(False)
        rv_err.set("")
        name = d["screen"].get("fileName") or f["name"]
        # Build the screening and structure tables again from the rows saved in the session, with standard headers.
        sf = io.ScreenFile.from_sheets({"Session": [["Raw Data", "Plate Number", "384 well ID", "Batch ID"]] + d["screen"]["rows"]}, name)
        rows_df = set_screen(sf, f"{name} (from session)")
        smap = None
        if d.get("structures") and d["structures"].get("rows"):
            stf = io.StructFile.from_sheets({"s": [["Batch id", "Compound id", "SMILES"]] + d["structures"]["rows"]}, d["structures"].get("fileName") or "structures", set(rows_df["key"]))
            set_struct(stf, f"{d['structures'].get('fileName') or 'structures'} (from session)")
            smap = stf.mapping()
        else:
            rv_struct.set(None)

        def upd(fn, iid, key, cast=lambda x: x, **kw):
            """Update input ``iid`` from session setting ``key`` if present (``arg`` names the update argument)."""
            if st.get(key) not in (None, ""):
                fn(iid, **{kw.get("arg", "selected"): cast(st[key])})
        # Put every saved setting back in its place in the sidebar.
        ui.update_selectize("pos", selected=st.get("pos"))
        ui.update_selectize("neg", selected=st.get("neg"))
        upd(ui.update_radio_buttons, "scope", "scope")
        upd(ui.update_radio_buttons, "stat", "stat")
        upd(ui.update_select, "cls_method", "clsMethod")
        upd(ui.update_radio_buttons, "k", "k", str)
        upd(ui.update_numeric, "thr_a", "thrA", float, arg="value")
        upd(ui.update_numeric, "thr_w", "thrW", float, arg="value")
        upd(ui.update_numeric, "rob_n", "robN", float, arg="value")
        upd(ui.update_radio_buttons, "chemo_method", "chemoMethod")
        upd(ui.update_slider, "sim", "sim", float, arg="value")
        upd(ui.update_numeric, "min_size", "minSize", int, arg="value")
        if "inclSingle" in st:
            ui.update_checkbox("incl_single", value=bool(st["inclSingle"]))
        if "showCtrl" in st:
            ui.update_checkbox("show_ctrl", value=bool(st["showCtrl"]))
        upd(ui.update_select, "flag_mode", "flagMode")
        if "showZ" in st:
            ui.update_checkbox("show_z", value=bool(st["showZ"]))
        # Colours: use them here, and also send them to the browser (see _push_session_colors).
        if isinstance(st.get("colors"), dict):
            c = dict(plots.DEFAULT_COLORS)
            c.update({k: v for k, v in st["colors"].items() if k in c and isinstance(v, str) and HEX.match(v)})
            rv_colors.set(c)
            session_colors["pending"] = c
        upd(ui.update_slider, "cliff_sim", "cliffSim", float, arg="value")
        upd(ui.update_numeric, "cliff_delta", "cliffDelta", float, arg="value")
        upd(ui.update_numeric, "rep_t", "repT", float, arg="value")
        upd(ui.update_radio_buttons, "chemo_tab", "tab")
        campaign = st.get("campaign") or (d.get("triage") or {}).get("campaign") or re.sub(r"\.[^.]+$", "", name)
        camp.update(auto=False, name=campaign, last_auto="")
        ui.update_text("campaign", value=campaign)
        # Run the analysis again with the saved settings. _restore_selection then puts the selection back.
        if d.get("analyzed", True):
            do_analyze(rows_df, st.get("pos", ""), st.get("neg", ""), st.get("scope", "plate"), st.get("stat", "mean"), smap, name)
            rv_plate.set(int(st.get("plateSel") or 0))
            rv_pending.set(d.get("selection") or {})
        # The decisions from the session are added to what is already on the board.
        # For each compound the newer decision wins.
        was_empty = not len(tboard)
        added = tboard.merge((d.get("triage") or {}).get("decisions") or [])
        if was_empty:
            tboard.mark_saved()
        bump(rv_bump)
        rv_note.set(f"{added} decision{'s' if added != 1 else ''} loaded from the session." if added else "")
        ui.notification_show("Session loaded", duration=3)

    session_colors: dict = {}

    @reactive.effect
    async def _push_session_colors():
        """Send colours restored from a session to the browser."""
        rv_colors.get()
        c = session_colors.pop("pending", None)
        if c:
            await push_colors(c)

    @render.download(filename=lambda: f"{base_name()}_plot-data.csv")
    def plot_csv():
        """Download every point of the % response plot, controls included."""
        an = rv_an.get()
        req(an)
        c, ch = cl(), chemo()
        flags()
        zmap = {p["label"]: p["z"] for p in an.plates}
        kind = {"sample": "sample", "pos": "positive control", "neg": "negative control"}
        buf = _io.StringIO()
        w = csv.writer(buf)
        w.writerow(["compound_number", "plot_x", "plate", "well", "well_type", "batch_id", "compound_name", "raw", "pct_response",
                    "compound_mean_pct", "n_wells", "activity_class", "chemotype", "plate_z_prime", "pains_alerts", "covalent_alerts", "smiles"])
        # Rows in the order they appear on the plot. Compound wells get the compound columns filled in;
        # for control wells those columns stay empty.
        W = an.wells.sort_values("x", kind="stable")
        for wi, r in W.iterrows():
            z = zmap.get(r["plate"], float("nan"))
            zs = f"{z:.3f}" if math.isfinite(z) else ""
            if r["type"] == "sample":
                ci = int(r["cpd"])
                cp = an.cpds.loc[ci]
                w.writerow([int(r["x"]), f"{r['x']:g}", r["plate"], r["well"], kind[r["type"]], r["bid"], cp["name"], f"{r['raw']:g}", f"{r['pct']:.2f}",
                            f"{cp['mean']:.2f}", cp["n"], c.labels[c.cls[ci]][0], (ch.cpd_chemo.get(ci, "") if ch else ""), zs, cp["pains"], cp["covalent"], cp["smiles"]])
            else:
                w.writerow(["", f"{r['x']:g}", r["plate"], r["well"], kind[r["type"]], r["bid"], "", f"{r['raw']:g}", f"{r['pct']:.2f}", "", "", "", "", zs, "", "", ""])
        yield buf.getvalue()

    @reactive.effect
    def _restore_selection():
        """After a session is opened and analysed, re-select its pinned compound and chemotype."""
        p = rv_pending.get()
        if not p:
            return
        an = rv_an.get()
        if an is None:
            return
        k2i = key_to_idx()
        with reactive.isolate():
            rv_pending.set(None)
            if p.get("chemoRepOf") in k2i:
                rv_anchor.set(k2i[p["chemoRepOf"]])
            if p.get("cpd") in k2i:
                rv_cpd.set(k2i[p["cpd"]])


# ================================================================ app
# Serve /www (CSS and JavaScript) and /plotly (plotly.min.js from the plotly package) as static files.
app = App(app_ui, server, static_assets={"/www": str(WWW), "/plotly": str(PLOTLY_JS)})
