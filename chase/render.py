"""HTML snippets for the CHASE interface (rendered server-side, sent as ui.HTML).

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import math
from datetime import datetime
from html import escape as _esc

import numpy as np

from .core.chem import svg
from .core.triage import STATUS, STATUSES
from .plots import C, col, heat_rgb
from .core.plates import plate_grid


def esc(s) -> str:
    """HTML-escape any value (None becomes an empty string)."""
    return _esc("" if s is None else str(s), quote=True)


def fmt(x, d=1) -> str:
    """Format a number with ``d`` decimals; '-' for missing or non-finite values."""
    try:
        return f"{x:.{d}f}" if x is not None and math.isfinite(x) else "-"
    except (TypeError, ValueError):
        return "-"


def fint(x) -> str:
    """Format an integer with thousands separators."""
    return f"{int(x):,}"


def fmt_p(p) -> str:
    """Format a p- or q-value: three decimals, or scientific notation (e.g. 1.2e-05) below 0.001."""
    if p is None or not math.isfinite(p):
        return "-"
    if p >= 0.001:
        return f"{p:.3f}"
    return f"{p:.1e}"


def flags_of(cp, i) -> tuple[str, str]:
    """Return the PAINS and covalent alert names of compound ``i`` (empty strings if not flagged or not yet checked)."""
    # The 'pains' and 'covalent' columns hold the names of the alerts, separated by '; '.
    # They only exist once the alerts have been calculated.
    if i is None or "pains" not in cp.columns:
        return "", ""
    return cp.at[i, "pains"] or "", cp.at[i, "covalent"] or ""


def flag_chips(cp, i) -> str:
    """Return the red PAINS and violet Covalent chips for compound ``i`` (alert names in the tooltip)."""
    p, c = flags_of(cp, i)
    return ((f'<span class="chip bad" title="PAINS: {esc(p)}">PAINS</span>' if p else "")
            + (f'<span class="chip cov" title="Covalent warhead: {esc(c)}">Covalent</span>' if c else ""))


def tile_hl(smiles: str) -> str:
    """Structure with PAINS atoms (red) and covalent-warhead atoms (violet) highlighted."""
    from .core.alerts import alert_atoms
    pa, ca = alert_atoms(smiles or "")
    return f'<div class="tilebox">{svg(smiles or "", 280, 210, pa, ca)}</div>'


def tile(smiles: str) -> str:
    """Return a structure tile (SVG in a light box) for ``smiles``."""
    return f'<div class="tilebox">{svg(smiles or "")}</div>'


def cls_chip(cl, ci: int) -> str:
    """Return a chip with the colour dot and name of activity class ``ci``."""
    label, tok = cl.labels[ci]
    return f'<span class="chip"><i style="background:{col(cl, tok)}"></i>{label}</span>'


def tri_chip(d) -> str:
    """Return a chip showing the triage decision ``d``, or an empty string."""
    if not d or d["status"] not in STATUS:
        return ""
    _, label, cls = STATUS[d["status"]]
    return f'<span class="chip {cls}">{label}</span>'


def when(t) -> str:
    """Format a millisecond timestamp as '05 Nov 14:32'."""
    try:
        return datetime.fromtimestamp(float(t) / 1000).strftime("%d %b %H:%M")
    except (TypeError, ValueError, OSError):
        return ""


def meta(d) -> str:
    """Return 'by <name> | <time>' for a triage decision."""
    if not d:
        return ""
    who = f"by {esc(d.get('by'))} | " if d.get("by") else ""
    return f'<span class="tmeta">{who}{when(d.get("at"))}</span>'


# ---------------------------------------------------------------- stats strip
def stats_strip(an, cl, chemo) -> str:
    """Return the row of summary tiles (compounds, plates, structures, actives, any activity, chemotypes, flagged hits)."""
    cp = an.cpds
    n_s = int((cp["smiles"] != "").sum())
    # Median Z' over the plates where Z' could be calculated.
    zs = [p["z"] for p in an.plates if math.isfinite(p["z"])]
    zmed = float(np.median(zs)) if zs else math.nan
    K = len(cl.labels)
    hits = int((cl.cls > 0).sum())
    top = cl.counts[-1]
    n_w = int((an.wells["type"] == "sample").sum())
    # For each tile: (title, big number, small text underneath).
    items = [
        ("Compounds", fint(len(cp)), f"{n_w:,} wells"),
        ("Plates", fint(len(an.plates)), f"median Z' {zmed:.2f}" if math.isfinite(zmed) else "no Z'"),
        ("Structures", fint(n_s), ("100" if n_s == len(cp) else fmt(100 * n_s / max(1, len(cp)), 1)) + "% matched"),
        (cl.labels[-1][0], fint(top), f"{100 * top / max(1, len(cp)):.2f}% of compounds"),
        ("Any activity", fint(hits), " | ".join(f"{cl.counts[j]} {l.lower()}" for j, (l, _) in enumerate(cl.labels) if j > 0) if K > 2 else "non-inactive"),
        ("Chemotypes", fint(len(chemo.chemos)) if chemo else ("..." if n_s else "-"),
         f"{sum(1 for c in chemo.chemos if c.n_hit > 0):,} with hits" if chemo else ("clustering" if n_s else "needs structures")),
    ]
    # Flagged hits: compounds above Inactive with a PAINS or covalent alert. A compound with both counts once.
    if "pains" in cp.columns:
        hit = cl.cls > 0
        hp, hc = int(((cp["pains"] != "") & hit).sum()), int(((cp["covalent"] != "") & hit).sum())
        hf = int((((cp["pains"] != "") | (cp["covalent"] != "")) & hit).sum())
        items.append(("Flagged hits", fint(hf), f"{hp} PAINS | {hc} covalent"))
    return '<div class="stats">' + "".join(f'<div class="stat"><span class="k">{esc(k)}</span><span class="v">{v}</span><span class="s">{esc(s)}</span></div>' for k, v, s in items) + "</div>"


# ---------------------------------------------------------------- chemotype cards
def strip(ch, means, cl) -> str:
    """Return a small SVG strip with one dot per chemotype member at its % response (-20 to 120, ticks at 0 and 100)."""
    w, h, lo, hi = 180, 18, -20, 120
    # Turn a % response into a position along the strip, keeping it inside the strip.
    # At most 200 dots are drawn.
    X = lambda v: 4 + ((max(lo, min(hi, v)) - lo) / (hi - lo)) * (w - 8)
    dots = "".join(f'<circle cx="{X(means[i]):.1f}" cy="9" r="3.2" fill="{col(cl, cl.labels[cl.cls[i]][1])}" fill-opacity=".85" stroke="#fff" stroke-width=".8"/>' for i in ch.members[:200])
    return (f'<svg class="strip" viewBox="0 0 {w} {h}" aria-label="Member % response from -20 to 120">'
            f'<line x1="4" x2="{w - 4}" y1="9" y2="9" stroke="{C["line"]}"/>'
            f'<line x1="{X(0):.1f}" x2="{X(0):.1f}" y1="3" y2="15" stroke="{C["ink3"]}" stroke-width=".8"/>'
            f'<line x1="{X(100):.1f}" x2="{X(100):.1f}" y1="3" y2="15" stroke="{C["ink3"]}" stroke-width=".8"/>{dots}</svg>')


def chemo_cards(lst, tab: str, sel_id, means, cl, decisions, keys, cp=None) -> str:
    """Return the chemotype cards for the Active or Inactive tab.

    Args:
        lst: Chemotypes to show, already filtered and sorted.
        tab: ``"act"`` or ``"ina"`` (active cards also show EF, p and the enriched chip).
        sel_id: ID of the selected chemotype.
        means: Mean % response per compound.
        cl: Activity classes.
        decisions: Triage decisions (for the 'n triaged' chip).
        keys: Compound keys by index.
        cp: Compound table, for PAINS/covalent counts.
    """
    out = []
    for ch in lst:
        # Cards in the 'Active' tab show hits, hit rate and maximum, and for chemotypes of two or more
        # compounds also the enrichment (EF and p, plus an 'enriched' chip when q < 0.05).
        # Cards in the 'Inactive' tab show the mean and maximum response.
        if tab == "act":
            l2 = f'<span><b>{ch.n_hit}</b>/{ch.n} hits</span><span><b>{ch.rate * 100:.0f}</b>%</span><span>max <b>{ch.max:.0f}</b>%</span>'
            enr = (f'<span class="l2 enr"><span>EF <b>{fmt(ch.ef, 1)}</b></span><span>p <b>{fmt_p(ch.p)}</b></span>'
                   + ('<span class="chip good">enriched</span>' if ch.q < 0.05 else "") + "</span>") if ch.n > 1 else ""
        else:
            l2 = f'<span>mean <b>{ch.mean:.1f}</b>%</span><span>max <b>{ch.max:.0f}</b>%</span>'
            enr = ""
        # Chips on the card: how many members already have a triage decision, and how many have PAINS or covalent alerts.
        tri = sum(1 for i in ch.members if keys[i] in decisions)
        fl = ""
        if cp is not None and "pains" in cp.columns:
            npn = sum(1 for i in ch.members if cp.at[i, "pains"])
            ncv = sum(1 for i in ch.members if cp.at[i, "covalent"])
            fl = (f'<span class="chip bad">{npn} PAINS</span>' if npn else "") + (f'<span class="chip cov">{ncv} covalent</span>' if ncv else "")
        # Each card is a button. chase.js reads data-chase and data-id and tells the server it was clicked.
        out.append(
            f'<button class="chemo" type="button" data-chase="sel_chemo" data-id="{ch.id}" aria-pressed="{str(sel_id == ch.id).lower()}">'
            f'{tile(ch.rep)}<span class="l1"><span class="id">{ch.id}</span><span class="n">{"singleton" if ch.n == 1 else f"{ch.n} compounds"}</span></span>'
            f'<span class="l2">{l2}</span>{enr}' + (f'<span class="l2">' + (f'<span class="chip">{tri} triaged</span>' if tri else "") + fl + "</span>" if (tri or fl) else "")
            + f"{strip(ch, means, cl)}</button>")
    return "".join(out)


# ---------------------------------------------------------------- detail panel
def triage_controls(key: str, d) -> str:
    """Return the decision buttons, note box and 'by ... at ...' line for compound ``key``."""
    # One button per decision, with the current decision pressed. The note box only becomes active once
    # there is a decision, because the note is saved together with it.
    btns = "".join(f'<button type="button" data-chase="tri" data-key="{esc(key)}" data-status="{sid}" aria-pressed="{str(bool(d and d["status"] == sid)).lower()}">{label}</button>' for sid, label, _ in STATUSES)
    return (f'<div class="trictl"><span class="lbl">Triage</span><div class="seg tri">{btns}</div>'
            f'<textarea id="tri_note" data-key="{esc(key)}" rows="2" placeholder="{"Add a note for the team" if d else "Pick a decision, then add a note"}" {"" if d else "disabled"}>{esc(d["note"] if d else "")}</textarea>'
            f'<div class="trirow">{meta(d) if d else "<span class=hint>No decision yet.</span>"}'
            + (f'<button class="btn btn-sm btn-outline-secondary" type="button" data-chase="tri" data-key="{esc(key)}" data-status="clear">Clear decision</button>' if d else "")
            + "</div></div>")


def detail(an, cl, chemo, sel_cpd, sel_chemo, sel_pair, decisions) -> str:
    """Return the detail panel: the selected cliff pair, the pinned compound and the selected chemotype with its members."""
    cp = an.cpds
    means = cp["mean"].to_numpy()
    parts = []
    # 1. The selected activity-cliff pair, side by side.
    if sel_pair:
        hi, lo = sel_pair
        # Find the pair in the cliff list to get its Tanimoto similarity.
        p = next((x for x in (chemo.cliffs or []) if {x[0], x[1]} == {hi, lo}), None) if chemo else None
        same = same_cpd(cp, hi, lo)

        def half(i):
            """One side of the activity-cliff pair view (structure, IDs, % response, class)."""
            return (f'<button class="pairside" type="button" data-chase="sel_cpd" data-cpd="{i}">{tile(cp.at[i, "smiles"])}'
                    f'<span class="who"><b>{esc(cp.at[i, "bid"])}</b><span>{esc(cp.at[i, "name"])}</span></span><span class="mono">{fmt(means[i])}%</span>{cls_chip(cl, cl.cls[i])}</button>')
        d = abs(means[hi] - means[lo])
        parts.append('<div class="dsec"><h4>Activity-cliff pair<span class="sp"></span><button class="btn btn-sm btn-outline-secondary" data-chase="act" data-act="clear_pair">Clear</button></h4>'
                     f'<div class="pair">{half(hi)}{half(lo)}</div>')
        if p:
            parts.append(f'<dl class="kv"><dt>Tanimoto</dt><dd>{p[2]:.2f}</dd><dt>Delta % response</dt><dd>{d:.1f} pts</dd><dt>SALI</dt><dd>{d / max(0.01, 1 - p[2]):.0f}</dd></dl>')
        if same:
            parts.append('<p class="hint"><b>Same compound, different batch.</b> A difference this large points to a batch problem (purity, identity or stock solution) rather than SAR.</p>')
        parts.append('<p class="hint">Click either structure to pin it and set a triage decision.</p></div>')
    # 2. The pinned compound: structure (with alert atoms highlighted if it is flagged), alerts, triage buttons, its wells and its data.
    if sel_cpd is not None:
        c = cp.loc[sel_cpd]
        ws = an.wells.loc[an.cpd_wells[sel_cpd]]
        wl = "<br>".join(f"P{esc(w.plate)}:{esc(w.well)} {fmt(w.pct)}%" for w in ws.itertuples())
        chid = chemo.cpd_chemo.get(sel_cpd) if chemo else None
        d = decisions.get(c["key"])
        parts.append(
            f'<div class="dsec big"><h4><span class="mono">{esc(c["bid"])}</span>{cls_chip(cl, cl.cls[sel_cpd])}<span class="sp"></span>'
            f'<button class="btn btn-sm btn-outline-secondary" data-chase="act" data-act="clear_cpd">Unpin</button></h4>{tile_hl(c["smiles"]) if any(flags_of(cp, sel_cpd)) else tile(c["smiles"])}'
            + alerts_block(cp, sel_cpd)
            + f'{triage_controls(c["key"], d)}<dl class="kv">'
            + (f'<dt>Name</dt><dd>{esc(c["name"])}</dd>' if c["name"] else "")
            + f'<dt>% response</dt><dd>{fmt(c["mean"])}%' + (f' +/- {fmt(c["sd"])} (n={c["n"]})' if c["n"] > 1 else "") + "</dd>"
            f"<dt>Wells</dt><dd>{wl}</dd>"
            f'<dt>Chemotype</dt><dd>' + (f'<a href="#" data-chase="act" data-act="open_chemo" data-id="{chid}">{chid}</a>' if chid else ("..." if c["smiles"] else "no structure")) + "</dd>"
            + (f'<dt>SMILES</dt><dd>{esc(c["smiles"])}</dd>' if c["smiles"] else "") + "</dl></div>")
    # 3. The selected chemotype: counts per class, statistics, buttons to set all members at once, and the list of members (the 120 most active).
    ch = chemo.by_id.get(sel_chemo) if (chemo and sel_chemo) else None
    if ch:
        chips = "".join(f'<span class="chip"><i style="background:{col(cl, t)}"></i>{ch.counts[j]} {l.lower()}</span>' for j, (l, t) in enumerate(cl.labels) if ch.counts[j])
        # The text on the 'set all' buttons tells the user whether existing decisions will be overwritten.
        und = sum(1 for i in ch.members if cp.at[i, "key"] not in decisions)
        lab = f"Set all {ch.n} members" if und == ch.n else (f"Set the {und} undecided member{'s' if und != 1 else ''}" if und else f"Set all {ch.n} members (overwrites their decisions)")
        rows = "".join(
            f'<button class="mrow" type="button" data-chase="sel_cpd" data-cpd="{i}" aria-pressed="{str(sel_cpd == i).lower()}">{tile(cp.at[i, "smiles"])}'
            f'<span class="who"><b>{esc(cp.at[i, "bid"])}</b><span>{esc(cp.at[i, "name"])}</span><span class="chips">{tri_chip(decisions.get(cp.at[i, "key"]))}{flag_chips(cp, i)}</span></span>'
            f'<span class="pct">{fmt(means[i])}%{cls_chip(cl, cl.cls[i])}</span></button>' for i in ch.members[:120])
        parts.append(
            f'<div class="dsec"><h4><span class="mono">{ch.id}</span>{chips}' + ('<span class="chip good">enriched</span>' if ch.q < 0.05 and ch.n > 1 else "")
            + '<span class="sp"></span><button class="btn btn-sm btn-outline-secondary" data-chase="act" data-act="clear_chemo">Clear</button></h4>'
            f'<dl class="kv"><dt>{esc(ch.rep_kind)}</dt><dd>{esc(ch.rep)}</dd><dt>Compounds</dt><dd>{ch.n}</dd>'
            f'<dt>Hit rate</dt><dd>{ch.rate * 100:.0f}% (library {chemo.base * 100:.1f}%)</dd>'
            f'<dt>Enrichment</dt><dd>EF {fmt(ch.ef, 1)} | p {fmt_p(ch.p)} | q {fmt_p(ch.q)}</dd><dt>Mean / max</dt><dd>{fmt(ch.mean)}% / {fmt(ch.max)}%</dd></dl>'
            f'<div class="trictl"><span class="lbl">{lab}</span><div class="seg tri">'
            + "".join(f'<button type="button" data-chase="triall" data-status="{sid}">{label}</button>' for sid, label, _ in STATUSES) + "</div></div>"
            f'<div class="mlist">{rows}' + (f'<p class="hint">Showing the 120 most active of {ch.n} members. Export the CSV for the full list.</p>' if ch.n > 120 else "") + "</div></div>")
    if not parts:
        return '<div class="empty">Click a chemotype card to list its compounds, or click a point on the plot to pin a compound and record a triage decision.</div>'
    return "".join(parts)


def alerts_block(cp, i) -> str:
    """Return the list of PAINS and covalent alerts matched by compound ``i`` (empty if none)."""
    p, c = flags_of(cp, i)
    if not (p or c):
        return ""
    rows = "".join(f'<li><span class="chip bad">PAINS</span> {esc(x)}</li>' for x in p.split("; ") if x)
    rows += "".join(f'<li><span class="chip cov">Covalent</span> {esc(x)}</li>' for x in c.split("; ") if x)
    return f'<div class="alerts"><span class="lbl">Structural alerts</span><ul>{rows}</ul><span class="hint">Matching atoms are highlighted: red for PAINS, violet for covalent warheads.</span></div>'


def same_cpd(cp, a, b) -> bool:
    """Return True if compounds ``a`` and ``b`` are the same substance (same name or identical SMILES), i.e. two batches."""
    A, B = cp.loc[a], cp.loc[b]
    if not A["smiles"] or not B["smiles"]:
        return False
    # Same compound name (for example the same CBK number) or the same SMILES: for example, these are two batches of the same compound, not an activity cliff.
    if A["name"] and A["name"] == B["name"]:
        return True
    return A["smiles"] == B["smiles"]


# ---------------------------------------------------------------- tooltip
def tooltip(an, cl, chemo, wi: int, decisions) -> str:
    """Return the hover card for well ``wi``: structure, IDs, class, flags and plate position (controls: no structure)."""
    w = an.wells.loc[wi]
    # Compound well: structure, Batch ID, chips for class, triage and alerts, the % response of this well
    # and of the compound overall, and where the well is.
    if w["type"] == "sample":
        ci = int(w["cpd"])
        c = an.cpds.loc[ci]
        chid = chemo.cpd_chemo.get(ci) if chemo else None
        d = decisions.get(c["key"])
        fc = flag_chips(an.cpds, ci)
        return (f'{tile(c["smiles"])}<div class="tt"><b>{esc(c["bid"])}</b>{cls_chip(cl, cl.cls[ci])}</div>'
                + (f'<div class="chips">{tri_chip(d)}{fc}</div>' if (d or fc) else "")
                + f'<div class="tt"><span class="pv">{fmt(w["pct"])}%</span><span class="sm">' + (f'mean {fmt(c["mean"])}% (n={c["n"]})' if c["n"] > 1 else "") + "</span></div>"
                f'<div class="sm">' + (f'<span>{esc(c["name"])}</span>' if c["name"] else "") + f'<span>P{esc(w["plate"])} | {esc(w["well"])}</span><span>raw {w["raw"]:,.0f}</span>'
                + (f"<span>{chid}</span>" if chid else "") + "</div>")
    # Control well: there is no structure to show.
    kind = "Positive control" if w["type"] == "pos" else "Negative control"
    return (f'<div class="tt"><b>{esc(w["bid"])}</b><span class="chip">{kind}</span></div><div class="tt"><span class="pv">{fmt(w["pct"])}%</span></div>'
            f'<div class="sm"><span>P{esc(w["plate"])} | {esc(w["well"])}</span><span>raw {w["raw"]:,.0f}</span></div>')


# ---------------------------------------------------------------- triage board
def board(decisions: dict, an, key_to_idx) -> str:
    """Return the triage board: one column per decision, most active compound first."""
    items = list(decisions.values())
    cols = []
    for sid, label, cls in STATUSES:
        # Most active first, decisions without a % response at the end. At most 200 cards per column.
        lst = sorted([d for d in items if d["status"] == sid], key=lambda d: -(d["pct"] if d["pct"] is not None else -1e9))
        body = []
        for d in lst[:200]:
            # A decision for a compound that is not in the data loaded now (for example imported from another
            # screen) is shown in grey.
            ci = key_to_idx.get(d["key"]) if key_to_idx else None
            smi = an.cpds.at[ci, "smiles"] if (an is not None and ci is not None) else ""
            body.append(f'<button class="titem{"" if ci is not None else " gone"}" type="button" data-chase="board" data-key="{esc(d["key"])}">'
                        + (tile(smi) if smi else "")
                        + f'<span class="tl"><b>{esc(d["bid"])}</b><span class="mono">{fmt(d["pct"]) + "%" if d["pct"] is not None else ""}</span></span>'
                        + (f'<span class="tn">{esc(d["name"])}</span>' if d.get("name") else "")
                        + (f'<span class="chips">{flag_chips(an.cpds, ci)}</span>' if (an is not None and ci is not None and any(flags_of(an.cpds, ci))) else "")
                        + (f'<span class="tnote">{esc(d["note"])}</span>' if d.get("note") else "") + meta(d) + "</button>")
        cols.append(f'<div class="bcol"><h5><span class="chip {cls}">{label}</span><span class="c">{len(lst)}</span></h5>'
                    + ("".join(body) if body else '<p class="hint">Nothing here yet.</p>') + "</div>")
    return '<div class="board">' + "".join(cols) + "</div>"


# ---------------------------------------------------------------- cliffs
def cliffs(lst, an, cl, chemo, sel_pair, shown: int) -> str:
    """Return the activity-cliff list (first ``shown`` pairs)."""
    cp = an.cpds
    means = cp["mean"].to_numpy()

    def side(i):
        """One compound of a cliff pair: structure, IDs, % response, class and flags."""
        return (f'<span class="cside">{tile(cp.at[i, "smiles"])}<span class="who"><b>{esc(cp.at[i, "bid"])}</b><span>{esc(cp.at[i, "name"])}</span></span>'
                f'<span class="pct mono">{fmt(means[i])}%</span><span class="chips">{cls_chip(cl, cl.cls[i])}{flag_chips(cp, i)}</span></span>')
    out = []
    for p in lst[:shown]:
        # Pairs that are two batches of the same compound, or that have identical fingerprints (Tanimoto 1, where SALI is capped), get a label saying so.
        on = bool(sel_pair) and tuple(sel_pair) == (p["hi"], p["lo"])
        same = same_cpd(cp, p["hi"], p["lo"])
        a, b = chemo.cpd_chemo.get(p["hi"]), chemo.cpd_chemo.get(p["lo"])
        flag = '<span class="chip warn">same compound, other batch</span>' if same else ('<span class="chip">identical fingerprint</span>' if p["sim"] >= 0.995 else "")
        out.append(f'<button class="cliff" type="button" data-chase="sel_pair" data-hi="{p["hi"]}" data-lo="{p["lo"]}" aria-pressed="{str(on).lower()}">{side(p["hi"])}'
                   f'<span class="cmid"><span class="mono big">Delta {p["d"]:.0f}</span><span class="mono">Tc {p["sim"]:.2f}</span>'
                   f'<span class="mono sm">SALI {">= " if p["sim"] >= 0.995 else ""}{p["sali"]:.0f}</span>{flag}' + (f'<span class="sm">{a}</span>' if a and a == b else "")
                   + f"</span>{side(p['lo'])}</button>")
    return "".join(out)


# ---------------------------------------------------------------- plates
def mini_maps(an, dims, sel: int) -> str:
    """Return one small SVG plate map per plate (coloured by % response) with its Z' chip."""
    R, Cc = dims
    # Size of one well in the drawing: smaller for 1536-well plates. Empty wells are drawn as outlines only.
    cell = 3 if Cc > 24 else 6
    out = []
    for p in an.plates:
        z, _ = plate_grid(an.wells, p["idx"], dims)
        rects = []
        for r in range(R):
            for c in range(Cc):
                v = z[r][c]
                if v is None:
                    rects.append(f'<rect x="{c * cell + .75}" y="{r * cell + .75}" width="{cell - 1.5}" height="{cell - 1.5}" fill="none" stroke="{C["line"]}" stroke-width=".5"/>')
                else:
                    rgb = heat_rgb(v)
                    g = 1 if cell > 4 else 0
                    rects.append(f'<rect x="{c * cell}" y="{r * cell}" width="{cell - g}" height="{cell - g}" fill="rgb({rgb[0]},{rgb[1]},{rgb[2]})"/>')
        zc = "" if not math.isfinite(p["z"]) else f'<span class="chip {"good" if p["z"] >= 0.5 else "warn" if p["z"] > 0 else "bad"}">Z\' {p["z"]:.2f}</span>'
        out.append(f'<button class="mini" type="button" data-chase="plate" data-pi="{p["idx"]}" aria-pressed="{str(p["idx"] == sel).lower()}">'
                   f'<svg viewBox="0 0 {Cc * cell} {R * cell}" shape-rendering="crispEdges" aria-hidden="true">{"".join(rects)}</svg>'
                   f'<span class="ml"><b>P{esc(p["label"])}</b>{zc}</span></button>')
    return '<div class="minigrid">' + "".join(out) + "</div>"


def qc_table(an, cl) -> str:
    """Return the plate QC table: controls, CVs, S/B, Z' and the number of active wells per plate."""
    # Count, for each plate, the wells of compounds that are not inactive.
    hits: dict[str, int] = {}
    for ci in np.nonzero(cl.cls > 0)[0]:
        for wi in an.cpd_wells[ci]:
            pl = an.wells.at[wi, "plate"]
            hits[pl] = hits.get(pl, 0) + 1

    def zchip(z):
        """Return a coloured chip rating Z' (>= 0.5 excellent, 0-0.5 marginal, < 0 fail)."""
        if not math.isfinite(z):
            return '<span class="chip">n/a</span>'
        return f'<span class="chip good">{z:.2f} | excellent</span>' if z >= 0.5 else (f'<span class="chip warn">{z:.2f} | marginal</span>' if z > 0 else f'<span class="chip bad">{z:.2f} | fail</span>')
    rows = "".join(
        f'<tr><td>P{esc(p["label"])}' + (' <span class="chip warn">uses all-plate controls</span>' if p["borrowed"] else "") + f'</td><td>{p["nSample"]}</td>'
        f'<td>{fmt(p["mN"], 0)} +/- {fmt(p["sN"], 0)}</td><td>{fmt(p["mP"], 0)} +/- {fmt(p["sP"], 0)}</td><td>{fmt(p["cvN"])}%</td><td>{fmt(p["cvP"])}%</td>'
        f'<td>{fmt(p["sb"])}</td><td>{zchip(p["z"])}</td><td>{hits.get(p["label"], 0)}</td></tr>' for p in an.plates)
    return ('<div class="tablewrap"><table class="qc"><thead><tr><th>Plate</th><th>Samples</th><th>Neg mean +/- SD</th><th>Pos mean +/- SD</th><th>CV neg</th><th>CV pos</th>'
            f'<th>S/B</th><th>Z\'</th><th>Wells with activity</th></tr></thead><tbody>{rows}</tbody></table></div>')
