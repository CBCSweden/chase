"""Plotly figure specs (rendered by plotly.js in the browser).

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import math

import numpy as np

from .core.plates import plate_grid, row_label

# Colours used in the plots. They match the ones in chase.css (background, greys for text, lines,
# the accent colour and the default class colours).
C = dict(surface="#ffffff", ink="#131715", ink2="#48514d", ink3="#737d78", line="#d9dfdc", grid="#edf0ee",
         band="rgba(19,23,21,0.035)", accent="#0e6a5e", active="#eb6834", moderate="#1baf7a", weak="#2a78d6",
         inactive="#a3aba7", control="#48514d", hm_neg="#6250d6", hm_zero="#eef1ef")
FONT = "Public Sans, Segoe UI, Helvetica, Arial, sans-serif"
# Settings for plotly.js. CHASE has its own JPEG button, so plotly's camera button is removed.
CONFIG = dict(responsive=True, displaylogo=False, scrollZoom=False,
              modeBarButtonsToRemove=["toImage", "lasso2d", "select2d", "autoScale2d"])
# Range of the heatmap colours, in % response.
HM_MIN, HM_MAX = -50.0, 120.0
# The colours the user can change under 'Colors'.
DEFAULT_COLORS = dict(active=C["active"], moderate=C["moderate"], weak=C["weak"], inactive=C["inactive"], pos=C["control"], neg=C["control"],
                      pains="#c0392b", covalent="#6250d6")
# Colours of the Z' bars: green for 0.5 and above, amber between 0 and 0.5, red for 0 and below.
Z_COLORS = dict(good="#2e8b57", warn="#d39b1a", bad="#c0392b")


def col(cl, key: str) -> str:
    """Return the colour for ``key`` (class token, 'pos', 'neg', 'pains' or 'covalent'): the user's choice or the default."""
    return (getattr(cl, "colors", None) or {}).get(key) or DEFAULT_COLORS.get(key) or C.get(key, "#888888")


def _f(a):
    """Convert values to JSON-friendly floats, with None for NaN."""
    return [None if (v is None or (isinstance(v, float) and math.isnan(v))) else float(v) for v in a]


def scatter(an, cl, show_ctrl: bool, sel_chemo, sel_cpd, sel_pair, title: str | None = None, show_z: bool = True, flag_mode: str = "hits") -> dict:
    """Build the '% response by compound' figure.

    Traces: control wells, one trace per activity class, PAINS/covalent marks, rings for the selected chemotype,
    cliff pair and pinned compound, a marginal histogram of compound means (log count axis) and, above the plot,
    one Z' bar per plate. Plates are shaded alternately and class borders are drawn as dashed lines.

    Args:
        an: The :class:`Analysis`.
        cl: Activity classes (with the user's colours).
        show_ctrl: Draw control wells.
        sel_chemo: Selected :class:`Chemotype` or None.
        sel_cpd: Index of the pinned compound or None.
        sel_pair: (high, low) compound indices of the selected cliff pair, or None.
        title: Optional title (used for exports).
        show_z: Draw the Z' strip.
        flag_mode: ``"hits"``, ``"all"`` or ``"none"`` - which flagged compounds get PAINS/covalent marks.

    Returns:
        A plotly.js figure as ``{"data": [...], "layout": {...}}``.
    """
    W = an.wells
    samples = W[W["type"] == "sample"]
    cls_w = cl.cls[samples["cpd"].to_numpy()]
    # For big screens use WebGL, which stays fast with many points. Small screens use normal SVG,
    # which looks sharper.
    typ = "scattergl" if len(samples) > 8000 else "scatter"
    traces = []
    # Control wells are drawn first, as open symbols, so the compounds appear on top of them.
    # customdata holds the well number, so hovering and clicking know which well it is.
    if show_ctrl:
        for t, sym, nm in (("neg", "diamond-open", f"Negative control | {an.neg_name}"), ("pos", "square-open", f"Positive control | {an.pos_name}")):
            ws = W[W["type"] == t]
            traces.append(dict(type=typ, mode="markers", name=nm, x=_f(ws["x"]), y=_f(ws["pct"]), customdata=ws.index.tolist(), hoverinfo="none",
                               marker=dict(symbol=sym, size=6, color=col(cl, t), opacity=0.6, line=dict(width=1))))
    for ci, (label, tok) in enumerate(cl.labels):
        # One set of points per activity class. Inactive compounds are smaller and paler, so the hits stand out.
        ws = samples[cls_w == ci]
        traces.append(dict(type=typ, mode="markers", name=f"{label} ({cl.counts[ci]:,})", x=_f(ws["x"]), y=_f(ws["pct"]), customdata=ws.index.tolist(),
                           hoverinfo="none", marker=dict(color=col(cl, tok), size=5 if ci == 0 else 7.5, opacity=0.6 if ci == 0 else 0.95,
                                                         line=dict(width=0 if ci == 0 else 0.8, color=C["surface"]))))

    # Mark the wells of flagged compounds: an x for PAINS and an open triangle for covalent alerts.
    # With 'hits' only compounds above Inactive are marked. The columns only exist once the alerts
    # have been calculated.
    if flag_mode != "none" and "pains" in an.cpds.columns:
        cp = an.cpds
        for colname, label, sym, size, key in (("pains", "PAINS", "x-thin", 11, "pains"), ("covalent", "Covalent warhead", "triangle-up-open", 13, "covalent")):
            flagged = (cp[colname] != "")
            if flag_mode == "hits":
                flagged &= cl.cls > 0
            idx = [w for c in np.nonzero(flagged.to_numpy())[0] for w in an.cpd_wells[c]]
            if idx:
                ws = W.loc[idx]
                traces.append(dict(type=typ, mode="markers", name=f"{label} ({int(flagged.sum()):,})", x=_f(ws["x"]), y=_f(ws["pct"]), customdata=idx,
                                   hoverinfo="none", marker=dict(symbol=sym, size=size, color=col(cl, key), line=dict(width=1.8, color=col(cl, key)))))

    def ring(cpd_list, name, sym, size, color, width, legend=True):
        """Add a trace of open markers around all wells of the given compounds (selection highlight)."""
        idx = [w for c in cpd_list for w in an.cpd_wells[c]]
        ws = W.loc[idx]
        traces.append(dict(type=typ, mode="markers", name=name, showlegend=legend, x=_f(ws["x"]), y=_f(ws["pct"]), customdata=idx, hoverinfo="none",
                           marker=dict(symbol=sym, size=size, color=color, line=dict(width=width))))

    # Rings around the selection, drawn last so they are on top: chemotype, cliff pair, pinned compound.
    if sel_chemo is not None:
        ring(sel_chemo.members, f"{sel_chemo.id} ({sel_chemo.n} compound{'s' if sel_chemo.n != 1 else ''})", "circle-open", 15, C["accent"], 2.2)
    if sel_pair:
        ring(list(sel_pair), "Activity-cliff pair", "diamond-open", 18, C["ink"], 2.2)
    if sel_cpd is not None:
        ring([sel_cpd], an.cpds.at[sel_cpd, "bid"], "circle-open", 22, C["ink"], 2.4, legend=False)
    # A histogram of the compound means on the right-hand side (the last 11% of the width).
    # It shares the % response axis with the main plot.
    traces.append(dict(type="histogram", y=_f(an.cpds["mean"]), xaxis="x2", yaxis="y", nbinsy=70, marker=dict(color=C["inactive"], line=dict(width=0)),
                       opacity=0.9, showlegend=False, hovertemplate="%{y}%: %{x} compounds<extra></extra>", name="Distribution"))
    shapes, ann = [], []
    # Every other plate gets a shaded background. Plate labels are only shown when there are few enough
    # plates to read them. With the Z' strip switched on, the bars carry the plate labels instead.
    show_labels = len(an.plates) <= 36
    # With the Z' strip on, the main plot uses the bottom 78% of the height and the strip the top 17%.
    main_top = 0.78 if show_z else 1.0
    zs = []
    for i, p in enumerate(an.plates):
        if not math.isfinite(p["x0"]):
            continue
        if i % 2 == 1:
            shapes.append(dict(type="rect", xref="x", yref="paper", x0=p["x0"] - 0.5, x1=p["x1"] + 0.5, y0=0, y1=1, fillcolor=C["band"], line=dict(width=0), layer="below"))
        if show_labels and not show_z:
            ann.append(dict(xref="x", yref="paper", x=(p["x0"] + p["x1"]) / 2, y=1, yanchor="top", text=f"P{p['label']}", showarrow=False, font=dict(size=10, color=C["ink3"])))
        if show_z and math.isfinite(p["z"]):
            zs.append(p)
    if show_z and zs:
        # The Z' strip: one bar per plate, as wide as that plate's compounds. Values below -1 are drawn at -1
        # (the hover text still shows the real value), so one failed plate doesn't squeeze all the others.
        lab = len(an.plates) <= 24
        status = lambda z: "good" if z >= 0.5 else ("warn" if z > 0 else "bad")
        traces.append(dict(type="bar", xaxis="x", yaxis="y3", x=[(p["x0"] + p["x1"]) / 2 for p in zs], y=[max(-1.0, p["z"]) for p in zs],
                           width=[max(1.0, (p["x1"] - p["x0"] + 1) * 0.86) for p in zs], marker=dict(color=[Z_COLORS[status(p["z"])] for p in zs], line=dict(width=0)),
                           opacity=0.85, showlegend=False, name="Z'",
                           text=[f"P{p['label']} | {p['z']:.2f}" if lab else "" for p in zs], textposition="inside", insidetextanchor="end",
                           textfont=dict(size=10, color="#ffffff"), cliponaxis=False,
                           hovertext=[f"Plate {p['label']}<br>Z' {p['z']:.2f}{' (shown as -1)' if p['z'] < -1 else ''}<br>S/B {p['sb']:.1f}" for p in zs],
                           hovertemplate="%{hovertext}<extra></extra>"))
        # The Z' axis runs from the lowest bar (but not below -1) up to 1. The dotted line marks 0.5.
        zmin = min(0.0, max(-1.0, min(p["z"] for p in zs)))
        shapes.append(dict(type="line", xref="paper", yref="y3", x0=0, x1=0.86, y0=0.5, y1=0.5, line=dict(color=C["ink3"], width=1, dash="dot")))
        ann.append(dict(xref="paper", yref="y3", x=0.86, y=0.5, xanchor="left", yanchor="middle", text="0.5", showarrow=False, font=dict(size=10, color=C["ink3"])))
    # A dashed line with a label at every border between two classes.
    for j, b in enumerate(cl.bounds):
        shapes.append(dict(type="line", xref="paper", yref="y", x0=0, x1=0.86, y0=b, y1=b, line=dict(color=C["ink3"], width=1, dash="dash")))
        ann.append(dict(xref="paper", yref="y", x=0.86, y=b, xanchor="right", yanchor="bottom", text=f"{cl.labels[j + 1][0]} >= {b:.1f}%", showarrow=False,
                        font=dict(size=11, color=C["ink2"]), bgcolor=C["surface"], borderpad=1))
    layout = dict(
        paper_bgcolor=C["surface"], plot_bgcolor=C["surface"], font=dict(family=FONT, size=12, color=C["ink2"]),
        # uirevision keeps the user's zoom when the plot is redrawn, for example after selecting something.
        margin=dict(l=60, r=14, t=70 if title else 34, b=50), uirevision="main",
        xaxis=dict(domain=[0, 0.86], title=dict(text="Compound number", standoff=8), zeroline=False, gridcolor=C["grid"], linecolor=C["line"], tickcolor=C["line"], ticks="outside", showline=True, rangemode="tozero"),
        yaxis=dict(domain=[0, main_top], title=dict(text="% response", standoff=6), zeroline=True, zerolinecolor=C["line"], zerolinewidth=1.5, gridcolor=C["grid"], linecolor=C["line"], tickcolor=C["line"], ticks="outside", showline=True),
        # The count axis is logarithmic, so a handful of actives stays visible next to thousands of inactives.
        xaxis2=dict(domain=[0.885, 1], anchor="y", type="log", showgrid=False, zeroline=False, title=dict(text="Compounds (log)", font=dict(size=11)), tickfont=dict(size=10), linecolor=C["line"], showline=True, dtick=1),
        legend=dict(orientation="h", x=0, xanchor="left", y=1.01, yanchor="bottom", font=dict(size=12, color=C["ink"]), bgcolor="rgba(0,0,0,0)"),
        shapes=shapes, annotations=ann, hovermode="closest", hoverdistance=14, dragmode="zoom", bargap=0.05)
    if show_z and zs:
        layout["yaxis3"] = dict(domain=[0.83, 1.0], anchor="x", title=dict(text="Z'", standoff=6), range=[zmin - 0.05, 1.0], fixedrange=True,
                                gridcolor=C["grid"], zeroline=True, zerolinecolor=C["line"], linecolor=C["line"], showline=True, tickfont=dict(size=10), nticks=4)
        layout["legend"]["y"] = 1.02
    if title:
        layout["title"] = dict(text=title, x=0.01, xanchor="left", y=0.985, font=dict(size=15, color=C["ink"]))
    return dict(data=traces, layout=layout)


def _hex(h):
    """Convert ``#rrggbb`` to an [r, g, b] list."""
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) for i in (0, 2, 4)]


def heat_rgb(v: float) -> tuple[int, int, int]:
    """Colour for a % response on the plate maps: neutral at 0, towards orange at 100, towards violet at -50."""
    z, pos, neg = _hex(C["hm_zero"]), _hex(C["active"]), _hex(C["hm_neg"])
    # Mix the neutral colour gradually towards red or blue. The colour stops changing at 100% and at -50%.
    if v >= 0:
        t, pole = min(1.0, v / 100), pos
    else:
        t, pole = min(1.0, -v / 50), neg
    return tuple(round(a + (b - a) * t) for a, b in zip(z, pole))


def heatmap(an, plate_idx: int, dims, title: str | None = None) -> dict:
    """Build the large heatmap of one plate (% response per well); hover data carry the well index."""
    R, Cc = dims
    z, cd = plate_grid(an.wells, plate_idx, dims)
    # Build the plotly colour scale with the same colour function as the small plate maps, so they match.
    stops = [[(v - HM_MIN) / (HM_MAX - HM_MIN), "rgb({},{},{})".format(*heat_rgb(v))] for v in (HM_MIN, -25, 0, 25, 50, 75, 100, HM_MAX)]
    data = [dict(type="heatmap", z=z, customdata=cd, x=[str(c + 1) for c in range(Cc)], y=[row_label(r) for r in range(R)], zmin=HM_MIN, zmax=HM_MAX,
                 colorscale=stops, xgap=2, ygap=2, hoverinfo="none", hoverongaps=False,
                 colorbar=dict(title=dict(text="% resp.", side="right", font=dict(size=11, color=C["ink2"])), thickness=10, len=0.9, tickfont=dict(size=10, color=C["ink2"]), outlinewidth=0))]
    layout = dict(paper_bgcolor=C["surface"], plot_bgcolor=C["surface"], font=dict(family=FONT, size=11, color=C["ink2"]),
                  margin=dict(l=34, r=10, t=56 if title else 40, b=10), uirevision=f"heat{plate_idx}",
                  xaxis=dict(side="top", type="category", tickfont=dict(size=10), showgrid=False, zeroline=False, ticks=""),
                  # Row A at the top and the column numbers above the plate, the way a real plate looks.
                  yaxis=dict(autorange="reversed", type="category", tickfont=dict(size=10), showgrid=False, zeroline=False, ticks=""))
    if title:
        layout["title"] = dict(text=title, x=0.01, xanchor="left", font=dict(size=14, color=C["ink"]))
    return dict(data=data, layout=layout)


def replicate_plot(an, cl, rep: dict) -> dict:
    """Build the replicate-1 versus replicate-2 scatter plot with the identity line and +/- threshold lines."""
    pts, T = rep["pts"], rep["T"]
    traces = []
    for ci, (label, tok) in enumerate(cl.labels):
        ps = [p for p in pts if cl.cls[p["ci"]] == ci]
        traces.append(dict(type="scatter", mode="markers", name=label, x=[p["a"] for p in ps], y=[p["b"] for p in ps], customdata=[p["ci"] for p in ps], hoverinfo="none",
                           marker=dict(color=col(cl, tok), size=5.5 if ci == 0 else 7.5, opacity=0.6 if ci == 0 else 0.95, line=dict(width=0 if ci == 0 else 0.8, color=C["surface"]))))
    if rep["disc"]:
        d = rep["disc"]
        traces.append(dict(type="scatter", mode="markers", name=f"|Delta| > {T:.0f} pts", x=[p["a"] for p in d], y=[p["b"] for p in d], customdata=[p["ci"] for p in d], hoverinfo="none",
                           marker=dict(symbol="circle-open", size=13, color=C["ink"], line=dict(width=1.6))))
    allv = [v for p in pts for v in (p["a"], p["b"])]
    # Use square axes from at least -10 to 110%. Draw the diagonal (both replicates equal) and dotted lines
    # at plus and minus the threshold for replicates that disagree.
    lo, hi = min([-10.0] + allv) - 5, max([110.0] + allv) + 5
    shapes = [dict(type="line", x0=lo, y0=lo, x1=hi, y1=hi, line=dict(color=C["ink3"], width=1)),
              dict(type="line", x0=lo, y0=lo + T, x1=hi - T, y1=hi, line=dict(color=C["ink3"], width=1, dash="dot")),
              dict(type="line", x0=lo + T, y0=lo, x1=hi, y1=hi - T, line=dict(color=C["ink3"], width=1, dash="dot"))]
    layout = dict(paper_bgcolor=C["surface"], plot_bgcolor=C["surface"], font=dict(family=FONT, size=12, color=C["ink2"]), margin=dict(l=56, r=14, t=30, b=48),
                  shapes=shapes, hovermode="closest", dragmode="zoom", uirevision="rep",
                  xaxis=dict(title=dict(text="Replicate 1 (% response)"), range=[lo, hi], gridcolor=C["grid"], zeroline=False, linecolor=C["line"], showline=True),
                  yaxis=dict(title=dict(text="Replicate 2 (% response)"), range=[lo, hi], gridcolor=C["grid"], zeroline=False, linecolor=C["line"], showline=True, scaleanchor="x", scaleratio=1),
                  legend=dict(orientation="h", x=0, y=1.02, yanchor="bottom", font=dict(size=12, color=C["ink"])))
    return dict(data=traces, layout=layout)
