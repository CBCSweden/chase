"""Plate layouts, heatmap grids, edge effects and replicate agreement.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# A well ID is one or two row letters (A-Z, and AA-AF on 1536-well plates) followed by a column number,
# for example A01, P24 or AF48.
_WELL = re.compile(r"^([A-Za-z]{1,2})\s*0*(\d{1,3})$")


def well_rc(well: str) -> tuple[int, int] | None:
    """Convert a well ID such as ``A01``, ``A1`` or ``AF48`` to zero-based (row, column); None if unreadable."""
    m = _WELL.match(str(well or "").strip())
    if not m:
        return None
    L = m.group(1).upper()
    # Row letters to numbers: A is 0 and Z is 25. With two letters the count goes on: AA is 26, AB is 27, ...
    r = ord(L) - 65 if len(L) == 1 else (ord(L[0]) - 64) * 26 + (ord(L[1]) - 65)
    return r, int(m.group(2)) - 1


def row_label(r: int) -> str:
    """Convert a zero-based row index to its plate letter(s): 0 -> A, 26 -> AA."""
    return chr(65 + r) if r < 26 else chr(64 + r // 26) + chr(65 + r % 26)


def plate_dims(wells: pd.DataFrame) -> tuple[int, int]:
    """Return the plate layout (rows, columns) seen in the data, snapped to 8 x 12, 16 x 24 or 32 x 48 wells."""
    R = C = 0
    for w in wells["well"]:
        rc = well_rc(w)
        if rc:
            R, C = max(R, rc[0] + 1), max(C, rc[1] + 1)
    # Pick the smallest standard plate (96, 384 or 1536 wells) that has room for every well we have seen.
    for r, c in ((8, 12), (16, 24), (32, 48)):
        if R <= r and C <= c:
            return r, c
    return R, C


def plate_grid(wells: pd.DataFrame, plate_idx: int, dims: tuple[int, int]):
    """Return two row x column grids for one plate: % response (None for empty wells) and well indices (-1)."""
    R, C = dims
    # z holds the % response of each well for the heatmap. cd holds the well number at each position,
    # so we know which well was clicked.
    z = [[None] * C for _ in range(R)]
    cd = [[-1] * C for _ in range(R)]
    sub = wells[wells["plate_idx"] == plate_idx]
    for wi, w, p in zip(sub.index, sub["well"], sub["pct"]):
        rc = well_rc(w)
        if rc and rc[0] < R and rc[1] < C:
            z[rc[0]][rc[1]] = float(p)
            cd[rc[0]][rc[1]] = int(wi)
    return z, cd


def edge_stats(wells: pd.DataFrame, plate_idx: int, dims: tuple[int, int]) -> dict:
    """Mean % response of sample wells on the outer ring of the plate versus the interior.

    Returns:
        ``{"edge", "inner", "nE", "nI"}`` - the two means and the number of wells behind each.
    """
    R, C = dims
    e, m = [], []
    sub = wells[(wells["plate_idx"] == plate_idx) & (wells["type"] == "sample")]
    for w, p in zip(sub["well"], sub["pct"]):
        rc = well_rc(w)
        if not rc:
            continue
        # Wells in the first or last row or column are edge wells. All the others are interior wells.
        (e if rc[0] in (0, R - 1) or rc[1] in (0, C - 1) else m).append(p)
    return dict(edge=float(np.mean(e)) if e else math.nan, inner=float(np.mean(m)) if m else math.nan, nE=len(e), nI=len(m))


def replicates(an, threshold: float) -> dict:
    """Compare the first two wells of every compound tested more than once.

    Args:
        an: The :class:`Analysis`.
        threshold: Difference in % points above which replicates are called discordant.

    Returns:
        ``{"pts", "T", "disc"}`` and, with at least three compounds, Pearson ``r``, Spearman ``rho`` and the median difference ``med``.
    """
    pts = []
    for ci, ws in enumerate(an.cpd_wells):
        if len(ws) >= 2:
            # Replicate 1 and replicate 2 are the first two wells of the compound, in plate and file order.
            a, b = an.wells.loc[ws[:2], "pct"].tolist()
            pts.append(dict(ci=ci, a=a, b=b, d=abs(a - b)))
    out = dict(pts=pts, T=threshold, disc=sorted([p for p in pts if p["d"] > threshold], key=lambda p: -p["d"]))
    # Correlations don't mean much with fewer than three compounds, so we only compute them from three on.
    if len(pts) >= 3:
        A = [p["a"] for p in pts]
        B = [p["b"] for p in pts]
        out["r"] = float(pearsonr(A, B).statistic)
        out["rho"] = float(spearmanr(A, B).statistic)
        out["med"] = float(np.median([p["d"] for p in pts]))
    return out
