"""Control-based normalization, plate QC and per-compound aggregation.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .io import alt_key


def nat_key(s: str):
    """Sort key that orders plate labels naturally (P2 before P10)."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


def _sd(a: np.ndarray) -> float:
    """Sample standard deviation (ddof = 1); NaN for fewer than two values."""
    return float(np.std(a, ddof=1)) if len(a) > 1 else math.nan


@dataclass
class Analysis:
    """Result of normalization: wells, compounds, plate QC and the control names used."""
    wells: pd.DataFrame          # one row for each well that is used (compounds and controls)
    cpds: pd.DataFrame           # one row for each compound
    plates: list[dict]           # quality numbers for each plate, in plate order
    pos_name: str
    neg_name: str
    scope: str
    stat: str
    n_empty: int
    cpd_wells: list[list[int]] = field(default_factory=list)

    @property
    def key_to_idx(self) -> dict[str, int]:
        """Map each compound key to its row index in :attr:`cpds`."""
        return {k: i for i, k in enumerate(self.cpds["key"])}


class NormalizationError(ValueError):
    """Raised when the controls cannot be used (missing, identical or not found in the data)."""
    pass


def analyze(rows: pd.DataFrame, pos: str, neg: str, scope: str = "plate", stat: str = "mean",
            struct_map: dict[str, dict] | None = None) -> Analysis:
    """Normalize raw readouts to % response and aggregate wells into compounds.

    % response = 100 * (raw - neg) / (pos - neg), where pos and neg are the mean (or median) of the control wells
    on the same plate (``scope="plate"``) or on all plates (``scope="global"``). A plate missing a control falls back
    to the all-plate value and is flagged as ``borrowed``.

    Args:
        rows: Output of :meth:`ScreenFile.extract`.
        pos: Batch ID of the positive control (defines 100%).
        neg: Batch ID of the negative control (defines 0%).
        scope: ``"plate"`` or ``"global"``.
        stat: ``"mean"`` or ``"median"`` of the control wells.
        struct_map: Output of :meth:`StructFile.mapping`, used to attach SMILES and names.

    Returns:
        An :class:`Analysis`.

    Raises:
        NormalizationError: if a control ID is empty, duplicated or absent.
    """
    # Match the controls the same way as the wells: on the Batch ID in capitals.
    pos_k, neg_k = pos.strip().upper(), neg.strip().upper()
    if not pos_k or not neg_k:
        raise NormalizationError("Enter both a positive and a negative control Batch ID in step 2.")
    if pos_k == neg_k:
        raise NormalizationError("The positive and negative controls must be different Batch IDs.")
    if not (rows["key"] == pos_k).any():
        raise NormalizationError(f"No wells have the Batch ID '{pos}'. Check the positive control ID.")
    if not (rows["key"] == neg_k).any():
        raise NormalizationError(f"No wells have the Batch ID '{neg}'. Check the negative control ID.")
    fstat = np.median if stat == "median" else np.mean
    # Control values over all plates together. They are used for 'All plates' normalization, and for any
    # plate that is missing one of its controls.
    g_pos = fstat(rows.loc[rows["key"] == pos_k, "raw"].to_numpy())
    g_neg = fstat(rows.loc[rows["key"] == neg_k, "raw"].to_numpy())
    plates_order = sorted(rows["plate"].unique(), key=nat_key)
    recs, plates = [], []
    x = 0
    n_empty_total = 0
    for pi, pl in enumerate(plates_order):
        # The wells of this plate in file order, and the raw values of its positive and negative controls.
        pr = rows[rows["plate"] == pl].sort_values("ord")
        P = pr.loc[pr["key"] == pos_k, "raw"].to_numpy()
        N = pr.loc[pr["key"] == neg_k, "raw"].to_numpy()
        local = scope == "plate"
        # The control values that mean 100% (p_s) and 0% (n_s) on this plate.
        p_s = fstat(P) if (local and len(P)) else g_pos
        n_s = fstat(N) if (local and len(N)) else g_neg
        # Mean and standard deviation of this plate's own controls, for the QC table (NaN if a control is missing).
        mP, mN = (float(np.mean(P)) if len(P) else math.nan), (float(np.mean(N)) if len(N) else math.nan)
        sP, sN = _sd(P), _sd(N)
        # Plate quality always uses this plate's own controls:
        # Z' = 1 - 3 * (SD pos + SD neg) / |mean pos - mean neg|
        z = 1 - 3 * (sP + sN) / abs(mP - mN) if len(P) > 1 and len(N) > 1 else math.nan
        # Everything we report about the plate. S/B is mean positive / mean negative, CV is in %.
        # 'borrowed' is True when the plate had to use the all-plate control values because one of its own
        # controls is missing. x0 and x1 are where the plate's wells start and end on the plot.
        qc = dict(label=pl, idx=pi, nPos=len(P), nNeg=len(N), mP=mP, mN=mN, sP=sP, sN=sN, z=z,
                  sb=mP / mN if mN else math.nan, cvP=100 * sP / mP if mP else math.nan, cvN=100 * sN / mN if mN else math.nan,
                  borrowed=local and (not len(P) or not len(N)), nSample=0, nEmpty=0, x0=math.inf, x1=-math.inf)
        span = p_s - n_s
        for r in pr.itertuples(index=False):
            # Wells without a Batch ID are counted, but otherwise left out.
            if not r.key:
                qc["nEmpty"] += 1
                continue
            typ = "pos" if r.key == pos_k else "neg" if r.key == neg_k else "sample"
            # Position on the x axis of the plot: compound wells are numbered 1, 2, 3, ... in plate and file order.
            # Control wells are placed half-way between the compound wells next to them, so they appear inside their plate.
            if typ == "sample":
                x += 1
                wx = float(x)
                qc["nSample"] += 1
            else:
                wx = x + 0.5
            qc["x0"], qc["x1"] = min(qc["x0"], wx), max(qc["x1"], wx)
            recs.append((r.raw, pl, pi, r.well, r.bid, r.key, typ, 100 * (r.raw - n_s) / span, wx))
        n_empty_total += qc["nEmpty"]
        plates.append(qc)
    wells = pd.DataFrame(recs, columns=["raw", "plate", "plate_idx", "well", "bid", "key", "type", "pct", "x"])
    # Group the compound wells into compounds: one compound per Batch ID, in the order they first appear.
    wells["cpd"] = -1
    samples = wells[wells["type"] == "sample"]
    order = list(dict.fromkeys(samples["key"]))
    idx_of = {k: i for i, k in enumerate(order)}
    wells.loc[samples.index, "cpd"] = samples["key"].map(idx_of).to_numpy()
    cpd_wells: list[list[int]] = [[] for _ in order]
    for wi, ci in zip(samples.index, wells.loc[samples.index, "cpd"]):
        cpd_wells[int(ci)].append(int(wi))
    sm = struct_map or {}
    # One record per compound: mean, SD and maximum of its replicate wells, and its structure.
    # The structure is looked up by Batch ID, and if that fails by the Batch ID without leading zeros.
    crec = []
    for ci, k in enumerate(order):
        v = wells.loc[cpd_wells[ci], "pct"].to_numpy()
        st = sm.get(k) or sm.get(alt_key(k))
        crec.append(dict(idx=ci, key=k, bid=wells.at[cpd_wells[ci][0], "bid"], n=len(v), mean=float(v.mean()),
                         sd=_sd(v), max=float(v.max()), smiles=st["smiles"] if st else "", name=st["name"] if st else ""))
    cpds = pd.DataFrame(crec, columns=["idx", "key", "bid", "n", "mean", "sd", "max", "smiles", "name"])
    # The control names as they are written in the file (not in capitals), for the labels.
    pos_name = rows.loc[rows["key"] == pos_k, "bid"].iloc[0]
    neg_name = rows.loc[rows["key"] == neg_k, "bid"].iloc[0]
    return Analysis(wells=wells, cpds=cpds, plates=plates, pos_name=pos_name, neg_name=neg_name, scope=scope, stat=stat,
                    n_empty=n_empty_total, cpd_wells=cpd_wells)


def relink(an: Analysis, struct_map: dict[str, dict] | None) -> None:
    """Re-attach SMILES and names to the compounds of ``an`` after the structures file or its columns changed."""
    sm = struct_map or {}
    smiles, names = [], []
    for k in an.cpds["key"]:
        st = sm.get(k) or sm.get(alt_key(k))
        smiles.append(st["smiles"] if st else "")
        names.append(st["name"] if st else "")
    an.cpds["smiles"] = smiles
    an.cpds["name"] = names
