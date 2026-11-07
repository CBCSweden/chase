"""Chemotype statistics (hit rate, enrichment, BH q-values) and activity cliffs.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import hypergeom

from .chem import Pairs, Prepared, butina

# The lowest Tanimoto similarity the activity-cliff search can use. Similar pairs are kept down to this value.
CLIFF_MIN = 0.6
# With more compounds than this, comparing every pair would take too long and need too much memory
# (the work grows with n^2), so CHASE groups compounds by scaffold instead.
PAIR_LIMIT = 20000


@dataclass
class Chemotype:
    """One chemotype (Butina cluster or Murcko scaffold group) and its activity statistics."""
    id: str
    members: list[int]            # compound numbers, most active first once stats() has run
    rep: str                      # SMILES shown for the chemotype (its centre compound or its scaffold)
    rep_kind: str
    n: int = 0
    counts: list[int] = field(default_factory=list)
    n_hit: int = 0
    n_top: int = 0
    rate: float = 0.0
    mean: float = 0.0
    max: float = 0.0
    ef: float = float("nan")
    p: float = 1.0
    q: float = 1.0


@dataclass
class ChemoResult:
    """All chemotypes of a data set, plus the similarity pairs used for activity cliffs and the library hit rate."""
    chemos: list[Chemotype]
    by_id: dict[str, Chemotype]
    cpd_chemo: dict[int, str]     # for each compound number, the ID of its chemotype
    method: str
    note: str = ""
    cliffs: list[tuple[int, int, float]] | None = None   # similar pairs (compound a, compound b, similarity) with similarity >= CLIFF_MIN
    cliff_note: str = ""
    base: float = 0.0
    N: int = 0
    K: int = 0


def build_chemotypes(prep: Prepared, smiles_by_cpd: list[str], method: str, sim_cut: float, pairs: Pairs | None) -> ChemoResult:
    """Group compounds into chemotypes.

    Args:
        prep: Fingerprints and scaffolds from :func:`chem.prepare`.
        smiles_by_cpd: SMILES per compound (for representative structures).
        method: ``"butina"`` (Tanimoto clustering) or ``"murcko"`` (shared scaffold).
        sim_cut: Tanimoto cut-off for Butina clustering.
        pairs: Pre-computed similarity pairs; None if the data set is too large, which forces scaffold grouping.

    Returns:
        A :class:`ChemoResult`, chemotypes numbered CT-0001... from largest to smallest.
    """
    note = ""
    # 'ok' turns a fingerprint's position back into the compound number.
    ok = prep.idx
    if method == "butina" and pairs is None:
        method = "murcko"
        note = f"Similarity clustering is limited to {PAIR_LIMIT:,} compounds; {len(ok):,} compounds were grouped by scaffold instead."
    groups: list[tuple[list[int], str, str]] = []
    if method == "butina":
        # Butina gives clusters as fingerprint positions. The first member of each cluster is its centre.
        for m in butina(pairs, sim_cut):
            groups.append(([ok[i] for i in m], smiles_by_cpd[ok[m[0]]], "Cluster centroid"))
    else:
        # Grouping by scaffold: compounds with the same Murcko scaffold form one chemotype.
        # Compounds without rings have no scaffold, so each of them is a chemotype on its own.
        by: dict[str, list[int]] = {}
        for pos, ci in enumerate(ok):
            sc = prep.scaffold[pos]
            if not sc:
                groups.append(([ci], smiles_by_cpd[ci], "Acyclic compound"))
            else:
                by.setdefault(sc, []).append(ci)
        groups += [(v, k, "Murcko scaffold") for k, v in by.items()]
    # Number the chemotypes from the largest to the smallest. Chemotypes of the same size keep their order.
    groups.sort(key=lambda g: -len(g[0]))
    chemos = [Chemotype(id=f"CT-{i + 1:04d}", members=g[0], rep=g[1], rep_kind=g[2]) for i, g in enumerate(groups)]
    cliffs = None
    cliff_note = ""
    if pairs is not None:
        # Keep the similar pairs for the activity-cliff table, written with compound numbers.
        keep = pairs.s >= CLIFF_MIN - 1e-9
        cliffs = [(ok[a], ok[b], float(s)) for a, b, s in zip(pairs.i[keep].tolist(), pairs.j[keep].tolist(), pairs.s[keep].tolist())]
    else:
        cliff_note = f"Activity cliffs are limited to {PAIR_LIMIT:,} compounds with structures."
    return ChemoResult(chemos, {c.id: c for c in chemos}, {ci: c.id for c in chemos for ci in c.members}, method, note, cliffs, cliff_note)


def chemo_stats(res: ChemoResult, means: np.ndarray, cls: np.ndarray, n_classes: int) -> None:
    """Fill in counts, hit rate, mean/max response, enrichment factor, p and q for every chemotype.

    "Hits" are compounds in any class above Inactive. p is the hypergeometric probability of at least that many hits
    in a chemotype of that size; q applies the Benjamini-Hochberg correction across all chemotypes.
    """
    # Totals for the whole library: N compounds with a structure, K of them hits.
    # 'base' is the hit rate of the whole library.
    N = sum(len(c.members) for c in res.chemos)
    K = int(sum(int(cls[i] > 0) for c in res.chemos for i in c.members))
    base = K / N if N else 0.0
    for c in res.chemos:
        m = np.array(c.members)
        c.n = len(m)
        cc = cls[m]
        c.counts = [int((cc == j).sum()) for j in range(n_classes)]
        c.n_hit = c.n - c.counts[0]
        c.n_top = c.counts[-1]
        c.rate = c.n_hit / c.n
        c.mean = float(means[m].mean())
        c.max = float(means[m].max())
        c.members = sorted(c.members, key=lambda i: -means[i])
        # Enrichment factor: how much higher the chemotype's hit rate is than the library's hit rate.
        c.ef = c.rate / base if base > 0 else float("nan")
        # The p-value is the chance of finding at least n_hit hits in a random group of the same size
        # (hypergeometric distribution). scipy's sf(x) gives the chance of more than x, so we pass n_hit - 1.
        c.p = float(hypergeom.sf(c.n_hit - 1, N, K, c.n)) if c.n_hit > 0 else 1.0
    # Benjamini-Hochberg correction for testing many chemotypes at once. Sort the p-values, scale each by
    # m / rank, and go from the largest p down so that a q-value is never smaller than the one before it.
    order = sorted(res.chemos, key=lambda c: c.p)
    m = len(order)
    q = 1.0
    for r in range(m, 0, -1):
        c = order[r - 1]
        q = min(q, c.p * m / r)
        c.q = q
    res.base, res.N, res.K = base, N, K


def cliff_list(res: ChemoResult, means: np.ndarray, sim_min: float, delta_min: float, only: set[int] | None = None) -> list[dict]:
    """Return activity-cliff pairs at or above ``sim_min`` Tanimoto whose % response differs by at least ``delta_min``.

    Args:
        only: If given, keep only pairs with at least one compound in this set (e.g. the selected chemotype).

    Returns:
        Dicts with ``hi`` and ``lo`` (more and less active compound), ``sim``, ``d`` (Delta %) and ``sali``, sorted by SALI.
    """
    if not res.cliffs:
        return []
    out = []
    for a, b, s in res.cliffs:
        if s < sim_min - 1e-9:
            continue
        d = abs(means[a] - means[b])
        if d < delta_min:
            continue
        if only is not None and not (a in only or b in only):
            continue
        hi, lo = (a, b) if means[a] >= means[b] else (b, a)
        # SALI = difference in % response / (1 - Tanimoto). Identical fingerprints would divide by zero,
        # so the bottom of the fraction is never allowed below 0.01.
        out.append(dict(hi=hi, lo=lo, sim=s, d=float(d), sali=float(d / max(0.01, 1 - s))))
    out.sort(key=lambda x: -x["sali"])
    return out
