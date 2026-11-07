"""RDKit chemistry: fingerprints, Murcko scaffolds, similarity pairs, Butina clustering, depictions.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold

# Turn off RDKit's own warnings. When a SMILES can't be read we report it for that compound instead.
RDLogger.DisableLog("rdApp.*")
FP_RADIUS, FP_BITS = 3, 2048   # Morgan fingerprints with radius 3 (similar to ECFP6), 2048 bits
# Creating the fingerprint generator takes some time, so we make one and reuse it everywhere.
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


@dataclass
class Prepared:
    """Per-compound chemistry for the compounds that have a SMILES (indexes into Analysis.cpds)."""
    idx: list[int] = field(default_factory=list)         # compound numbers of the compounds whose SMILES could be read
    fps: list = field(default_factory=list)
    scaffold: list[str] = field(default_factory=list)    # Murcko scaffold as canonical SMILES ('' when the molecule has no rings)
    errors: dict[int, str] = field(default_factory=dict)


def prepare(smiles_by_cpd: list[str], progress: Callable[[float, str], None] | None = None) -> Prepared:
    """Parse SMILES and compute fingerprints and Murcko scaffolds for every compound that has a structure.

    Args:
        smiles_by_cpd: One SMILES per compound (``""`` when the compound has no structure).
        progress: Optional callback ``progress(fraction, message)``.

    Returns:
        A :class:`Prepared` object; ``Prepared.idx[k]`` is the compound index of the k-th fingerprint.
    """
    out = Prepared()
    n = len(smiles_by_cpd)
    for i, smi in enumerate(smiles_by_cpd):
        # Update the progress bar every 250 structures. More often would only slow the loop down.
        if progress and i % 250 == 0:
            progress(i / max(n, 1), f"Reading structures {i:,} / {n:,}")
        if not smi:
            continue
        # If a SMILES can't be read, the compound gets no fingerprint and is left out of the chemotypes.
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out.errors[i] = "SMILES could not be parsed"
            continue
        out.idx.append(i)
        out.fps.append(_GEN.GetFingerprint(m))
        try:
            # The Murcko scaffold keeps the ring systems and the linkers between them and drops the side chains.
            # We write it without stereochemistry, so stereoisomers end up with the same scaffold.
            # Molecules without rings get an empty scaffold.
            sc = MurckoScaffold.GetScaffoldForMol(m)
            out.scaffold.append(Chem.MolToSmiles(sc, isomericSmiles=False) if sc.GetNumAtoms() else "")
        except Exception:  # noqa: BLE001 - scaffold failures fall back to acyclic
            out.scaffold.append("")
    return out


@dataclass
class Pairs:
    """All compound pairs above a Tanimoto threshold, stored as parallel arrays (i < j)."""
    n: int
    min_sim: float
    i: np.ndarray
    j: np.ndarray
    s: np.ndarray


def sim_pairs(fps: list, min_sim: float, progress: Callable[[float, str], None] | None = None) -> Pairs:
    """All pairs (i < j) with Tanimoto >= min_sim."""
    n = len(fps)
    I, J, Sv = [], [], []
    # Total number of comparisons, n*(n-1)/2, used only for the progress bar. 'or 1' avoids dividing by zero.
    total = n * (n - 1) / 2 or 1
    done = 0
    for i in range(n - 1):
        # Compare fingerprint i with all fingerprints after it in one fast RDKit call, and keep only the pairs
        # above the threshold. Memory then grows with the number of similar pairs instead of with n^2.
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]), dtype=np.float32)
        # The tiny tolerance makes sure a pair that sits exactly on the threshold is not lost to rounding.
        hit = np.nonzero(sims >= min_sim - 1e-9)[0]
        if len(hit):
            I.append(np.full(len(hit), i, dtype=np.int32))
            # 'hit' counts positions inside fps[i+1:], so we add i+1 to get the real fingerprint numbers.
            J.append((hit + i + 1).astype(np.int32))
            Sv.append(sims[hit])
        done += n - i - 1
        if progress and i % 200 == 0:
            progress(done / total, f"Comparing {n:,} fingerprints")
    # Join the pieces collected above. If no pair reached the threshold, return empty arrays.
    cat = lambda a, t: np.concatenate(a) if a else np.zeros(0, dtype=t)
    return Pairs(n, min_sim, cat(I, np.int32), cat(J, np.int32), cat(Sv, np.float32))


def butina(pairs: Pairs, sim_cut: float) -> list[list[int]]:
    """Taylor-Butina clustering (as in rdkit.ML.Cluster.Butina); the first member of each cluster is its centroid."""
    n = pairs.n
    nbrs: list[list[int]] = [[] for _ in range(n)]
    # The pair list may have been made with a lower threshold than the one asked for here, because it is
    # reused when the user raises the cut-off. So I first keep only the pairs at or above this cut-off,
    # and then build each compound's list of neighbours from them.
    keep = pairs.s >= sim_cut - 1e-9
    for a, b in zip(pairs.i[keep].tolist(), pairs.j[keep].tolist()):
        nbrs[a].append(b)
        nbrs[b].append(a)
    # How Butina clustering works: go through the compounds starting with the one that has the most neighbours.
    # If a compound is not in a cluster yet, it becomes the centre of a new cluster, and all its neighbours
    # that are still free join it. When two compounds have the same number of neighbours, the one with the
    # higher number goes first. From my understanding RDKit does the same, so our clusters should be identical to RDKit's.
    # As far as I know rdkit.ML.Cluster.Butina.ClusterData needs the distance for every pair of compounds, 
    # so I had to make this function to take only the pairs above the similarity cut-off. Worth checking if possible with RDkit in the next version? 
    order = sorted(range(n), key=lambda x: (len(nbrs[x]), x), reverse=True)
    assigned = np.zeros(n, dtype=bool)
    clusters = []
    for c in order:
        if assigned[c]:
            continue
        members = [c]
        assigned[c] = True
        for j in nbrs[c]:
            if not assigned[j]:
                assigned[j] = True
                members.append(j)
        clusters.append(members)
    return clusters


def svg(smiles: str, w: int = 280, h: int = 210, hl_pains: tuple = (), hl_cov: tuple = ()) -> str:
    """SVG depiction, cached for the current session only."""
    from .cache import session_cache
    c = session_cache("svg")
    key = (smiles, w, h, hl_pains, hl_cov)
    # Outside the app (tests, scripts) there is no cache, so the drawing is made each time.
    if c is not None and key in c:
        return c[key]
    out = _svg(smiles, w, h, hl_pains, hl_cov)
    if c is not None:
        c[key] = out
    return out


def _svg(smiles: str, w: int, h: int, hl_pains: tuple, hl_cov: tuple) -> str:
    """Draw ``smiles`` as an SVG string with RDKit, optionally highlighting PAINS (red) and covalent (violet) atoms."""
    if not smiles:
        return '<span class="na">No structure</span>'
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        # Some SMILES fail RDKit's checks (e.g., wrong valence) but can still be drawn if we skip them.
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is None:
            return '<span class="na">Could not draw</span>'
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:  # noqa: BLE001
            return '<span class="na">Could not draw</span>'
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    o = d.drawOptions()
    # No background, a little padding and no stereo labels, so the drawing blends into the card it sits on.
    o.clearBackground = False
    o.padding = 0.06
    o.bondLineWidth = 1
    o.minFontSize = 11
    o.addStereoAnnotation = False
    hl_atoms, atom_cols, hl_bonds, bond_cols = [], {}, [], {}
    if hl_pains or hl_cov:
        # Highlight colours as RGB values between 0 and 1: light red for PAINS atoms, light violet for covalent
        # warheads. Covalent atoms are coloured second, so they win where the two alerts overlap.
        pc, cc = (1.0, 0.72, 0.68), (0.80, 0.76, 1.0)
        for a in hl_pains:
            atom_cols[a] = pc
        for a in hl_cov:
            atom_cols[a] = cc
        # A bond is highlighted when both of its atoms are. It is violet only if both atoms belong to a covalent
        # warhead, otherwise it is red.
        hl_atoms = list(atom_cols)
        for b in m.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i in atom_cols and j in atom_cols:
                hl_bonds.append(b.GetIdx())
                bond_cols[b.GetIdx()] = cc if (i in hl_cov and j in hl_cov) else pc
    try:
        if hl_atoms:
            rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=hl_atoms, highlightAtomColors=atom_cols, highlightBonds=hl_bonds, highlightBondColors=bond_cols)
        else:
            rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    except Exception:  # noqa: BLE001
        return '<span class="na">Could not draw</span>'
    d.FinishDrawing()
    s = d.GetDrawingText()
    # RDKit writes a complete SVG file. To put it inside a web page need to remove the XML header line and scale the picture to the size of the box it is placed in.
    s = s.replace("<?xml version='1.0' encoding='iso-8859-1'?>", "")
    s = s.replace(f"width='{w}px' height='{h}px'", "width='100%' height='100%'", 1)
    return s
