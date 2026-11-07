"""Structural alerts: PAINS (RDKit FilterCatalog, Baell & Holloway 2010) and covalent warheads (editable SMARTS list).

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# The list of covalent warheads that comes with CHASE. It is a tab-separated text file with
# three columns: name, SMARTS pattern and a short note.
DEFAULT_COVALENT = Path(__file__).with_name("covalent_alerts.tsv")


def covalent_file() -> Path:
    """Return the covalent-alert file to use: ``$CHASE_COVALENT``, else ``~/.chase/covalent_alerts.tsv``, else the bundled list."""
    env = os.environ.get("CHASE_COVALENT")
    if env and Path(env).expanduser().is_file():
        return Path(env).expanduser()
    user = Path.home() / ".chase" / "covalent_alerts.tsv"
    return user if user.is_file() else DEFAULT_COVALENT


@dataclass(frozen=True)
class Alert:
    """One covalent-warhead alert: display name, SMARTS, a short note and the compiled RDKit query."""
    name: str
    smarts: str
    note: str
    pattern: Chem.Mol


@lru_cache(maxsize=1)
def covalent_alerts() -> tuple[tuple[Alert, ...], tuple[str, ...], str]:
    """(alerts, problems, source path)."""
    path = covalent_file()
    out, bad = [], []
    for ln in path.read_text(encoding="utf-8").splitlines():
        # Blank lines and lines starting with # are ignored. A line with a missing or broken SMARTS is
        # reported to the user, but it does not stop the rest of the file from loading.
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            bad.append(ln)
            continue
        patt = Chem.MolFromSmarts(parts[1].strip())
        if patt is None:
            bad.append(parts[0])
            continue
        out.append(Alert(parts[0].strip(), parts[1].strip(), parts[2].strip() if len(parts) > 2 else "", patt))
    return tuple(out), tuple(bad), str(path)


@lru_cache(maxsize=1)
def pains_catalog() -> FilterCatalog:
    """Return RDKit's PAINS filter catalogue (families A, B and C, 480 patterns), built once per process."""
    p = FilterCatalogParams()
    for c in (FilterCatalogParams.FilterCatalogs.PAINS_A, FilterCatalogParams.FilterCatalogs.PAINS_B, FilterCatalogParams.FilterCatalogs.PAINS_C):
        p.AddCatalog(c)
    return FilterCatalog(p)


def flags(mol: Chem.Mol) -> tuple[list[str], list[str]]:
    """Names of matching PAINS families and covalent warheads."""
    # One PAINS family can match the same molecule several times. We list each family only once.
    pains = sorted({e.GetDescription() for e in pains_catalog().GetMatches(mol)})
    cov = [a.name for a in covalent_alerts()[0] if mol.HasSubstructMatch(a.pattern)]
    return pains, cov


def alert_atoms(smiles: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Atom indices matched by PAINS and by covalent alerts (for highlighting); cached for the current session only."""
    from .cache import session_cache
    c = session_cache("alert_atoms")
    if c is not None and smiles in c:
        return c[smiles]
    out = _alert_atoms(smiles)
    if c is not None:
        c[smiles] = out
    return out


def _alert_atoms(smiles: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Compute the atom indices matched by PAINS and covalent alerts (uncached helper of :func:`alert_atoms`)."""
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return (), ()
    pa: set[int] = set()
    for e in pains_catalog().GetMatches(m):
        for fm in e.GetFilterMatches(m):
            # atomPairs gives (atom in the pattern, atom in the molecule). We only need the molecule atoms.
            pa.update(a for _, a in fm.atomPairs)
    ca: set[int] = set()
    for a in covalent_alerts()[0]:
        for match in m.GetSubstructMatches(a.pattern):
            ca.update(match)
    return tuple(sorted(pa)), tuple(sorted(ca))


def flag_many(smiles: list[str], progress=None, chunk: int = 400) -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    """PAINS and covalent flags for many SMILES (multithreaded PAINS screen), cached for the current session only."""
    from rdkit.Chem.FilterCatalog import RunFilterCatalog
    from .cache import session_cache

    _CACHE = session_cache("flags")
    # When CHASE is used outside the app (in the tests or a script) there is no session to store results
    # in, so we use a temporary dictionary instead.
    if _CACHE is None:
        _CACHE = {}
    # Every distinct SMILES is checked only once, and SMILES already checked earlier in this session are skipped.
    todo = [s for s in dict.fromkeys(smiles) if s and s not in _CACHE]
    cat = pains_catalog()
    cov = covalent_alerts()[0]
    for k in range(0, len(todo), chunk):
        part = todo[k:k + chunk]
        if progress:
            progress(k / max(1, len(todo)), f"Checking PAINS and covalent alerts {k:,} / {len(todo):,}")
        # Run the PAINS filters on this chunk using all CPU cores (numThreads=0 means 'all').
        # We work in chunks only so the progress bar can move while this runs.
        res = RunFilterCatalog(cat, part, numThreads=0)
        for s, entries in zip(part, res):
            m = Chem.MolFromSmiles(s)
            if m is None:
                _CACHE[s] = ((), ())
                continue
            p = tuple(sorted({e.GetDescription() for e in entries}))
            # Covalent alerts are simple substructure searches with the SMARTS patterns.
            c = tuple(a.name for a in cov if m.HasSubstructMatch(a.pattern))
            _CACHE[s] = (p, c)
    return [_CACHE.get(s, ((), ())) if s else ((), ()) for s in smiles]
