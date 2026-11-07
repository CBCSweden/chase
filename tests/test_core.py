"""Self-checks for the CHASE analysis core.

    python -m tests.test_core                         # built-in example data
    python -m tests.test_core SCREEN.xlsx STRUCT.csv  # plus your own files
Also collected by pytest if you have it installed.
Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# Let Python find the 'chase' package when the tests are run from the project folder
# without installing CHASE.
sys.path.insert(0, str(ROOT))

from chase.core import chem, classes, demo, io, normalize, plates, session, stats  # noqa: E402
from chase.core import triage  # noqa: E402


def _run(screen_sheets, screen_name, struct_sheets, struct_name):
    """Run the whole analysis (controls, normalization, classes, fingerprints, chemotypes) on one data set."""
    sf = io.ScreenFile.from_sheets(screen_sheets, screen_name)
    rows = sf.extract(sf.guess())
    _, pos, neg = io.control_candidates(rows)
    stf = io.StructFile.from_sheets(struct_sheets, struct_name, set(rows["key"]))
    an = normalize.analyze(rows, pos, neg, "plate", "mean", stf.mapping())
    cl = classes.assign(an.cpds["mean"].to_numpy(), "kmeans", 3)
    prep = chem.prepare(an.cpds["smiles"].tolist())
    pairs = chem.sim_pairs(prep.fps, 0.5)
    res = stats.build_chemotypes(prep, an.cpds["smiles"].tolist(), "butina", 0.5, pairs)
    stats.chemo_stats(res, an.cpds["mean"].to_numpy(), cl.cls, len(cl.labels))
    return sf, rows, pos, neg, stf, an, cl, prep, pairs, res


def test_example_data():
    """The built-in example gives the known controls, class counts, Z', enrichment and replicate correlation."""
    s, st = demo.make_demo()
    sf, rows, pos, neg, stf, an, cl, prep, pairs, res = _run({"Example": s}, "example.xlsx", {"s": st}, "s.csv")
    # The example data are made with a fixed seed, so we know exactly what to expect:
    # two 384-well plates with every compound tested twice.
    assert (pos, neg) == ("CTRL", "DMSO")
    assert len(an.cpds) == 352 and len(an.plates) == 2
    assert cl.counts == [236, 89, 27], cl.counts
    z = np.nanmedian([p["z"] for p in an.plates])
    assert abs(z - 0.68) < 0.01, z
    assert all(n > 0 for n in [len(res.chemos), len(res.cliffs or [])])
    # The chemotype with the lowest p-value has to be clearly enriched in hits.
    top = min((c for c in res.chemos if c.n > 1 and c.n_hit), key=lambda c: c.p)
    assert top.q < 1e-10 and top.ef > 1
    rep = plates.replicates(an, 20)
    assert len(rep["pts"]) == 352 and abs(rep["r"] - 0.965) < 0.01
    assert plates.plate_dims(an.wells) == (16, 24)
    assert chem.svg(an.cpds.at[0, "smiles"]).lstrip().startswith("<svg")


def test_butina_matches_rdkit():
    """CHASE's Butina clustering gives exactly the same clusters as rdkit.ML.Cluster.Butina."""
    from rdkit.ML.Cluster import Butina
    from rdkit import DataStructs
    s, st = demo.make_demo()
    *_, prep, pairs, res = _run({"Example": s}, "example.xlsx", {"s": st}, "s.csv")
    fps = prep.fps
    dist = []
    # RDKit's own Butina function wants the distances (1 - Tanimoto) as a flat lower-triangle list.
    for i in range(1, len(fps)):
        dist.extend(1 - x for x in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]))
    ref = Butina.ClusterData(dist, len(fps), 0.5, isDistData=True)
    ours = chem.butina(pairs, 0.5)
    # Compare the clusters as sets, so the order of the clusters and of their members doesn't matter.
    key = lambda cl: sorted(tuple(sorted(c)) for c in cl)
    assert key(ref) == key(ours)


def test_session_roundtrip_and_triage():
    """Session files round-trip; the triage board, its CSV export/import and merge rules work."""
    s, st = demo.make_demo()
    sf, rows, *_ = _run({"Example": s}, "example.xlsx", {"s": st}, "s.csv")
    txt = session.build("test", "example.xlsx", rows[["raw", "plate", "well", "bid"]].values.tolist(), "s.csv", st[1:], {"pos": "CTRL", "neg": "DMSO"},
                        True, {}, "Example", [])
    d = session.parse(txt)
    assert d["app"] == "CHASE" and len(d["screen"]["rows"]) == len(rows)
    b = triage.TriageBoard()
    b.set(dict(key="EX100001", bid="EX100001", status="confirm", note='ok, "quoted", comma', by="me", pct=12.5))
    b.set(dict(key="EX100008", bid="EX100008", status="discard", note="", by="you"))
    assert b.unsaved == 2 and len(b) == 2
    text = triage.to_csv(b.decisions(), "Example")
    ds, bad, camp = triage.from_csv(text)
    assert not bad and camp == "Example" and {d["key"]: d["status"] for d in ds} == {"EX100001": "confirm", "EX100008": "discard"}
    assert next(d for d in ds if d["key"] == "EX100001")["note"] == 'ok, "quoted", comma'
    b2 = triage.TriageBoard()
    assert b2.merge(ds) == 2
    # An older decision must not overwrite a newer one.
    older = dict(ds[0], status="resynth", at=1.0)
    assert b2.merge([older]) == 0 and b2.decisions()[older["key"]]["status"] == ds[0]["status"]
    b.mark_saved()
    b.delete("EX100001")
    assert b.unsaved == 1 and len(b) == 1


def test_alerts():
    """PAINS and covalent alerts flag reference compounds correctly and the batched check agrees."""
    from rdkit import Chem
    from chase.core import alerts
    al, bad, _ = alerts.covalent_alerts()
    assert len(al) >= 30 and not bad
    # For each test SMILES: the PAINS families and covalent alerts we expect to find.
    # Imatinib should trigger nothing.
    cases = {
        "C=CC(=O)N1CCC[C@@H](C1)n1nc(-c2ccc(Oc3ccccc3)cc2)c2c(N)ncnc21": ([], ["Acrylamide-type (alpha,beta-unsaturated amide)"]),   # ibrutinib
        "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O": ([], ["Boronic acid / ester"]),                       # bortezomib
        "O=C1N(CC(=O)O)C(=S)S/C1=C/c1ccccc1": (["ene_rhod_A(235)"], []),                                                # rhodanine
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1": ([], []),                                         # imatinib
    }
    for smi, (p, c) in cases.items():
        got = alerts.flags(Chem.MolFromSmiles(smi))
        assert got == (p, c), (smi, got)
    many = alerts.flag_many(list(cases))
    assert [list(x[1]) for x in many] == [c for _, c in cases.values()]


def check_files(screen: str, struct: str) -> None:
    """Run the analysis on the user's own files and print a short summary."""
    sf, rows, pos, neg, stf, an, cl, prep, pairs, res = _run(io.read_table(screen), Path(screen).name, io.read_table(struct), Path(struct).name)
    print(f"  {Path(screen).name}: {len(rows):,} wells, {len(an.plates)} plates, controls {pos} / {neg}, median Z' {np.nanmedian([p['z'] for p in an.plates]):.2f}")
    print(f"  {len(an.cpds):,} compounds, {(an.cpds['smiles'] != '').sum():,} with structures; classes {dict(zip([l for l, _ in cl.labels], cl.counts))}")
    print(f"  {len(res.chemos):,} chemotypes, {sum(c.n_hit > 0 for c in res.chemos):,} with hits; {len(res.cliffs or []):,} pairs with Tanimoto >= 0.6")


if __name__ == "__main__":
    # Run every function whose name starts with test_, so pytest is not needed.
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    if len(sys.argv) == 3:
        print("Your files:")
        check_files(sys.argv[1], sys.argv[2])
    print("All checks passed.")
