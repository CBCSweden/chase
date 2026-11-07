"""Generated example screen (two replicate 384-well plates, five analogue series plus reference drugs).

Same generator and random sequence as the HTML version of CHASE, so both show identical example data.
Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import math

# Used to keep numbers within 32 bits, which is what JavaScript does with bit operations.
# This way the Python and the HTML version produce exactly the same random numbers.
M32 = 0xFFFFFFFF


def _imul(a: int, b: int) -> int:
    """32-bit integer multiplication, as JavaScript's ``Math.imul``."""
    return (a * b) & M32


class _Rng:
    """Mulberry32 random generator, identical to the one in the HTML version so both build the same example data."""
    def __init__(self, seed: int):
        """Start the generator from ``seed``."""
        self.s = seed & M32

    def rnd(self) -> float:
        """Return the next uniform random number in [0, 1)."""
        self.s = (self.s + 0x6D2B79F5) & M32
        s = self.s
        t = _imul(s ^ (s >> 15), 1 | s)
        t = ((t + _imul(t ^ (t >> 7), 61 | t)) & M32) ^ t
        return ((t ^ (t >> 14)) & M32) / 4294967296

    def gauss(self) -> float:
        """Return a standard normal random number (Box-Muller transform)."""
        u = v = 0.0
        while not u:
            u = self.rnd()
        while not v:
            v = self.rnd()
        return math.sqrt(-2 * math.log(u)) * math.cos(2 * math.pi * v)


# The example library: five series of 8 x 8 analogues and 32 known drugs with a made-up activity.
# Each series has a core structure, a list of R2 groups and a rule that gives the 'true' % response
# for R1 group i and R2 group j. Series 1 and 5 contain active chemotypes with clear SAR, series 2 is
# weakly active, and series 3 and 4 are inactive.
R1 = ["F", "Cl", "Br", "C", "OC", "C(F)(F)F", "C#N", "N(C)C"]
SERIES = [
    (lambda a, b: f"{b}c1ccc(cc1)S(=O)(=O)Nc1nc(cs1)-c1ccc({a})cc1", ["C", "Cl", "F", "Br", "CC", "N#C", "CO", "CC(C)"],
     lambda i, j: max(0, [78, 85, 88, 70, 45, 92, 60, 35][i] - j * 11)),
    (lambda a, b: f"COc1cc2ncnc(Nc3cccc({a})c3)c2cc1{b}", ["OC", "OCC", "OCCOC", "OCCN1CCOCC1", "Cl", "C", "F", "OC(C)C"],
     lambda i, j: 18 + ((i * 7 + j * 5) % 22)),
    (lambda a, b: f"O=C(c1ccc({a})cc1)N1CCN(CC1)c1ncccc1{b}", ["C", "F", "Cl", "OC", "C#N", "C(F)(F)F", "CC", "Br"], lambda i, j: 0),
    (lambda a, b: f"Cn1c(CNC(=O)c2ccc({a})cc2)nc2cc({b})ccc12", ["C", "F", "Cl", "OC", "C#N", "C(F)(F)F", "CC", "Br"], lambda i, j: 1),
    (lambda a, b: f"O=c1cc(-c2ccc({a})cc2)oc2cc({b})ccc12", ["O", "OC", "F", "Cl", "C", "N(C)C", "OCC", "Br"],
     lambda i, j: 72 - i * 4 if j == 0 else 2),
]
DRUGS = [
    ("Metformin", "CN(C)C(=N)NC(=N)N", 0), ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", 0), ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", 0),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 0), ("Naproxen", "COc1ccc2cc(ccc2c1)C(C)C(=O)O", 0), ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12", 0),
    ("Lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C", 0), ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", 0), ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", 0),
    ("Fluoxetine", "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1", 0), ("Imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1", 64),
    ("Gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", 38), ("Celecoxib", "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F", 0),
    ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", 0), ("Benzocaine", "CCOC(=O)c1ccc(N)cc1", 0), ("Carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc12", 0),
    ("Chlorpromazine", "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc12", 5), ("Haloperidol", "O=C(CCCN1CCC(O)(CC1)c1ccc(Cl)cc1)c1ccc(F)cc1", 0),
    ("Indomethacin", "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1", 0), ("Diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", 0),
    ("Omeprazole", "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C", 0), ("Losartan", "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1", 0),
    ("Tamoxifen", "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1", 6), ("Resveratrol", "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1", 46),
    ("Quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", 91), ("Sulfamethoxazole", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", 0),
    ("Trimethoprim", "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", 0), ("Ciprofloxacin", "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", 0),
    ("Phenytoin", "O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1", 0), ("Nifedipine", "COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC", 0),
    ("Curcumin", "COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O", 77), ("Genistein", "O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12", 29),
]


def make_demo() -> tuple[list[list], list[list]]:
    """Return (screen rows incl. header, structure rows incl. header)."""
    # A fixed seed, so the example data are the same every time and in both versions of CHASE.
    rng = _Rng(20261104)
    cpds = []
    for t, r2, act in SERIES:
        for i, a in enumerate(R1):
            for j, b in enumerate(r2):
                cpds.append(dict(smiles=t(a, b), act=act(i, j), name=None))
    for n, smi, act in DRUGS:
        cpds.append(dict(smiles=smi, act=act, name=n))
    for i, c in enumerate(cpds):
        # Batch IDs EX100001, EX100008, ... and names EXC-0001, EXC-0002, ... for the analogues.
        c["bid"] = "EX" + str(100001 + i * 7).zfill(6)
        c["name"] = c["name"] or f"EXC-{i + 1:04d}"
    rows = []
    letters = "ABCDEFGHIJKLMNOP"
    for plate in (1, 2):
        # Every plate contains every compound once, each plate in a different random order (Fisher-Yates shuffle).
        order = list(range(len(cpds)))
        for i in range(len(order) - 1, 0, -1):
            j = math.floor(rng.rnd() * (i + 1))
            order[i], order[j] = order[j], order[i]
        k = 0
        # The raw signal of the controls shifts a little from plate to plate, as in a real screen.
        # Columns 23 and 24 hold the controls (DMSO and CTRL).
        neg, pos = 340 + plate * 25, 5900 - plate * 180
        for r in range(16):
            for col in range(1, 25):
                well = f"{letters[r]}{col:02d}"
                if col == 23:
                    bid, raw = "DMSO", neg + rng.gauss() * 70
                elif col == 24:
                    bid, raw = "CTRL", pos + rng.gauss() * 380
                else:
                    c = cpds[order[k]]
                    k += 1
                    bid = c["bid"]
                    # Measured % = true activity plus some noise (more noise for active compounds). We then turn it back
                    # into a raw signal and add a bit of reading noise on top.
                    p = max(-5, c["act"] + rng.gauss() * (6 if c["act"] > 10 else 2.2))
                    raw = neg + (p / 100) * (pos - neg) + rng.gauss() * 55
                rows.append([math.floor(max(20, raw) + 0.5), plate, well, bid])
    screen = [["Raw Data", "Plate Number", "384 well ID", "Batch ID"]] + rows
    struct = [["Batch id", "Compound id", "mol"]] + [[c["bid"], c["name"], c["smiles"]] for c in cpds]
    return screen, struct
