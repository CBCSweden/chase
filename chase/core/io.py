"""Reading screening and structure files (xlsx/xls/csv/tsv/txt/smi) and guessing their columns.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import csv
import io
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Patterns used to recognise the four columns of a screening file by their header names (case is ignored).
RX = {
    "raw": re.compile(r"raw|signal|readout|value|rfu|count|intensity|lum|fluor|response", re.I),
    "plate": re.compile(r"plate", re.I),
    "well": re.compile(r"well", re.I),
    "bid": re.compile(r"batch|compound|sample|substance|\bid\b", re.I),
}


def id_str(v) -> str:
    """Batch IDs as text. Excel turns numeric-looking IDs into numbers; bring them back to plain integers."""
    if v is None:
        return ""
    # Excel stores IDs that look like numbers (e.g. 117983001) as floating-point numbers.
    # We turn 117983001.0 back into '117983001'.
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if v.is_integer():
            return str(int(v))
        return repr(v)
    # True and False count as integers in Python, so they have to be handled before integers.
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    return str(v).strip()


def norm_key(v) -> str:
    """Return the matching key of a Batch ID: text form, trimmed, upper case."""
    return id_str(v).strip().upper()


def alt_key(k: str) -> str:
    """Return ``k`` without leading zeros, so that ``00123`` and ``123`` can still be matched."""
    return re.sub(r"^0+(?=\d)", "", k)


def to_num(v) -> float:
    """Convert a cell to float; empty or non-numeric cells give NaN."""
    if v is None:
        return math.nan
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v)
    s = str(v).strip().replace(" ", "")
    if not s:
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def parse_delimited(text: str) -> list[list[str]]:
    """Parse CSV/TSV/semicolon or whitespace-separated text into rows, guessing the delimiter from the first line."""
    # Excel puts an invisible byte-order mark at the start of UTF-8 files. Remove it.
    text = text.lstrip("\ufeff")
    first = text.splitlines()[0] if text else ""
    best, sep = -1, ","
    # Guess the separator: whichever of comma, tab or semicolon splits the first line into the most fields.
    for c in (",", "\t", ";"):
        n = len(first.split(c))
        if n > best:
            best, sep = n, c
    # If none of them splits the line but it contains spaces, the columns are separated by spaces
    # (as in .smi files).
    if best <= 1 and re.search(r"\s", first.strip()):
        return [ln.split() for ln in text.splitlines() if ln.strip()]
    rows = list(csv.reader(io.StringIO(text), delimiter=sep))
    return [r for r in rows if any(c != "" for c in r)]


def read_table(path: str | Path, name: str | None = None) -> dict[str, list[list]]:
    """Return {sheet name: rows (lists of cells)} for any supported file."""
    path = Path(path)
    name = name or path.name
    if re.search(r"\.(csv|tsv|txt|smi)$", name, re.I):
        raw = path.read_bytes()
        # Try UTF-8 first. If that fails, use latin-1, which can read any file. Older Windows programs often
        # save files that way.
        for enc in ("utf-8-sig", "latin-1"):
            try:
                text = raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        return {name: parse_delimited(text)}
    # Excel file: read every sheet as raw cell values, so that IDs and numbers stay exactly as stored.
    engine = "xlrd" if name.lower().endswith(".xls") else "openpyxl"
    sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=object, engine=engine)
    out = {}
    for sname, df in sheets.items():
        rows = []
        for rec in df.itertuples(index=False):
            # pandas gives empty cells as NaN. We use None instead, and drop rows that are completely empty.
            row = [None if (isinstance(x, float) and math.isnan(x)) else x for x in rec]
            if any(x is not None and x != "" for x in row):
                rows.append(row)
        out[str(sname)] = rows
    return out


def find_header(rows: list[list]) -> int:
    """Return the index of the first row (within the first 30) that looks like a header, or -1."""
    for r, row in enumerate(rows[:30]):
        txt = [c for c in row if isinstance(c, str)]
        hits = sum(1 for rx in RX.values() if any(rx.search(t) for t in txt))
        # We call a row the header when at least two of the four column names are recognised in it.
        if hits >= 2:
            return r
    return -1


def guess_screen_cols(header: list[str]) -> dict[str, int]:
    """Guess which columns hold raw signal, plate, well and Batch ID from the header names.

    Returns:
        ``{"raw": i, "plate": j, "well": k, "bid": l}``; columns not recognised by name fall back to the order Raw, Plate, Well, Batch ID.
    """
    used: set[int] = set()

    def pick(rx, excl=None):
        """Return the first unused column whose header matches ``rx`` (and not ``excl``), marking it used; -1 if none."""
        for i, h in enumerate(header):
            h = str(h or "")
            if i not in used and re.search(rx, h, re.I) and not (excl and re.search(excl, h, re.I)):
                used.add(i)
                return i
        return -1

    # The order of these lines matters. 'Plate' must not pick up a column called '384 well ID', and 'Well'
    # has to be found before the general 'ID' pattern, which would otherwise take '384 well ID' as the Batch ID.
    c = {"plate": pick(r"plate", r"well"), "well": pick(r"well")}
    c["bid"] = pick(r"batch")
    if c["bid"] < 0:
        c["bid"] = pick(r"compound|sample|substance|\bid\b")
    c["raw"] = pick(RX["raw"].pattern)
    # Any column that was not recognised by its name is taken from the usual order: Raw, Plate, Well, Batch ID.
    for i, k in enumerate(["raw", "plate", "well", "bid"]):
        if c[k] < 0 and i < len(header) and i not in used:
            c[k] = i
            used.add(i)
    return c


@dataclass
class ScreenFile:
    """A loaded screening workbook or text file: all sheets, the chosen sheet and its header row."""
    file_name: str
    sheets: dict[str, list[list]]
    sheet: str = ""
    hdr: int = -1
    header: list[str] = field(default_factory=list)

    @classmethod
    def from_sheets(cls, sheets: dict[str, list[list]], file_name: str) -> "ScreenFile":
        """Create a ScreenFile, choosing the first sheet that has a recognisable header."""
        names = list(sheets)
        # Use the first sheet that has a recognisable header. Workbooks often also contain notes or summary sheets.
        sheet = next((n for n in names if len(sheets[n]) > 1 and find_header(sheets[n]) >= 0), names[0])
        sf = cls(file_name=file_name, sheets=sheets)
        sf.set_sheet(sheet)
        return sf

    def set_sheet(self, sheet: str) -> None:
        """Switch to ``sheet`` and detect its header row (generic column names if there is none)."""
        self.sheet = sheet
        rows = self.sheets[sheet]
        self.hdr = find_header(rows)
        # No header row found: call the columns 'Column 1', 'Column 2', ... and read data from the very first row.
        if self.hdr < 0:
            width = max([len(r) for r in rows[:20]] + [4])
            self.header = [f"Column {i + 1}" for i in range(width)]
        else:
            self.header = [str(h) if h not in (None, "") else f"Column {i + 1}" for i, h in enumerate(rows[self.hdr])]

    def guess(self) -> dict[str, int]:
        """Return the guessed column mapping for the current sheet (see :func:`guess_screen_cols`)."""
        return guess_screen_cols(self.header)

    def extract(self, cols: dict[str, int]) -> pd.DataFrame:
        """One row per well with a numeric readout and a plate: raw, plate, well, bid, key, ord."""
        out = []
        rows = self.sheets[self.sheet]
        get = lambda r, i: r[i] if 0 <= i < len(r) else None
        for i in range(self.hdr + 1, len(rows)):
            r = rows[i]
            raw = to_num(get(r, cols["raw"]))
            plate = get(r, cols["plate"])
            # Skip rows without a number or without a plate. These are blank lines and summary rows like 'Average'.
            if math.isnan(raw) or plate is None or str(plate).strip() == "":
                continue
            bid = id_str(get(r, cols["bid"]))
            # 'key' is the Batch ID in capitals and is used for matching. 'ord' remembers the order of the wells
            # in the file.
            out.append((raw, id_str(plate), id_str(get(r, cols["well"])), bid, bid.upper(), i))
        return pd.DataFrame(out, columns=["raw", "plate", "well", "bid", "key", "ord"])


def control_candidates(rows: pd.DataFrame) -> tuple[list[tuple[str, int]], str, str]:
    """Batch IDs repeated on most plates, plus a guess for the positive and negative control."""
    n_plates = rows["plate"].nunique()
    sub = rows[rows["key"] != ""]
    g = sub.groupby("key").agg(n=("raw", "size"), m=("raw", "mean"), bid=("bid", "first"))
    # Controls are on every plate, usually in several wells, so their IDs appear much more often than
    # any compound's ID.
    thr = max(3, 2 * n_plates)
    cands = g[g["n"] >= thr].sort_values("n", ascending=False).head(30)
    keys = list(cands.index)
    # First try to recognise the controls by their names...
    pos = next((k for k in keys if re.search(r"CTRL|POS|REF|HIGH|MAX|CONTROL|STAURO|INHIB", k) and "NEG" not in k), None)
    neg = next((k for k in keys if re.search(r"DMSO|NEG|VEH|BLANK|LOW|MIN|BUFFER", k)), None)
    # ...and if that doesn't work, take the two most frequent IDs and call the one with the higher signal
    # the positive control. The user checks this guess and can change it before normalizing.
    if (pos is None or neg is None) and len(keys) >= 2:
        a, b = sorted(keys[:2], key=lambda k: -cands.loc[k, "m"])
        if pos is None and a != neg:
            pos = a
        if neg is None and b != pos:
            neg = b
    disp = lambda k: cands.loc[k, "bid"] if k in cands.index else ""
    return [(disp(k), int(cands.loc[k, "n"])) for k in keys], (disp(pos) if pos else ""), (disp(neg) if neg else "")


@dataclass
class StructFile:
    """A loaded structures file: header, body rows and the indices of the ID, SMILES and display-name columns."""
    file_name: str
    header: list[str]
    body: list[list]
    smi: int
    id: int = -1
    name: int = -1

    @classmethod
    def from_sheets(cls, sheets: dict[str, list[list]], file_name: str, screen_keys: set[str] | None = None) -> "StructFile":
        """Read a structures table, detecting the header, the SMILES column and the ID column.

        Args:
            sheets: Output of :func:`read_table`; only the first sheet is used.
            file_name: Original file name (``.smi`` files are read without a header).
            screen_keys: Batch-ID keys of the screen; when given, the ID column is the one that matches most of them.
        """
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        rows = next(iter(sheets.values()))
        if not rows:
            raise ValueError("The structures file is empty.")
        header = [str(h if h is not None else "") for h in rows[0]]

        def looks_smiles(v: str) -> bool:
            """Return True if ``v`` parses as SMILES and does not look like a plain compound ID."""
            v = str(v or "").strip()
            if len(v) < 2 or not re.search(r"[A-Za-z]", v) or re.fullmatch(r"[A-Z]{2,}\d+[A-Z]?", v):
                return False
            return Chem.MolFromSmiles(v, sanitize=False) is not None

        # If the first row already contains a SMILES, there is no header row.
        has_header = not any(len(h) > 3 and re.search(r"[\[\]=()#@]|c1|C1", h) and looks_smiles(h) for h in header)
        if file_name.lower().endswith(".smi") and "smiles" not in " ".join(header).lower():
            has_header = False
        if not has_header:
            header = [f"Column {i + 1}" for i in range(len(header))]
        body = rows[1:] if has_header else rows
        # Find the SMILES column by its header name. If no header says SMILES, take the column whose first
        # 40 values are most often valid SMILES.
        smi = next((i for i, h in enumerate(header) if re.fullmatch(r"\s*(smiles|mol|structure|smi)\s*", h, re.I) or re.search(r"smiles", h, re.I)), -1)
        if smi < 0:
            best, best_n = -1, 0
            for c in range(len(header)):
                ok = sum(1 for r in body[:40] if c < len(r) and looks_smiles(r[c]))
                if ok > best_n:
                    best, best_n = c, ok
            smi = best
        if smi < 0:
            raise ValueError("Could not find a SMILES column in the structures file.")
        sf = cls(file_name=file_name, header=header, body=body, smi=smi)
        sf.auto_id(screen_keys)
        # An optional column with a compound name to show (for example 'Compound id' with CBK numbers).
        sf.name = next((i for i, h in enumerate(header) if i not in (smi, sf.id) and re.search(r"compound|name|title|cbk|label", h, re.I)), -1)
        return sf

    def auto_id(self, screen_keys: set[str] | None) -> None:
        """Choose the ID column: the one matching most screen Batch IDs, else a header containing 'batch' or 'id'."""
        best, pick = -1, -1
        for c, h in enumerate(self.header):
            if c == self.smi:
                continue
            # The best clue: the column whose values match the most Batch IDs of the screening file already loaded.
            if screen_keys:
                score = 0
                for r in self.body:
                    k = norm_key(r[c] if c < len(r) else "")
                    if k in screen_keys or alt_key(k) in screen_keys:
                        score += 1
            else:
                score = 2 if re.search("batch", h, re.I) else 1 if re.search("id", h, re.I) else 0
            if score > best:
                best, pick = score, c
        self.id = pick

    def mapping(self) -> dict[str, dict]:
        """Return ``{key: {"id", "smiles", "name"}}`` for every row, indexed by key and by key without leading zeros."""
        m: dict[str, dict] = {}
        get = lambda r, i: r[i] if 0 <= i < len(r) else None
        for r in self.body:
            bid = id_str(get(r, self.id))
            smiles = str(get(r, self.smi) or "").strip()
            if not bid or not smiles:
                continue
            rec = {"id": bid, "smiles": smiles, "name": id_str(get(r, self.name)) if self.name >= 0 else ""}
            k = bid.upper()
            # Store each structure under its ID and also under the ID without leading zeros. If an ID appears
            # twice, the first one is kept.
            m.setdefault(k, rec)
            a = alt_key(k)
            if a != k:
                m.setdefault(a, rec)
        return m

    def rows(self) -> list[list[str]]:
        """Return ``[id, name, smiles]`` for every row with an ID and a SMILES (the format stored in session files)."""
        get = lambda r, i: r[i] if 0 <= i < len(r) else None
        out = []
        for r in self.body:
            bid, smi = id_str(get(r, self.id)), str(get(r, self.smi) or "").strip()
            if bid and smi:
                out.append([bid, id_str(get(r, self.name)) if self.name >= 0 else "", smi])
        return out
