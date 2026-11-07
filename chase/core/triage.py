"""Triage decisions, kept in memory for one browser session only.

Nothing is written to disk by CHASE. Decisions leave the app only when the user saves a session file
or exports the triage CSV, and come back by opening that session or importing that CSV.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import csv
import io
import time
from datetime import datetime

# The triage decisions. Each entry is (name saved in files, button label, colour name used in the CSS).
STATUSES = [
    ("confirm", "Confirm", "good"),
    ("resynth", "Resynthesize", "warn"),
    ("analogues", "Buy analogues", "acc"),
    ("discard", "Discard", "bad"),
]
STATUS = {s[0]: s for s in STATUSES}
# When reading a CSV we accept both the label ('Buy analogues') and the short name ('analogues').
_BY_LABEL = {s[1].lower(): s[0] for s in STATUSES} | {s[0]: s[0] for s in STATUSES}

CSV_COLUMNS = ["campaign", "batch_id", "compound_name", "decision", "note", "pct_response", "activity_class", "chemotype", "decided_by", "decided_at"]


class TriageBoard:
    """In-memory decisions keyed by compound key (upper-case Batch ID)."""

    def __init__(self) -> None:
        """Start with an empty board and nothing unsaved."""
        self._d: dict[str, dict] = {}
        self.unsaved = 0          # number of changes since the board was last saved, exported or loaded

    def decisions(self) -> dict[str, dict]:
        """Return a copy of all decisions, keyed by compound key."""
        return {k: dict(v) for k, v in self._d.items()}

    def __len__(self) -> int:
        """Number of decided compounds."""
        return len(self._d)

    def set(self, d: dict) -> None:
        """Add or replace the decision for ``d["key"]`` (timestamped now if ``at`` is missing)."""
        d = dict(d)
        # Times are stored in milliseconds since 1970, like JavaScript and the HTML version's session files do.
        d.setdefault("at", time.time() * 1000)
        self._d[d["key"]] = d
        self.unsaved += 1

    def delete(self, key: str) -> None:
        """Remove the decision for ``key``, if any."""
        if self._d.pop(key, None) is not None:
            self.unsaved += 1

    def clear(self) -> None:
        """Empty the board (the user confirmed it), leaving nothing unsaved."""
        self._d.clear()
        self.unsaved = 0

    def mark_saved(self) -> None:
        """Record that the board has just been saved or loaded, so nothing is unsaved."""
        self.unsaved = 0

    def merge(self, ds: list[dict]) -> int:
        """Add decisions; for a compound already on the board the newer decision wins. Returns how many were taken."""
        n = 0
        for d in ds:
            # Skip broken entries instead of rejecting the whole import.
            if not isinstance(d, dict) or d.get("status") not in STATUS or not d.get("key"):
                continue
            cur = self._d.get(d["key"])
            # The newer decision wins. If both have the same time, the imported one is taken.
            if cur is None or (cur.get("at") or 0) <= (d.get("at") or 0):
                self._d[d["key"]] = dict(d)
                n += 1
        if n:
            self.unsaved += n
        return n


def to_csv(decisions: dict[str, dict], campaign: str) -> str:
    """Write the board as a triage CSV (columns in :data:`CSV_COLUMNS`), grouped by decision."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(CSV_COLUMNS)
    order = {s[0]: i for i, s in enumerate(STATUSES)}
    # Rows are grouped by decision (in the same order as the buttons), with the most active compound first
    # in each group. Times are written in local time, in ISO format.
    for d in sorted(decisions.values(), key=lambda d: (order.get(d["status"], 9), -(d.get("pct") or 0))):
        when = datetime.fromtimestamp(d["at"] / 1000).isoformat(timespec="seconds") if d.get("at") else ""
        w.writerow([campaign, d.get("bid", d["key"]), d.get("name", ""), STATUS.get(d["status"], (0, d["status"]))[1], d.get("note", ""),
                    "" if d.get("pct") is None else d["pct"], d.get("cls", ""), d.get("chemo", ""), d.get("by", ""), when])
    return buf.getvalue()


def from_csv(text: str) -> tuple[list[dict], list[str], str]:
    """Parse a triage CSV (as written by to_csv). Returns (decisions, problems, campaign)."""
    # Excel adds an invisible byte-order mark at the start when it saves a CSV as UTF-8. Remove it.
    text = text.lstrip("\ufeff")
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows or "batch_id" not in rows[0] or "decision" not in rows[0]:
        raise ValueError("This is not a CHASE triage file: it needs at least the columns batch_id and decision.")
    out, bad = [], []
    campaign = ""
    # Line numbers in the problem messages count the header as line 1, so data start at line 2.
    for i, r in enumerate(rows, start=2):
        bid = (r.get("batch_id") or "").strip()
        status = _BY_LABEL.get((r.get("decision") or "").strip().lower())
        if not bid or not status:
            bad.append(f"line {i}")
            continue
        # The campaign name is taken from the first row that has one.
        campaign = campaign or (r.get("campaign") or "").strip()
        at = None
        if (r.get("decided_at") or "").strip():
            try:
                # If the time is missing or can't be read, we use 'now', so the imported decision wins the merge.
                at = datetime.fromisoformat(r["decided_at"].strip()).timestamp() * 1000
            except ValueError:
                at = None
        try:
            pct = float(r["pct_response"]) if (r.get("pct_response") or "").strip() else None
        except ValueError:
            pct = None
        out.append(dict(key=bid.upper(), bid=bid, name=(r.get("compound_name") or "").strip(), status=status, note=r.get("note") or "", pct=pct,
                        cls=(r.get("activity_class") or "").strip(), chemo=(r.get("chemotype") or "").strip(), by=(r.get("decided_by") or "").strip(),
                        at=at if at is not None else time.time() * 1000))
    return out, bad, campaign
