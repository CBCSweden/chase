"""Session files (.json), interchangeable with the HTML version of CHASE (format 1).

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone

# Version of the session file format. The HTML version of CHASE reads and writes the same format.
FORMAT = 1


def build(version: str, screen_file: str, screen_rows: list[list], struct_file: str | None, struct_rows: list[list] | None,
          settings: dict, analyzed: bool, selection: dict, campaign: str, decisions: list[dict]) -> str:
    """Serialise a CHASE session to JSON (format 1, shared with the HTML version).

    Args:
        version: CHASE version that wrote the file.
        screen_file: Name of the screening file.
        screen_rows: ``[raw, plate, well, batch_id]`` for every well.
        struct_file: Name of the structures file, or None.
        struct_rows: ``[batch_id, name, smiles]`` rows, or None.
        settings: All user settings (controls, normalization, classes, chemotypes, filters, colours ...).
        analyzed: Whether normalization had been run.
        selection: Keys of the pinned compound and of a member of the selected chemotype.
        campaign: Campaign name.
        decisions: Triage decisions.

    Returns:
        The session as a JSON string.
    """
    doc = {
        "app": "CHASE", "format": FORMAT, "version": version, "savedAt": datetime.now(timezone.utc).isoformat(),
        "screen": {"fileName": screen_file, "rows": screen_rows},
        "structures": {"fileName": struct_file, "rows": struct_rows} if struct_rows else None,
        "settings": settings, "analyzed": analyzed, "selection": selection,
        "triage": {"campaign": campaign, "decisions": decisions},
    }
    # numpy numbers are turned into normal Python numbers with .item(). Anything else unusual is saved as text.
    return json.dumps(doc, default=lambda o: o.item() if hasattr(o, "item") else str(o))


class SessionError(ValueError):
    """Raised when a file is not a readable CHASE session."""
    pass


def parse(text: str) -> dict:
    """Parse and minimally validate a session file; raise :class:`SessionError` if it is not one."""
    try:
        d = json.loads(text)
    except json.JSONDecodeError as e:
        raise SessionError("That file is not a CHASE session (it could not be read as JSON).") from e
    # A quick check that this is a CHASE file and that it contains the screening data. The settings
    # are checked one by one when they are applied.
    if not isinstance(d, dict) or d.get("app") != "CHASE" or not isinstance(d.get("screen", {}).get("rows"), list):
        raise SessionError("That file is not a CHASE session.")
    return d
