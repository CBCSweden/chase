"""Per-session cache. Nothing is shared between users and everything is removed when the session ends.

Outside a Shiny session (tests, scripts) the helpers return None and callers simply don't cache.
Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

# The layout is {session id: {cache name: {key: value}}}. Everything stays in memory, and a session's
# entry is thrown away when that browser tab is closed.
_STORE: dict[str, dict[str, dict]] = {}


def session_cache(name: str) -> dict | None:
    """Return the cache dict called ``name`` for the current Shiny session, or None outside a session.

    Args:
        name: Cache name, e.g. ``"svg"`` or ``"flags"``.

    Returns:
        A dict private to the current session, or None when there is no session (then callers do not cache).
    """
    try:
        from shiny.session import get_current_session
        s = get_current_session()
    except Exception:  # noqa: BLE001 - no Shiny, no session
        return None
    sid = getattr(s, "id", None) if s is not None else None
    if not sid:
        return None
    return _STORE.setdefault(sid, {}).setdefault(name, {})


def drop(session_id: str) -> None:
    """Forget every cache belonging to ``session_id`` (called when the session ends)."""
    _STORE.pop(session_id, None)


def n_sessions() -> int:
    """Return the number of sessions that currently hold caches (used in tests)."""
    return len(_STORE)
