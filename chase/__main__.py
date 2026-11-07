"""Start CHASE:  python -m chase [--host 127.0.0.1] [--port 8765] [--no-browser] [--tmp DIR]

CHASE keeps all user data in memory and stores nothing on the server:
uploaded files are read and deleted at once, and Shiny's temporary upload folder is placed in
RAM (/dev/shm) on Linux, so uploads never reach the disk.
Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import argparse
import atexit
import os
import shutil
import sys
import tempfile
from pathlib import Path


def _pid_alive(pid: int) -> bool:
    """Return True if a process with this PID exists (used to find upload folders left by crashed runs)."""
    try:
        # Sending signal 0 does nothing to the process. It only tells us whether a process with this PID still exists.
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # The process is there, it just belongs to another user, so it is still alive.
    except PermissionError:
        return True
    return True


def _memory_tmp(port: int, requested: str | None) -> tuple[Path | None, str]:
    """A private upload folder for this process only; folders of dead CHASE processes should be removed."""
    # Decide where the temporary upload folder goes. The --tmp folder is used if the user gave one,
    # otherwise RAM (/dev/shm) on Linux, or whatever the system uses for temporary files.
    if requested:
        base = Path(requested).expanduser()
        where = str(base)
    elif sys.platform.startswith("linux") and os.path.isdir("/dev/shm") and os.access("/dev/shm", os.W_OK):
        base = Path("/dev/shm")
        where = "RAM"
    else:
        return None, "temporary uploads in the system temp folder (deleted as soon as they are read)"
    # Each running copy of CHASE gets its own folder, named chase-<user id>-<port>-<process id>.
    # If CHASE crashed earlier, its old folder may still be lying around. We delete those leftovers here,
    # but only when the process that made them is gone, so we never touch a copy that is still running.
    prefix = f"chase-{os.getuid() if hasattr(os, 'getuid') else 0}-{port}-"
    for old in base.glob(prefix + "*"):
        try:
            if not _pid_alive(int(old.name.rsplit("-", 1)[1])):
                shutil.rmtree(old, ignore_errors=True)
        except ValueError:
            pass
    root = base / f"{prefix}{os.getpid()}"
    # Only the current user may read this folder. Python and Shiny are told to put their temporary files
    # here, and the folder is deleted again when CHASE stops.
    root.mkdir(parents=True, mode=0o700, exist_ok=True)
    tempfile.tempdir = str(root)
    os.environ["TMPDIR"] = str(root)
    atexit.register(shutil.rmtree, root, True)
    return root, f"temporary uploads in {where} ({root})"


def _port_free(host: str, port: int) -> bool:
    """Return True if ``host:port`` can be bound, i.e. no other program (or CHASE instance) is using it."""
    import socket
    fam = socket.AF_INET6 if ":" in host else socket.AF_INET
    with socket.socket(fam, socket.SOCK_STREAM) as s:
        # SO_REUSEADDR lets us ignore connections left over from a CHASE that just stopped. That way the port
        # only counts as busy if another program is really listening on it.
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind((host, port))
            return True
        except OSError:
            return False


def main() -> None:
    """Parse the command line, prepare the in-memory upload folder and start the Shiny app."""
    ap = argparse.ArgumentParser(prog="chase", description="CHASE - Chemotype-Hit Analysis of Screening Experiments")
    ap.add_argument("--host", default="127.0.0.1", help="interface to listen on (0.0.0.0 serves your network)")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    ap.add_argument("--tmp", help="folder for Shiny's short-lived upload copies (default on Linux: /dev/shm, i.e. RAM)")
    a = ap.parse_args()
    # Don't start a second copy of CHASE on a port that is already in use.
    if not _port_free(a.host, a.port):
        sys.exit(f"Port {a.port} is already in use (is CHASE already running?). Choose another with --port.")
    root, note = _memory_tmp(a.port, a.tmp)
    if root is not None:
        # Ctrl+C is handled by the 'finally' below. Here we also clean up if CHASE is stopped with SIGTERM
        # (for example by 'kill' or by a service manager).
        import signal

        def _stop(signum, frame):
            """On SIGTERM remove the upload folder before exiting."""
            shutil.rmtree(root, ignore_errors=True)
            raise SystemExit(0)
        signal.signal(signal.SIGTERM, _stop)
    # Shiny is only imported here, so that 'python -m chase --help' and the port check answer quickly.
    from shiny import run_app
    from chase import __version__
    shown = "localhost" if a.host in ("127.0.0.1", "0.0.0.0") else a.host
    print(f"CHASE {__version__} - http://{shown}:{a.port}  ({note}; nothing is written to disk). Ctrl+C to stop.")
    try:
        # access_log=False: CHASE does not keep a record of the requests it receives. Only warnings are printed.
        run_app("chase.app:app", host=a.host, port=a.port, launch_browser=not a.no_browser,
                access_log=False, log_level="warning")
    finally:
        if root is not None:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
