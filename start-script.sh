#!/bin/bash
# Start CHASE on SciLifeLab Serve
exec /opt/conda/bin/python -m chase --host 0.0.0.0 --port 8765 --no-browser --tmp /tmp
