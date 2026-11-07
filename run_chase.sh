#!/usr/bin/env bash
# Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
# Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
# for valuable scientific input and feedback during the development of CHASE.
# Start CHASE from its conda environment. Extra options are passed on, e.g. ./run_chase.sh --port 9000
set -e
cd "$(dirname "$0")"
if ! command -v conda >/dev/null 2>&1; then
  for d in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
    [ -f "$d/etc/profile.d/conda.sh" ] && . "$d/etc/profile.d/conda.sh" && break
  done
fi
if ! command -v conda >/dev/null 2>&1; then
  echo "conda was not found. Install Miniforge (https://conda-forge.org/download/), then run: conda env create -f environment.yml"; exit 1
fi
eval "$(conda shell.bash hook)"
if ! conda env list | grep -qE '^chase\s'; then
  echo "Creating the 'chase' environment (first run only)..."
  conda env create -f environment.yml
fi
conda activate chase
exec python -m chase "$@"
