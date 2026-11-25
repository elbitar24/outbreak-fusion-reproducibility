$ErrorActionPreference = "Stop"
python analysis/run_uncertainty.py
python tests/test_reproducibility.py
