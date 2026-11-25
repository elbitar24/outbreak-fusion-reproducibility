#!/usr/bin/env bash
set -euo pipefail
python analysis/run_uncertainty.py
python tests/test_reproducibility.py
