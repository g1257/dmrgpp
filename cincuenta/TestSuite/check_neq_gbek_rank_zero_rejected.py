#!/usr/bin/env python3
"""Verify production GBEK rejects NeqBathRank=0."""

import subprocess
import sys


if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} CINCUENTA INPUT_DIR")

result = subprocess.run(
    [sys.argv[1], "-I", sys.argv[2], "-f", "inputNeqGBEKRankZeroRejected.ain"],
    capture_output=True,
    text=True,
)
output = result.stdout + result.stderr
expected = "Non-equilibrium GBEK requires NeqBathRank>0"
if result.returncode == 0:
    raise SystemExit("rank-zero GBEK input unexpectedly succeeded")
if expected not in output:
    raise SystemExit(f"rank-zero rejection did not report {expected!r}:\n{output}")
