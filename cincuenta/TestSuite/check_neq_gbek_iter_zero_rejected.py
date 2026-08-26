#!/usr/bin/env python3
"""Verify production GBEK rejects NeqDmftIter=0."""

import subprocess
import sys


if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} CINCUENTA INPUT_DIR")

result = subprocess.run(
    [sys.argv[1], "-I", sys.argv[2], "-f", "inputNeqGBEKIterZeroRejected.ain"],
    capture_output=True,
    text=True,
)
output = result.stdout + result.stderr
expected = "Non-equilibrium GBEK requires NeqDmftIter>0"
if result.returncode == 0:
    raise SystemExit("zero-iteration GBEK input unexpectedly succeeded")
if expected not in output:
    raise SystemExit(f"zero-iteration rejection did not report {expected!r}:\n{output}")
