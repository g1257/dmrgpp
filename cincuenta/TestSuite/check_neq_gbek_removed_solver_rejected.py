#!/usr/bin/env python3
"""Verify production input parsing rejects removed NeqSolver=."""

import subprocess
import sys


if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} CINCUENTA INPUT_DIR")

result = subprocess.run(
    [sys.argv[1], "-I", sys.argv[2], "-f", "inputNeqGBEKRemovedSolverRejected.ain"],
    capture_output=True,
    text=True,
)
output = result.stdout + result.stderr
expected = "NeqSolver= was removed"
if result.returncode == 0:
    raise SystemExit("removed NeqSolver input unexpectedly succeeded")
if expected not in output:
    raise SystemExit(f"removed NeqSolver rejection did not report {expected!r}:\n{output}")
