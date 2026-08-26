#!/usr/bin/env python3
"""Check the fast non-atomic, positive-rank GBEK end-to-end regression."""

import argparse
import math
import pathlib
import sys


EPS = 1.0e-12


def fail(message):
    raise RuntimeError(message)


def read_four_column_file(path):
    """Return {(first, second): complex(third, fourth)} after strict validation."""
    data = {}
    try:
        lines = path.read_text().splitlines()
    except OSError as exc:
        fail(f"cannot read {path}: {exc}")

    if not lines:
        fail(f"{path} is empty")
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if len(fields) != 4:
            fail(f"{path}:{line_number}: expected four columns, got {len(fields)}")
        try:
            first, second, real, imag = (float(field) for field in fields)
        except ValueError as exc:
            fail(f"{path}:{line_number}: non-numeric row: {exc}")
        if not all(math.isfinite(value) for value in (first, second, real, imag)):
            fail(f"{path}:{line_number}: non-finite value")
        key = (round(first, 8), round(second, 8))
        if key in data:
            fail(f"{path}:{line_number}: duplicate row for {key}")
        data[key] = complex(real, imag)
    return data


def require_full_time_grid(label, data, nt, dt, triangular):
    expected = set()
    for n in range(nt + 1):
        js = range(n + 1) if triangular else range(nt + 1)
        for j in js:
            expected.add((round(n * dt, 8), round(j * dt, 8)))
    actual = set(data)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        fail(f"{label}: incomplete grid; missing={missing}, unexpected={unexpected}")


def read_matsubara_file(path, expected_rows):
    try:
        lines = path.read_text().splitlines()
    except OSError as exc:
        fail(f"cannot read {path}: {exc}")

    if len(lines) != expected_rows:
        fail(f"{path}: expected {expected_rows} rows, got {len(lines)}")
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if len(fields) != 3:
            fail(f"{path}:{line_number}: expected three columns, got {len(fields)}")
        try:
            values = [float(field) for field in fields]
        except ValueError as exc:
            fail(f"{path}:{line_number}: non-numeric row: {exc}")
        if not all(math.isfinite(value) for value in values):
            fail(f"{path}:{line_number}: non-finite value")


def read_cholesky(path, nt, rank):
    data = read_four_column_file(path)
    expected = {(float(n), float(p)) for n in range(nt + 1) for p in range(rank)}
    if set(data) != expected:
        fail(f"{path}: expected Cholesky rows {sorted(expected)}, got {sorted(data)}")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=pathlib.Path)
    parser.add_argument("--prefix", default="neq-gbek-fast")
    parser.add_argument("--nt", type=int, default=3)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--rank", type=int, default=2)
    parser.add_argument("--matsubaras", type=int, default=16)
    args = parser.parse_args()

    directory = args.run_directory
    prefix = directory / args.prefix
    g_retarded = read_four_column_file(pathlib.Path(f"{prefix}-green-retarded"))
    g_lesser = read_four_column_file(pathlib.Path(f"{prefix}-green-lesser"))
    lambda_lesser = read_four_column_file(pathlib.Path(f"{prefix}-lambda-lesser"))
    plus_lesser = read_four_column_file(pathlib.Path(f"{prefix}-plus-bath-lesser"))
    read_matsubara_file(
        pathlib.Path(f"{prefix}-equilibrium-gimp-matsubara"), 2 * args.matsubaras)
    cholesky = read_cholesky(pathlib.Path(f"{prefix}-cholesky-V"), args.nt, args.rank)

    # Complete, finite output grids demonstrate the real-time run reached its last step.
    require_full_time_grid("Gimp retarded", g_retarded, args.nt, args.dt, triangular=True)
    require_full_time_grid("Gimp lesser", g_lesser, args.nt, args.dt, triangular=False)
    require_full_time_grid("total lesser hybridization", lambda_lesser, args.nt, args.dt, triangular=False)
    require_full_time_grid("reconstructed plus-bath lesser", plus_lesser, args.nt, args.dt, triangular=False)

    gr00 = g_retarded[(0.0, 0.0)]
    if abs(gr00.real) > 1.0e-8 or abs(gr00.imag + 1.0) > 1.0e-8:
        fail(f"G^R(0,0)={gr00} is not -i within 1e-8")

    v10 = cholesky[(1.0, 0.0)]
    v21 = cholesky[(2.0, 1.0)]
    if abs(v10) <= EPS:
        fail(f"Cholesky V(1,0) is zero: {v10}")
    if abs(v21) <= EPS:
        fail(f"Cholesky V(2,1) is zero: {v21}")

    # dumpPlusBath uses Lambda^+ directly; the real-time lesser dump is Delta^<,
    # so multiply it by -i before comparing the two hybridizations.
    total00 = -1j * lambda_lesser[(0.0, 0.0)]
    plus00 = plus_lesser[(0.0, 0.0)]
    if abs(total00) <= EPS:
        fail(f"total Lambda^<(0,0) is zero: {total00}")
    if abs(plus00) > EPS:
        fail(f"reconstructed Lambda^+<(0,0) is not zero: {plus00}")

    print("PASS: fast non-atomic rank-2 GBEK output is complete and finite")
    print(f"G^R(0,0)={gr00}; V(1,0)={v10}; V(2,1)={v21}")
    print(f"Lambda^<(0,0)={total00}; Lambda^+<(0,0)={plus00}")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
