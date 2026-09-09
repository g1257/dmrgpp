#!/usr/bin/env python3
"""Compare rank-count-independent outputs from two cincuenta DMFT runs."""

import argparse
import math
import re
from pathlib import Path


DEFAULT_FILES = ("gimp_dmrg.txt", "latticeG_dmrg.txt", "gimp_dmrg_real.txt")


def read_table(path: Path):
    rows = []
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            try:
                rows.append(tuple(float(field) for field in fields))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-numeric row") from error
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    return rows


def compare_file(reference: Path, candidate: Path, tolerance: float):
    expected = read_table(reference)
    actual = read_table(candidate)
    if len(expected) != len(actual):
        raise ValueError(
            f"{candidate}: row count {len(actual)} != reference {len(expected)}"
        )

    maximum = 0.0
    location = None
    for row_number, (expected_row, actual_row) in enumerate(
        zip(expected, actual), 1
    ):
        if len(expected_row) != len(actual_row):
            raise ValueError(
                f"{candidate}:{row_number}: column count {len(actual_row)} "
                f"!= reference {len(expected_row)}"
            )
        for column, (expected_value, actual_value) in enumerate(
            zip(expected_row, actual_row), 1
        ):
            difference = abs(expected_value - actual_value)
            if not math.isfinite(difference):
                raise ValueError(
                    f"{candidate}:{row_number}:{column}: non-finite difference"
                )
            if difference > maximum:
                maximum = difference
                location = (row_number, column, expected_value, actual_value)

    if maximum > tolerance:
        row, column, expected_value, actual_value = location
        raise ValueError(
            f"{candidate}:{row}:{column}: |{actual_value} - {expected_value}| "
            f"= {maximum:.6g} exceeds tolerance {tolerance:.6g}"
        )
    return maximum


def verify_iterations(directory: Path, log_name: str, iterations: int):
    log_path = directory / log_name
    matches = [
        int(value)
        for value in re.findall(r"SelfConsistLoop iter= ([0-9]+)", log_path.read_text())
    ]
    expected = list(range(iterations))
    if matches != expected:
        raise ValueError(f"{log_path}: iterations {matches} != expected {expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("candidate_dir", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    parser.add_argument("--files", nargs="+", default=DEFAULT_FILES)
    parser.add_argument("--log")
    parser.add_argument("--iterations", type=int)
    args = parser.parse_args()

    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        parser.error("--tolerance must be finite and nonnegative")
    if (args.log is None) != (args.iterations is None):
        parser.error("--log and --iterations must be given together")
    if args.iterations is not None:
        if args.iterations < 1:
            parser.error("--iterations must be positive")
        verify_iterations(args.reference_dir, args.log, args.iterations)
        verify_iterations(args.candidate_dir, args.log, args.iterations)

    for filename in args.files:
        maximum = compare_file(
            args.reference_dir / filename,
            args.candidate_dir / filename,
            args.tolerance,
        )
        print(f"{filename}: max absolute difference = {maximum:.6g}")


if __name__ == "__main__":
    main()
