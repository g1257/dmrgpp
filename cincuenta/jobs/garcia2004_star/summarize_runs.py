#!/usr/bin/env python3
"""Summarize completed Garcia-2004 star-DMFT runs using only the Python stdlib."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


def read_table(path: Path) -> list[tuple[float, float, float]]:
    rows = []
    if not path.exists():
        return rows
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: expected three columns")
        row = tuple(map(float, fields))
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"{path}:{line_number}: non-finite value")
        rows.append(row)
    return rows


def input_value(text: str, key: str, default: str = "-") -> str:
    match = re.search(
        rf"^\s*(?:\w+\s+)?{re.escape(key)}\s*=\s*([^;]+);", text, re.MULTILINE
    )
    return match.group(1).strip().strip('"') if match else default


def manifest_value(text: str, key: str, default: str = "-") -> str:
    match = re.search(rf"^{re.escape(key)}=(.*)$", text, re.MULTILINE)
    return match.group(1).strip() if match else default


def summarize(run_dir: Path) -> tuple[list[str], bool]:
    input_text = (run_dir / "input.ain").read_text() if (run_dir / "input.ain").exists() else ""
    manifest = (run_dir / "manifest.txt").read_text() if (run_dir / "manifest.txt").exists() else ""
    log = (run_dir / "cincuenta.log").read_text(errors="replace") if (run_dir / "cincuenta.log").exists() else ""
    mats = read_table(run_dir / "gimp_dmrg.txt")
    lattice = read_table(run_dir / "latticeG_dmrg.txt")
    real = read_table(run_dir / "gimp_dmrg_real.txt")

    expected_mats = 2 * int(input_value(input_text, "Matsubaras"))
    expected_real = int(input_value(input_text, "OmegaTotal"))
    exit_status = manifest_value(manifest, "exit_status")
    complete = (
        len(mats) == expected_mats
        and len(lattice) == expected_mats
        and len(real) == expected_real
        and exit_status == "0"
    )

    positive = [row for row in mats if row[0] > 0]
    low_mats = min(positive, key=lambda row: row[0]) if positive else None
    zero_real = min(real, key=lambda row: abs(row[0])) if real else None
    iterations = re.findall(r"SelfConsistLoop iter=\s*([0-9]+)", log)
    errors = re.findall(r"SelfConsistLoop error=([-+0-9.eE]+)", log)
    converged = "yes" if "Converged after" in log else ("no" if log else "-")

    fields = [
        str(run_dir),
        input_value(input_text, "NumberOfBathPoints"),
        input_value(input_text, "HubbardU"),
        manifest_value(manifest, "ranks"),
        str(len(iterations)),
        converged,
        errors[-1] if errors else "-",
        f"{low_mats[0]:.12g}" if low_mats else "-",
        f"{low_mats[2]:.12g}" if low_mats else "-",
        f"{zero_real[0]:.12g}" if zero_real else "-",
        f"{-zero_real[2] / math.pi:.12g}" if zero_real else "-",
        manifest_value(manifest, "elapsed_seconds"),
        exit_status,
        "yes" if complete else "no",
    ]
    return fields, complete


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path, help="run directories or parents containing runs")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="return success even if a selected run lacks complete output tables",
    )
    args = parser.parse_args()

    run_dirs: list[Path] = []
    for path in args.paths:
        if (path / "input.ain").exists():
            run_dirs.append(path)
        elif path.is_dir():
            run_dirs.extend(sorted(child for child in path.iterdir() if (child / "input.ain").exists()))
        else:
            parser.error(f"not a run directory or parent: {path}")

    headings = [
        "run", "Nb", "U", "ranks", "iterations", "converged", "final_error",
        "lowest_wn", "ImG_lowest_wn", "nearest_real_w", "A_nearest_zero",
        "elapsed_s", "status", "complete",
    ]
    print("\t".join(headings))
    all_complete = True
    for run_dir in run_dirs:
        fields, complete = summarize(run_dir)
        print("\t".join(fields))
        all_complete = all_complete and complete
    return 0 if all_complete or args.allow_incomplete else 1


if __name__ == "__main__":
    raise SystemExit(main())
