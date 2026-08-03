"""
Provenance tracking for generated plots/data artifacts in this directory.

This code changes frequently while debugging (multiple real bugs found and
fixed in a single session -- see README.md), and work is often not
committed as it goes. A bare git commit hash is not enough to reconstruct
which exact code produced a given historical plot if the tree was dirty at
the time -- which is the common case here. write_provenance() records the
commit hash, the dirty-file list, and a content hash of each relevant
source file, so a later bisection (by hand or by re-checking out old file
content) is actually possible.

Usage: call write_provenance(output_path, extra_files=[...]) right after
writing any plot or data artifact this package produces. Writes
"<output_path>.provenance.txt" alongside it.
"""
import hashlib
import subprocess
from pathlib import Path

# Source files that affect the physics/algorithm of anything this package
# generates. Extend this list as new files are added.
DEFAULT_TRACKED_FILES = [
    "gbek_ed.py",
    "gbek_dynamics.py",
    "gbek_cholesky.py",
    "gbek_selfconsistency.py",
    "../../src/NeqBathDecomposition.h",
    "../../src/ImpuritySolverNeqGBEK.h",
]


def _run(cmd, cwd):
    try:
        return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                               check=False).stdout.strip()
    except Exception as e:
        return f"<error: {e}>"


def _sha256(path):
    p = Path(path)
    if not p.exists():
        return "<missing>"
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def write_provenance(output_path, extra_files=None, notes=""):
    here = Path(__file__).resolve().parent
    repo_root = here
    for _ in range(6):
        if (repo_root / ".git").exists():
            break
        repo_root = repo_root.parent

    commit = _run(["git", "rev-parse", "HEAD"], cwd=repo_root)
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo_root)
    dirty_files = _run(["git", "status", "--porcelain", "--untracked-files=no"],
                        cwd=repo_root)

    files = list(DEFAULT_TRACKED_FILES) + list(extra_files or [])
    lines = [
        f"generated: {output_path}",
        f"git commit: {commit}",
        f"git branch: {branch}",
        f"working tree dirty: {'yes' if dirty_files else 'no'}",
    ]
    if dirty_files:
        lines.append("dirty files (git status --porcelain):")
        lines.append(dirty_files)
    lines.append("")
    lines.append("content hashes (sha256, first 16 hex chars) of files relevant "
                  "to this computation:")
    for f in files:
        full = (here / f).resolve()
        lines.append(f"  {_sha256(full)}  {f}")
    if notes:
        lines.append("")
        lines.append(f"notes: {notes}")

    prov_path = f"{output_path}.provenance.txt"
    with open(prov_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return prov_path
