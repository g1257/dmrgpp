"""Repository-relative paths shared by GBEK reference diagnostics."""

import os
from pathlib import Path


def default_build_dir() -> Path:
    """Return DMRGPP_BUILD_DIR, or this checkout's build/ directory."""
    return Path(os.environ.get("DMRGPP_BUILD_DIR", Path(__file__).resolve().parents[3] / "build"))
