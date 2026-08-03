#!/usr/bin/env bash
# Run every module's built-in self-test (its `if __name__ == "__main__":`
# block) and fail loudly if any of them fails.
#
# These are smoke tests, not a full regression suite: each core library
# module validates itself against an independent analytic/exact result
# (see README.md's "Validation status" section for what each one checks),
# but nothing here runs them automatically -- this script is that missing
# "run everything and fail on the first problem" entry point. It is NOT
# part of the cmake build; run it by hand after changing any of these
# modules, or before relying on a fresh checkout.
#
# Usage:
#   ./run_self_tests.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

UV_NUMPY="uv run --with numpy --with scipy python3"

MODULES=(
    gbek_ed.py
    gbek_dynamics.py
    gbek_cholesky.py
    atomic_limit_reference.py
)

fail=0
for m in "${MODULES[@]}"; do
    echo "=== $m ==="
    if $UV_NUMPY "$m"; then
        echo "--- $m: OK ---"
    else
        echo "--- $m: FAILED ---" >&2
        fail=1
    fi
    echo
done

if [ "$fail" -ne 0 ]; then
    echo "run_self_tests.sh: one or more self-tests FAILED" >&2
    exit 1
fi
echo "run_self_tests.sh: all self-tests passed"
