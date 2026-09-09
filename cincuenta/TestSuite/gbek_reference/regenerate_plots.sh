#!/usr/bin/env bash
# Regenerate every plot this directory's tooling can produce.
#
# None of these .png files are committed to the repo (see the "Regenerating
# plots" section of README.md for why) -- this script is the documented,
# reproducible way to get them back locally. Run from anywhere; it cds into
# this directory itself.
#
# Two families of plots:
#   (A) Pure Python, no C++ build needed: Fig. 7/8/9/10 double-occupation and
#       energy-conservation reproduction (fig7_docc.png, fig8_docc.png,
#       fig9_energy.png, fig10_docc.png) and their prerequisite .npz data.
#   (B) Compare the independent Python atomic-limit reference with an actual
#       atomic-limit cincuenta run. This script configures/builds cincuenta
#       and runs the prerequisite .ain input when its dumps are absent or
#       stale. Fitted-first-bath stress tests are deliberately not part of
#       plot or report generation.
#
# Usage:
#   ./regenerate_plots.sh            # regenerate everything (groups A+B)
#   ./regenerate_plots.sh --group-a  # just the Fig. 7/8 tooling (fast)
#   ./regenerate_plots.sh --group-b  # just the C++-dependent plots (slower,
#                                     # builds cincuenta and runs it if needed)
#   ./regenerate_plots.sh --report   # build report.pdf (needs A+B's plots
#                                     # already present, plus network access
#                                     # for fetch_arxiv_figures.sh and a TeX
#                                     # install); NOT run by default, and not
#                                     # part of the cmake build.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

REPO_ROOT="$(cd ../../.. && pwd)"
# Keep Python diagnostics, generated C++ dumps, and the runner on the same
# configurable build tree. An explicit environment setting wins.
BUILD_DIR="${DMRGPP_BUILD_DIR:-$REPO_ROOT/build}"
export DMRGPP_BUILD_DIR="$BUILD_DIR"
UV_NUMPY="uv run --with numpy --with scipy --with matplotlib python3"

GROUP_A=1
GROUP_B=1
GROUP_REPORT=0
if [ "${1:-}" = "--group-a" ]; then GROUP_B=0; fi
if [ "${1:-}" = "--group-b" ]; then GROUP_A=0; fi
if [ "${1:-}" = "--report" ]; then GROUP_A=0; GROUP_B=0; GROUP_REPORT=1; fi

# ---------------------------------------------------------------------------
# Group A: Fig. 7/8 double-occupation reproduction (pure Python)
# ---------------------------------------------------------------------------
if [ "$GROUP_A" = "1" ]; then
	echo "=== Group A: Fig. 7/8 reproduction ==="

	if [ ! -f fig7_docc_L2_cholesky.npz ]; then
		$UV_NUMPY run_fig7_scan.py --L 2
	fi
	if [ ! -f fig7_docc_L4_cholesky.npz ]; then
		$UV_NUMPY run_fig7_scan.py --L 4
	fi
	$UV_NUMPY plot_docc_scan.py --figure 7 --L 2,4

	if [ ! -f fig8_docc_A_cholesky.npz ]; then
		$UV_NUMPY run_fig8_scan.py
	fi
	if [ ! -f fig8_docc_E_cholesky.npz ]; then
		$UV_NUMPY investigate_L4_tmax4.py
	fi
	$UV_NUMPY plot_docc_scan.py --figure 8

	if [ ! -f fig9_energy_U2_L2.npz ]; then
		$UV_NUMPY run_fig9_scan.py
	fi
	$UV_NUMPY plot_energy_scan.py

	if [ ! -f fig10_docc_U0_L2.npz ]; then
		$UV_NUMPY run_fig10_scan.py
	fi
	$UV_NUMPY plot_docc_scan.py --figure 10

	echo "Group A done: fig7_docc.png, fig8_docc.png, fig9_energy.png, fig10_docc.png"
fi

# ---------------------------------------------------------------------------
# Group B: plots depending on actual cincuenta C++ dumps
# ---------------------------------------------------------------------------
if [ "$GROUP_B" = "1" ]; then
	echo "=== Group B: C++-dependent validation plots ==="

	AL_INPUT="$REPO_ROOT/cincuenta/TestSuite/inputs/inputNeqAtomicLimitGBEKL3.ain"
	AL_PREFIX="$BUILD_DIR/atomic-limit-gbek-L3"

	if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
		echo "Configuring cincuenta build tree in $BUILD_DIR..."
		cmake -S "$REPO_ROOT" -B "$BUILD_DIR" -DBUILD_TESTING=OFF
	fi

	echo "Building cincuenta in $BUILD_DIR..."
	cmake --build "$BUILD_DIR" --target cincuenta -j4
	CINCUENTA="$BUILD_DIR/cincuenta/src/cincuenta"

	need_atomic_run() {
		local suffix
		for suffix in lambda-lesser plus-bath-lesser cholesky-V docc-energy; do
			[ -s "${AL_PREFIX}-${suffix}" ] || return 0
			[ "${AL_PREFIX}-${suffix}" -nt "$CINCUENTA" ] || return 0
			[ "${AL_PREFIX}-${suffix}" -nt "$AL_INPUT" ] || return 0
		done
		return 1
	}

	if need_atomic_run; then
		echo "Running inputNeqAtomicLimitGBEKL3.ain to produce ${AL_PREFIX}-*..."
		rm -f "${AL_PREFIX}-"{green-retarded,green-lesser,lambda-retarded,lambda-lesser,plus-bath-lesser,cholesky-V,docc-energy,equilibrium-gimp-matsubara}
		( cd "$BUILD_DIR" && "$CINCUENTA" -f "$AL_INPUT" )
	fi

	if [ ! -f gbek-atomic-limit-exact-lesser ]; then
		echo "Generating the pure-Python exact atomic-limit reference..."
		$UV_NUMPY gbek_selfconsistency.py \
			--L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser
	fi

	# Report plots are all based on the atomic-limit target. Diagnostics use
	# DMRGPP_BUILD_DIR (exported above) to find the matching C++ dump.
	$UV_NUMPY plot_fig3_errstep.py --target gbek-atomic-limit-exact-lesser
	$UV_NUMPY plot_fig4_hybridization.py --target gbek-atomic-limit-exact-lesser
	$UV_NUMPY plot_atomic_limit_2d.py

	# The default labels are correct because arg1 is the independent,
	# undecomposed Python target and arg2 is cincuenta's atomic-limit
	# rank-L Cholesky reconstruction.
	$UV_NUMPY compare_reference.py \
		gbek-atomic-limit-exact-lesser "${AL_PREFIX}-plus-bath-lesser" \
		--tmax 4.0 --out gbek_reference_comparison.png

	echo "Group B done: fig3_errstep.png, fig4_hybridization.png,"
	echo "  atomic_limit_2d_rank_comparison.png, gbek_reference_comparison.png"
fi

# ---------------------------------------------------------------------------
# Group report: build report.pdf (LaTeX write-up, groups A+B's plots plus
# the paper's own figures fetched fresh from arXiv)
# ---------------------------------------------------------------------------
if [ "$GROUP_REPORT" = "1" ]; then
	echo "=== Group report: building report.pdf and report-summary.pdf ==="

	if [ ! -d arxiv_figures ]; then
		./fetch_arxiv_figures.sh
	fi

	latexmk -pdf -interaction=nonstopmode report.tex
	# report-summary.pdf: same content with the "bug found and fixed"
	# development-history callouts omitted (see report.tex's \PIVERSION
	# header comment) -- a status/capabilities-only copy for sharing with
	# people who aren't interested in development-process detail.
	latexmk -pdf -interaction=nonstopmode report-summary.tex

	echo "Group report done: report.pdf, report-summary.pdf"
fi
