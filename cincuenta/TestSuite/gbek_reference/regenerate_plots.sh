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
#   (B) Depend on actual cincuenta C++ runs having produced dump files in
#       build/ first (the GBEK Fig. 3 causal-Cholesky-vs-eigenvector
#       validation plots from the earlier seeding-timing investigation).
#       This script builds cincuenta and runs the two prerequisite .ain
#       inputs if their dumps aren't already present in build/.
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
BUILD_DIR="$REPO_ROOT/build"
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

	AL_PREFIX="$BUILD_DIR/atomic-limit-gbek-L3"       # inputNeqAtomicLimitGBEKL3.ain
	FIG3_PREFIX="$BUILD_DIR/gebk-fig3-L3"             # inputNeqGBEKFig3L3.ain

	need_cincuenta_run() {
		local prefix="$1"
		[ ! -f "${prefix}-weiss-delta-lesser" ] || [ ! -f "${prefix}-plus-bath-lesser" ]
	}

	if need_cincuenta_run "$AL_PREFIX" || need_cincuenta_run "$FIG3_PREFIX"; then
		echo "Missing prerequisite cincuenta dumps -- building cincuenta..."
		cmake --build "$BUILD_DIR" --target cincuenta -j4
	fi

	if need_cincuenta_run "$AL_PREFIX"; then
		echo "Running inputNeqAtomicLimitGBEKL3.ain to produce ${AL_PREFIX}-*..."
		( cd "$BUILD_DIR" && ./cincuenta/src/cincuenta \
			-f "$REPO_ROOT/cincuenta/TestSuite/inputs/inputNeqAtomicLimitGBEKL3.ain" )
	fi
	if need_cincuenta_run "$FIG3_PREFIX"; then
		echo "Running inputNeqGBEKFig3L3.ain to produce ${FIG3_PREFIX}-*..."
		( cd "$BUILD_DIR" && ./cincuenta/src/cincuenta \
			-f "$REPO_ROOT/cincuenta/TestSuite/inputs/inputNeqGBEKFig3L3.ain" )
	fi

	if [ ! -f gbek-atomic-limit-exact-lesser ]; then
		echo "Generating the pure-Python exact atomic-limit reference..."
		$UV_NUMPY gbek_selfconsistency.py \
			--L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser
	fi

	# Each of these has its input paths hardcoded as module-level constants
	# pointing at $BUILD_DIR/atomic-limit-gbek-L3-* or
	# $BUILD_DIR/gebk-fig3-L3-* -- see each script's own header if you need
	# to point it elsewhere.
	$UV_NUMPY plot_atomic_limit_2d.py
	$UV_NUMPY plot_collapse_evidence_summary.py
	$UV_NUMPY plot_errstep_t3scan.py
	$UV_NUMPY plot_fig3l3_post_fix.py
	$UV_NUMPY scan_t3_activation.py
	# These two calls use compare_reference.py's DEFAULT labels ("Exact
	# reference" / "cincuenta rank-L Cholesky approx"), which is correct
	# ONLY because arg1 here is genuinely the independent, undecomposed
	# Python target and arg2 is genuinely cincuenta's reconstruction of it.
	# gbek_reference_comparison.png is embedded in the GBEK progress report
	# artifact -- if you ever repoint either argument at something else
	# (e.g. comparing two cincuenta runs against each other), you MUST add
	# --ref-label/--approx-label explicitly, or the report will silently
	# ship a mislabeled plot (this happened once already -- see
	# small_bath_vs_atomic_limit_L2.png's generation for the fix).
	$UV_NUMPY compare_reference.py \
		gbek-atomic-limit-exact-lesser "${AL_PREFIX}-plus-bath-lesser" \
		--tmax 4.0 --out atomic_limit_true_comparison.png
	$UV_NUMPY compare_reference.py \
		gbek-atomic-limit-exact-lesser "${AL_PREFIX}-plus-bath-lesser" \
		--tmax 4.0 --out gbek_reference_comparison.png
	$UV_NUMPY quantify_lambda_minus_leak.py "$FIG3_PREFIX"

	echo "Group B done: atomic_limit_2d_rank_comparison.png,"
	echo "  atomic_limit_true_comparison.png, collapse_evidence_summary.png,"
	echo "  lambda_minus_leak_check.png, errstep_t3scan_comparison.png,"
	echo "  fig3L3_near_atomic_post_fix.png, gbek_reference_comparison.png,"
	echo "  t3_activation_scan.png"
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
