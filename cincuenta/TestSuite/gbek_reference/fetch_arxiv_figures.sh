#!/usr/bin/env bash
# Fetch the GBEK paper's own figures directly from its arXiv e-print source
# (arXiv:1306.6315), for embedding in report.tex.
#
# We deliberately do NOT check these figures into the repo (they're the
# authors'/publisher's copyrighted material, and they're trivially
# re-fetchable) -- this script is the reproducible way to get them back,
# same spirit as regenerate_plots.sh for our own plots. Output goes to
# arxiv_figures/, which is gitignored.
#
# Requires: curl, epstopdf (TeX Live), pdftoppm (poppler), gs (ghostscript;
# epstopdf shells out to it).
#
# Usage:
#   ./fetch_arxiv_figures.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

ARXIV_ID="1306.6315"
OUT_DIR="arxiv_figures"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

mkdir -p "$OUT_DIR"

echo "Fetching arXiv:$ARXIV_ID e-print source..."
curl -sL -o "$WORK_DIR/src.tar.gz" "https://arxiv.org/e-print/$ARXIV_ID"
tar xzf "$WORK_DIR/src.tar.gz" -C "$WORK_DIR"

# Mapping confirmed against \label{fig:...}/\caption in the paper's own
# siambbl.tex -- figure numbering in the compiled PRB paper matches the
# \includegraphics order in the source exactly (contour, dmft_scheme,
# lambdaappr=Fig.3, hybrid=Fig.4, initialstate, selfcon,
# causalityA=Fig.7, causalityB=Fig.8, energy_cnsrv=Fig.9, double_occ=Fig.10).
# (Plain "N:source" pairs, not an associative array -- macOS ships bash 3.2,
# which doesn't have those.)
FIG_SOURCE_LIST="3:fig_cholesky_eigenvalue.eps 4:fig_hopping.eps 7:fig_causality_A.eps 8:fig_causality_B.eps 9:fig_energy_cnsrv.eps 10:fig_double_occ.eps"

for pair in $FIG_SOURCE_LIST; do
	n="${pair%%:*}"
	src="${pair#*:}"
	if [ ! -f "$WORK_DIR/$src" ]; then
		echo "ERROR: expected $src in arXiv source, not found" >&2
		exit 1
	fi
	( cd "$WORK_DIR" && epstopdf "$src" )
	pdf="$WORK_DIR/${src%.eps}.pdf"
	pdftoppm -png -r 200 -singlefile "$pdf" "$OUT_DIR/paper_fig${n}"
	echo "Wrote $OUT_DIR/paper_fig${n}.png (from $src)"
done

cat > "$OUT_DIR/PROVENANCE.txt" <<EOF
Source: https://arxiv.org/e-print/$ARXIV_ID (Gramsch, Balzer, Eckstein &
Kollar, Phys. Rev. B 88, 235106 (2013); arXiv:1306.6315)
Fetched: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
Figure mapping verified against siambbl.tex \\label{}/\\caption{} order.
Regenerate with: ./fetch_arxiv_figures.sh
EOF

echo "Done. Figures in $OUT_DIR/."
