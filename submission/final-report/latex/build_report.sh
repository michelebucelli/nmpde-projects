#!/usr/bin/env bash
set -euo pipefail

REPORT=report.tex
OUTPUT=WaveFEM_Solver_Report.pdf

if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$REPORT"
elif command -v pdflatex >/dev/null 2>&1; then
  pdflatex -interaction=nonstopmode -halt-on-error "$REPORT"
  pdflatex -interaction=nonstopmode -halt-on-error "$REPORT"
else
  echo "No LaTeX compiler found. Install latexmk or pdflatex."
  exit 1
fi

if [ -f report.pdf ]; then
  mv report.pdf "$OUTPUT"
fi

echo "Wrote $OUTPUT"
