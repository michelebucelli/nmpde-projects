# LaTeX Report

This folder contains the LaTeX source for the final report.

Build command:

```bash
cd submission/final-report/latex
./build_report.sh
```

Main file:

```text
report.tex
```

Structure:

```text
Frontespizio.tex
sections/
assets/
```

Expected output:

```text
WaveFEM_Solver_Report.pdf
```

The Polimi logo used in the front page is stored in:

```text
assets/Logo_Politecnico_Milano.png
```

The report explicitly includes the four existing screenshots in
`submission/assets/paraview-screenshots/`; there is no automatic image discovery.
Generate its numerical tables and plots with `python scripts/review_validation.py`
from the repository root in the deal.II environment before compiling. Generate
the empirical phase table and figure with `python scripts/dispersion_check.py`.
