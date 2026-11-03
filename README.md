# 2DOR reproducibility materials, version 1.0.1

This release contains frozen programme definitions, analysis-code provenance, group-level results, derived sample-level scores and quality-control data, and the final figure versions for the 2DOR cross-cohort transcriptomic analysis. The authors approved public release of the included sample identifiers and derived data. Raw sequencing data remain available through GEO.

## Contents

- `group_level_results/`: frozen node definitions and group-level estimates.
- `sample_level_results/Additional_file_2.csv`: full-precision long-format supporting tables, including sample identifiers, scores and quality-control measures. Old Fig3 boxplot display jitter was excluded because the final Fig3 uses a raincloud layout.
- `figure_reproduction/`: plot scripts and the source tables used to draw the figures. Run `plot_JOR_figures.py`, then `plot_fig3_raincloud.py` to replace the earlier Fig3 layout; run `plot_JOR_supplementary.py` and `plot_JOR_audit_figure.py` for the supplementary figures. Plots use Python, NumPy, pandas, SciPy, seaborn, matplotlib and Pillow. Horizontal jitter is deterministic and used solely for display.
- `audit_inputs/`: the row-disposition and sample-retention audit tables used by the audit-figure script.
- `figures_vector/`: PDF and SVG final versions of Figures 1–4, Supplementary Figure S1, and Audit Figure A1. Fig3, S1, and A1 display sample-level derived values.
- `code_provenance/`: original analysis scripts preserved for provenance. The scripts expect the original project layout and source GEO inputs; they are not a standalone analysis library.

## Data and citation

The underlying expression datasets are available in GEO under GSE286918, GSE193136, GSE232306, GSE274832 and GSE191322. Cite this version's Zenodo DOI once registered. The original v1.0.0 archive remains available at https://doi.org/10.5281/zenodo.23885987.

The analysis code is released under the MIT License. Tables and figures are released under CC BY 4.0; see `DATA_LICENSE.md`. The manuscript and credentials are not included.
