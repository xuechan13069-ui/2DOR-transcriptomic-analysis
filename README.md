# 2DOR reproducibility materials

This repository contains the release-approved code provenance, frozen programme definitions,
group-level source tables, and vector figures for the 2DOR cross-cohort transcriptomic analysis.

## Contents

- `group_level_results/`: frozen node definitions and group-level estimates used in the manuscript figures and tables.
- `code_provenance/`: analysis scripts retained for provenance. They expect the original project layout and the public GEO inputs listed below; they are not packaged as a standalone software library.
- `figures_vector/`: PDF and SVG versions of Figures 1, 2, and 4, which contain group-level information.

## Data boundary

Derived sample-level identifiers, quality-control values, programme scores, and display coordinates are not included in this public release. The corresponding vector figures (Figure 3, Supplementary Figure S1, and Audit Figure A1) are also excluded because exact point coordinates may be recoverable from vector files. Group-level estimates and confidence intervals are included.

The underlying expression datasets are available from NCBI GEO under GSE286918, GSE193136, GSE232306, GSE274832, and GSE191322. Source-study conditions and eligibility definitions should be checked against the cited publications.

## Licence

The code is released under the MIT License. Group-level tables and figures are released under
Creative Commons Attribution 4.0 International (CC BY 4.0); see `DATA_LICENSE.md`.

## Citation

Citation metadata are provided in `CITATION.cff` and `.zenodo.json`. Cite the archived Zenodo
release once its DOI is registered.

No manuscript file, clinical source material, credential, or access token is included.
