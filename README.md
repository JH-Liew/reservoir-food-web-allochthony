# Terrestrial resource use in tropical reservoir food webs

This repository contains the data and R code supporting the manuscript “Catchment urbanisation is associated with greater terrestrial resource use in tropical reservoir food webs”.

These analyses build on the posterior-median prey–consumer matrices reported and archived by Wilkinson et al. (2022). The original stable-isotope mixing models used to construct these matrices are not refitted here.

## Files

- `reservoir_allochthony_data.xlsx`: reservoir and taxon records, mixing-model inputs, fish-biomass summaries and model outputs.
- `code/Liew_et_al_analysis_archive.R`: R code for the allochthony calculations, statistical analyses, diagnostics and figure components.

## Precursor food-web matrices

Obtain the 12 matrix files from the Dryad release [Empirical food webs of 12 tropical reservoirs in Singapore](https://doi.org/10.5061/dryad.jsxksn088). Download `Bottom-up_predation_matrices_Res_1.csv` through `Bottom-up_predation_matrices_Res_12.csv`, place them in a local directory and set the `ALLOCHTHONY_MATRIX_DIR` environment variable to that directory.

## Running the analyses

1. Download the repository and open R in its top-level directory.
2. Set `ALLOCHTHONY_MATRIX_DIR` to the directory containing the Wilkinson et al. matrices.
3. Run `code/Liew_et_al_analysis_archive.R`.

By default, the script recalculates the matrix-based summaries and reads the computationally intensive model results from the workbook. Set the relevant `run_...` option at the beginning of the script to `TRUE` to refit an analysis. Recalculated results and a record of the R session are written to a newly created `outputs` directory.

The script requires `readxl`, `ggplot2`, `NetIndices`, `igraph`, `vegan`, `lme4`, `rjags`, `coda` and `loo`. A working JAGS installation is required. The default run reads stored model results, while refit options may take several hours.

If the manuscript is accepted, a fixed release of the repository will be permanently archived.

## References

Wilkinson, C., R. B. H. Lim, J. H. Liew, J. T. B. Kwik, C. L. Y. Tan, H. H. Tan, and D. C. J. Yeo. 2022. Empirical food webs of 12 tropical reservoirs in Singapore. *Biodiversity Data Journal* 10: e86192. [Publication record](https://doi.org/10.3897/BDJ.10.e86192)

Wilkinson, C. L., R. B. H. Lim, J. H. Liew, J. T. Kwik, D. C. J. Yeo, C. L. Y. Tan, and H. H. Tan. 2022. Empirical food webs of 12 tropical reservoirs in Singapore. Dryad Digital Repository. [Dataset](https://doi.org/10.5061/dryad.jsxksn088)
