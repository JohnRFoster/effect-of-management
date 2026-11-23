# effect-of-management

This repository contains the model-fitting and analysis code for the manuscript “The effect of long-term management on wild pig (*Sus scrofa x domesticus*) populations across the southeastern United States.” 
The model-fitting script is `R/states.R`. 
The manuscript analysis uses the `targets` pipeline in `_targets.R`.

## Model Fitting

Run `R/states.R` in the configured high-performance computing (HPC) environment.
The script fits the state-by-state [removal model](https://github.com/JohnRFoster/boaR/blob/08e885212707c4a92b368d2084ce5f60c27a499e/R/nimble_removal_model.R#L8) for one state and for Texas by ecoregion. 
It uses NIMBLE for model fitting via the `boaR` [package](https://github.com/JohnRFoster/boaR/tree/main).

Set these environment variables before you run the script:

- `data_store`: path to the data store.
- `STATENAME`: uppercase state name to fit.
- `SLURM_ARRAY_TASK_ID`: ecoregion index from 1 to 5 for Texas.
- `SLURM_CPUS_PER_TASK`: number of MCMC chains.

The script reads processed MIS event data and `county_ecoregion.csv` from `data_store`. 
The script reads `data/original_ids.rds` from the project directory, which contains the original property IDs used to fit the original model. 
It uses the `prod` profile from the external `config` setup.

The script writes model output under `out/states2/states2-2026-03-25/<STATE>/`. 
It also writes MCMC diagnostics under the `analysis` subdirectory.

`R/ecoregions.R` creates the county-to-ecoregion table. 
It needs the 2024 TIGER county data and the `NorthAmericaEcoregionsLevel2/NA_CEC_Eco_Level2.shp` shapefile in the data store. 
The fitting script also requires project helper functions and configuration from the HPC environment. 
An example bash script is below.

## Manuscript Analysis

The `targets` pipeline reads fitted-model files.
The pipeline requires these files:

- `model_data.rds`
- `all_density.rds`
- `all_params_samples.rds`
- `all_params.rds`
- `priors_df.rds`
- `all_abundance_samples.rds`

From the project directory, run:

```r
targets::tar_make()
```

The pipeline prepares the model data. It calculates realized growth rates,
method-specific removal proportions, management bootstrap simulations,
effort-population relationship (EPR) data, and density intervals. It creates
the manuscript figures. The figure functions write TIFF files to `plots/`.
The targets store is `_targets/`.

## Report

The Quarto document `results.qmd` compiles the results text, analysis, and figures. 

## Example bash script

```
#!/bin/env bash
#SBATCH --job-name=fs-illinois
#SBATCH --output=outfiles/%x-%j-%a.txt 
#SBATCH --cpus-per-task=3
#SBATCH --mem=128g
#SBATCH --ntasks=1
#SBATCH --nodes=1

echo "Loading modules"
module load R/4.4.1
module load openblas/dynamic/0.3.18

# path to R library
echo "Set R library path"
MY_R_LIBRARY_PATH="path/to/R/library"

# export to R
echo "Export library path to R"
export R_LIBS_USER="$MY_R_LIBRARY_PATH"

# check statename
set -euo pipefail
echo "STATENAME in job: '${STATENAME-UNSET}'"
env | grep '^STATENAME=' || true

# run
echo "Running R script..."
Rscript R/states.R
echo "Complete"
```