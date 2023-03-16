# Shielding individuals at high risk of COVID-19: a micro-simulation study. Model and analyses.

This repository contains the code used for all analyses described in our manuscript:

Van Zandvoort K, Favas, C, Checchi, F, *Shielding individuals at high risk of COVID-19: a micro-simulation study. Model and analyses.*. Available at <https://doi.org/10.1101/2022.01.03.22268675>.

## Requirements

C++ libraries: GSL, openmp
R packages: pacman, Rcpp, RcppGSL, data.table, distcrete, magrittr, socialmixr, qs, ggplot2)

All model runs and analyses were conducted using Ubuntu 20.04.3 LTS and R4.1.2.

### Generating model output

*Note that we specified 5292 scenarios that were each ran 1000 times for the analyses in our main manuscript. If ran with the same configuration (default), the model generates a lot of output (22.1GB compressed data files) and would take a long time to run requiring parallelization of individual model runs.*

The following files are used to run the model:

- `./create_scenarios.R`
  - Generates a list with all scenarios that will be ran by the model.
  - We used 5292 different scenarios in our analysis.
- `./create_prevalence_by_age.R`
    - Calculates the average prevalence by age across model runs for the unmitigated scenarios.
    - Used in scenarios where prior immunity is assumed (through setPriorImmunity).
- `./single_run.R`
  - Randomly generates a new population from the source data
  - Runs all scenarios for a single model iteration
- `./scripts/functions.R`
  - Helper functions used by the model. Used in `single_run.R`.
- `./scripts/read_data.R`
  - Downloads Somaliland IDP demographic and contact data from [GitHub](https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019) [1]
  - Prepares non-intervention parameters and distributions used in the model
- `./data`
  - Other datasets used in the model
  - Stores backup of Somaliland contact and demographic datasets used in the analysis
- `./model/covidIBM_shielding.cpp`
  - Code for model. Is compiled in R.
  - Requires the GSL library to run (though RcppGSL).

We used a SLURM scheduler to parallelize our model runs. The files we used are provided here for reference, but may require adapting to local systems.
- `./single_run.sh`
  - Shell script to call `./single_run.R` for a single iteration of the model
- `./submit_all_runs.sbatch`
  - Calls `./single_run.sh` for each submitted task (model iteration)

### Processing model output

We saved our model output as a single file per model iteration with output for all scenarios. This ensured that random numbers drawn and used in preparing the modelled populations and running the stochastic model (up to the point where an intervention changes the course of the outbreak) are equivalent between model iterations. However, analysis of the output requires analysis of all iterations per scenario, not analysis of all scenarios per model iteration.

We first split our data for each model iteration and save an individual file for each scenario per iteration. Again, we used a SLURM scheduler to speed up this process (*Nb. This could already have been done in generating the data in `single_run.R`*).

- `split_run.R`
- `split_run.sh`
- `submit_split_all_runs.sbatch`

We then recombine data from all model iterations per scenario.

- `combine_single_scen.R`
- `combine_single_scen.sh`
- `submit_combine_all_scens.sbatch`

### Data analysis

All analysis scripts can be found in the `./analysis` folder. We use one script for each analysis, which correspond to the different figures presented in our manuscript.

Model output requires different processing in each scenario. To process and summarize output of sometimes a large number of model output files, we again utilized SLURM to speed up the process.

We specify a different `singleFunction` function in each analysis file. These specify what scenarios are needed in this analysis, and how model output files should be summarized.
The `prepareSlurmBatchRuns` function creates a folder with all files needed for a SLURM `sbatch` command, akin to the `.R`, `.sh`, and `.sbatch` files in the previous sections.
Upon completion of the jobs, the `combineSlurmResults` function combines all processed results one final time, before the specific analysis is performed on this data object.

The SLURM related functions are specified in `./analyis/helper_functions.R`

The different analyses are included as:

- `./analysis/breached.R`
- `./analysis/breached_mitigation.R`
- `./analysis/breached_time.R`
- `./analysis/cumulative_risk.R`
- `./analysis/risk_bytime_full.R`
- `./analysis/risk_bytime_full2.R`

[1] Van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman MS, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C, Mulholland K, Egeh MM, Hassan MM, Hergeeye MA, Eggo RM, Checchi F, Flasche S, Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Available at https://doi.org/10.1016/j.epidem.2022.100625.
