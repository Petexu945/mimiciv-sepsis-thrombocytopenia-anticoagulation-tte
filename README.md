# mimiciv-sepsis-thrombocytopenia-anticoagulation-tte
Code for a target trial emulation of early prophylactic anticoagulation in severe sepsis-associated thrombocytopenia using MIMIC-IV.
# MIMIC-IV Sepsis Thrombocytopenia Anticoagulation Target Trial Emulation

This repository contains the SQL and R code used for the study:

**Early prophylactic anticoagulation in patients with severe sepsis-associated thrombocytopenia: a target trial emulation**

The study emulates a target trial comparing early pharmacological venous thromboembolism prophylaxis versus no early pharmacological prophylaxis among critically ill adults with severe sepsis-associated thrombocytopenia, defined by platelet counts of 30,000–50,000/µL.

## Repository contents

```text
sql/
  01_build_primary_cohort.sql
  02_revision_supplementary_analyses.sql

R/
  01_main_analysis.R
  02_make_tables_and_figures.R

docs/
  target_trial_protocol_summary.md
  variable_definitions.md
  itemid_and_term_lists.md

config/
  .env.example
```

## Data availability

This repository contains code only. It does not contain patient-level MIMIC-IV data, derived datasets, intermediate tables, analysis-ready datasets, or any database extracts.

The study uses MIMIC-IV data, which are available from PhysioNet to credentialed users who have completed the required data use training and have signed the data use agreement. Users who wish to reproduce the analysis must obtain independent authorized access to MIMIC-IV and run the code within their own credentialed environment.

## Software requirements

The SQL scripts were written for a PostgreSQL implementation of MIMIC-IV.

The R scripts were developed using R 4.4.3. Main R packages include:

```text
tidyverse
DBI
RPostgres
survival
survminer
cobalt
WeightIt
MatchIt
tableone
survey
broom
patchwork
forestploter
```

Package versions can be recorded using `sessionInfo()` after running the analysis.

## Database configuration

Database credentials should not be written directly into the analysis scripts. Instead, set the following environment variables before running the R code:

```text
MIMIC_DBNAME
MIMIC_HOST
MIMIC_PORT
MIMIC_USER
MIMIC_PASSWORD
```

An example configuration file is provided in `config/.env.example`. The actual `.env` file should not be committed to GitHub.

## Analysis workflow

Run the scripts in the following order:

1. `sql/01_build_primary_cohort.sql`

   This script constructs the primary target trial emulation cohort, including adult first ICU stays, Sepsis-3 eligibility, platelet-based time zero, exclusion criteria, treatment assignment during the 24-hour grace period, baseline covariates, outcome definitions, and the final analysis table.

2. `sql/02_revision_supplementary_analyses.sql`

   This script creates reviewer-requested supplementary and descriptive tables, including bleeding endpoint components, VTE-related imaging intensity, pharmacological prophylaxis agent distribution, available documentation of mechanical prophylaxis, calendar-era summaries, and post-landmark treatment switching summaries.

3. `R/01_main_analysis.R`

   This script performs the main statistical analyses, including inverse probability weighting, covariate balance diagnostics, primary weighted Cox models, entropy balancing, propensity score matching, secondary outcome analyses, net clinical benefit estimation, subgroup analyses, and sensitivity analyses.

4. `R/02_make_tables_and_figures.R`

   This script generates publication-ready tables and figures from the outputs of the main analysis script.

## Main outputs

The analysis scripts generate summary outputs such as:

```text
Primary Cox model results
Entropy balancing sensitivity analysis
Propensity score matching sensitivity analysis
Baseline balance tables
Secondary outcome tables
Net clinical benefit estimates
Subgroup analysis results
Kaplan–Meier curves
Forest plots
Reviewer-requested supplementary tables
```

These outputs are generated locally and are not included in this repository because they may be derived from restricted-access MIMIC-IV data.

## Important notes

No patient-level data are included in this repository.

No intermediate derived datasets are included in this repository.

No database credentials are included in this repository.

Users must have authorized MIMIC-IV access and must configure their own local PostgreSQL database environment before running the scripts.

## License

This code is provided for academic transparency and reproducibility. Please cite the corresponding manuscript if using or adapting this code.
