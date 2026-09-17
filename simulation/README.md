Simulation code (reproducibility snapshot)
==========================================

This folder is the frozen January 2026 pipeline that generated the
simulation results used in

> Panaro, R., Röver, C., & Friede, T.
> *Subgroup comparisons within and across studies in meta-analysis.*
> Accepted for publication in *Research Synthesis Methods*.

The published simulation figures are produced from

```
results/result_sim_2026-01-13 23:09:47.033524.rds
```

in the parent repository (also copied here as
`../results/result_sim_2026-01-13 23:09:47.033524.rds` if you have the
local data files). That RDS is about 850 MB and is **not** stored in
Git. Please request it from the authors or from the journal/OSF data
deposit.

What this snapshot is
---------------------

| File | Role |
|------|------|
| `fastR/main.R` | Experiment grid, seeds, `run_batch()`, save RDS |
| `fastR/simulate_trial.R` | Data-generating process (IPD → mean-difference AgD) |
| `fastR/computations.R` | The 11 estimators (`fit_all`) |
| `fastR/cluster_engine.R` | L'Ecuyer streams and `clustermq` runner |
| `fastR/load_packages.R` | R packages |
| `R/wt_framework.R` | Within-trial method (sourced by workers) |
| `fastR/clustermq_slurm.tmpl`, `fastR/r_job.sh` | Cluster launch |
| `fastR/main_2026-01-09_snapshot.R` | 9 Jan 2026 driver (same grid; 5,000 replications) |

The experiment grid is 120 cells: 5 prevalence scenarios × 2 aggregation
settings (`a_interaction` = 2 or −1) × \(k \in \{10, 15, 20\}\) ×
\(\tau \in \{0, 0.1, 0.2, 0.5\}\), with base seed `123` and
**5,000 replications** per cell.

How to run a local check
------------------------

From this `simulation/` directory, with the packages listed below
installed:

```sh
Rscript fastR/main.R
```

This is `SIM_MODE=smoke`: 2 replications of 2 experiments on local
cores. It writes `results/result_sim_<timestamp>.rds`. It is only a
check that the code runs, not the paper results.

How to run the published simulation
-----------------------------------

The published simulation uses **5,000 replications** per experiment ×
method (~600,000 jobs). That is a cluster run.

```sh
# from simulation/
export SIM_MODE=paper
# optional, HPC library path:
# export LIB_PATH=~/R/x86_64-pc-linux-gnu-library/4.4
sbatch fastR/r_job.sh
```

Or, if `clustermq` is already configured for Slurm:

```sh
SIM_MODE=paper Rscript fastR/main.R
```

Output path: `results/result_sim_<Sys.time()>.rds`. The published
figures then come from `R/sim_viz_preamble.r` and `R/figure6.r` …
`R/figure23.r` in the parent repository, with

```sh
export SIM_RESULTS_RDS=/path/to/that.rds
```

Packages
--------

R 4.4.x (cluster template: `module load r/4.4.1`) and

- clustermq, tidyverse, assertthat, evaluate
- MASS, metafor, mixmeta, magrittr, nloptr, EnvStats

Caveats (please read)
---------------------

1. **Numeric identity.** Monte Carlo summaries should match within Monte
   Carlo error. Bitwise identical draws are not expected across BLAS / R
   patch versions.

2. **`evaluate` in `load_packages.R`.** `cluster_engine.R` calls
   `evaluate::try_capture_stack()`. The original `load_packages.R` did
   not attach `evaluate`; it is attached here so a clean session runs.

Applied examples
----------------

IL-6 and corticosteroid 2×2 counts are in `data/`. Forest plots and
method-comparison figures are the parent scripts `R/figure1.r`,
`R/figure2.r`, `R/figure9.r`–`R/figure12.r`.

File provenance
---------------

- `simulate_trial.R`, `cluster_engine.R`, `load_packages.R`, Slurm
  template: `archive_outside_inst_2026-04-16/fastR.zip` (pipeline as of
  January 2026).
- `main.R` experiment grid and `replications <- 5000`.
- `computations.R`: `R/computations.r` from the parent repository
  (method names match the plotted RDS, e.g. `"Equal weights SWADA"`).
- `R/wt_framework.R`: `R/wt_framework.r` from the parent repository.
