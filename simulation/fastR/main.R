# Simulation driver for Panaro, Röver & Friede
# "Subgroup comparisons within and across studies in meta-analysis"
# Accepted for publication in Research Synthesis Methods.
#
# Frozen snapshot of the January 2026 pipeline that produced
#   result_sim_2026-01-13 23:09:47.033524.rds
# (the file used by the simulation figures).
#
# Usage (from the simulation/ directory):
#   Rscript fastR/main.R                  # small local check (default)
#   SIM_MODE=paper Rscript fastR/main.R   # published 5,000-replication cluster run
#
# The 9 Jan 2026 copy of the driver is fastR/main_2026-01-09_snapshot.R.

find_simulation_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    script_dir <- dirname(script_path)
    if (identical(basename(script_dir), "fastR")) {
      return(dirname(script_dir))
    }
    return(script_dir)
  }
  if (file.exists("fastR/main.R")) {
    return(normalizePath("."))
  }
  if (file.exists("main.R") && dir.exists("../R")) {
    return(normalizePath(".."))
  }
  getwd()
}

root <- find_simulation_root()
setwd(root)
relative_path <- file.path(root, "fastR")

# One BLAS thread per R process (avoids oversubscription with clustermq).
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

lib_path <- Sys.getenv("LIB_PATH", unset = paste(.libPaths(), collapse = .Platform$path.sep))
log_file <- Sys.getenv("LOG_FILE", unset = file.path(root, "logs", "clustermq.log"))

source(file.path("fastR", "load_packages.R"))
source(file.path("fastR", "simulate_trial.R"))
source(file.path("fastR", "cluster_engine.R"))
source(file.path("fastR", "computations.R"))

# -----------------------------------------------------------------------------
# Run mode
#   smoke (default): 2 replications of 2 experiments, local multiprocess
#   paper:           5,000 replications of all 120 experiments (cluster)
#
# The plotted RDS has replication_index 1..5000.
# -----------------------------------------------------------------------------
sim_mode <- Sys.getenv("SIM_MODE", unset = "smoke")
if (!sim_mode %in% c("smoke", "paper")) {
  stop("SIM_MODE must be 'smoke' or 'paper', not: ", sim_mode, call. = FALSE)
}

tau_values <- c(0.0, 0.1, 0.2, 0.5)

experiments <- tribble(
  ~scenario, ~computation, ~seed, ~n_patients, ~phi, ~a_interaction, ~w_interaction, ~prevalence, ~V1, ~c, ~sigma,
  "Identical/balanced (1:1)",            fit_all, 123, exp(4.4), 2,  2, -1, 0.5, 0.0, 0.5, 2,
  "Identical/imbalanced (1:3)",          fit_all, 123, exp(4.4), 2,  2, -1, 1/3, 0.0, 1/3, 2,
  "Less variation",                      fit_all, 123, exp(4.4), 2,  2, -1, 0.5, 0.2, 0.5, 2,
  "More variation",                      fit_all, 123, exp(4.4), 2,  2, -1, 0.5, 0.4, 0.5, 2,
  "Skewed variation",                    fit_all, 123, exp(4.4), 2,  2, -1, 1/3, 0.2, 1/3, 2,
  # No aggregation bias (a_interaction == w_interaction)
  "Identical/balanced (1:1)",            fit_all, 123, exp(4.4), 2, -1, -1, 0.5, 0.0, 0.5, 2,
  "Identical/imbalanced (1:3)",          fit_all, 123, exp(4.4), 2, -1, -1, 1/3, 0.0, 1/3, 2,
  "Less variation",                      fit_all, 123, exp(4.4), 2, -1, -1, 0.5, 0.2, 0.5, 2,
  "More variation",                      fit_all, 123, exp(4.4), 2, -1, -1, 0.5, 0.4, 0.5, 2,
  "Skewed variation",                    fit_all, 123, exp(4.4), 2, -1, -1, 1/3, 0.2, 1/3, 2
)

experiment_args <- experiments |>
  expand_grid(n_trials = c(10, 15, 20)) |>
  expand_grid(size = c("A: All similar")) |>
  expand_grid(tau1 = tau_values) |>
  expand_grid(a = c(0.5)) |>
  mutate(tau2 = a * tau1) |>
  mutate(experiment = paste0(
    "Collapsible =", a_interaction - w_interaction == 0,
    ", scenario =", scenario,
    ", V1 =", round(V1, 4),
    ", n_trials =", n_trials,
    ", size = ", size,
    ", mode = ", c,
    ", prevalence = ", prevalence,
    ", a =", round(a, 4),
    ", tau1 =", round(tau1, 4),
    ", tau2 =", round(tau2, 4)
  ))

if (any(duplicated(experiments))) {
  stop("At least one experiment is duplicated in the tibble.")
}

if (identical(sim_mode, "paper")) {
  replications <- 5000L
  hpc_parallelize <- TRUE
  n_jobs_local <- max(1L, parallel::detectCores() - 2L)
} else {
  replications <- as.integer(Sys.getenv("SIM_REPLICATIONS", unset = "2"))
  hpc_parallelize <- FALSE
  n_jobs_local <- max(1L, parallel::detectCores() - 2L)
  experiment_args <- experiment_args |>
    dplyr::filter(
      n_trials == 10,
      tau1 == 0,
      scenario == "Identical/balanced (1:1)"
    )
  message(
    "SIM_MODE=smoke: ", nrow(experiment_args), " experiment(s) x ",
    replications, " replication(s). Set SIM_MODE=paper for the published run."
  )
}

dir.create(file.path(root, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "logs"), showWarnings = FALSE, recursive = TRUE)

# Targets used in the paper: within-trial interaction = w_interaction = -1
target_gamma <- experiment_args$w_interaction
target_b1 <- with(experiment_args, phi + prevalence * a_interaction + w_interaction * (0 - prevalence))
target_b2 <- with(experiment_args, phi + prevalence * a_interaction + w_interaction * (1 - prevalence))
stopifnot(isTRUE(all.equal(as.numeric(target_b2 - target_b1), as.numeric(target_gamma))))

if (hpc_parallelize) {
  options(
    clustermq.scheduler = "slurm",
    clustermq.template  = file.path(root, "fastR", "clustermq_slurm.tmpl"),
    clustermq.ssh.log = log_file,
    clustermq.ssh.timeout = 30
  )
} else {
  options(
    clustermq.scheduler = "multiprocess",
    clustermq.ssh.log = log_file,
    clustermq.ssh.timeout = 30
  )
}

begin_time <- Sys.time()
message("Starting ", sim_mode, " run at ", begin_time)

result_sim <- run_batch(
  experiments = experiment_args,
  batch_replications = replications,
  job_function = simulate_trial,
  worker_setup_function = default_worker_setup(
    lib_paths = strsplit(lib_path, .Platform$path.sep, fixed = TRUE)[[1]],
    working_dir = relative_path,
    additional_source_files = c("simulate_trial.R", "../R/wt_framework.R", "computations.R")
  ),
  job_resources = list(
    walltime = "2-00:00:00",
    cores = 1
  ),
  log_worker = FALSE,
  test_single_job_externally = FALSE,
  timeout = 360000,
  do_gc = FALSE,
  n_jobs = ifelse(hpc_parallelize, 1000, n_jobs_local)
)

reduced_sim <- reduce_to_tibble(result_sim)
out_file <- file.path("results", paste0("result_sim_", begin_time, ".rds"))
saveRDS(reduced_sim |> dplyr::select(-seed), file = out_file)
message("Wrote ", out_file)

# The published figures load this file (or a copy of it) via:
#   results/result_sim_2026-01-13 23:09:47.033524.rds
# Then run the scripts R/figure6.r ... R/figure23.r in the parent repository.
