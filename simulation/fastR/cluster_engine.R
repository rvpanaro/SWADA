run_batch <- function(
    experiments, # Which experiments to run (tibble)
    batch_replications, # How many replications per experiment?
    job_function, # This function is run on the workers after replication setup
    worker_setup_function = default_worker_setup(), # This function is run on the workers on worker startup
    job_wrapper = default_job_wrapper, # This function is run on the workers to set up each replication
    test_single_job_index = integer(), # Use this to test single jobs (e.g. if job_id 11 crashes, set this to 11)
    test_single_job_per_experiment = FALSE, # Use to test 1 replication per experiment
    test_single_job_externally = TRUE, # By FALSE test locally, setting this to TRUE tests on a worker
    job_resources = list( # Each worker will request these resources
      walltime = 60, # Wall time in minutes
      memory = 3072, # Memory in megabytes
      cores = 1 # Number of CPU cores per worker (R is single-threaded by default, so default to 1)
    ),
    do_gc = FALSE, # Garbage collect between jobs?
    n_jobs = 200, # The maximum simultaneous jobs to start
    timeout = 300, # Job timeout in seconds
    log_worker = FALSE, # Set this to true to create logs
    max_calls_worker = Inf # This option can limit the number of jobs run per worker (useful for debugging random segfaults)
    ) {
  if (test_single_job_per_experiment) {
    assert_that(length(test_single_job_index) == 0)
  }

  options(clustermq.data.warning = 999)

  cat(paste0("Clustermq run - enumerating clustermq experiments...\n"))

  clustermq_args <- list()

  replication_index <- integer(nrow(experiments) * batch_replications) * NA
  experiment_vec <- character()
  seed_vec <- list()

  max_job_id <- 0L
  for (experiment_index in seq_len(nrow(experiments)))
  {
    experiment_seeds <- setup_lecuyer_seeds(seed = experiments[["seed"]][experiment_index], num = batch_replications)

    if (test_single_job_per_experiment) {
      test_single_job_index <- c(test_single_job_index, max_job_id + 1)
    }

    replication_index[max_job_id + seq_len(batch_replications)] <- seq_len(batch_replications)
    max_job_id <- max_job_id + batch_replications
    experiment_vec <- c(experiment_vec, rep(experiments[["experiment"]][experiment_index], batch_replications))

    seed_vec <- c(seed_vec, experiment_seeds)
  }

  clustermq_jobs <- tibble(
    job_id = seq_len(max_job_id),
    replication_index = replication_index,
    experiment = experiment_vec,
    seed = seed_vec
  )

  Q_const_arg <- list(
    "job_function"     = job_function,
    "setup_function"   = worker_setup_function,
    "num_cores"        = job_resources$cores,
    "do_gc"            = do_gc,
    "experiments"      = experiments
  )

  cat(paste0("Clustermq run - submitting and running experiments...\n"))
  cat(paste0("Clustermq run - submitting jobs at ", Sys.time(), "\n"))

  if (length(test_single_job_index) > 0) { # Set this for local testing
    assert_that(
      all(test_single_job_index > 0 & test_single_job_index <= nrow(clustermq_jobs)),
      msg = "Test job indices out of range"
    )

    clustermq_jobs <- clustermq_jobs[test_single_job_index, ]

    if (!test_single_job_externally) {
      options("backup_scheduler" = getOption("clustermq.scheduler"))
      options("backup_error" = getOption("error"))
      options("clustermq.scheduler" = "local")
      options("worker_setup_required" = TRUE)
      options(error = \(){
        options("clustermq.scheduler" = getOption("backup_scheduler"))
        options(error = getOption("backup_error"))
        recover()
      })
    }
  }

  n_jobs_final <- max(min(nrow(clustermq_jobs), n_jobs), ceiling(nrow(clustermq_jobs) / max_calls_worker))

  outputs <- Q_rows(
    clustermq_jobs,
    fun = job_wrapper,
    const = Q_const_arg,
    template = job_resources,
    n_jobs = n_jobs_final,
    max_calls_worker = max_calls_worker,
    fail_on_error = TRUE,
    timeout = timeout,
    log_worker = log_worker
  )

  errored_jobs <- sapply(outputs, \(x) {
    is(x, "simpleError")
  })

  if (any(errored_jobs)) {
    errored_job_indices <- which(errored_jobs)
    for (job_index in errored_job_indices) {
      print(outputs[[job_index]])
      message("Error in job_id ", clustermq_jobs[job_index, ]$job_id, ": ", outputs[[job_index]]$message)
      message("Call stack:")
      for (i in seq_along(outputs[[job_index]]$calls)) {
        message(i, ": ", outputs[[job_index]]$calls[i])
      }
      cat("\n")
    }

    stop("At least one job failed.")
  }

  cat(paste0("Clustermq run - finished waiting for jobs at ", Sys.time(), "\n"))

  return(list(clustermq_jobs = clustermq_jobs, outputs = outputs))
}

default_worker_setup <- function(
    working_dir = getwd(), # Working directory
    additional_source_files = character(0L), # Source these additional files
    lib_paths = .libPaths(), # Use these .libPaths() - by default same as on master node
    transfer_options = TRUE, # Transfer options() from master node to workers
    source_files = \(working_dir) {
      file.path(working_dir, c("load_packages.R", "cluster_engine.R"))
    },
    additional_setup = function(working_dir) {}) {
  # Ensure arguments are evaluated so the function that this function returns has these arguments available
  working_dir
  lib_paths
  transfer_options
  source_files
  additional_source_files
  additional_setup

  if (transfer_options) {
    master_options <- options()
  }

  function() {
    # This is the function that will get run on worker startup
    if (transfer_options) {
      options(master_options)
    }
    setwd(working_dir)

    if (!is.null(lib_paths)) {
      .libPaths(lib_paths)
    }
    all_source_files <- c(source_files(working_dir), additional_source_files)

    for (file_name in all_source_files) {
      source(file_name)
    }

    if (transfer_options) {
      options(master_options)
    }

    additional_setup(working_dir = working_dir)
  }
}

default_job_wrapper <- function(
    setup_function,
    job_function,
    job_id,
    replication_index,
    experiment,
    experiments,
    seed,
    do_gc,
    num_cores) {
  cat("Current job: ", job_id, "\n")

  if (getOption("worker_setup_required", default = TRUE)) {
    cat("Calling setup function:\n")
    setup_function()
    options("worker_setup_required" = FALSE)
  } else {
    if (do_gc) {
      cat("Calling Garbage Collector:\n")
      print(gc())
    }
  }

  RNGkind("L'Ecuyer-CMRG")
  .Random.seed <<- seed

  # Use as many cores as specified
  options(mc.cores = num_cores)

  evaluate::try_capture_stack(
    {
      job_function(
        experiment = experiment,
        experiment_args = experiments[experiments$experiment == experiment, ],
        job_id = job_id,
        replication_index = replication_index
      )
    },
    env = environment()
  )
}


setup_lecuyer_seeds <- function(seed = NULL, lecuyer_seed = NULL, num) {
  RNGkind("L'Ecuyer-CMRG")

  assert_that(is.null(lecuyer_seed) + is.null(seed) == 1)

  if (!is.null(seed)) {
    set.seed(seed)
    lecuyer_seed <- .Random.seed
  }

  job_seeds <- list()
  job_seeds[[1]] <- parallel::nextRNGStream(lecuyer_seed)
  i <- 2
  while (i < num + 1) {
    job_seeds[[i]] <- parallel::nextRNGStream(job_seeds[[i - 1]])
    i <- i + 1
  }
  job_seeds
}

reduce_to_tibble <- function(result, reduce_what = "output_tibble") {
  cat(paste0("Clustermq run - starting reduce at ", Sys.time(), "\n"))

  outputs_tibble <- result[["clustermq_jobs"]] %>% inner_join(
    bind_rows(lapply(result[["outputs"]], function(x) x[[reduce_what]])),
    by = "job_id"
  )

  cat(paste0("Clustermq run - finished reduce at ", Sys.time(), "\n"))
  cat(paste0("Clustermq run - reducing results... done\n"))

  outputs_tibble
}
