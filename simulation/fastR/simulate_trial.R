simulate_trial <- function(experiment, experiment_args, job_id, replication_index) {
  
  # Extract parameters from experiment argumentz
  # Generate random tau2 values
  tau1_values <- experiment_args$tau1
  tau2_values <- experiment_args$tau2
  
  # Define the size of the studies
  size <- experiment_args$size
  k <-  experiment_args$n_trials
  
  # Define study sizes based on scenario
  if (size == "A: All similar") { # A: All similar studies
    n <- rep(1, k)
  } else if (size == "B: One Large Study") { # B: one large study
    n <- c(rep(1, k - 1), 10)
  } else if (size == "C: 10% Large Studies") { # C: 10% large studies
    n <- c(rep(1, k * 0.9), rep(10, k * 0.1))
  } else { # D: 25% large studies
    n <- c(rep(1, k * 0.25), rep(10, k * 0.75))
  }
  
  # Correlation between sample sizes
  ss_cor <- .75
  
  # Covariance between sample size values
  mcov <- diag(1 - ss_cor, k) + matrix(ss_cor, k, k)
  
  # Function to generate a population for each meta-analysis
  population_gen <- function(meta) {
    
    # Heterogeneity values
    tau1 <- tau1_values[meta]
    tau2 <- tau2_values[meta]
    
    # Generate study sizes
    trial_sizes <- n * exp(mvrnorm(n = 1, mu = rep(log(experiment_args$n_patients), k), Sigma = mcov)) |> ceiling()
    
    # Force each sample size to be divisible by the twice the prevalence since prevalence = 2/integer:
    divisor <- 2 / experiment_args$prevalence
    trial_sizes <- ceiling(trial_sizes / divisor) * divisor
    
    trial_sizes[trial_sizes < 12] <- 12
    
    # This creates a generic patient id.
    patient_ids <- sapply(trial_sizes, function(x) 1:x) |> unlist()
    
    #This creates the trial ids
    trial_ids <- sapply(1:k, function(x) rep(x, times = trial_sizes[x])) |> unlist()
    
    # This creates, for each trial, a vector with half 0's and half 1's (shuffled randomly).
    treatment_ids <- lapply(trial_sizes, function(x) {
      rep(c(0, 1), times = x / 2)
    }) |> unlist()
    
    # Generated data on IPD level indexes
    population <- experiment_args |>
      tibble(
        patient_id = patient_ids,
        trial_id = trial_ids,
        treatment_id = treatment_ids) 
    
    
    # For each trial, determine the prevalence
    population <- population |> 
      group_by(trial_id)  |> 
      mutate(
        subgroup_prevalence = switch(first(scenario),
                                     "Skewed (triangular [0.1, 0.5, c= 0.3])" = rtri(n = 1, min = prevalence - V1, max = prevalence + V1, mode = c),
                                     runif(n = 1, min = prevalence - V1, max = prevalence + V1))
      ) |> 
      ungroup()
    
    # For each trial, determine the number of subgroup 1's and randomly assign them
    population <- population |> 
      group_by(trial_id) |> 
      mutate(
        n_trial   = n(), 
        n_success = round(n_trial * subgroup_prevalence),
        subgroup_id = c(rep(x = 1, times = first(n_success)),
                        rep(x = 0, times = first(n_trial) - first(n_success)))
      )   |>
      ungroup() |>
      # Compute the trial-specific subgroup proportion for every observation
      mutate(
        subgroup_proportion = ave(subgroup_id, trial_id, FUN = mean))
    
    population <- population |> 
      group_by(trial_id) |> 
      mutate(
        alphai    = runif(1, 0, 0.5),
        beta1i    = runif(1, 0, 0.5),
        a_interaction_id = treatment_id * subgroup_proportion,
        w_interaction_id = treatment_id * (subgroup_id - subgroup_proportion)) |> 
      mutate(
        beta2i    = (phi + a_interaction * subgroup_proportion) + rnorm(n(), mean = 0, sd = tau1),
        gammai    = w_interaction + rnorm(n(), mean = 0, sd = tau2),
        predictor = alphai + beta1i * subgroup_id + beta2i * treatment_id + gammai * (treatment_id * (subgroup_id - subgroup_proportion)),
        measurement = predictor + rnorm(n(), mean = 0 , sd = experiment_args$sigma),
        interaction_id = factor(treatment_id * subgroup_id)
      )
    
    # Aggregate the data using the native pipe
    aggregated <- population |>
      mutate(
        subgroup_id = subgroup_id,
        trial_id = trial_id,
        treatment_id = treatment_id
      ) |>
      group_by(trial_id, treatment_id, subgroup_id) |>
      summarise(
        mi = mean(measurement),
        ni = n(),
        sdi = experiment_args$sigma, # sigma assumed
        .groups = "drop"
      ) 
    
    rm(population)
    
    # Pivot to wide format so each row has both treatment groups’ stats
    aggregated_wide <- aggregated |>
      pivot_wider(
        names_from = treatment_id,
        values_from = c(mi, sdi, ni),
        names_prefix = "treatment_"
      )
    
    # Compute effect sizes using metafor::escalc; use an anonymous function
    effect_size <- aggregated_wide |>
      (\(df) escalc(
        measure = "MD",
        m1i = mi_treatment_1,   m2i = mi_treatment_0,
        sd1i = sdi_treatment_1, sd2i = sdi_treatment_0,
        n1i = ni_treatment_1,   n2i = ni_treatment_0,
        data = df
      ))()
    
    # Optionally, adjust the total sample size per subgroup if needed
    effect_size$ni <- aggregated_wide$ni_treatment_1 + aggregated_wide$ni_treatment_0
    
    effect_size$sample_prevalence <- rep(effect_size$ni[seq(2, nrow(effect_size), 2)]/
                                           (effect_size$ni[seq(1, nrow(effect_size), 2)] + effect_size$ni[seq(2, nrow(effect_size), 2)]), 
                                         each = 2)
    
    effect_size$sample_sigma <- sqrt(effect_size$ni * effect_size$vi)
    
    rm(aggregated)
    
    return(list(# ipd = population, 
      agd = effect_size))
  }
  
  # Generate population data for all meta-analyses
  population <- population_gen(1)
  
  population$agd <- population$agd |> 
    dplyr::rename(trial = trial_id, subgroup = subgroup_id) |>  
    as.data.frame() |>
    arrange(trial, desc(subgroup))
  
  # Define the joint data (IPD and AgD)
  dat <- tribble(
    # ~ipd, 
    ~agd, 
    # population$ipd,
    population$agd)
  
  # Simulation targets
  target_gamma <- experiment_args$w_interaction
  target_A <- with(experiment_args, phi + prevalence * (a_interaction - w_interaction))
  target_B <- with(experiment_args, phi + prevalence * (a_interaction - w_interaction)+  w_interaction)
  target_tau <- experiment_args$tau1
  target_tauW <- experiment_args$tau2
  prevalence <- experiment_args$prevalence
  
  ## Fit the models and handle errors
  computation <- experiment_args$computation[[1]](dat)
  
  # Pre-computing a tibble with stats for each simulation study
  # on the workers allows for rapid result aggregation on the master node
  
  regional <- as.logical(apply(cbind(computation$subgroupA_b,  
                                     computation$subgroupB_b,
                                     computation$variance_A,
                                     computation$variance_B,
                                     computation$covariance_B_A
  ), 1, function(x){
    V <- matrix(c(x[1+2],x[3+2],x[3+2],x[1+2]), byrow = T, nrow = 2)
    result <- try({(c((x[1:2]- c(target_A, target_B)) %*% MASS::ginv(V) %*%
                        c(x[1:2]- c(target_A, target_B))) <= qchisq(.95, df = 2))})
    
    if(class(result) == "try-error"){
      return(NA)
    } else{
      return(result)        
    } 
  }))
  
  result <- list(output_tibble = tibble(
    "job_id" = job_id,
    "method" = computation$analysis,
    
    # Interaction estimates metric
    "interaction_b" = c(computation$interaction_b),
    "interaction_bias" = c(computation$interaction_b) - target_gamma,
    "interaction_se" = c(computation$interaction_se),
    "interaction_tau" = c(computation$interaction_tau),
    "interaction_tau_bias" = c(computation$interaction_tau) - target_tauW,
    "interaction_tau_zero" = round(c(computation$interaction_tau), 4) == 0,
    "interaction_width" = 1.959964 * c(computation$interaction_se),
    "interaction_lb" = interaction_b-interaction_width,
    "interaction_ub" = interaction_b+interaction_width,
    "interaction_contained" = (interaction_lb < target_gamma) & (interaction_ub > target_gamma),
    
    # Subgroup A estimates metric
    "subgroupA_b" = c(computation$subgroupA_b),
    "subgroupA_bias" = c(computation$subgroupA_b) - target_A,
    "subgroupA_se" = c(computation$subgroupA_se),
    "subgroupA_tau" = c(computation$subgroupA_tau),
    "subgroupA_tau_bias" = c(computation$subgroupA_tau) - target_tau,
    "subgroupA_tau_zero" = round(c(computation$subgroupA_tau), 4) == 0,
    "subgroupA_width" = 1.959964 * c(computation$subgroupA_se),
    "subgroupA_lb" = subgroupA_b - subgroupA_width,
    "subgroupA_ub" = subgroupA_b + subgroupA_width,
    "subgroupA_contained" = (subgroupA_lb < target_A) & (subgroupA_ub > target_A),
    
    # Subgroup B estimates metric
    "subgroupB_b" = c(computation$subgroupB_b),
    "subgroupB_bias" = c(computation$subgroupB_b) - target_B,
    "subgroupB_se" = c(computation$subgroupB_se),
    "subgroupB_tau" = c(computation$subgroupB_tau),
    "subgroupB_tau_bias" = c(computation$subgroupB_tau) - sqrt(target_tau^2 + target_tauW^2),
    "subgroupB_tau_zero" = round(c(computation$subgroupB_tau), 4) == 0,
    "subgroupB_width" = 1.959964 * c(computation$subgroupB_se),
    "subgroupB_lb" = subgroupB_b - subgroupB_width,
    "subgroupB_ub" = subgroupB_b + subgroupB_width,
    "subgroupB_contained" = (subgroupB_lb < target_B) & (subgroupB_ub > target_B),
    
    # Target values
    "tau" = target_tau,
    "tauW" = target_tauW,
    "V1" = experiment_args$V1,
    "a_interaction" = experiment_args$a_interaction,
    "delta" = c((computation$subgroupB_b -computation$subgroupA_b) - computation$interaction_b),
    "converge" = c(!is.na(computation$interaction_b)),
    "k" = experiment_args$n_trials,
    "scenario" = experiment_args$scenario,
    "size" = experiment_args$size,
    "collapsible" = experiment_args$a_interaction == 0,
    "a" =  experiment_args$a, 
    "true_hetero" = computation$true_hetero,
    "target_gamma" = target_gamma,
    "target_A" = target_A,
    "target_B" = target_B,
    "target_tau" = target_tau,
    "target_tauW" = target_tauW,
    "regional_coverage" = c(regional)
  ))
  
  return(result)
}

simulate_error <- function() {
  stop("This is a simulated error")
}
