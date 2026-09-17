# wt_framework.R
#
#    Copyright (C) 2024  Renato Panaro
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# R code to allow for the calculations shown in:
#
#   Panaro, R., Röver, C., & Friede, T.
#   Subgroup comparisons within and across studies in meta-analysis.
#   Accepted for publication in Research Synthesis Methods.

# within_meta_analysis Function: Perform within-trial interaction and subgroup-specific meta-analysis
# ------------------------------------------------------------------------------
# This function conducts a multivariate meta-analysis to estimate within-trial 
# interactions and subgroup-specific effects. It takes effect sizes, subgroup 
# information, and other parameters as input, then calculates the interaction 
# and subgroup-specific effects using specified heterogeneity structure and
# estimation method.
# 
# Args:
#   - effect_sizes: A data frame or escalc object containing effect sizes and 
# variances.
#   - reference_group: The reference group for interactions.
#   - use_last_reference: If TRUE, the last subgroup is taken as the reference 
# group.
#   - trial_labels: Name of the column containing trial identifiers.
#   - subgroup_moderator: Name of the column containing subgroup information.
#   - heterogeneity_structure: Structure of the heterogeneity matrix ("UN" or 
#                                                                     "CS").
#   - estimation_method: Estimation method ("FE" or "REML").
# 
# Returns:
#   A list containing interaction and subgroup-specific estimates, their 
# confidence intervals,
#   standard deviations, and hat values.
# ------------------------------------------------------------------------------

wt_framework <- 
  function(effect_sizes, reference_group, use_last_reference = TRUE, trial_labels, subgroup_moderator, heterogeneity_structure = "UN", estimation_method, tau_2 = NA) {
    # Determine the number of subgroups (k) and the total number of obs (n)
    k <- length(table(effect_sizes[, subgroup_moderator]))
    n <- nrow(effect_sizes)
    
    
    effect_sizes <- effect_sizes[order(effect_sizes[[trial_labels]]), ]
    
    # Input Error Handling (auxiliary_functions.R)
    validate_input(effect_sizes, heterogeneity_structure, estimation_method)
    
    # Step 1: Estimate the Within-Trial Interaction(s) (auxiliary_functions.R)
    # ---------------------------------------------
    step1_model <- 
      calculate_interaction_model(num_subgroups = k, effect_sizes, reference_group,
                                  subgroup_col = subgroup_moderator, labels_col = trial_labels,
                                  heterogeneity_structure,estimation_method, tau_2 = tau_2)
    
    # Step 2: Estimate the Floating Subgroup-Specific Treatment Effects
    # ------------------------------------------------------------------------------
    
    # Define the final subgroup as the reference subgroup
    final_reference <- ifelse(use_last_reference, k, reference_group)
    
    # Define a k by k matrix of ones
    J <- matrix(1, ncol = k, nrow = k)
    
    # Overwrite Z and the Contrast Matrix (M) based on the final subgroup
    Z <- as.matrix(diag(k)[, -final_reference])
    M <- diag(k)
    M[, final_reference] <- rep(-1, k)
    
    # Calculate heterogeneity estimates (tau^2 and rho) based on estimation method
    if (estimation_method == "REML") {
      if (k > 2) {
        # Calculate tau^2 for multiple subgroups
        tau2 <- c(diag(step1_model$G), NA)
        
        # Calculate rho for multiple subgroups
        rho <- c(mixmeta::vechMat(cov2cor(step1_model$G), diag = FALSE))
        rho <- c(rho, rep(NA, choose(k, 2) - length(rho)))
      } else {
        # Calculate tau^2 for two subgroups
        tau2 <- c(step1_model$tau2, NA)
        rho <- NA
      }
      
      if (heterogeneity_structure == "CS") {
        # Set missing rho values to the first tau^2 value for a compound symmetry structure
        rho[is.na(rho)] <- tau2[1]
      }
    } else {
      # Set tau^2 and rho to NA for Fixed-Effect (FE) method
      tau2 <- NA
      rho <- NA
    }
    
    # Permutation Matrix and Adjustments
    if (estimation_method == "REML") {
      # Define the permutation order
      x <- append(1:k, final_reference, after = k)[-final_reference]
      
      # Create permutation indices
      L <- diag(1:k)
      L[lower.tri(diag(1:k))] <- (k + 1):(k + choose(k, 2))
      L[upper.tri(diag(1:k))] <- t(L)[upper.tri(diag(1:k))]
      
      # Permute tau^2 and rho
      tau2 <- tau2[x]
      rho <- rho[L[x, x][lower.tri(L)] - k]
      
      # Create the Permutation Matrix P
      P <- diag(length(x))[, x]
      
      # Adjust the R matrix
      R <- rbind(cbind(diag(k - 1), -1), c(rep(0, k - 1), 1))
      R <- diag(1, n / k) %x% (P %*% R %*% t(P))
      
      # Initialize ek vector
      ek <- rep(x = 0, times = k)
      
      # Adjust the ek vector
      ek <- rep(append(ek, values = 1, after = final_reference)[-final_reference], n / k)
      
      
    } else {
      # Initialize R matrix for Fixed-Effects method
      R <- diag(1, n / k) %x% diag(k)
      
      # Initialize ek vector for Fixed-Effects method
      ek <- rep(1, n)
    }
    
    # Fit the Second-Stage Model
    step1_final <- 
      calculate_interaction_model(num_subgroups = k, effect_sizes, final_reference, subgroup_col = subgroup_moderator, labels_col = trial_labels, heterogeneity_structure, estimation_method, tau_2 = tau_2)
    
    # Fit the Third-Stage Model
    step2_model <-
      metafor::rma.mv(yi = R %*% c(effect_sizes$yi - rep(Z %*% c(step1_final$b), times = n / k)),
                      V = R %*% diag(effect_sizes$vi) %*% t(R),
                      method = estimation_method,
                      struct = "UN", tau2 = tau2, rho = rho,
                      intercept = FALSE, mods = ek,
                      random = as.formula(paste0("~ -1 + factor(", subgroup_moderator, ") | ", trial_labels)),
                      data = effect_sizes,
                      cvvc = ifelse(estimation_method == "REML", TRUE, FALSE)) %>%
      suppressWarnings()
    
    # Step 3: Correct the Variances of the Floating Subgroup-Specific Treatment Effects
    # ------------------------------------------------------------------------------
    
    # Conditional block to calculate Sigma_b based on the selected estimation method
    if (estimation_method == "REML") {
      R_inv <- MASS::ginv(R[1:k, 1:k])
      Sigma_b <- R_inv %*% step2_model$G %*% t(R_inv)
    } else {
      Sigma_b <- diag(k) * 0
    }
    
    # Calculate A matrix
    A <- A_comp(Sigma_b)$A
    
    # Post-process the results
    post_process_result <- 
      post_process(step1_model, step1_final, step2_model, effect_sizes, reference_group, subgroup_moderator, k, estimation_method, final_reference)
    
    return(post_process_result)
  }

# Function to Print Weighted Multivariate Meta-Analysis Results
print.wt <- function(wt, digits = 4, ...) {
  k <- which(rownames(wt$interaction) == "tau1^2")
  
  if (length(k) == 0) {
    k <- nrow(wt$interaction) + 1
  }
  
  interaction <- wt$interaction[!duplicated(rownames(wt$interaction)), ]
  subgroups <- wt$subgroups[!duplicated(rownames(wt$subgroups)), ]
  
  cat("\nInteraction estimates \n")
  cat("\n")
  
  if (k > 2) {
    print(interaction[1:(k - 1), ], digits = digits, ...)
    
    if (length(interaction[-c(1:(k - 1)), ]) > 0) {
      print(interaction[-c(1:(k - 1)), ], digits = digits, ...)
    }
  } else {
    print(interaction, digits = digits, ...)
  }
  
  cat("\n")
  cat("---\n")
  cat("\nSubgroup-specific estimates \n")
  cat("\n")
  print(subgroups[1:k, ], digits = digits, ...)
  cat("\n")
  
  if (length(subgroups[-c(1:k), ]) > 0) {
    print(subgroups[-c(1:k), ], digits = digits, ...)
  }
  
  cat("---\n")
}


# meta_analysis_auxiliary_functions.R

# Function to Calculate the Interaction Model
# ----------------------------------------------------------------------------
# This function calculates the interaction model for effect sizes using the specified
# heterogeneity structure and estimation method. It constructs the Z and M matrices,
# performs the Kronecker product, and fits the interaction model using metafor::rma.mv.
#
# Args:
#   - num_subgroups: Number of subgroups.
#   - effect_sizes: A data frame or escalc object containing effect sizes and variances.
#   - reference_group: The reference group of the interactions.
#   - subgroup_col: Name of the column containing subgroup information.
#   - labels_col: Name of the column containing trial identifiers.
#   - heterogeneity_structure: Structure of the heterogeneity matrix ("UN," "CS," or "FE").
#   - estimation_method: Estimation method ("FE" or "REML").
#
# Returns:
#   An object representing the fitted interaction model.
# ----------------------------------------------------------------------------

# Function to Calculate the Interaction Model
calculate_interaction_model <- 
  function(num_subgroups, effect_sizes, reference_group, subgroup_col,
           labels_col, heterogeneity_structure,
           estimation_method, tau_2 = NA) {
    # Create the Z matrix
    Z_matrix <- as.matrix(diag(num_subgroups)[, -reference_group])
    
    # Create the M matrix
    M_matrix <- diag(num_subgroups)
    M_matrix[, reference_group] <- rep(-1, num_subgroups)
    M_matrix[reference_group, ] <- rep(0, num_subgroups)
    
    # Calculate the Kronecker product of the M matrix
    kronecker_M <- diag(nrow(effect_sizes) / num_subgroups) %x% M_matrix
    
    # Reference-group sequence
    sequence <- seq(reference_group, nrow(effect_sizes), num_subgroups)
    
    # Calculate the effect sizes for the interaction model
    effect_sizes_interacted <- c(kronecker_M %*% matrix(effect_sizes$yi))
    effect_sizes_interacted <- effect_sizes_interacted[-sequence]
    
    # Calculate the V matrix (variance-covariance matrix)
    V_matrix <- (kronecker_M %*% diag(effect_sizes$vi) %*% t(kronecker_M))
    V_matrix <- V_matrix[-sequence, -sequence]
    
    # Define models for fixed and random effects for more than 2 subgroups
    if (num_subgroups > 2) {
      fixed_effects_model <-
        as.formula(paste0(" ~ -1 + factor(", subgroup_col, ")"))
      random_effects_model <-
        as.formula(paste0(" ~ -1 + factor(", subgroup_col, ") | factor(", labels_col, ")"))
    } else {
      fixed_effects_model <- ~1
      random_effects_model <-
        as.formula(paste0("~ 1 | factor(", labels_col, ")"))
    }
    
    # Fit the interaction model using metafor::rma.mv
    interaction_model <- metafor::rma.mv(
      yi = effect_sizes_interacted,
      V = V_matrix,
      mods = fixed_effects_model,
      random = random_effects_model,
      test = "z",
      struct = heterogeneity_structure,
      sigma = tau_2,
      rho = ifelse((heterogeneity_structure == "CS") & (num_subgroups - 1 > 1), 1 / 2, NA),
      method = estimation_method,
      data = effect_sizes[-sequence, ],
      cvvc = ifelse(estimation_method == "REML", TRUE, FALSE)
    ) %>%
      suppressWarnings()
    
    return(interaction_model)
  }


# Function to Validate and Prepare Input for the Interaction Model
# ----------------------------------------------------------------------------
# This function validates and prepares the input variables required for the interaction model.
# It checks the validity of effect size and variance data, the chosen heterogeneity structure,
# and the format of the input data. If any input is invalid, it raises an error.
#
# Args:
#   - effect_sizes: A data frame or escalc object containing effect sizes and variances.
#   - heterogeneity_structure: Structure of the heterogeneity matrix ("UN," "CS," or "FE").
#   - estimation_method: Estimation method ("FE" or "REML").
#
# Returns:
#   This function doesn't return a value but raises an error if any input is invalid.
# ----------------------------------------------------------------------------

validate_input <- function(effect_sizes, heterogeneity_structure, estimation_method) {
  # Validate and prepare the input variables
  
  # Ensure that the chosen structure is one of the allowed options
  heterogeneity_structure <<- match.arg(heterogeneity_structure, choices = c(NULL, "UN", "CS"))    
  
  # Ensure that the chosen estimation is one of the allowed options
  estimation_method <<- match.arg(estimation_method, choices = c("FE", "ML", "REML"))    
  
  
  # Check if effect size and variance data are provided and are valid
  if (any(is.na(effect_sizes$yi)) || any(is.na(effect_sizes$vi)) ||
      !all(effect_sizes$vi > 0)) {
    stop("Invalid effect size or variance data.")
  }
  
  # Check if the input data is in the correct format
  if (all(!(class(effect_sizes) %in% c("escalc", "data.frame")))) {
    stop("Effect size data should be in 'escalc' or 'data.frame' format.")
  }
  
  # Check if estimation_method is valid
  if (!estimation_method %in% c("FE", "ML", "REML")) {
    stop("Invalid estimation_method. Choose from 'FE' or 'REML'.")
  }
}

# Calculate Blocksum Matrix
# ----------------------------------------------------------------------------
# This function takes an input matrix 'M' and divides it into smaller block
# matrices based on specified row and column block sizes. It then computes
# the blocksum matrix by summing the individual block matrices.
#
# Args:
#   - M: Input matrix to be divided and summed.
#   - row_block_size: The number of rows in each block.
#   - col_block_size: The number of columns in each block.
#
# Returns:
#   A blocksum matrix resulting from the summation of block matrices.
# ----------------------------------------------------------------------------

blocksum_matrix <- function(M, row_block_size, col_block_size) {
  # Calculate the row group and column group indices based on the block sizes
  row_group <- (row(M) - 1) %/% row_block_size + 1
  col_group <- (col(M) - 1) %/% col_block_size + 1
  
  # Calculate the group index for each element of the original matrix
  group_index <- (row_group - 1) * max(col_group) + col_group
  
  # Calculate the total number of groups
  num_groups <- prod(dim(M)) / row_block_size / col_block_size
  
  # Split the matrix into individual group matrices
  group_matrices <- lapply(1:num_groups, function(group) {
    # Create a boolean mask for the current group
    group_mask <- group_index == group
    
    # Extract elements belonging to the current group and create a matrix
    group_matrix <- matrix(M[group_mask, , drop = FALSE], row_block_size, col_block_size)
  })
  
  # Sum up the group matrices to obtain the blocksum matrix
  blocksum_matrix <- Reduce(`+`, group_matrices)
  
  return(blocksum_matrix)
}

# Check and Prepare Variance-Covariance Matrix (vcc)
# ----------------------------------------------------------------------------
# This function checks and prepares the variance-covariance matrix 'vcc' used
# for estimating variances of subgroup-specific treatment effects. It ensures
# that the matrix is square and has no NA values. If needed, it replaces NA
# values with Inf and creates a square matrix. It also subsets 'vcc' based on
# the provided index matrix 'L'.
#
# Args:
#   - vcc: Variance-covariance matrix for subgroup-specific treatment effects.
#   - m: Number of subgroups.
#
# Returns:
#   A square variance-covariance matrix with any NA values replaced with Inf,
#   and subset based on the index matrix 'L'.
# ----------------------------------------------------------------------------

check_vcc <- function(vcc, m) {
  # Create an index matrix L for subsetting vcc
  L <- diag(1:m)
  L[lower.tri(diag(1:m))] <- (m + 1):(m + choose(m, 2))
  L[upper.tri(diag(1:m))] <- t(L)[upper.tri(diag(1:m))]
  
  # Check if vcc contains any NA values
  if (any(is.na(vcc))) {
    # Replace NA values with Inf and create a square matrix
    vcc <- diag(Inf, (m + choose(m, 2)), (m + choose(m, 2)))
    colnames(vcc) <- rownames(vcc) <- 1:(m + choose(m, 2))
  } else if (nrow(vcc) < (m + 1)) {
    # If vcc has fewer rows than expected, create a square matrix with Inf
    vcc <- diag(Inf, (m + choose(m, 2)), (m + choose(m, 2)))
    colnames(vcc) <- rownames(vcc) <- 1:(m + choose(m, 2))
  }
  
  # Return the subset of vcc matrix based on the index matrix L
  return(as.matrix(vcc[c(L), c(L)]))
}

# Function to compute confidence intervals for Sigma
confint_tau <- function(Sigma_b, Sigma_vb) {
  m <- ncol(Sigma_b)
  
  # Create an index matrix L for subsetting Sigma_b
  L <- diag(1:m)
  L[lower.tri(diag(1:m))] <- (m + 1):(m + choose(m, 2))
  L[upper.tri(diag(1:m))] <- t(L)[upper.tri(diag(1:m))]
  
  if (m > 1) {
    # Ensure the diagonals of Sigma_b contain correlations (convert to correlation matrix)
    Sigma_b[row(Sigma_b) != col(Sigma_b)] <- 
      cov2cor(Sigma_b)[row(Sigma_b) != col(Sigma_b)]
  }
  
  # Compute confidence intervals for Sigma
  Sigma_confint <- cbind(
    estimate = c(Sigma_b),
    ci.lb = c(Sigma_b + qnorm(0.025) * sqrt(diag(Sigma_vb))),
    ci.ub = c(Sigma_b + qnorm(0.975) * sqrt(diag(Sigma_vb)))
  )
  
  # Define row names for the confidence intervals
  rownames(Sigma_confint) <-   c(
    paste0("tau", 1:m, "^2"),
    paste0("rho", apply(combn(m, min(c(m, 2))), 2, paste, collapse = ""))
  )[c(L)]
  
  
  return(Sigma_confint)
}

# Function to compute A matrix based on the estimated heterogeneity matrix
A_comp <- function(Sigma) {
  
  # Retrieve variables from parent environment
  J <- get(x = "J", envir = parent.frame())
  Z <- get(x = "Z", envir = parent.frame())
  k <- get(x = "k", envir = parent.frame())
  n <- get(x = "n", envir = parent.frame())
  es <- get(x = "effect_sizes", envir = parent.frame())
  
  # Initialize variables to store weights
  wei_total <- 0
  wei <- matrix(nrow = n / k, ncol = k)
  
  # Calculate weights for each subgroup
  for (i in seq(1, n, k)) {
    
    # Calculate total weight for the subgroup
    wei_total <- wei_total + J[, 1] %*% 
      MASS::ginv(diag(es$vi)[i:(i + k - 1), i:(i + k - 1)] + Sigma) %*% J[, 1]
    
    # Store weights for each subgroup
    wei[(i + k - 1) / k, ] <- J[, 1] %*%
      MASS::ginv(diag(es$vi)[i:(i + k - 1), i:(i + k - 1)] + Sigma)
  }
  
  # Calculate the A matrix based on weights and Z
  A <- -as.matrix(J[, 1]) %*% colSums(wei / c(wei_total)) %*% Z + Z
  
  # Return A matrix and weights as a list
  return(list(A = A, weights = wei / c(wei_total)))
}

# post_process Function: Post-processes the results from step 1 and step 2 analyses
# ----------------------------------------------------------------------------
# This function performs post-processing on the results obtained from the step 1 
# and step 2 analyses. It computes heterogeneity confidence intervals and 
# combines the subgroup-specific estimates to provide a comprehensive output 
# containing both interaction and subgroup-specific information.
#
# Args:
#   - step1_model: Results from step 1 analysis.
#   - step2_model: Results from step 2 analysis.
#   - es: Effect size (usually an escalc object).
#   - reference: The reference group of the interactions.
#   - subgroup: Subgroup moderator.
#   - k: Number of subgroups.
#   - estimation_method: Estimation method ("FE" or "REML").
#   - reference_last: Whether the last subgroup is the reference group.
#
# Returns:
#   A list containing interaction and subgroup-specific estimates, their 
# confidence intervals, standard deviations, and hat values.
# ----------------------------------------------------------------------------

post_process <- function(step1_model, step1_final, step2_model, es, reference, 
                         subgroup, k, estimation_method, reference_last) {
  
  # Retrieve the J, A, Z, Sigma_b matrices from the parent environment
  J <- get(x = "J", envir = parent.frame())
  A <- get(x = "A", envir = parent.frame())
  Z <- get(x = "Z", envir = parent.frame())
  Sigma_b <- get(x = "Sigma_b", envir = parent.frame())
  
  if(estimation_method == "REML"){
    R_inv   <- get(x = "R_inv", envir = parent.frame())
  }
  
  # Calculate subgroup-specific estimates
  beta_b <- append(c(step2_model$b) + 
                     c(step1_final$b),
                   c(step2_model$b), 
                   after = reference_last - 1)
  
  beta_vb <- c(step2_model$vb) * J + A %*% 
    step1_final$vb %*% t(A)
  
  # Confidence intervals for subgroup-specific estimates
  beta_confint <- cbind(
    beta_b - sqrt(diag(beta_vb)) * qnorm(1 - 0.05 / 2),
    beta_b + sqrt(diag(beta_vb)) * qnorm(1 - 0.05 / 2)
  )
  
  if (estimation_method == "REML") {
    
    ## Define Subgroup Heterogeneity Estimates
    
    Sigma_vb <- MASS::ginv(R_inv %x% R_inv) %*% 
      check_vcc(vcc = step2_model$vvc, m = k)
    
    diag(Sigma_vb)[diag(Sigma_vb) < 0] <- 0
    Sigma_confint <- confint_tau(Sigma_b = Sigma_b, Sigma_vb = Sigma_vb)
    
    ## Define Interaction Heterogeneity Estimates
    if ((k - 1) == 1) {
      Sigma.gamma_b <- as.matrix(step1_model$sigma2)
    } else {
      Sigma.gamma_b <- step1_model$G
    }
    
    Sigma.gamma_vb <- check_vcc(step1_model$vvc, m = k - 1)
    Sigma.gamma_confint <-
      confint_tau(Sigma_b = Sigma.gamma_b, Sigma_vb = Sigma.gamma_vb)
  } else {
    Sigma_confint <- Sigma_vb <- Sigma.gamma_b <- Sigma.gamma_confint <-
      Sigma.gamma_vb <- NULL
  }
  # browser()
  output <- list(
    interaction = cbind(
      estimate = c((step1_model$b), Sigma.gamma_confint[, 1]),
      ci.lb = c((step1_model$ci.lb), Sigma.gamma_confint[, 2]),
      ci.ub = c((step1_model$ci.ub), Sigma.gamma_confint[, 3]),
      se = sqrt(c(diag(step1_model$vb), diag(Sigma.gamma_vb)))
    ),
    subgroups = cbind(
      estimate = c((beta_b), Sigma_confint[, 1]),
      ci.lb = c((beta_confint[, 1]), Sigma_confint[, 2]),
      ci.ub = c((beta_confint[, 2]), Sigma_confint[, 3]),
      se = sqrt(c(diag(beta_vb), diag(Sigma_vb)))
    ), 
    Sigma_subgroups = beta_vb,
    hatvalues = hatvalues(step2_model)
  )
  
  # browser()
  
  rownames(output$subgroups) <- c(unique(es[, subgroup]), rownames(Sigma_confint))
  rownames(output$interaction) <- c(unique(es[, subgroup])[-reference], rownames(Sigma.gamma_confint))
  
  class(output) <- "wt"
  return(output)
}


# Load necessary libraries
optimized_weights <- function(vi1, vi2, G, off_diagonal = TRUE) {
  
  # Define the determinant function
  determinant_func <- function(params, vi1, vi2, off_diagonal) {
    # Extract weights
    w1i <- params[1:length(vi1)]
    w2i <- params[(length(w1i) + 1):(2 * length(w1i))]
    wi <- params[(2 * length(w1i) + 1):(3 * length(w1i))]
    
    if (off_diagonal) {
      # Construct P and S matrices
      P <- matrix(c(c(rbind(w1i + wi, w2i)), c(rbind(w1i, w2i + wi))), nrow = 2, byrow = TRUE)
    } else {
      P <- matrix(c(c(rbind(wi, 0)), c(rbind(0, wi))), nrow = 2, byrow = TRUE)
    }
    
    S <- diag(c(rbind(vi1, vi2)))
    
    # Compute P S P'
    WSWT <- P %*% (S+ diag(1, length(w1i)) %x% G) %*% t(P)
    
    # Compute the determinant
    det_WSWT <- det(WSWT)
    
    return(det_WSWT)
  }
  
  # Define constraints function
  constraint_func <- function(params, off_diagonal) {
    # Extract weights
    w1i <- params[1:length(vi1)]
    w2i <- params[(length(w1i) + 1):(2 * length(w1i))]
    wi <- params[(2 * length(w1i) + 1):(3 * length(w1i))]
    
    # Constraints
    constraint1 <- sum(wi) - 1 # Sum of weights should be 1
    constraint2 <- sum(w1i + w2i) - 1 # Sum of w1i and w2i should be 1
    
    if(off_diagonal){
      return(c(constraint1, constraint2))      
    }
    else{
      return(constraint1)
    }
  }
  
  # Define the optimization problem
  optimize_determinant <- function(initial_values, vi1, vi2) {
    
    # Options for the auglag function
    opts <- list(
      # algorithm = "NLOPT_LD_AUGLAG",
      xtol_rel = 1e-10,  # Adjust tolerance
      xtol_abs = 1e-10,  # Adjust tolerance
      maxeval = 2000,   # Further increase the maximum number of iterations
      local_opts = list(
        algorithm = "NLOPT_LD_LBFGS",
        xtol_rel = 1e-10,
        xtol_abs = 1e-10
      )
    )
    
    # Optimization using auglag
    result <- auglag(
      x0 = initial_values,
      fn = function(params) determinant_func(params, vi1, vi2, off_diagonal),
      heq = function(params) constraint_func(params, off_diagonal),
      lower = c(rep(c(0, -1), each = N), rep(0, N)),
      upper = rep(c(1, 1, 1), each = N),
      localsolver = "LBFGS",
      control = opts
    )
    return(result)
  }
  
  # Example usage
  # Initial guess for weights
  N <- length(vi1)
  initial_weights <- rep(c(1 / (2 * N), 1 / (2 * N), 1 / N), each = N)
  
  # Perform optimization
  opt_result <- optimize_determinant(initial_values = initial_weights, vi1, vi2)
  
  w1i <- opt_result$par[1:length(vi1)]
  w2i <- opt_result$par[(length(w1i) + 1):(2 * length(w1i))]
  wi <- opt_result$par[(2 * length(w1i) + 1):(3 * length(w1i))]
  
  if (off_diagonal) {
    # Construct P and S matrices
    P <- matrix(c(c(rbind(w1i, w2i)), c(rbind(w1i - wi, w2i + wi))), nrow = 2, byrow = TRUE)
  } else {
    P <- matrix(c(c(rbind(wi, 0)), c(rbind(0, wi))), nrow = 2, byrow = TRUE)
  }
  
  S <- diag(c(rbind(vi1, vi2)))
  
  return(list(w1i = w1i, w2i = w2i, wi = wi, P = P, vb = P %*% (S+ diag(1, length(w1i)) %x% G) %*% t(P), opt_result = opt_result))
}
