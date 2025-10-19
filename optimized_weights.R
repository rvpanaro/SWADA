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