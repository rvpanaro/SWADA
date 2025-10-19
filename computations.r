# Define IPD meta-analytical model computations given an AgD dataset
fit_all <- function(x, prev_mean) {

  agd <- x$agd[[1]]

  agd <- agd |> 
    mutate(
      subgroup12 = as.numeric(subgroup) -0.5,
      prevalence = rep(ni[subgroup == 0] /(ni[subgroup == 0] + ni[subgroup == 1]), each = 2),
      prevalence2 = (abs(as.numeric(subgroup)- 1)  -prevalence)
    )
  
  # Define results tibbles
  interaction_results <- tribble(
    ~analysis, ~interaction_b, ~interaction_se, ~interaction_tau
  )
  
  subgroupB_results <- tribble(
    ~analysis, ~subgroupB_b, ~subgroupB_se, ~subgroupB_tau
  )
  
  subgroupA_results <- tribble(
    ~analysis, ~subgroupA_b, ~subgroupA_se, ~subgroupA_tau
  )
  
  covariance_results <- tribble(
    ~analysis, ~variance_B, ~variance_A, ~covariance_B_A
  )
  
  # Helper: add a row of NA values for a failed analysis
  include_na <- function(analysis) {
    interaction_results <<- interaction_results |>
      add_row(analysis = analysis, interaction_b = NA, interaction_se = NA, interaction_tau = NA)
    
    subgroupB_results <<- subgroupB_results |>
      add_row(analysis = analysis, subgroupB_b = NA, subgroupB_se = NA, subgroupB_tau = NA)
    
    subgroupA_results <<- subgroupA_results |>
      add_row(analysis = analysis, subgroupA_b = NA, subgroupA_se = NA, subgroupA_tau = NA)
    
    covariance_results <<- covariance_results |>
      add_row(analysis = analysis, variance_B = NA, variance_A = NA, covariance_B_A = NA)
  }
  
  ## "Prevalence-adjusted DA"
  fit4 <- try(
    rma.mv(
      yi = yi, V = vi,
      mods = ~ 1 + prevalence + factor(subgroup, levels = c(0, 1)),
      random = list( ~ -1 + prevalence2  | trial,  ~ 1 | trial), 
      struct = c("GEN", "GEN"),
      method = "REML",
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit4, "try-error")) {
    if(missing(prev_mean)){
      prev_mean <- mean(unique(agd$prevalence))      
    }
    
    # Contrast matrices
    if(length(fit4$b) == 2){
      C_interaction <- c(0, 1)
      C_subgroupB <- c(1, 1)
      C_subgroupA <- c(1, 0)
    } else{
      C_interaction <- c(0, 0, 1)
      C_subgroupB <- c(1, prev_mean, 1)
      C_subgroupA <- c(1, prev_mean, 0)
    }
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Prevalence-adjusted DA",
        interaction_b = C_interaction %*% fit4$b,
        interaction_se = sqrt(C_interaction %*% vcov(fit4) %*% C_interaction),
        interaction_tau = sqrt(fit4$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Prevalence-adjusted DA",
        subgroupB_b = C_subgroupB %*% fit4$b,
        subgroupB_se = sqrt(C_subgroupB %*% vcov(fit4) %*% C_subgroupB),
        subgroupB_tau = sqrt(fit4$sigma2 + prev_mean^2 * fit4$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Prevalence-adjusted DA",
        subgroupA_b = C_subgroupA %*% fit4$b,
        subgroupA_se = sqrt(C_subgroupA %*% vcov(fit4) %*% C_subgroupA),
        subgroupA_tau = sqrt(fit4$sigma2 + (1-prev_mean)^2 * fit4$tau2)
      )
    
    vb <- matrix(
      c(
        C_subgroupB,
        C_subgroupA
      ),
      byrow = TRUE, nrow = 2
    ) %*% fit4$vb %*% t(matrix(
      c(
        C_subgroupB,
        C_subgroupA
      ),
      byrow = TRUE, nrow = 2
    ))
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Prevalence-adjusted DA",
        variance_B = vb[1, 1],
        variance_A = vb[2, 2],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Prevalence-adjusted DA")
  }
  
  ## "Within-trial framework" (WT)
  fit3 <- try(
    wt_framework(
      reference_group = 2,
      effect_sizes = agd,
      trial_labels = "trial",
      subgroup_moderator = "subgroup",
      estimation_method = "REML",
      heterogeneity_structure = "CS"
    ),
    silent = TRUE
  )
  
  if (!inherits(fit3, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Within-trial framework",
        interaction_b = fit3$interaction[1, "estimate"],
        interaction_se = fit3$interaction[1, "se"],
        interaction_tau = sqrt(fit3$interaction[2, "estimate"])
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Within-trial framework",
        subgroupB_b = fit3$subgroups[1, "estimate"],
        subgroupB_se = fit3$subgroups[1, "se"],
        subgroupB_tau = sqrt(fit3$subgroups[3, "estimate"])
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Within-trial framework",
        subgroupA_b = fit3$subgroups[2, "estimate"],
        subgroupA_se = fit3$subgroups[2, "se"],
        subgroupA_tau = sqrt(fit3$subgroups[6, "estimate"])
      )
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Within-trial framework",
        variance_B = fit3$Sigma_subgroup[1, 1],
        variance_A = fit3$Sigma_subgroup[2, 2],
        covariance_B_A = fit3$Sigma_subgroup[1, 2]
      )
    
  } else {
    include_na(analysis = "Within-trial framework")
  }
  
  ## "Difference of Averages" (DA)
  fit5 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      method = "REML",
      struct = "UN",
      mods = ~ 0 + factor(subgroup, levels = c(0, 1)),
      random = ~ 0 + subgroup | trial,
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit5, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "DA",
        interaction_b = predict(fit5, newmods = c(-1, 1))$pred,
        interaction_se = predict(fit5, newmods = c(-1, 1))$se,
        interaction_tau = sqrt(c(1, -1) %*% fit5$G %*% c(1, -1))
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "DA",
        subgroupB_b = predict(fit5, newmods = c(0, 1))$pred,
        subgroupB_se = predict(fit5, newmods = c(0, 1))$se,
        subgroupB_tau = sqrt(c(0, 1) %*% fit5$G %*% c(0, 1))
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "DA",
        subgroupA_b = predict(fit5, newmods = c(1, 0))$pred,
        subgroupA_se = predict(fit5, newmods = c(1, 0))$se,
        subgroupA_tau = sqrt(c(1, 0) %*% fit5$G %*% c(1, 0))
      )
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "DA",
        variance_B = fit5$vb[2, 2],
        variance_A = fit5$vb[1, 1],
        covariance_B_A = fit5$vb[1, 2]
      )
    
  } else {
    include_na("DA")
  }
  
  ## "Separate AgD-MAs"
  uni <- try(metafor::rma.uni(
    yi = yi[subgroup == 1] - yi[subgroup == 0],
    vi = vi[subgroup == 1] + vi[subgroup == 0],
    method = "REML",
    data = agd
  ), silent = TRUE)
  
  uni1 <- try(metafor::rma.uni(
    yi = yi,
    vi = vi,
    method = "REML",
    data = agd[agd$subgroup == 0, ]
  ), silent = TRUE)
  
  uni2 <- try(metafor::rma.uni(
    yi = yi,
    vi = vi,
    method = "REML",
    data = agd[agd$subgroup == 1, ]
  ), silent = TRUE)
  
  if (!inherits(uni, "try-error") &&
      !inherits(uni1, "try-error") &&
      !inherits(uni2, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Separate MAs (DA and AD)",
        interaction_b = uni$b,
        interaction_se = uni$se,
        interaction_tau = sqrt(uni$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Separate MAs (DA and AD)",
        subgroupB_b = uni2$b,
        subgroupB_se = uni2$se,
        subgroupB_tau = sqrt(uni2$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Separate MAs (DA and AD)",
        subgroupA_b = uni1$b,
        subgroupA_se = uni1$se,
        subgroupA_tau = sqrt(uni1$tau2)
      )
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Separate MAs (DA and AD)",
        variance_B = uni2$se^2,
        variance_A = uni1$se^2,
        covariance_B_A = 0
      )
    
  } else {
    include_na("Separate MAs (DA and AD)")
  }
  
  ## "Separate AgD-MAs CE"
  uni_CE <- try(metafor::rma.uni(
    yi = yi[subgroup == 1] - yi[subgroup == 0],
    vi = vi[subgroup == 1] + vi[subgroup == 0],
    method = "CE",
    data = agd
  ), silent = TRUE)
  
  uni1_CE <- try(metafor::rma.uni(
    yi = yi,
    vi = vi,
    method = "CE",
    data = agd[agd$subgroup == 0, ]
  ), silent = TRUE)
  
  uni2_CE <- try(metafor::rma.uni(
    yi = yi,
    vi = vi,
    method = "CE",
    data = agd[agd$subgroup == 1, ]
  ), silent = TRUE)
  
  if (!inherits(uni, "try-error") &&
      !inherits(uni1, "try-error") &&
      !inherits(uni2, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Separate MAs CE (DA and AD)",
        interaction_b = uni_CE$b,
        interaction_se = uni_CE$se,
        interaction_tau = sqrt(uni_CE$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Separate MAs CE (DA and AD)",
        subgroupB_b = uni2_CE$b,
        subgroupB_se = uni2_CE$se,
        subgroupB_tau = sqrt(uni2_CE$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Separate MAs CE (DA and AD)",
        subgroupA_b = uni1_CE$b,
        subgroupA_se = uni1_CE$se,
        subgroupA_tau = sqrt(uni1_CE$tau2)
      )
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Separate MAs CE (DA and AD)",
        variance_B = uni2_CE$se^2,
        variance_A = uni1_CE$se^2,
        covariance_B_A = 0
      )
    
  } else {
    include_na("Separate MAs CE (DA and AD)")
  }
  
  ## "Equal weights SWADA"
  wi <- try(rep(2 / length(agd$vi), times = length(agd$vi)), silent = TRUE)
  
  fit7 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      W = diag(wi),
      method = "REML",
      mods = ~ 1 + subgroup12,
      random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("ID", "CS"),  # first term shared, second independent  
      tau2 = uni$tau2/2,
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit7, "try-error") && !inherits(wi, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Equal weights SWADA",
        interaction_b = c(0, 1) %*% fit7$b,
        interaction_se = sqrt(c(0, 1) %*% vcov(fit7) %*% c(0, 1)),
        interaction_tau = sqrt(uni$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Equal weights SWADA",
        subgroupB_b = c(1, 0.5) %*% fit7$b,
        subgroupB_se = sqrt(c(1, 0.5) %*% vcov(fit7) %*% c(1, 0.5)),
        subgroupB_tau = sqrt(fit7$sigma2 + 0.25* uni$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Equal weights SWADA",
        subgroupA_b = c(1, -0.5) %*% fit7$b,
        subgroupA_se = sqrt(c(1, -0.5) %*% vcov(fit7) %*% c(1, -0.5)),
        subgroupA_tau = sqrt(fit7$sigma2 + 0.25* uni$tau2)
      )
    
    vb <- (rbind(c(1, -0.5), c(1, 0.5)) %*% fit7$vb %*% cbind(c(1, -0.5), c(1, 0.5)))
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Equal weights SWADA",
        variance_B = vb[2, 2],
        variance_A = vb[1, 1],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Equal weights SWADA")
  }
  
  ## "Interaction weights SWADA"
  wi <- try(diag(weights(uni, type = "matrix")), silent = TRUE)
  fit8 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      W = diag(rep(wi, each = 2)),
      method = "REML",
      mods = ~ 1 + subgroup12,
      random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("CS", "CS"),  # first term shared, second independent
      tau2 = uni$tau2/4,
      rho  = -1, 
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit8, "try-error") && !inherits(wi, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Interaction weights SWADA",
        interaction_b = c(0, 1) %*% fit8$b,
        interaction_se = sqrt(c(0, 1) %*% vcov(fit8) %*% c(0, 1)),
        interaction_tau = sqrt(uni$tau2)
      )

    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Interaction weights SWADA",
        subgroupB_b = c(1, 0.5)  %*% fit8$b,
        subgroupB_se = sqrt(c(1, 0.5)  %*% fit8$vb  %*% c(1, 0.5)),
        subgroupB_tau = sqrt(fit8$sigma2 + 0.25* uni$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Interaction weights SWADA",
        subgroupA_b = c(1, -0.5)  %*% fit8$b,
        subgroupA_se = sqrt(c(1, -0.5)  %*% fit8$vb  %*% c(1, -0.5)),
        subgroupA_tau = sqrt(fit8$sigma2 + 0.25* uni$tau2)
      )
    
    vb <- (rbind(c(1, -0.5), c(1, 0.5)) %*% fit8$vb %*% cbind(c(1, -0.5), c(1, 0.5)))
      
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Interaction weights SWADA",
        variance_B = vb[2, 2],
        variance_A = vb[1, 1],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Interaction weights SWADA")
  }

  ## "Study size weights SWADA"
  wi <- try(agd$ni[agd$subgroup == 0] + agd$ni[agd$subgroup == 1], silent = TRUE)
  
  fit9 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      W = diag(rep(wi, each = 2)),
      method = "REML",
      mods = ~ 1 + subgroup12,
      random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("CS", "CS"),  # first term shared, second independent
      tau2 = uni$tau2/4,
      rho  = -1, 
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit9, "try-error") && !inherits(wi, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Study size weights SWADA",
        interaction_b = c(0, 1) %*% fit9$b,
        interaction_se = sqrt(c(0, 1) %*% vcov(fit9) %*% c(0, 1)),
        interaction_tau = sqrt(uni$tau2)
      )

    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Study size weights SWADA",
        subgroupB_b =  c(1, 0.5) %*% fit9$b,
        subgroupB_se = sqrt(c(1, 0.5) %*% fit9$vb %*% c(1, 0.5)),
        subgroupB_tau = sqrt(fit9$sigma2 + 0.25* uni$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Study size weights SWADA",
        subgroupA_b =  c(1, -0.5) %*% fit9$b,
        subgroupA_se = sqrt(c(1, -0.5) %*% fit9$vb %*% c(1, -0.5)),
        subgroupA_tau = sqrt(fit9$sigma2 + 0.25* uni$tau2)
      )
    
    vb <- (rbind(c(1, -0.5), c(1, 0.5)) %*% fit9$vb %*% cbind(c(1, -0.5), c(1, 0.5)))
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Study size weights SWADA",
        variance_B = vb[2, 2],
        variance_A = vb[1, 1],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Study size weights SWADA")
  }
  
  ## "Smaller subgroup weights SWADA"
  wi <- try(
    apply(cbind(
      (agd$ni[agd$subgroup == 0]),
      (agd$ni[agd$subgroup == 1])
    ), 1, min),
    silent = TRUE
  )
  
  fit10 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      W = diag(rep(wi, each = 2)),
      method = "REML",
      mods = ~ 1 + subgroup12,
      random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("CS", "CS"),  # first term shared, second independent
      tau2 = uni$tau2/4,
      rho  = -1, 
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit10, "try-error") && !inherits(wi, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Smaller subgroup weights SWADA",
        interaction_b = c(0, 1) %*% fit10$b,
        interaction_se = sqrt(c(0, 1) %*% vcov(fit10) %*% c(0, 1)),
        interaction_tau = sqrt(uni$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Smaller subgroup weights SWADA",
        subgroupB_b = c(1, 0.5) %*% fit10$b,
        subgroupB_se = sqrt(c(1, 0.5) %*% fit10$vb %*% c(1, 0.5)),
        subgroupB_tau = sqrt(fit10$sigma2 + 0.25* uni$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Smaller subgroup weights SWADA",
        subgroupA_b = c(1, -0.5) %*% fit10$b,
        subgroupA_se = sqrt(c(1, -0.5) %*% fit10$vb %*% c(1, -0.5)),
        subgroupA_tau = sqrt(fit10$sigma2 + 0.25* uni$tau2)
      )
    
    vb <- (rbind(c(1, -0.5), c(1, 0.5)) %*% fit10$vb %*% cbind(c(1, -0.5), c(1, 0.5)))
  
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Smaller subgroup weights SWADA",
        variance_B = vb[2, 2],
        variance_A = vb[1, 1],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Smaller subgroup weights SWADA")
  }
  
  ## "Minimum of three weights SWADA"
  wi <- try(
    apply(cbind(
      diag(weights(uni, type = "matrix")),
      diag(weights(uni1, type = "matrix")),
      diag(weights(uni2, type = "matrix"))
    ), 1, min),
    silent = TRUE
  )
  
  fit11 <- try(
    rma.mv(
      yi = yi,
      V = diag(vi),
      W = diag(rep(wi, each = 2)),
      method = "REML",
      mods = ~ 1 + subgroup12,
      random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("CS", "CS"),  # first term shared, second independent
      tau2 = uni$tau2/4,
      rho  = -1, 
      data = agd
    ),
    silent = TRUE
  )
  
  if (!inherits(fit11, "try-error") && !inherits(wi, "try-error")) {
    
    interaction_results <- interaction_results |>
      add_row(
        analysis = "Minimum of three weights SWADA",
        interaction_b = c(0, 1) %*% fit11$b,
        interaction_se = sqrt(c(0, 1) %*% vcov(fit11) %*% c(0, 1)),
        interaction_tau = sqrt(uni$tau2)
      )
    
    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Minimum of three weights SWADA",
        subgroupB_b = c(1, 0.5) %*% fit11$b,
        subgroupB_se = sqrt(c(1, 0.5) %*% fit11$vb %*% c(1, 0.5)),
        subgroupB_tau = sqrt(fit11$sigma2 + 0.25* uni$tau2)
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Minimum of three weights SWADA",
        subgroupA_b = c(1, -0.5) %*% fit11$b,
        subgroupA_se = sqrt(c(1, -0.5) %*% fit11$vb %*% c(1, -0.5)),
        subgroupA_tau = sqrt(fit11$sigma2 + 0.25* uni$tau2)
      )
    
    vb <- (rbind(c(1, -0.5), c(1, 0.5)) %*% fit11$vb %*% cbind(c(1, -0.5), c(1, 0.5)))
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Minimum of three weights SWADA",
        variance_B = vb[2, 2],
        variance_A = vb[1, 1],
        covariance_B_A = vb[1, 2]
      )
    
  } else {
    include_na("Minimum of three weights SWADA")
  }
  
  opt3 <- try(
    optimized_weights(
      vi1 = agd$vi[agd$subgroup == 1],
      vi2 = agd$vi[agd$subgroup == 0],
      G = fit8$G,
      off_diagonal = FALSE
    ),
    silent = TRUE
  )
  
  if (!inherits(opt3, "try-error")) {
    P <- opt3$P
    vb <- opt3$vb
    b_hat <- P %*% agd$yi

    interaction_results <- interaction_results |>
      add_row(
        analysis = "Minimum of total variance",
        interaction_b = c(1, -1) %*% b_hat,
        interaction_se = sqrt(c(1, -1) %*% vb %*% c(1, -1)),
        interaction_tau = sqrt(c(1, -1) %*% vb %*% c(1, -1))
      )

    subgroupB_results <- subgroupB_results |>
      add_row(
        analysis = "Minimum of total variance",
        subgroupB_b = b_hat[1, ],
        subgroupB_se = sqrt(vb[1, 1]),
        subgroupB_tau = sqrt(c(1, 0) %*% fit8$G %*% c(1, 0))
      )
    
    subgroupA_results <- subgroupA_results |>
      add_row(
        analysis = "Minimum of total variance",
        subgroupA_b = b_hat[2, ],
        subgroupA_se = sqrt(vb[2, 2]),
        subgroupA_tau = sqrt(c(0, 1) %*% fit8$G %*% c(0, 1))
      )
    
    covariance_results <- covariance_results |>
      add_row(
        analysis = "Minimum of total variance",
        variance_B = vb[1, 1],
        variance_A = vb[2, 2],
        covariance_B_A = vb[1, 2]
      )
  } else {
    include_na("Minimum of total variance")
  }
  
  # compute confidence intervals
  interaction_results <- interaction_results |> 
    mutate(interaction_ci.lb = interaction_b - qnorm(.975) * interaction_se) |>
    mutate(interaction_ci.ub = interaction_b + qnorm(.975) * interaction_se) 

  subgroupB_results <- subgroupB_results |> 
    mutate(subgroupB_ci.lb = subgroupB_b - qnorm(.975) * subgroupB_se) |>
    mutate(subgroupB_ci.ub = subgroupB_b + qnorm(.975) * subgroupB_se) 
  
  subgroupA_results <- subgroupA_results |> 
    mutate(subgroupA_ci.lb = subgroupA_b - qnorm(.975) * subgroupA_se) |>
    mutate(subgroupA_ci.ub = subgroupA_b + qnorm(.975) * subgroupA_se) 
  
  
  covariance_results <- covariance_results |> 
    mutate(correlation_B_A = cov2cor(matrix(c(variance_B, covariance_B_A, 
                                            variance_A, covariance_B_A), 
                                            nrow =  2,
                                            ncol =  2, byrow = TRUE))[1, 2]) 
  
  return(c(interaction_results, subgroupB_results, subgroupA_results, covariance_results))
}
