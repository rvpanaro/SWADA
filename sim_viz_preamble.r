
# 1. Set Working Directory
# (no need) setwd("/run/user/1000/gvfs/sftp:host=hpc-gwdg/home/uni08/panaro/Documents/interaction-MA/")

# 2. Load Required Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(dplyr)

# 3. Read Simulation Data & Assign Factors
dat <- readRDS("result_sim_2025-10-18 21_10_12.178239.rds")

dat$size <- factor(
  dat$size,
  levels = unique(dat$size),
  labels = rep("All similar", length(unique(dat$size)))
)

dat$rho <- factor(
  dat$a,
  levels = c(0.5),
  labels = c(expression(rho == 0.5))
)

dat$a_interaction <- factor(
  dat$a_interaction,
  levels = c(2, -1),
  labels = c("Aggregation bias", "No aggregation bias")
)

# 4. Consistency checks --------------------------------------------------------
# Check for flipped confidence intervals
# (The lower bound should never be greater than the upper bound.)
dat %>%
  summarise(
    bad_interaction = sum(interaction_lb > interaction_ub, na.rm = TRUE),
    bad_subgroupA   = sum(subgroupA_lb > subgroupA_ub, na.rm = TRUE),
    bad_subgroupB   = sum(subgroupB_lb > subgroupB_ub, na.rm = TRUE)
  )
# EXPECTATION:
# → All counts should be 0.
# If >0, it means some CIs were constructed incorrectly (bounds reversed).

# Check if the CI width matches ±1.96 × SE (within rounding error)
dat %>%
  summarise(
    int_width_diff  = mean(abs((interaction_ub - interaction_lb)/2 - interaction_width), na.rm = TRUE),
    subA_width_diff = mean(abs((subgroupA_ub - subgroupA_lb)/2 - subgroupA_width), na.rm = TRUE),
    subB_width_diff = mean(abs((subgroupB_ub - subgroupB_lb)/2 - subgroupB_width), na.rm = TRUE)
  )

dat <- dat |>
  mutate(
    subgroupA_width = 1.959964 * subgroupA_se,
    subgroupA_lb = subgroupA_b - subgroupA_width,
    subgroupA_ub = subgroupA_b + subgroupA_width,
    subgroupA_contained = (subgroupA_lb < target_A) & (subgroupA_ub > target_A),
    
    subgroupB_width = 1.959964 * subgroupB_se,
    subgroupB_lb = subgroupB_b - subgroupB_width,
    subgroupB_ub = subgroupB_b + subgroupB_width,
    subgroupB_contained = (subgroupB_lb < target_B) & (subgroupB_ub > target_B)
  )  

# Check if the CI width matches ±1.96 × SE (within rounding error)
dat %>%
  summarise(
    int_width_diff  = mean(abs((interaction_ub - interaction_lb)/2 - interaction_width), na.rm = TRUE),
    subA_width_diff = mean(abs((subgroupA_ub - subgroupA_lb)/2 - subgroupA_width), na.rm = TRUE),
    subB_width_diff = mean(abs((subgroupB_ub - subgroupB_lb)/2 - subgroupB_width), na.rm = TRUE)
  )

# Check correlation between SE and CI width
# Width should be proportional to SE, so correlation ≈ 1.
dat %>%
  summarise(
    cor_int  = cor(interaction_width, interaction_se, use = "complete.obs"),
    cor_subA = cor(subgroupA_width, subgroupA_se, use = "complete.obs"),
    cor_subB = cor(subgroupB_width, subgroupB_se, use = "complete.obs")
  )

# Check coverage — proportion of times the true value was contained in CI
dat %>%
  summarise(
    coverage_interaction = mean(interaction_contained, na.rm = TRUE),
    coverage_subA = mean(subgroupA_contained, na.rm = TRUE),
    coverage_subB = mean(subgroupB_contained, na.rm = TRUE)
  )

# Sanity check for impossible or extreme values
dat %>%
  summarise(
    neg_se = sum(interaction_se < 0 | subgroupA_se < 0 | subgroupB_se < 0, na.rm = TRUE),
    huge_width = sum(interaction_width > 10 | subgroupA_width > 10 | subgroupB_width > 10, na.rm = TRUE)
  )

# 5. Calculate Reference SEs & Relative Differences
dat <- dat |>
  group_by(job_id) |>
  mutate(
    se_AD = interaction_se[method == "Separate MAs (DA and AD)"],
    se_DA = subgroupA_se[method == "DA"]
  ) |>
  ungroup() |>
  mutate(
    se_interaction_rel_AD = round(interaction_se / se_AD, 2),
    se_subgroupA_rel_DA   = round(subgroupA_se / se_DA,   2)
  )

# 6. Summarise Simulation Results
summarised <- dat |>
  group_by(experiment, method) |>
  summarise(
    experiment            = first(experiment),
    k                     = mean(k, na.rm = TRUE), 
    tau                   = mean(tau, na.rm = TRUE), 
    tauW                   = mean(tauW, na.rm = TRUE), 
    a_interaction         = first(a_interaction), 
    scenario              = first(scenario), 
    convergence           = mean(is.na(interaction_contained)),
    cov_interaction       = mean(interaction_contained, na.rm = TRUE),
    target_subgroup1      = mean(target_A, na.rm = TRUE),
    mean_subgroup1        = mean(subgroupA_b, na.rm = TRUE),
    target_subgroup2      = mean(target_B, na.rm = TRUE),
    mean_subgroup2        = mean(subgroupB_b, na.rm = TRUE),
    cov_subgroup1         = mean(subgroupA_contained,   na.rm = TRUE),
    cov_subgroup2         = mean(subgroupB_contained,   na.rm = TRUE),
    rmse_interaction      = mean(interaction_bias^2,    na.rm = TRUE),
    rmse_subgroup1        = mean(subgroupA_bias^2,      na.rm = TRUE),
    rmse_subgroup2        = mean(subgroupB_bias^2,      na.rm = TRUE),
    tau_bias_interaction  = mean(interaction_tau_bias,   na.rm = TRUE),
    tau_bias_subgroup1    = mean(subgroupA_tau_bias,     na.rm = TRUE),
    tau_bias_subgroup2    = mean(subgroupB_tau_bias,     na.rm = TRUE),
    se_interaction        = mean(interaction_se,         na.rm = TRUE),
    se_subgroup1          = mean(subgroupA_se,           na.rm = TRUE),
    se_subgroup2          = mean(subgroupB_se,           na.rm = TRUE),
    se_interaction_rel_AD = median(se_interaction_rel_AD,   na.rm = TRUE),
    se_subgroup1_rel_DA   = median(se_subgroupA_rel_DA,     na.rm = TRUE)
  ) |>
  mutate(
    type = if_else(
      method %in% c(
        "Collapsible estimator",
        "Prevalence-adjusted DA (location-scale)",
        "Prevalence-adjusted DA",
        "Separate MAs (DA and AD)",
        "Separate MAs CE (DA and AD)",
        "Within-trial framework",
        "DA"
      ),
      "Model-based methods", "Weigthing methods"
    )
  ) |>
  ungroup()

# 4. Define ggplot2 Theme & Shading Function
defined_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title           = element_text(size = 20),
    panel.background     = element_rect(fill = "white"),
    legend.position      = c(.99, .99),
    legend.text          = element_text(size = 18),
    legend.justification = c("right", "top"),
    legend.box.just      = "right",
    legend.margin        = margin(6, 6, 6, 6),
    legend.key.size      = unit(1, "lines"),
    legend.title         = element_blank(),
    axis.title.x         = element_text(size = 20),
    axis.title.y         = element_text(size = 20),
    axis.text            = element_text(size = 18),
    strip.text.x         = element_text(size = 18),
    strip.background     = element_blank(),
    strip.placement      = "outside",
    strip.text.y.left    = element_text(size = 14, angle = 90)
  )

defined_shade <- function(intrcpt = 0.95, replications = 1000) {
  annotate(
    geom = "rect",
    ymin = 100 * (intrcpt - qnorm(.025) * sqrt(intrcpt * (1 - intrcpt) / replications)),
    ymax = 100 * (intrcpt + qnorm(.025) * sqrt(intrcpt * (1 - intrcpt) / replications)),
    xmin = -Inf, xmax = Inf,
    alpha = 0.10
  )
}

