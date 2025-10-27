# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 14. Check & Print Coverage Table for Underperforming Methods
# =============================================================================

factor = 2.5 

summ_under <- summarised |>
  filter(
    (cov_interaction < 0.9364919) | (cov_subgroup1 < 0.9364919)
  ) |>
  select(method, tauW, cov_interaction, cov_subgroup1, cov_subgroup2) |>
  mutate(
    cov_interaction = cov_interaction * 100,
    cov_subgroup1   = cov_subgroup1   * 100,
    cov_subgroup2   = cov_subgroup2   * 100
  )

print(summ_under)

summ2 <- summarised |>
  select(method, cov_interaction, cov_subgroup1, cov_subgroup2) |>
  mutate(
    cov_interaction = cov_interaction * 100,
    cov_subgroup1   = cov_subgroup1   * 100,
    cov_subgroup2   = cov_subgroup2   * 100
  )

print(summ2)

# =============================================================================
# 10. Plot: Interaction Coverage Proportions
# =============================================================================
df1 <- summarised |> filter(a_interaction == "Aggregation bias")
df2 <- summarised |> filter(a_interaction != "Aggregation bias")

df1$undercover1 <- if_else(
  df1$method %in% unique(df1[df1$cov_interaction < 0.9364919, "method"]),
  "dashed", "solid"
)
df2$undercover2 <- if_else(
  df2$method %in% unique(df2[df2$cov_interaction < 0.9364919, "method"]),
  "dashed", "solid"
)

p1 <- df1 |>
  ggplot(aes(
    x        = as.factor(round(tauW, 2)),
    y        = 100 * cov_interaction,
    color    = method,
    shape    = method,
    group    = method
  )) +
  geom_abline(intercept = 95, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11)+
  coord_cartesian(ylim = c(92, 98)) +
  labs(
    x = expression("Interaction heterogeneity" ~ (tau[W])),
    y = "",
    title = "Noncollapsibility with confounding \n (aggregation bias)"
  ) +
  scale_y_continuous(breaks = seq(92, 98, 2), labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_shade() +
  defined_theme +
  theme(
    text            = element_text(size = 11 * factor),
    axis.text.x     = element_text(size = 9 * factor),
    axis.text.y     = element_text(size = 9 * factor),
    axis.title.x    = element_text(size = 10 * factor),
    axis.title.y    = element_text(size = 10 * factor),
    strip.text.x    = element_text(size = 10 * factor),
    strip.text.y    = element_text(size = 10 * factor),
    legend.justification = "center",
    legend.position = "bottom",
    plot.title      = element_text(size = 12 * factor),
    plot.subtitle   = element_text(size = 10 * factor),
    plot.caption    = element_text(size = 8 * factor)
  )

p2 <- df2 |>
  ggplot(aes(
    x        = as.factor(round(tauW, 2)),
    y        = 100 * cov_interaction,
    color    = method,
    shape    = method,
    group    = method
  )) +
  geom_abline(intercept = 95, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11)+
  coord_cartesian(ylim = c(92, 98)) +
  labs(
    x = expression("Interaction heterogeneity" ~ (tau[W])),
    y = "",
    title = "Noncollapsibility without confounding \n (no aggregation bias)"
  ) +
  scale_y_continuous(breaks = seq(92, 98, 2), labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_shade() +
  defined_theme +
  theme(
    text            = element_text(size = 11 * factor),
    axis.text.x     = element_text(size = 9 * factor),
    axis.text.y     = element_text(size = 9 * factor),
    axis.title.x    = element_text(size = 10 * factor),
    axis.title.y    = element_text(size = 10 * factor),
    strip.text.x    = element_text(size = 10 * factor),
    strip.text.y    = element_text(size = 10 * factor),
    legend.position = "none",
    plot.title      = element_text(size = 12 * factor),
    plot.subtitle   = element_text(size = 10 * factor),
    plot.caption    = element_text(size = 8 * factor)
  )

p <- wrap_elements(p2) / wrap_elements(p1)
ggsavepdf("img/viz_coverage_gamma.png", plot = p, width = 22, height = 22)

