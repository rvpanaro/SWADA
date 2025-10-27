# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 15. Plot: Treatment Effect Coverage (Subgroup 1)
# =============================================================================
factor = 2.5 

p1 <- summarised %>%
  filter(a_interaction == "Aggregation bias") %>%
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * cov_subgroup1,
    color = method, shape = method, group = method
  )) +
  geom_abline(intercept = 95, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11) +
  scale_y_continuous(breaks = c(seq(30, 90, 10)), labels = ~ paste0(.x, "%")) +
  labs(
    x = expression("Treatment heterogeneity" ~ (tau)),
    y = "",
    title = "Noncollapsibility with confounding \n (aggregation bias)"
  ) +
  # scale_y_continuous(breaks = seq(90, 100, 2), labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_theme + 
  defined_shade() +
  theme(
    text            = element_text(size = 11 * factor),
    axis.text.x     = element_text(size = 9 * factor),
    axis.text.y     = element_text(size = 9 * factor),
    axis.title.x    = element_text(size = 10 * factor),
    axis.title.y    = element_text(size = 10 * factor),
    strip.text.x    = element_text(size = 10 * factor),
    strip.text.y    = element_text(size = 10 * factor),
    legend.position = "bottom",
    legend.justification = "bottom",
    plot.title      = element_text(size = 12 * factor),
    plot.subtitle   = element_text(size = 10 * factor),
    plot.caption    = element_text(size = 8 * factor)
  ) 

p2 <- summarised %>%
  filter(a_interaction != "Aggregation bias") %>%
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * cov_subgroup1,
    color = method, shape = method, group = method
  )) +
  geom_abline(intercept = 95, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11)+
  scale_y_continuous(breaks = seq(92, 98, 2), labels = ~ paste0(.x, "%")) +
  labs(
    x = expression("Treatment heterogeneity" ~ (tau)),
    y = "",
    title = "Noncollapsibility without confounding \n (no aggregation bias)"
  ) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_theme + 
  defined_shade() +
  theme(
    text            = element_text(size = 11 * factor),
    axis.text.x     = element_text(size = 9 * factor),
    axis.text.y     = element_text(size = 9 * factor),
    axis.title.x    = element_text(size = 10 * factor),
    axis.title.y    = element_text(size = 10 * factor),
    strip.text.x    = element_text(size = 10 * factor),
    strip.text.y    = element_text(size = 10 * factor),
    legend.position = "bottom",
    legend.justification = "center",
    plot.title      = element_text(size = 12 * factor),
    plot.subtitle   = element_text(size = 10 * factor),
    plot.caption    = element_text(size = 8 * factor)
  )

p <- wrap_elements(p2) / wrap_elements(p1)
ggsavepdf("img/viz_coverage_beta1.png", plot = p, width = 22, height = 22)
