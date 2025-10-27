# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 12. Plot: Relative SE Length (Interaction)
# =============================================================================
factor <- 2.5 

df1 <- summarised |>
  group_by(method) |>
  filter(a_interaction == "Aggregation bias")

p1 <- df1 |>
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * se_interaction_rel_AD,
    color = method, shape = method, group = method
  )) +
  ylim(-110,115) + 
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  labs(
    x = expression("Interaction heterogeneity" ~ (tau[W])),
    y = "",
    title = "Noncollapsibility with confounding \n (aggregation bias)"
  ) +
  scale_shape_manual(values = 0:11) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
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

df2 <- summarised |>
  group_by(method) |>
  filter(a_interaction != "Aggregation bias")

p2 <- df2 |>
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * se_interaction_rel_AD,
    color = method, shape = method, group = method
  )) +
  ylim(-110,115) + 
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11)+
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  labs(
    x = expression("Interaction heterogeneity" ~ (tau[W])),
    y = "",
    title = "Noncollapsibility without confounding \n (no aggregation bias)"
  ) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
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

p <- wrap_elements(p2) / wrap_elements(p1) + plot_layout(guides = "collect") 

ggsavepdf("img/viz_se_gamma.png", plot = p, width = 22, height = 22)