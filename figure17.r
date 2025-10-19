# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")


# =============================================================================
# 13. Plot: Relative SE Length (Subgroup 1)
# =============================================================================
df1 <- summarised %>% filter(a_interaction == "Aggregation bias")

p1 <- df1 %>%
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * se_subgroup1_rel_DA,
    color = method, shape = method, group = method
  )) +
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  coord_cartesian(ylim = c(100-50, 100+20)) +
  scale_shape_manual(values = 0:11)+
  labs(
    title = "Noncollapsibility with confounding \n (aggregation bias)",
    x = expression("Treatment heterogeneity" ~ (tau)),
    y = ""
  ) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_theme +
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.box.just      = "center",
    legend.direction     = "horizontal",
    text                 = element_text(size = 16),
    strip.clip           = "on"
  )

df2 <- summarised %>% filter(a_interaction != "Aggregation bias")

p2 <- df2 %>%
  ggplot(aes(
    x     = as.factor(round(tau, 2)),
    y     = 100 * se_subgroup1_rel_DA,
    color = method,  shape = method, group = method
  )) +
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  coord_cartesian(ylim = c(100-50, 100+20)) +
  scale_shape_manual(values = 0:11)+
  labs(
    title = "Noncollapsibility without confounding \n (no aggregation bias)",
    x = expression("Treatment heterogeneity" ~ (tau[W])),
    y = ""
  ) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  facet_grid(k ~ scenario, switch = "y",
             labeller = labeller(k = ~ paste0("k = ", .x))) +
  defined_theme +
  theme(
    legend.position      = "none",
    legend.justification = "center",
    legend.box.just      = "center",
    legend.direction     = "horizontal",
    text                 = element_text(size = 16),
    strip.clip           = "on"
  )

p <- wrap_elements(p2) / wrap_elements(p1)
ggsavepdf("img/viz_se_subgroup1.png", plot = p, width = 1400/72, height = 1400/72, dpi = 72)
