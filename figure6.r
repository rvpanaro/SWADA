# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 6. Plot: Mismatch Between Interaction Estimates (Boxplots)
# =============================================================================
factor = 1.5 

## 6.1 With Aggregation Bias
df1 <- dat |>
  filter(
    method %in% c("Separate MAs CE (DA and AD)", "Separate MAs (DA and AD)"),
    a_interaction == "Aggregation bias",
    !is.na(interaction_tau)
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("Separate MAs CE (DA and AD)", "Separate MAs (DA and AD)"),
      labels = c("CE-MA", "RE-MA")
    )
  )

## 6.2 Without Aggregation Bias
df2 <- dat %>%
  filter(
    method %in% c("Separate MAs CE (DA and AD)", "Separate MAs (DA and AD)"),
    a_interaction == "No aggregation bias",
    !is.na(interaction_tau)
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Separate MAs CE (DA and AD)", "Separate MAs (DA and AD)"),
      labels = c("CE-MA", "RE-MA")
    )
  )

boxp2 <- df2 %>%
  ggplot(aes(x = factor(tau2), y = delta, fill = scenario)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_boxplot(linewidth = 1.1) +
  ylim(-0.5, 1.7) +
  labs(
    x     = expression("Interaction heterogeneity" ~ (tau[W])),
    y     = "",
    color = "Method",
    title = "Noncollapsibility without confounding \n (no aggregation bias)"
  ) +
  facet_grid(
    method ~ k, switch = "y",
    labeller = labeller(k = function(x) paste0("k = ", x))
  ) +
  defined_theme +
  theme(
    legend.position      = "none",
    legend.justification = "center",
    legend.box.just      = "center",
    legend.direction     = "horizontal",
    text                 = element_text(size = 18),
    strip.clip           = "on"
  )

# Small-panel versions for k = 20, tau = 0.1
df1_small <- df1 %>% filter(method == "CE-MA", k == 20, tau == 0.1)
df2_small <- df2 %>% filter(method == "CE-MA", k == 20, tau == 0.1)

boxp1_small <- ggplot(df1_small, aes(x = scenario, y = delta, fill = scenario)) +
  geom_boxplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  ylim(-0.5, 1.7) +
  labs(
    title = "Noncollapsibility with confounding \n (aggregation bias)",
    x     = "",
    y     = "Observed interaction mismatch"
  ) +
  defined_theme +
  theme(
    text          = element_text(size = 11 * factor),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text     = element_text(size = 9 * factor),
    axis.title    = element_text(size = 10 * factor),
    strip.text    = element_text(size = 10 * factor),
    plot.title    = element_text(size = 12 * factor),
    plot.subtitle = element_text(size = 10 * factor),
    plot.caption  = element_text(size = 8 * factor)
  )

boxp2_small <- ggplot(df2_small, aes(x = scenario, y = delta, fill = scenario)) +
  geom_boxplot() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  ylim(-0.5, 1.7) +
  labs(
    title = "Noncollapsibility without confounding \n (no aggregation bias)",
    x     = "Subgroup size variation",
    y     = "Observed interaction mismatch"
  ) +
  defined_theme +
  theme(
    text          = element_text(size = 11 * factor),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text     = element_text(size = 9 * factor),
    axis.title    = element_text(size = 10 * factor),
    strip.text    = element_text(size = 10 * factor),
    plot.title    = element_text(size = 12 * factor),
    plot.subtitle = element_text(size = 10 * factor),
    plot.caption  = element_text(size = 8 * factor)
  )

boxp_small <- wrap_elements(boxp2_small) + wrap_elements(boxp1_small)
ggsavepdf("img/viz_delta_small.png", plot = boxp_small, width = 14, height = 7)
