# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 13. Plot: Mismatch Between Interaction Estimates (Boxplots)
# =============================================================================

factor = 2.5

## 13.1 With Aggregation Bias
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

boxp1 <- df1 %>%
  ggplot(aes(x = factor(tauW), y = delta, fill = scenario)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_boxplot(linewidth = 1.1) +
  ylim(-0.5, 1.7) +
  labs(
    x     = expression("Interaction heterogeneity" ~ (tau[W])),
    y     = "",
    title = "Noncollapsibility with confounding \n (aggregation bias)",
    color = "Method"
  ) +
  facet_grid(
    method ~ k, switch = "y",
    labeller = labeller(k = function(x) paste0("k = ", x))
  ) +
  defined_theme +
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

## 13.2 Without Aggregation Bias
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
  ggplot(aes(x = factor(tauW), y = delta, fill = scenario)) +
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

## 8.3 Combine & Save
boxp     <- wrap_elements(boxp2) / wrap_elements(boxp1)
ggsavepdf("img/viz_delta.png", plot = boxp, width = 22, height = 22)
