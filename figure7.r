# Source simulation data
source("sim_viz_preamble.r")
source("ggsavepdf.r")
# system("module load imagemagick/7.1.1-39")

# =============================================================================
# 7. Plot: Relative SE Length (Interaction)
# =============================================================================
factor = 2

df1 <- summarised |>
  group_by(method) |>
  filter(a_interaction == "Aggregation bias")

df2 <- summarised |>
  group_by(method) |>
  filter(a_interaction != "Aggregation bias")

df1_small <- df1 |> filter(k == 20, scenario == "More variation")
df2_small <- df2 |> filter(k == 20, scenario == "More variation")

p1_small <- df1_small |> 
  ggplot(aes(
    x     = as.factor(round(tauW, 2)),
    y     = 100 * se_interaction_rel_AD,
    color = method, shape = method, group = method
  ))+
  coord_cartesian(ylim = c(100-15,100+15)) + 
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4)  +
  scale_shape_manual(values = 0:11) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  labs(
    title = "Noncollapsibility with confounding \n (aggregation bias)",
    x = " ",
    y     = "Average change on interaction interval width"
  )  +
  defined_theme +
  theme(
    text          = element_text(size = 11 * factor),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text     = element_text(size = 9 * factor),
    axis.title    = element_text(size = 10 * factor),
    strip.text    = element_text(size = 10 * factor),
    plot.title    = element_text(size = 12 * factor),
    plot.subtitle = element_text(size = 10 * factor),
    plot.caption  = element_text(size = 8 * factor)
  )

p2_small <- df2_small |> 
  ggplot(aes(
    x     = as.factor(round(tauW, 2)),
    y     = 100 * se_interaction_rel_AD,
    color = method, shape = method, group = method
  ))+
  coord_cartesian(ylim = c(100-15,100+15)) + 
  geom_abline(intercept = 100, slope = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 4) +
  scale_shape_manual(values = 0:11) +
  scale_y_continuous(labels = ~ paste0(.x, "%")) +
  labs(
    title = "Noncollapsibility without confounding \n (no aggregation bias)",
    x = expression("Interaction heterogeneity" ~ (tau[W])),
    y     = "Average change on interactin interval width"
  )  +
  defined_theme +
  theme(
    text          = element_text(size = 11 * factor),
    legend.position = "none",
    axis.text     = element_text(size = 9 * factor),
    axis.title    = element_text(size = 10 * factor),
    strip.text    = element_text(size = 10 * factor),
    plot.title    = element_text(size = 12 * factor),
    plot.subtitle = element_text(size = 10 * factor),
    plot.caption  = element_text(size = 8 * factor)
  )

# assume p1_small and p2_small share exactly the same scale/aesthetic
p_small <- ( 
  p2_small + p1_small + 
    plot_layout(guides = "collect")    # collect all individual legends
) & 
  theme(
    legend.position  = "bottom" ,      # one legend at bottom
     legend.direction = "horizontal"   # laid out in a row
  )

p_small   
ggsavepdf("img/viz_se_gamma_small.png", plot = p_small, width = 20, height = 10)
