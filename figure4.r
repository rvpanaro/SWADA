df <- read.csv("data/steroids.csv")
names(df)

library(metafor)
library(dplyr)
library(ggplot2)
library(ggrepel)

es <- escalc("OR", ai = treat.events, n1i = treat.total,
             ci = cont.events, n2i = cont.total,
             data = df) |> 
  mutate(vi = if_else(vi == 8, 100^2, vi)) # single-subgroup trials

es <- es |>
  group_by(trial) |> 
  mutate(total = sum(treat.events + treat.total)) |>
  group_by(trial, ventilation) |> 
  mutate(sub.total = sum(treat.events + treat.total))

es$ref_prevalence = es$sub.total[rep(seq(1,14,2), each = 2)]/es$total

# ------------------------------------------------------------------------------
# Prepare segment data: only draw lines for studies with >0 IV prevalence
# ------------------------------------------------------------------------------
segments <- es |>
  group_by(trial) |>
  filter(n() == 2) |>  # both subgroups present
  summarize(
    x = ref_prevalence[ventilation == "IV"],
    xend = ref_prevalence[ventilation == "NIV"],
    y = yi[ventilation == "IV"],
    yend = yi[ventilation == "NIV"],
    connect = all(ref_prevalence > 0),
    .groups = "drop"
  ) |>
  filter(connect)

# ------------------------------------------------------------------------------
# Compute study-specific pooled effects (crosses)
# ------------------------------------------------------------------------------
study_means <- es |>
  group_by(trial) |>
  summarize(
    overall_yi = sum(yi / vi) / sum(1 / vi),
    overall_vi = 1 / sum(1 / vi),
    ref_prevalence = unique(ref_prevalence),
    .groups = "drop"
  )

# Prepare filtered NIV line data (exclude prevalence = 0)
niv_line_data <- es |>
  filter(ventilation == "NIV", ref_prevalence<1) |>
  arrange(ref_prevalence)

iv_line_data <- es |>
  filter(ventilation == "IV", ref_prevalence > 0) |>
  arrange(ref_prevalence)


# ------------------------------------------------------------------------------
# Final polished plot with all requested modifications
# ------------------------------------------------------------------------------
final_plot <- ggplot(es, aes(x = ref_prevalence, y = yi, color = ventilation)) +
  
  # Connect IV and NIV within same trial (only if prevalence > 0)
  geom_segment(
    data = segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    color = "grey45", size = 1, linetype = "solid"
  ) +
  
  # Connect all IV and all NIV points across trials
  # geom_line(aes(group = ventilation), size = 1.2, linetype = "solid") +
  
  # Partial NIV line (omit 0 prevalence segments)
  # Partial IV line (excluding prevalence = 0)
  geom_path(
    data = iv_line_data,
    aes(x = ref_prevalence, y = yi),
    inherit.aes = FALSE,
    color = "#1f78b4", size = 1.2, linetype = "solid"
  ) + 
  geom_path(
    data = niv_line_data,
    aes(x = ref_prevalence, y = yi),
    inherit.aes = FALSE,
    color = "#e31a1c", size = 1.2, linetype = "solid"
  ) + 
  
  # Plot subgroup points
  geom_point(aes(size = 1/vi), alpha = 0.8) +
  
  # Add study-level pooled effect as black crosses
  geom_point(
    data = study_means,
    aes(x = ref_prevalence, y = overall_yi),
    shape = 4, size = 4, stroke = 1.5, color = "black", inherit.aes = FALSE
  ) +
  
  # Dashed line tracing study-level pooled effects
  geom_line(
    data = study_means |> arrange(ref_prevalence),
    aes(x = ref_prevalence, y = overall_yi),
    linetype = "dashed", color = "black", linewidth = 1, inherit.aes = FALSE
  ) +
  
  # Add labels to trials
  geom_text_repel(
    data = es[!duplicated(es$trial), ],
    aes(label = trial),
    color = "black",
    size = 5,
    box.padding = 0.35,
    point.padding = 0.5,
    hjust = "left",
    segment.color = "grey60"
  ) +
  
  # Manual colors for subgroups
  scale_color_manual(values = c("IV" = "#1f78b4", "NIV" = "#e31a1c")) +
  
  # Labels
  labs(
    color = "Subgroup",
    size = "Weight (1/Variance)",
    x = "Proportion of invasively ventilated (IV) patients per trial",
    y = expression("Subgroup-specific log-ORs (" * y[Aj] * ", " * y[Bj] * ")"),
    title = "Association between treatment effect and the prevalence"
  ) +
  
  # Theme settings
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_size_continuous(range = c(1, 10))

# Save and display

path <- "img/viz_subgroup_proportion.png"

source("ggsavepdf.r")
ggsavepdf(path, plot = final_plot, width = 14, height = 7)
