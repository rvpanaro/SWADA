# Load necessary libraries
library(ggplot2)
library(patchwork)

# Define parameters
phi <- 0.15
gamma_W <- 4

# Define a function that creates a panel for a given gamma_Agg
plot_panel <- function(gamma_Agg) {
  
  # Create a data frame for the curves (p from 0 to 1)
  p_values <- seq(0, 1, length.out = 100)
  df <- data.frame(
    p = rep(p_values, 2),
    subgroup = rep(c("IV", "Non-IV"), each = length(p_values))
  )
  
  # Compute outcomes for each subgroup:
  # For Invasive ventilated (IV): subgroup0 = phi + (gamma_Agg - gamma_W)*p
  # For non-Invasive ventilated (Non-IV): subgroup1 = phi + gamma_W + (gamma_Agg - gamma_W)*p
  df$outcome <- ifelse(df$subgroup == "IV",
                       phi + (gamma_Agg - gamma_W) * df$p,
                       phi + gamma_W + (gamma_Agg - gamma_W) * df$p)
  
  # Compute annotation coordinates
  # Vertical segment at x = 0.1 on the IV curve:
  y0_01     <- phi + (gamma_Agg - gamma_W) * 0.7
  y0_01_top <- y0_01 + gamma_W
  
  # Endpoints for the curves:
  y0_0 <- phi  # subgroup0 at p = 0
  y1_1 <- phi + gamma_W + (gamma_Agg - gamma_W) * 1  # subgroup1 at p = 1
  
  # For the dashed segment joining (0, subgroup0(0)) to (1, subgroup0(1)+gamma_W)
  seg_y0 <- phi
  seg_y1 <- phi + (gamma_Agg - gamma_W) * 1 + gamma_W
  
  # Positions for annotation text
  text_y_gammaW   <- phi + (gamma_Agg - gamma_W) * 0.125 + gamma_W/2
  text_y_gammaAgg <- phi + (gamma_Agg - gamma_W) * 0.125
  
  # Build the ggplot
  gg <- ggplot(df, aes(x = p, y = outcome, color = subgroup)) +
    geom_line(size = 1.2) +
    scale_color_manual(
      values = c("IV" = "#e31a1c", "Non-IV" = "#1f78b4"),
      labels = c("Invasive ventilation (IV)", "Non-invasive ventilation (NIV)")
    ) +
    # Vertical dashed line at p = 0.7
    geom_segment(x = 0.7, xend = 0.7, y = y0_01, yend = y0_01_top,
                 color = "gray45", linetype = "solid") +
    # Points at the ends of the vertical segment
    geom_point(x = 0.7, y = y0_01, color = "black", size = 2) +
    geom_point(x = 0.7, y = y0_01_top, color = "black", size = 2) +
    geom_point(x = 0.7, y =  0.73 * gamma_W + 0.73 * ( gamma_Agg - gamma_W), color = "black", size = 4, shape = 4) + 
    # Annotation for gamma[W]
    annotate("text", x = 0.125, y = text_y_gammaW, label = "gamma[W]",
             parse = TRUE, hjust = 0) +
    # Points at (0, subgroup0(0)) and (1, subgroup1(1))
    geom_point(x = 0, y = y0_0, color = "black", size = 2) +
    geom_point(x = 1, y = y1_1, color = "black", size = 2) +
    # Dashed segment joining the two endpoints
    geom_segment(x = 0, xend = 1, y = seg_y0, yend = seg_y1,
                 color = "black", linetype = "dashed") +
    # Annotation for gamma[Agg]
    annotate("text", x = 0.5, y = text_y_gammaAgg + 1.5, label = "gamma[Agg]",
             parse = TRUE, hjust = 0) +
    labs(x = "Proportion of IV patients", y = "Subgroup-specific outcome",
         color = "Patient subgroup") +
    ylim(-4, 7) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    geom_text(
      x = 0.05, y = 6,
      label = as.character(
        as.expression(
          bquote(gamma[W] == .(gamma_W) ~ "," ~ 
                   gamma[Agg] == .(gamma_Agg) ~ "," ~ 
                   delta == .(gamma_Agg - gamma_W)))
    ), colour =  "black", 
      parse = TRUE, hjust = 0, size = 5
    )
  
  return(gg)
}

# Create individual panels:
# Main (large) panel with gamma_Agg = gamma_W = 4
main_plot <- plot_panel(gamma_W) + labs(title = "Illustration of Aggregation Bias in Ventilation Outcomes", 
                                        subtitle = "Comparing invasive (IV) vs. non-invasive (NIV) ventilation as the proportion of IV patients varies") 


# Three smaller panels below with different gamma_Agg values (remove only x & y axis labels)
plot_bottom1 <- plot_panel(gamma_W + 2) + labs(x = NULL, y = NULL)
plot_bottom2 <- plot_panel(gamma_W - 2) + labs(x = NULL, y = NULL)
plot_bottom3 <- plot_panel(0)             + labs(x = NULL, y = NULL)

# Combine panels using patchwork with merged legends
combined_plot <- main_plot / (plot_bottom1 | plot_bottom2 | plot_bottom3) +
  plot_layout(heights = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

path <- "img/viz_aggregation_bias_plot.png"

source("ggsavepdf.r")
ggsavepdf(path, plot = combined_plot, width = 14, height = 10)
