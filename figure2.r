df <- read.csv("data/steroids.csv")
names(df)

library("metafor")
library("dplyr")

es <- escalc("OR", ai = treat.events, n1i = treat.total,
             ci = cont.events, n2i = cont.total,
             data = df) |> 
  mutate(vi = if_else(vi == 8, 100^2, vi)) # single-subgroup trials

es <- es |>
  group_by(trial) |> 
  mutate(total = sum(treat.events + treat.total)) |>
  group_by(trial, ventilation) |> 
  mutate(sub.total = sum(treat.events + treat.total))

average_difference <- with(es,
                           rma.uni(
                             yi = yi[ventilation == "IV"] - yi[ventilation == "NIV"],
                             vi = vi[ventilation == "IV"] + vi[ventilation == "NIV"], 
                             method = "REML"
                           )
)

(int_ad <- predict(average_difference))
(sub_ad <- list(pred = rep(NA, 2), ci.lb = rep(NA, 2), ci.ub = rep(NA, 2)))

difference_of_averages <- rma.mv( 
  yi    = yi,
  V     = vi,
  mods   = ~ -1 + ventilation,           # coef for 'corticosteroidyes' = yes − no
  random = ~ - 1 + ventilation | trial,  # random intercept & slope per trial
  data   = es,
  method = "REML", 
  struct = "UN",
  rho = 0
)

(int_da <- predict(difference_of_averages, newmods = matrix(c(1, -1), nrow = 1)))
(sub_da <- predict(difference_of_averages, newmods = diag(c(1, 1))))

cols <- c("pred","ci.lb","ci.ub")

int_intervals <- rbind(print(int_ad)[cols], print(int_da)[cols])
sub_intervals <- rbind((sub_ad)[cols], print(sub_da)[cols]) |> as.matrix()
rownames(sub_intervals) <- rep(c("IV","NIV"), length.out = nrow(sub_intervals))

source("forestplot.R")

fp <- 
  interaction_forest_plot(yi = es$yi,
                              vi = es$vi, 
                              labels = es$trial, 
                              subgroups = es$ventilation,
                              interaction_intervals = int_intervals,
                              subgroup_intervals =  sub_intervals, 
                              model_names = c("Average Difference (AD)",
                                              "Difference of averages (DA)"),
                              reference_group = 2,
                              xlim = c(1/5, 5),
                              xticks = c(0.5, 1, 2),
                              col_right_heading = list(
                                paste0("Subgroup-specific \n(", 100 * 0.95, "% CI)"), 
                                paste0("Interaction \n(", 100 * 0.95, "% CI)")),
                              xlab = c("Odds ratio (OR)", "Ratio of odds ratios (RORs)"),
                              mid.space = unit(2, "mm"),
                              left.space = unit(30, "mm"),
                              plot.margin = margin(8, 8, 8, 8, "mm"), 
                              pointsize = 20, 
                              stroke = 3, 
                              scalepoints = TRUE,
                              base_size = 18,
                              treat.events = es$treat.events,
                              cont.events = es$cont.events,
                              treat.total = es$treat.total,
                              cont.total = es$cont.total, 
                              row.labels.space = c(2,2,2,2)
) 

path <- "img/viz_forest_plot_steroids.png"
source("ggsavepdf.r")
ggsavepdf(path, plot = fp$plot, width = 14, height = 10)
