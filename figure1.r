df <- read.csv("data/il6.csv")
head(df)

library("metafor")
library("dplyr")

es <- escalc("OR", 
             ai = treat.events, n1i = treat.total,
             ci = cont.events, n2i = cont.total, data = df) |> 
  mutate(vi = if_else(vi == 8, 100^2, vi)) # single-subgroup trials

es <- es |>
  group_by(trial) |> 
  mutate(total = sum(treat.events + treat.total)) |>
  group_by(trial, corticosteroid) |> 
  mutate(sub.total = sum(treat.events + treat.total))

average_difference <- with(es,
  rma.uni(
    yi = yi[corticosteroid == "yes"] - yi[corticosteroid == "no"],
    vi = vi[corticosteroid == "yes"] + vi[corticosteroid == "no"], 
    method = "REML"
    )
  )

(int_ad <- predict(average_difference))
(sub_ad <- list(pred = rep(NA, 2), ci.lb = rep(NA, 2), ci.ub = rep(NA, 2)))

difference_of_averages_ce <- rma.mv( 
  yi    = yi,
  V     = vi,
  mods   = ~ -1 + corticosteroid,           # coef for 'corticosteroidyes' = yes − no
  random = ~ - 1 + corticosteroid | trial,  # random intercept & slope per trial
  data   = es,
  method = "CE"
)

difference_of_averages <- rma.mv( 
  yi    = yi,
  V     = vi,
  mods   = ~ -1 + corticosteroid,           # coef for 'corticosteroidyes' = yes − no
  random = ~ - 1 + corticosteroid | trial,  # random intercept & slope per trial
  data   = es,
  method = "REML", 
  struct = "UN",
  rho = 0
)

(predict(difference_of_averages_ce, newmods = matrix(c(-1, 1), nrow = 1), transf = exp))
(predict(difference_of_averages, newmods = matrix(c(-1, 1), nrow = 1), transf = exp))

(int_da <- predict(difference_of_averages, newmods = matrix(c(-1, 1), nrow = 1)))
(sub_da <- predict(difference_of_averages, newmods = diag(c(1, 1))))

cols <- c("pred","ci.lb","ci.ub")

int_intervals <- rbind(print(int_ad)[cols], print(int_da)[cols])
sub_intervals <- rbind((sub_ad)[cols], print(sub_da)[2:1,cols]) |> as.matrix()

rownames(sub_intervals) <- rep(c("yes","no"), length.out = nrow(sub_intervals))

source("forestplot.R")

fp <- interaction_forest_plot(yi = es$yi,
                              vi = es$vi, 
                        labels = es$trial, 
                        subgroups = es$corticosteroid,
                        interaction_intervals = int_intervals,
                        subgroup_intervals =  sub_intervals, 
                        model_names = c("Average Difference (AD)", "Difference of averages (DA)"),
                        reference_group = 2,
                        xlim = c(.3, 3),
                        xticks = c(.4, 1, 1.6),
                        col_right_heading = list(paste0("Subgroup-specific \n(", 100 * 0.95, "% CI)"), 
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

source("ggsavepdf.r")
path <- "img/viz_forest_plot_il6.png"
ggsavepdf(path, plot = fp$plot, width = 14, height = 14)
