df <- read.csv("data/il6.csv")
head(df)

library("metafor")
library("dplyr")

es <- escalc("OR", 
             ai = treat.events, n1i = treat.total,
             ci = cont.events, n2i = cont.total, data = df) |> 
  mutate(vi = if_else(vi == 8, 100^2, vi)) # single-subgroup trials

# exclude single-subgroups
is_single <- es$trial %in% c("ARCHITECTS", "BACC-Bay", "ImmCoVA")
es <- es[!is_single, ]

es <- es |>
  group_by(trial) |> 
  mutate(total = sum(treat.events + treat.total)) |>
  group_by(trial, corticosteroid) |> 
  mutate(sub.total = sum(treat.events + treat.total))

es$subgroup <- factor(es$corticosteroid, levels = c("no", "yes")) |> 
  as.numeric() -1

es$ni <- es$treat.total + es$cont.total

library("metafor")
library("nloptr")
library("tibble")
library("dplyr")

source("computations.r")
source("wt_framework.r")
source("optimized_weights.r")

results <- fit_all(x = list(agd = list(data.frame(es)))) |>
  data.frame() |> 
  {\(.) {
    rownames(.) <- NULL              # remove any existing row names
    tibble::column_to_rownames(., var = "analysis")
  }}()

# rownames(results)[6:11] <- paste(rownames(results)[6:11], "SWADA")

cols <- c("b","ci.lb","ci.ub")

int_intervals <- results |>
  select(paste0("interaction_", cols)) |> 
  slice(rep(1:n(), each = 2)) |> 
  mutate(across(everything(), ~replace(.x, seq(2, n(), 2), NA)))  # NA every second row

sub1 <- results |> 
  select(paste0("subgroupB_", cols)) |>
  setNames(cols)

sub2 <- results |> 
  select(paste0("subgroupA_", cols))|>
  setNames(cols)

sub_intervals <- do.call(rbind, Map(rbind, split(sub1, seq_len(nrow(sub1))), split(sub2, seq_len(nrow(sub2))))) |> as.matrix()
rownames(sub_intervals) <- rep(c("yes", "no"), length.out = nrow(sub_intervals))

source("forestplot.R")

col_sub <- rep(c(NA, "white", "royalblue3"), times = c(nrow(es), 2*5, 2*6))
col_int <- rep(c(NA, "white", "chartreuse3"), times = c(nrow(es), 2*5, 2*6))

names(col_sub) <- c(es$trial, rep(rownames(results), each = 2))
names(col_int) <- c(es$trial, rep(rownames(results), each = 2))

col_int[names(col_int) %in% c( 
  "Within-trial framework", 
  "Separate MAs (DA and AD)", 
  "Interaction weights SWADA")] <- "gray30"

fill <- list(col_sub, col_int)

fp <- interaction_forest_plot(yi = es$yi,
                              vi = es$vi, 
                              labels = es$trial, 
                              subgroups = es$corticosteroid,
                              interaction_intervals = int_intervals,
                              subgroup_intervals =  sub_intervals, 
                              model_names = rownames(results),
                              reference_group = 2,
                              xlim = c(.3, 3),
                              xticks = c(.4, 1, 1.6),
                              col_right_heading = list(paste0("Subgroup-specific \n(", 100 * 0.95, "% CI)"), 
                                                       paste0("Interaction \n(", 100 * 0.95, "% CI)")),
                              xlab = c("Odds ratio (OR)", "Ratio of odds ratios (RORs)"),
                              mid.space = unit(2, "mm"),
                              left.space = unit(10, "mm"),
                              plot.margin = margin(8, 8, 8, 8, "mm"), 
                              pointsize = 4, 
                              stroke = 3, 
                              scalepoints = FALSE,
                              estimates_only = TRUE,
                              base_size = 18,
                              # grouping = rep(c(NA, "Existing methods", "SWADAs"), times = c(nrow(es), 2*5, 2*6)),
                              row.labels.space = c(2,2,2,2), 
                              fill = fill
)

# include the dashed reference line
bg_line <-  geom_vline(
  data        = data.frame(panel = "2"),
  aes(xintercept = 0.69),
  linetype    = "dashed",
  colour      = "gray30",
  linewidth = 1.2
)

fp$plot$layers <- c(list(bg_line), fp$plot$layers)

path <- "img/viz_comparison_IL6_single.png"
ggsave(path, plot = fp$plot, width = 14, height = 10)

library(magick)
img <- image_read(path)

out <- sub("\\.png$", ".pdf", path)
image_write(img, path = out, format = "pdf")

file.remove(path)

