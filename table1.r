# ==========================================================
# Side-by-side LaTeX coverage table with grouped headers
# ==========================================================
library(dplyr)
library(xtable)

# --- 1. Monte Carlo bounds for 1000 replications ---
lower_bound <- 0.95 - 1.96 * sqrt(0.95 * 0.05 / 1000)  # 0.9365
upper_bound <- 0.95 + 1.96 * sqrt(0.95 * 0.05 / 1000)  # 0.9635

# --- 2. Function to summarize coverage ---
summarise_coverage <- function(data) {
  data %>%
    group_by(method) %>%
    summarise(
      n_under = sum(cov_subgroup1 < lower_bound, na.rm = TRUE),
      n_nominal = sum(cov_subgroup1 >= lower_bound & cov_subgroup1 <= upper_bound, na.rm = TRUE),
      n_over = sum(cov_subgroup1 > upper_bound, na.rm = TRUE),
      total = n(),
      perc_under = 100 * n_under / total,
      perc_nominal = 100 * n_nominal / total,
      perc_over = 100 * n_over / total,
      .groups = "drop"
    )
}

# --- 3. Compute summaries ---
coverage_summary_ag <- summarised %>%
  filter(a_interaction == "Aggregation bias") %>%
  summarise_coverage()

coverage_summary_noag <- summarised %>%
  filter(a_interaction != "Aggregation bias") %>%
  summarise_coverage()

# --- 4. Order by increasing undercoverage in aggregation bias ---
ordered_methods <- coverage_summary_ag %>%
  arrange(perc_under) %>%
  pull(method)

coverage_summary_ag <- coverage_summary_ag %>%
  mutate(method = factor(method, levels = ordered_methods)) %>%
  arrange(method)

coverage_summary_noag <- coverage_summary_noag %>%
  mutate(method = factor(method, levels = ordered_methods)) %>%
  arrange(method)

# --- 5. Merge the two tables and format ---
merged <- coverage_summary_noag %>%
  select(method,
         perc_under_noag = perc_under,
         perc_nominal_noag = perc_nominal,
         perc_over_noag = perc_over) %>%
  left_join(
    coverage_summary_ag %>%
      select(method,
             perc_under_ag = perc_under,
             perc_nominal_ag = perc_nominal,
             perc_over_ag = perc_over),
    by = "method"
  ) %>%
  mutate(across(starts_with("perc_"), ~sprintf("%.1f", .)))

# --- 6. Create xtable with grouped headers ---
tab <- xtable(
  merged,
  caption = "Coverage percentages for Subgroup 1 without and with Aggregation Bias, ordered by increasing undercoverage under Aggregation Bias.",
  label = "tab:coverage_bias",
  align = c("l", "l", rep("r", 6))
)

# --- 7. Custom LaTeX header ---
header <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Coverage percentages for Subgroup 1 without and with Aggregation Bias, ordered by increasing undercoverage under Aggregation Bias.}",
  "\\begin{tabular}{l|rrr|rrr}",
  "\\hline",
  " & \\multicolumn{3}{c|}{\\textbf{Without Aggregation Bias}} & \\multicolumn{3}{c}{\\textbf{With Aggregation Bias}} \\\\",
  "\\cline{2-7}",
  "\\textbf{Method} & \\%Under & \\%Nominal & \\%Over & \\%Under & \\%Nominal & \\%Over \\\\",
  "\\hline",
  ""
)

footer <- c("\\hline", "\\end{tabular}", "\\end{table}")

# --- 8. Print LaTeX table with custom header ---
print(
  tab,
  include.rownames = FALSE,
  include.colnames = FALSE,
  hline.after = NULL,
  only.contents = TRUE,
  comment = FALSE,
  file = "coverage_side_by_side_grouped.tex"
)

# Append header and footer manually
lines <- c(header, readLines("coverage_side_by_side_grouped.tex"), footer)
writeLines(lines, "coverage_side_by_side_grouped.tex")
