SWADA
================

Source code and materials for ‘Subgroup comparisons within and across
studies in meta-analysis’ (Panaro et al., 2025)
<https://doi.org/10.48550/arXiv.2508.15531>

## Data Availability

Large data files generated during the analyses are **not included in
this repository** due to GitHub’s file size limitations.  
Researchers interested in accessing these data can **contact me** to request a copy.

All files directly related to the article (e.g., analysis scripts,
figures, and supplementary materials) are included in this repository.

## Same Weighting Across Different Analyses (SWADA) in Interaction Meta-Analysis

This README demonstrates how to apply **the same or user-chosen
weighting scheme across multiple meta-analytic analyses** (Difference of
Averages and Average Difference) using aggregate data (AgD).

The goal is to approximate IPD-level inference under a consistent
information structure while allowing flexible comparison of different
weighting options.

------------------------------------------------------------------------

**Here we begin by loading an example dataset** to illustrate the SWADA
concept with real numbers.

## Example Dataset: Steroids and Ventilation

``` r
steroids <- read.csv("data/steroids.csv")
print(steroids)
#>            trial               drug ventilation treat.events treat.total
#> 1     CAPE COVID     Hydrocortisone          IV           10          61
#> 2     CAPE COVID     Hydrocortisone         NIV            1          14
#> 3          CoDEX      Dexamethasone          IV           69         128
#> 4          CoDEX      Dexamethasone         NIV            0           0
#> 5  COVID STEROID     Hydrocortisone          IV            4           7
#> 6  COVID STEROID     Hydrocortisone         NIV            2           8
#> 7  DEXA-COVID 19      Dexamethasone          IV            2           7
#> 8  DEXA-COVID 19      Dexamethasone         NIV            0           0
#> 9       RECOVERY      Dexamethasone          IV           95         324
#> 10      RECOVERY      Dexamethasone         NIV            0           0
#> 11     REMAP-CAP     Hydrocortisone          IV           18          68
#> 12     REMAP-CAP     Hydrocortisone         NIV            8          37
#> 13 Steroids-SARI Methylprednisolone          IV           10          13
#> 14 Steroids-SARI Methylprednisolone         NIV            3          11
#>    cont.events cont.total
#> 1           17         59
#> 2            3         14
#> 3           76        128
#> 4            0          0
#> 5            0          6
#> 6            2          8
#> 7            2         12
#> 8            0          0
#> 9          283        683
#> 10           0          0
#> 11          10         49
#> 12          19         43
#> 13           9         14
#> 14           4          9
```

------------------------------------------------------------------------

**Next, we compute log-odds ratios and their variances** for each
subgroup comparison within trials, preparing the data for meta-analysis.

## Compute Effect Sizes and Variances

``` r
agd <- escalc(measure = "OR",
                ai = treat.events, bi = treat.total - treat.events,
                ci = cont.events, di = cont.total - cont.events, 
              data = steroids) |> 
  mutate(
    subgroup = ifelse(ventilation == "IV", 0, 1),
    ni = treat.total + cont.total,
    subgroup12 = as.numeric(subgroup) - 0.5
  )

print(agd)
#> 
#>            trial               drug ventilation treat.events treat.total 
#> 1     CAPE COVID     Hydrocortisone          IV           10          61 
#> 2     CAPE COVID     Hydrocortisone         NIV            1          14 
#> 3          CoDEX      Dexamethasone          IV           69         128 
#> 4          CoDEX      Dexamethasone         NIV            0           0 
#> 5  COVID STEROID     Hydrocortisone          IV            4           7 
#> 6  COVID STEROID     Hydrocortisone         NIV            2           8 
#> 7  DEXA-COVID 19      Dexamethasone          IV            2           7 
#> 8  DEXA-COVID 19      Dexamethasone         NIV            0           0 
#> 9       RECOVERY      Dexamethasone          IV           95         324 
#> 10      RECOVERY      Dexamethasone         NIV            0           0 
#> 11     REMAP-CAP     Hydrocortisone          IV           18          68 
#> 12     REMAP-CAP     Hydrocortisone         NIV            8          37 
#> 13 Steroids-SARI Methylprednisolone          IV           10          13 
#> 14 Steroids-SARI Methylprednisolone         NIV            3          11 
#>    cont.events cont.total      yi     vi subgroup   ni subgroup12 
#> 1           17         59 -0.7248 0.2022        0  120       -0.5 
#> 2            3         14 -1.2657 1.5012        1   28        0.5 
#> 3           76        128 -0.2229 0.0638        0  256       -0.5 
#> 4            0          0  0.0000 8.0000        1    0        0.5 
#> 5            0          6  2.8163 2.6618        0   13       -0.5 
#> 6            2          8  0.0000 1.3333        1   16        0.5 
#> 7            2         12  0.6931 1.3000        0   19       -0.5 
#> 8            0          0  0.0000 8.0000        1    0        0.5 
#> 9          283        683 -0.5338 0.0209        0 1007       -0.5 
#> 10           0          0  0.0000 8.0000        1    0        0.5 
#> 11          10         49  0.3393 0.2012        0  117       -0.5 
#> 12          19         43 -1.0542 0.2538        1   80        0.5 
#> 13           9         14  0.6162 0.7444        0   27       -0.5 
#> 14           4          9 -0.7577 0.9083        1   20        0.5
```

------------------------------------------------------------------------

**We now define a helper function** to choose the study weights
according to different strategies—equal, inverse-variance, by study
size, or by smaller subgroup size.

## Define Weighting Switch Function

``` r
choose_weights <- function(type, agd, uni = NULL) {
  switch(
    type,
    "equal" = rep(2 / length(agd$vi), length(agd$vi)),
    "inverse_variance" = 1 / agd$vi,
    "study_size" = rep(agd$ni[agd$subgroup == 0] + agd$ni[agd$subgroup == 1], each = 2),
    "smaller_subgroup" = rep(apply(cbind(agd$ni[agd$subgroup == 0],
                                         agd$ni[agd$subgroup == 1]), 1, min), each = 2),
    stop("Unknown weighting type.")
  )
}
```

------------------------------------------------------------------------

**We select one weighting scheme** (e.g., inverse-variance) and fit a
random-effects meta-analytic model using `rma.mv()` from **metafor**.

## Choose Weighting Type and Fit Model

``` r
weight_type <- "inverse_variance"  # try "equal", "study_size", "smaller_subgroup"

wi <- choose_weights(weight_type, agd)
W_common <- diag(wi)

fit_common <- rma.mv(
  yi = yi,
  V = diag(vi),
  W = W_common,
  random = list(~ -1 + subgroup | trial, ~ 1 | trial),
  struct = c("ID", "CS"),  
  method = "REML",
  data = agd
)

summary(fit_common)
#> 
#> Multivariate Meta-Analysis Model (k = 14; method: REML)
#> 
#>   logLik  Deviance       AIC       BIC      AICc   
#> -17.7302   35.4604   41.4604   43.1553   44.1271   
#> 
#> Variance Components:
#> 
#>             estim    sqrt  nlvls  fixed  factor 
#> sigma^2    0.0000  0.0000      7     no   trial 
#> 
#> outer factor: trial    (nlvls = 7)
#> inner factor: subgroup (nlvls = 2)
#> 
#>             estim    sqrt  fixed 
#> tau^2      0.0252  0.1587     no 
#> 
#> Test for Heterogeneity:
#> Q(df = 13) = 13.2662, p-val = 0.4275
#> 
#> Model Results:
#> 
#> estimate      se     zval    pval    ci.lb    ci.ub     
#>  -0.4150  0.1472  -2.8191  0.0048  -0.7034  -0.1265  ** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

------------------------------------------------------------------------

**Now we compare how each weighting approach affects the interaction
estimate**—showing the sensitivity of results to the assumed information
structure.

## Compare Weighting Schemes

``` r
weight_options <- c("equal", "inverse_variance", "study_size", "smaller_subgroup")

res_list <- lapply(weight_options, function(wt) {
  wi <- choose_weights(wt, agd)
  fit <- rma.mv(yi = yi, V = diag(vi), W = diag(wi),
                mods = ~ 1 + subgroup12,      
                random = list(~ -1 + subgroup | trial, ~ 1 | trial),
      struct = c("ID", "CS"),  
      method = "REML", data = agd)
  tibble(weight_type = wt, b1 = -fit$beta[2], se = sqrt(vcov(fit)[2,2]))
})

bind_rows(res_list) %>%
  ggplot(aes(weight_type, exp(b1))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = exp(b1 - 1.96*se), ymax = exp(b1 + 1.96*se)), width = 0.1) +
  coord_flip() +
  theme_minimal() +
  labs(
    y = "Interaction estimate (95% CI)",
    x = "Weighting scheme",
    title = "Impact of Weighting Scheme on Interaction Estimate"
  )
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

------------------------------------------------------------------------

**Finally, we interpret and summarize** how the choice of weighting
affects inference and what each option represents conceptually.

## Interpretation

| Scheme           | Description                                            |
|------------------|--------------------------------------------------------|
| Equal weights    | Treats all studies equally; ignores precision and size |
| Inverse-variance | Standard approach emphasizing interaction precision    |
| Study-size       | Reflects total study sample size; design-based         |
| Smaller subgroup | Rewards balanced subgroups; penalizes imbalance        |

## Citation

Panaro, R., Röver, C., & Friede, T. (2025). *Subgroup comparisons within
and across studies in meta-analysis.* arXiv preprint
**arXiv:2508.15531**.
