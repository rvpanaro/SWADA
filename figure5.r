# Load necessary libraries
library(dplyr)
library(tidyr)
library(reactable)
library(htmltools)
library(htmlwidgets)
library(webshot)
library(nloptr)  # required by optimized_weights

# Source the optimized_weights function
source("optimized_weights.R")

# Function to validate that each weight vector sums to 1 (within tolerance)
validate_probabilities <- function(probs) {
  sum_proba <- sapply(probs, sum)
  all(abs(sum_proba - 1) < 1e-6)
}

# ----------------------
# Case 1: 50+50 / 70+30 / 90+10 patients
# ----------------------
# Study 1: 50+50
n1A_case1 <- 50; n1B_case1 <- 50; total1_case1 <- n1A_case1 + n1B_case1  # 100
# Study 2: 70+30
n2A_case1 <- 70; n2B_case1 <- 30; total2_case1 <- n2A_case1 + n2B_case1  # 100
# Study 3: 90+10
n3A_case1 <- 90; n3B_case1 <- 10; total3_case1 <- n3A_case1 + n3B_case1  # 100

w_equal_case1 <- rep(1/3, 3)
w_int_case1 <- c( 1/(1/n1A_case1 + 1/n1B_case1),
                  1/(1/n2A_case1 + 1/n2B_case1),
                  1/(1/n3A_case1 + 1/n3B_case1) )
w_int_case1 <- w_int_case1 / sum(w_int_case1)
w_size_case1 <- c(total1_case1, total2_case1, total3_case1)
w_size_case1 <- w_size_case1 / sum(w_size_case1)
w_small_case1 <- c(min(n1A_case1, n1B_case1),
                   min(n2A_case1, n2B_case1),
                   min(n3A_case1, n3B_case1))
w_small_case1 <- w_small_case1 / sum(w_small_case1)
w_min_case1 <- c( min(n1A_case1, n1B_case1, 1/(1/n1A_case1 + 1/n1B_case1)),
                  min(n2A_case1, n2B_case1, 1/(1/n2A_case1 + 1/n2B_case1)),
                  min(n3A_case1, n3B_case1, 1/(1/n3A_case1 + 1/n3B_case1)) )
w_min_case1 <- w_min_case1 / sum(w_min_case1)

# For Scheme 6 in Case 1: compute subgroup variances (assuming UISD^2 = 1, so variance = 1/n)
vi1_case1 <- c(1/n1A_case1, 1/n2A_case1, 1/n3A_case1)
vi2_case1 <- c(1/n1B_case1, 1/n2B_case1, 1/n3B_case1)
opt_case1 <- optimized_weights(vi1_case1, vi2_case1, diag(0, 2), off_diagonal = FALSE)
w_totvar_case1 <- opt_case1$wi

# ----------------------
# Case 2: 50+50 / 50+100 / 100+100 patients
# ----------------------
# Study 1: 50+50
n1A_case2 <- 50; n1B_case2 <- 50; total1_case2 <- n1A_case2 + n1B_case2  # 100
# Study 2: 50+100
n2A_case2 <- 50; n2B_case2 <- 100; total2_case2 <- n2A_case2 + n2B_case2  # 150
# Study 3: 100+100
n3A_case2 <- 100; n3B_case2 <- 100; total3_case2 <- n3A_case2 + n3B_case2  # 200

w_equal_case2 <- rep(1/3, 3)
w_int_case2 <- c( 1/(1/n1A_case2 + 1/n1B_case2),
                  1/(1/n2A_case2 + 1/n2B_case2),
                  1/(1/n3A_case2 + 1/n3B_case2) )
w_int_case2 <- w_int_case2 / sum(w_int_case2)
w_size_case2 <- c(total1_case2, total2_case2, total3_case2)
w_size_case2 <- w_size_case2 / sum(w_size_case2)
w_small_case2 <- c(min(n1A_case2, n1B_case2),
                   min(n2A_case2, n2B_case2),
                   min(n3A_case2, n3B_case2))
w_small_case2 <- w_small_case2 / sum(w_small_case2)
w_min_case2 <- c( min(n1A_case2, n1B_case2, 1/(1/n1A_case2+1/n1B_case2)),
                  min(n2A_case2, n2B_case2, 1/(1/n2A_case2+1/n2B_case2)),
                  min(n3A_case2, n3B_case2, 1/(1/n3A_case2+1/n3B_case2)) )
w_min_case2 <- w_min_case2 / sum(w_min_case2)

vi1_case2 <- c(1/n1A_case2, 1/n2A_case2, 1/n3A_case2)
vi2_case2 <- c(1/n1B_case2, 1/n2B_case2, 1/n3B_case2)
opt_case2 <- optimized_weights(vi1_case2, vi2_case2, diag(0,2), off_diagonal = FALSE)
w_totvar_case2 <- opt_case2$wi

# ----------------------
# Case 3: 25+25 / 50+50 / 100+100 patients
# ----------------------
# Study 1: 25+25
n1A_case3 <- 25; n1B_case3 <- 25; total1_case3 <- n1A_case3 + n1B_case3  # 50
# Study 2: 50+50
n2A_case3 <- 50; n2B_case3 <- 50; total2_case3 <- n2A_case3 + n2B_case3  # 100
# Study 3: 100+100
n3A_case3 <- 100; n3B_case3 <- 100; total3_case3 <- n3A_case3 + n3B_case3  # 200

w_equal_case3 <- rep(1/3, 3)
w_int_case3 <- c( 1/(1/n1A_case3 + 1/n1B_case3),
                  1/(1/n2A_case3 + 1/n2B_case3),
                  1/(1/n3A_case3 + 1/n3B_case3) )
w_int_case3 <- w_int_case3 / sum(w_int_case3)
w_size_case3 <- c(total1_case3, total2_case3, total3_case3)
w_size_case3 <- w_size_case3 / sum(w_size_case3)
w_small_case3 <- c(min(n1A_case3, n1B_case3),
                   min(n2A_case3, n2B_case3),
                   min(n3A_case3, n3B_case3))
w_small_case3 <- w_small_case3 / sum(w_small_case3)
w_min_case3 <- c( min(n1A_case3, n1B_case3, 1/(1/n1A_case3+1/n1B_case3)),
                  min(n2A_case3, n2B_case3, 1/(1/n2A_case3+1/n2B_case3)),
                  min(n3A_case3, n3B_case3, 1/(1/n3A_case3+1/n3B_case3)) )
w_min_case3 <- w_min_case3 / sum(w_min_case3)

vi1_case3 <- c(1/n1A_case3, 1/n2A_case3, 1/n3A_case3)
vi2_case3 <- c(1/n1B_case3, 1/n2B_case3, 1/n3B_case3)
opt_case3 <- optimized_weights(vi1_case3, vi2_case3, diag(0,2), off_diagonal = FALSE)
w_totvar_case3 <- opt_case3$wi

# ----------------------
# Case 4: 25+50 / 50+100 / 100+200 patients
# ----------------------
# Study 1: 25+50
n1A_case4 <- 25; n1B_case4 <- 50; total1_case4 <- n1A_case4 + n1B_case4  # 75
# Study 2: 50+100
n2A_case4 <- 50; n2B_case4 <- 100; total2_case4 <- n2A_case4 + n2B_case4  # 150
# Study 3: 100+200
n3A_case4 <- 100; n3B_case4 <- 200; total3_case4 <- n3A_case4 + n3B_case4  # 300

w_equal_case4 <- rep(1/3, 3)
w_int_case4 <- c( 1/(1/n1A_case4 + 1/n1B_case4),
                  1/(1/n2A_case4 + 1/n2B_case4),
                  1/(1/n3A_case4 + 1/n3B_case4) )
w_int_case4 <- w_int_case4 / sum(w_int_case4)
w_size_case4 <- c(total1_case4, total2_case4, total3_case4)
w_size_case4 <- w_size_case4 / sum(w_size_case4)
w_small_case4 <- c(min(n1A_case4, n1B_case4),
                   min(n2A_case4, n2B_case4),
                   min(n3A_case4, n3B_case4))
w_small_case4 <- w_small_case4 / sum(w_small_case4)
w_min_case4 <- c( min(n1A_case4, n1B_case4, 1/(1/n1A_case4+1/n1B_case4)),
                  min(n2A_case4, n2B_case4, 1/(1/n2A_case4+1/n2B_case4)),
                  min(n3A_case4, n3B_case4, 1/(1/n3A_case4+1/n3B_case4)) )
w_min_case4 <- w_min_case4 / sum(w_min_case4)

vi1_case4 <- c(1/n1A_case4, 1/n2A_case4, 1/n3A_case4)
vi2_case4 <- c(1/n1B_case4, 1/n2B_case4, 1/n3B_case4)
opt_case4 <- optimized_weights(vi1_case4, vi2_case4, diag(0,2), off_diagonal = FALSE)
w_totvar_case4 <- opt_case4$wi

# ----------------------
# Create a tibble with the results for each case and each weighting method
# Each cell is a vector of 3 weights (one per study)
data <- tibble(
  Method = c("1. Equal Weights", 
             "2. Interaction Weights", 
             "3. Study Size Weights", 
             "4. Smaller Subgroup Weights", 
             "5. Minimum of Three",
             "6. Minimum Total Variance Weights"),
  `Case 1` = list(w_equal_case1, w_int_case1, w_size_case1, w_small_case1, w_min_case1, w_totvar_case1),
  `Case 2` = list(w_equal_case2, w_int_case2, w_size_case2, w_small_case2, w_min_case2, w_totvar_case2),
  `Case 3` = list(w_equal_case3, w_int_case3, w_size_case3, w_small_case3, w_min_case3, w_totvar_case3),
  `Case 4` = list(w_equal_case4, w_int_case4, w_size_case4, w_small_case4, w_min_case4, w_totvar_case4)
)

# Validate probabilities
is_valid <- sapply(names(data)[-1], function(case) {
  validate_probabilities(data[[case]])
})
print(is_valid)

# For global scaling in the histograms, find the maximum weight across all cases & methods
max_prob_global <- max(unlist(data %>% select(-Method)))

# Function to generate a vertical histogram (bar chart) for a given vector of weights
generate_vertical_histogram <- function(probs) {
  bars <- lapply(probs, function(prob) {
    div(
      style = list(
        background = "#0072B2",
        width = "30px",
        height = paste0((prob / max_prob_global) * 50, "px"),
        margin = "2px 6px",
        position = "relative"
      ),
      div(
        style = list(
          position = "absolute",
          top = "-20px",
          width = "100%",
          textAlign = "center",
          fontSize = "20px",
          fontFamily = "Times New Roman, serif",
          color = "#000000"
        ),
        paste0(round(prob * 100, 0), "%")
      )
    )
  })
  
  div(
    style = list(
      display = "flex",
      flexDirection = "row",
      alignItems = "flex-end",
      justifyContent = "center",
      height = "80px"
    ),
    bars
  )
}

# Create the reactable table with improved styling for publication
table_output <- reactable(
  data,
  columns = list(
    Method = colDef(
      name = "Weighting Scheme", 
      style = list(fontFamily = "Times New Roman, serif", fontSize = "20px")
    ),
    `Case 1` = colDef(
      name = "Case 1: 50+50 / 70+30 / 90+10",
      cell = function(value) generate_vertical_histogram(value)
    ),
    `Case 2` = colDef(
      name = "Case 2: 50+50 / 50+100 / 100+100",
      cell = function(value) generate_vertical_histogram(value)
    ),
    `Case 3` = colDef(
      name = "Case 3: 25+25 / 50+50 / 100+100",
      cell = function(value) generate_vertical_histogram(value)
    ),
    `Case 4` = colDef(
      name = "Case 4: 25+50 / 50+100 / 100+200",
      cell = function(value) generate_vertical_histogram(value)
    )
  ),
  defaultColDef = colDef(align = "center"),
  bordered = FALSE,
  striped = FALSE,
  highlight = FALSE,
  theme = reactableTheme(
    borderColor = "#e6e6e6",
    headerStyle = list(
      borderBottom = "1px solid #000000",
      fontFamily = "Times New Roman, serif",
      fontSize = "20px"
    ),
    cellStyle = list(
      padding = "4px 8px",
      fontFamily = "Times New Roman, serif",
      fontSize = "20px"
    )
  ),
  pagination = FALSE
)

# Save the table as an HTML file
path0 <- "img/viz_table_output.html" 
saveWidget(table_output, path0, selfcontained = TRUE)

# Convert the HTML to an image with improved quality
path <- "img/viz_table_output.png"
webshot2::webshot(path0, path, 
                  vwidth = 1400, vheight = 1400)

library(magick)
img <- image_read(path)

out <- sub("\\.png$", ".pdf", path)
image_write(img, path = out, format = "pdf")

file.remove(path0)
file.remove(path)
