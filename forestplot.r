# interaction_forest_plot.R
#
#    Copyright (C) 2024  Renato Panaro
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# R code to generate interaction forest plots for subgroup analyses, 
# visualizing effect size estimates and confidence intervals.
#
#   Authors: Renato Panaro
#   Last Updated: December 2024
#
#' Generate Interaction Forest Plots
#'
#' This function creates forest plots to visualize effect size estimates and confidence intervals
#' for subgroup analyses, including interaction effects.  
#' The function relies on the `ckbplotr` package for forest plot visualization and the `metafor`  
#' package for effect size and variance calculations.
#'
#' @param yi Numeric vector of effect size estimates.
#' @param vi Numeric vector of variances for the effect size estimates.
#' @param labels Character vector of labels for each subgroup.
#' @param subgroups Character vector indicating subgroup membership.
#' @param subgroup_intervals Data frame containing confidence intervals for subgroups.
#' @param interaction_intervals Data frame containing confidence intervals for interaction effects.
#' @param model_names Character vector of model names.
#' @param treat.events Numeric vector of events in the treatment group.
#' @param treat.total Numeric vector of total subjects in the treatment group.
#' @param cont.events Numeric vector of events in the control group.
#' @param cont.total Numeric vector of total subjects in the control group.
#' @param reference_group Integer specifying the reference group index.
#' @param group_labels Character vector for labeling groups.
#' @param xlab Character string for x-axis label.
#' @param exponentiate Logical; whether to exponentiate effect sizes.
#' @param col_right_heading Character vector for right column headings.
#' @param stroke Numeric; stroke width of points.
#' @param pointsize Numeric; size of points.
#' @param base_size Numeric; base font size.
#' @param prob Numeric; probability level for confidence intervals.
#' @param mid.space Unit; spacing for mid-section (from `grid::unit()`).
#' @param left.space Unit; spacing for left section (from `grid::unit()`).
#' @param plot.margin Margin; plot margin settings (from `ggplot2::margin()`).
#' @param color_palette Character vector of colors for subgroups.
#' @param xlim Numeric vector defining x-axis limits.
#' @param xticks Numeric vector of x-axis tick positions.
#' @param ... Additional arguments passed to `ckbplotr::forest_plot`.
#'
#' @return A forest plot object (from `ckbplotr`).
#'
#' @details
#' This function leverages:
#' - **`ckbplotr::forest_plot()`**: For rendering forest plots.
#' - **`metafor::escalc()`**: For computing effect sizes and variances.
#' - **`magrittr` pipes (`%>%`)**: For streamlined data manipulation.
#' - **`dplyr::group_by()` and `dplyr::do()`**: For grouped calculations.
#' - **`grid::unit()`**: To define spacing in plot layouts.
#' - **`ggplot2::margin()`**: For setting plot margins.
#'
#' @examples
#' # Example usage
#' interaction_forest_plot(
#'   yi = c(0.5, 1.2, 0.8),
#'   vi = c(0.04, 0.06, 0.05),
#'   labels = c("Subgroup 1", "Subgroup 2", "Subgroup 3"),
#'   subgroups = c("Group A", "Group B", "Group C"),
#'   subgroup_intervals = data.frame(
#'     lower = c(0.3, 1.0, 0.6),
#'     upper = c(0.7, 1.4, 1.0)
#'   ),
#'   interaction_intervals = data.frame(
#'     lower = c(0.2, 0.9, 0.5),
#'     upper = c(0.8, 1.5, 1.1)
#'   )
#' )
#'
#' @import ckbplotr
#' @import metafor
#' @import dplyr
#' @import magrittr
#' @importFrom grid unit
#' @importFrom ggplot2 margin
#' @export

interaction_forest_plot <-
  function(yi, vi, labels, subgroups,
           subgroup_intervals, interaction_intervals,
           model_names = "Model 1",
           treat.events = "",
           treat.total = "",
           cont.events = "",
           cont.total = "",
           subgroup.total = "",
           total = "",
           reference_group = 1,
           xlab, exponentiate = TRUE, col_right_heading, 
           stroke = 3,
           digits = 2,
           fill = "white",
           pointsize = 13, base_size = 20, prob = 0.95,
           left.space = unit(30, "mm"),
           xlim = c(0.125, 8),
           xticks = c(0.125, 1, 8),
           addaes = c("identity", "log"),            
           grouping,
           row.labels.space = c(0,1,0,0),
           estimates_only = FALSE, ...) {
    
    group_labels = rep("", length(labels))     # (deprecated)
    
    # Load required libraries
    library(ckbplotr)
    library(metafor)
    library(dplyr)
    library(magrittr)
    library(grid)
    library(ggplot2)
    library(ggtext)
    # Argument validation
    addaes <- match.arg(addaes)  # Garante que seja "identity" ou "log"
    
    # switch function application
    addaes <- switch(
      addaes,
      log = list(
        paste0(
          "stroke = stroke / 10 * ((2 * 1.96) / ",
          "(log(uci_transformed) - log(lci_transformed)))^(1/2)"
        )
      ),
      identity = list(
        paste0(
          "stroke = stroke / 10 * ((2 * 1.96) / ",
          "(uci_transformed - lci_transformed))^(1/2)"
        )
      )
    )
    
    # Ensure subgroup_intervals has row names
    stopifnot(!is.null(rownames(subgroup_intervals)))
    
    # Extract relevant properties for plotting
    num_subgroups <- length(unique(rownames(subgroup_intervals))) # Count unique subgroups
    k <- nrow(subgroup_intervals)                          # Number of rows in subgroups intervals
    
    # Prepare row labels data frame
    row_labels <- data.frame(
      subheading = labels,
      heading = group_labels,
      label = subgroups
    )
    
    # Extend row labels for model names and subgroups details
    n <- nrow(row_labels) # Get current number of rows in row_labels
    additional_rows <- num_subgroups * length(model_names) # Calculate additional rows needed
    
    row_labels[(n + 1):(n + additional_rows), ] <- cbind(
      rep(model_names, each = num_subgroups),                     # Repeat model names for each subgroups
      sub("\\|.*", "", rownames(subgroup_intervals)),         # Extract subgroup name before '|' as heading
      sub(".*\\|", "", rownames(subgroup_intervals))          # Extract subgroup detail after '|' as label
    )
    
    # Calculate effect size summaries for plotting
    resultsA <- metafor::escalc(yi = yi, vi = vi) %>%
      summary(level = 100 * prob, transf = identity) # Summarize with confidence intervals
    
    # Replace large confidence intervals with NA for better visualization
    resultsA[resultsA$ci.ub > 195, ] <- NA
    
    # Append subgroup intervals to results
    resultsA[(n + 1):(n + additional_rows), ] <- subgroup_intervals
    
    # Set colors and formatting for results
    resultsA$fill <- resultsA$colour <- NA
    resultsA$fill[row_labels$subheading %in% labels] <- "black"
    resultsA$colour[row_labels$subheading %in% labels] <- "black"
    resultsA$fill[!(row_labels$subheading %in% labels)] <- "white"
    resultsA$colour[!(row_labels$subheading %in% labels)] <- "royalblue3"
    resultsA$bold <- resultsA$colour != "black"
    resultsA$colour[resultsA$bold] <- "royalblue3"
    
    # Assign row identifiers for linking labels and results
    resultsA$variable <- row_labels$variable <- apply(row_labels, 1, paste0, collapse = "")
    resultsB <- resultsA
    
    # Calculate interaction effects based on subgroup endpoints
    compute_interactions <- function(grouped, reference = 1, point = TRUE) {
      if (!all(c("endpoints", "yi", "vi") %in% colnames(grouped))) {
        stop("Input data 'grouped' must contain 'endpoints', 'yi', and 'vi' columns.")
      }
      
      unique_subgroups <- unique(subgroups)
      
      if (reference < 1 || reference > length(unique_subgroups)) {
        stop("Reference group index is out of bounds.")
      }
      
      C <- diag(length(unique_subgroups))
      C[, reference] <- C[, reference] - 1
      
      idx <- match(grouped$endpoints, unique_subgroups)
      if (any(is.na(idx))) {
        stop("Some endpoints in 'grouped' do not match the subgroups.")
      }
      
      if (!point) {
        C <- abs(C)
        result <- C[idx, idx] %*% grouped$vi
      } else {
        result <- C[idx, idx] %*% grouped$yi
      }
      
      return(result)
    }
    
    # Step 1: Create the subset data frame
    data_subset <- data.frame(
      yi = yi,
      vi = vi,
      labels = labels,
      endpoints = subgroups,
      groups = group_labels
    )
    
    # Step 2: Group data by 'groups' and 'labels'
    grouped_data <- data_subset %>%
      group_by(groups, labels)
    
    # Step 3: Compute interactions for 'yi' and 'vi'
    interaction_results_yi <- grouped_data %>%
      do(data.frame(yi = compute_interactions(., reference = reference_group)))
    
    interaction_results_vi <- grouped_data %>%
      do(data.frame(vi = compute_interactions(., point = FALSE, reference = reference_group)))
    
    # Step 4: Match labels and extract corresponding 'yi' and 'vi' values
    matched_yi <- interaction_results_yi %>%
      .[pmatch(labels, .$labels), ] %>%
      .$yi
    
    matched_vi <- interaction_results_vi %>%
      .[pmatch(labels, .$labels), ] %>%
      .$vi
    
    # Step 5: Update 'resultsB$yi' and 'resultsB$vi' for entries where 'bold' is FALSE
    resultsB$yi[!resultsB$bold] <- matched_yi
    resultsB$vi[!resultsB$bold] <- matched_vi
    
    resultsB[, c("ci.lb", "ci.ub")] <-
      summary(escalc(yi = resultsB$yi, vi = resultsB$vi),
              level = 100 * prob)[, c("ci.lb", "ci.ub")]
    
    # Step 1: Identify rows to update based on subheading matching model names
    rows_to_update <- which(row_labels$subheading %in% model_names, c("yi", "ci.lb", "ci.ub"))
    
    # Step 2: Calculate indices to update
    include_indices <- seq(rows_to_update[1], rows_to_update[k], num_subgroups) + reference_group-1
    
    # Step 3: Update selected rows with non-missing interaction intervals
    resultsB[include_indices, c("yi", "ci.lb", "ci.ub")] <- na.omit(interaction_intervals)
    resultsB[resultsB$ci.ub > 195 & !is.na(resultsB$ci.ub), ] <- NA
    resultsB[which(resultsB$vi == 0), c("yi", "ci.lb", "ci.ub")] <- NA
    resultsB[is.na(resultsB$ci.lb), c("yi", "ci.lb", "ci.ub")] <- NA
    resultsB$colour[!is.na(resultsB$yi)] <- "black"
    resultsB$colour[resultsB$bold] <- "chartreuse3"
    
    if (missing(xlab)) {
      xlab <- c("", "")
    }
    if (missing(col_right_heading)) {
      col_right_heading <- list(paste0("Subgroup \n(", 100 * prob, "% interval)"), paste0("Interaction \n(", 100 * prob, "% interval)"))
    }
    
    row_labels$heading <- row_labels$subheading
    row_labels$subheading <- NA_character_
    
    # Define new shapes
    resultsA$shape <- 22
    resultsB$shape <- 21
    
    # Define diamond shape for summary estimates
    resultsA$shape[resultsA$bold] <- 23
    resultsB$shape[resultsB$bold] <- 23
    
    # Initialize outputs
    col.left <- NULL
    col.left.heading <- ""
    
    resultsA$treat <- resultsA$cont <- resultsA$prev <- ""
    
    # Check if totals are provided
    if(all(treat.total != "")) {
      # Totals for control and treatment groups
      resultsA$treat <- c(paste0(treat.events, "/", treat.total), 
                          rep("", nrow(resultsA) - length(treat.total)))
      resultsA$cont <- c(paste0(cont.events, "/", cont.total), 
                         rep("", nrow(resultsA) - length(treat.total)))
      col.left.heading_list1 <- list(c("Control \n (n/N)", "Treatment \n (n/N)"), c("", ""))
      
      resultsA$treat[resultsA$treat == "0/0"] <- ""
      resultsA$cont[resultsA$cont == "0/0"] <- ""
    }
    
    
    # Check if prevalence is provided
    if(all(total != "")) {
      # Prevalence formatting
      resultsA$prev <- c(paste0("       ", subgroup.total, " (",round(subgroup.total/total*100, digits =  0), "%)"), 
                         rep("", nrow(resultsA) - length(total)))
      col.left.heading_list2 <- list(paste0("Prevalence \n (n/N)"), "")
    }
    
    # Determine which columns to display
    if(all(treat.total != "") & !all(total != "")) {
      # Only totals
      col.left = c("cont", "treat")
      col.left.heading <- col.left.heading_list1
      
    } else if(!all(treat.total != "") & all(total != "")) {
      # Only prevalence
      col.left = c("prev")
      col.left.heading <- col.left.heading_list2
      
    } else if(all(treat.total != "") & all(total != "")) {
      # Both totals
      col.left = c("cont", "treat", "prev")
      col.left.heading <- c(col.left.heading_list1[[1]], col.left.heading_list2)
      
    } else {
      # Neither totals nor prevalence provided
      col.left <- NULL
      col.left.heading <- ""
    }
    
    resultsB$treat <- resultsB$cont <- resultsB$prev <- ""
    
    ## remove points with "no weight"
    idxA <- which(resultsA[,"ci.ub"] - resultsA[,"ci.lb"] > 39)
    idxB <- which(resultsB[,"ci.ub"] - resultsB[,"ci.lb"] > 39)
    
    resultsA[idxA, c("yi", "ci.lb", "ci.ub")] <- c(NA,NA,NA)
    resultsB[idxB, c("yi", "ci.lb", "ci.ub")] <- c(NA,NA,NA)
    
    if (is.list(fill)) {
      if (length(fill) != 2) {
        stop("`fill` must be a list of length 2 when provided as a list")
      }
      fill_list <- fill
      
    } else if (is.null(fill)) {
      # NULL → both white
      fill_list <- list("white", "white")
      
    } else {
      # anything else (a single string) → replicate to length 2
      fill_list <- rep(fill, 2)
    }
    
    ## define filling
    resultsA$fill <- fill_list[[1]]
    resultsB$fill <- fill_list[[2]]
    
    ## define grouping
    if(!missing(grouping)) {
      row_labels$subheading <- grouping
    }
    bold <- is.na(row_labels$heading)
    
    ## whether to plot only the estimates
    if(estimates_only){
      idx <- resultsB$bold
      idx[is.na(idx)] <- FALSE
      
      row_labels <- row_labels[idx, ]
      resultsA <- resultsA[idx, ]
      resultsB <- resultsB[idx, ]
      
      
      bold <- bold[resultsB$bold]
    }
    
    forestplot <-
      ckbplotr::forest_plot(
        base_size = base_size,
        panels = list(resultsA, resultsB),
        col.key = "variable", col.estimate = "yi",
        col.lci = "ci.lb", col.uci = "ci.ub",
        col.right.heading = col_right_heading,
        col.heading.space = 1,
        row.labels = row_labels, 
        row.labels.space = row.labels.space,
        xlab = xlab,
        showcode = FALSE,
        digits = digits,
        fill = "fill",
        shape = "shape",
        colour = "colour",
        col.left = col.left,
        col.left.heading = col.left.heading,
        bold.labels = bold,
        pointsize = pointsize,
        envir = environment(),
        left.space = left.space,
        exponentiate = exponentiate,
        xlim = xlim,
        xticks = xticks,
        addaes = list(point = paste0("stroke = stroke/10*((2*1.96)/(log(uci_transformed)-log(lci_transformed)))^(1/2)")), ...
      )
    
    attr(forestplot, "panels") <- list(resultsA, resultsB)
    
    return(forestplot)
  }
