source('basic_function.R')

library(movMF)
library(ggplot2)
library(dplyr)
library(patchwork)

# ================================================================
# ================================================================
## US energy composition
# ================================================================
# ================================================================

df <- read.csv('US_energy.csv', check.names = FALSE)
df[df == "--"] <- 0
period_list <- list()
for (i in 2:dim(df)[2]) {
  current_vec <- as.numeric(df[,i])
  compositional_vec <- rep(0,7)
  #coal
  compositional_vec[1] <- current_vec[1]/sum(current_vec)
  #petroleum
  compositional_vec[2] <- sum(current_vec[2:3])/sum(current_vec)
  #gas
  compositional_vec[3] <- sum(current_vec[4:5])/sum(current_vec)
  #nuclear
  compositional_vec[4] <- current_vec[6]/sum(current_vec)
  #Conventional hydroelectric
  compositional_vec[5] <- current_vec[7]/sum(current_vec)
  #Renewables (wind, geothermal, biomass (total))
  compositional_vec[6] <- sum(current_vec[8:10])/sum(current_vec)
  #solar
  compositional_vec[7] <- sum(current_vec[11:12])/sum(current_vec)
  
  period_list[[i-1]] <- sqrt(compositional_vec)
}


y_matrix <- matrix(0,nrow = length(period_list), ncol = length(period_list[[1]]))
for (i in 1:length(period_list)) {
  y_matrix[i,] <- period_list[[i]]
}

y_matrix <- y_matrix[1:120,] 
energy_dates <- seq.Date(
  from = as.Date("2001-01-01"),
  to = as.Date("2010-12-01"),
  by = "month"
)

remove_outlier <- FALSE
if (remove_outlier){
  y_matrix <- y_matrix[-c(60,61,74),]
  energy_dates <- energy_dates[-c(60,61,74)]
  final_residual <- get(load('energy_AR_residual_no_out.RData'))
} else{
  final_residual <- get(load('energy_AR_residual.RData'))
}


out <- residual_distance_plot(Y = y_matrix,fitted = NULL,space = "sphere", title = '(a) Distance plot for original data', mark_a_minus = FALSE, color_groups = FALSE)
p1 <- out$plot


out <- residual_distance_plot(Y = final_residual,fitted = NULL,space = "sphere", ylab = "", title = '(b) Distance plot for fitted residuals', mark_a_minus = FALSE, color_groups = FALSE)
p2 <- out$plot

p1 + p2

## test 1
reference <- sphere_intrinsic_mean(y_matrix)
n <- nrow(y_matrix)
bandwidth_grid_constants=exp(seq(log(0.1),log(1),length.out=15L))
bandwidth_exponent=1/6
h_grid <- bandwidth_grid_constants*n^(-bandwidth_exponent)
test_result <- residual_test(y_matrix,seq(1:n)/n,muY=reference,space="sphere",
                             method="nw",alpha=0.05,list(bw="CV",bwGrid=h_grid,kernel="gauss",tol=1e-6,maxit=500L),n_randomizations = 100)
test_result$p.value
test_result$reject

## test 2
reference <- sphere_intrinsic_mean(final_residual)
n <- nrow(final_residual)
bandwidth_grid_constants=exp(seq(log(0.1),log(1),length.out=15L))
bandwidth_exponent=1/6
h_grid <- bandwidth_grid_constants*n^(-bandwidth_exponent)
test_result <- residual_test(final_residual,seq(1:n)/n,muY=reference,space="sphere",
                             method="nw",alpha=0.05,list(bw="CV",bwGrid=h_grid,kernel="gauss",tol=1e-6,maxit=500L),n_randomizations = 100)
test_result$p.value
test_result$reject

## outlier detection
kappa <- 1.5
energy_boxplot <- distance_profile_boxplot(
  final_residual, space = "sphere", kappa = kappa, fence = "quartile"
)
energy_original <- residual_distance_data(y_matrix, space = "sphere")
energy_signed <- residual_distance_data(final_residual, space = "sphere")
energy_outlier_indices <- which(energy_boxplot$flagged)
energy_results <- data.frame(
  Index = energy_signed$data$Index,
  YearMonth = format(energy_dates, "%Y-%m"),
  Distance = energy_signed$data$Distance,
  SignedDistance = energy_signed$data$SignedDistance,
  DistanceProfileScore = energy_boxplot$scores,
  Fence = energy_boxplot$fence,
  Outlier = energy_boxplot$flagged
)
energy_outliers <- energy_results[energy_results$Outlier, , drop = FALSE]

energy_theme <- ggplot2::theme_classic(base_size = 10, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 11, face = "plain", hjust = 0.5,
                                       margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 10, face = "plain"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.4),
    axis.ticks = ggplot2::element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.width = grid::unit(1, "cm"),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.spacing.x = grid::unit(0.15, "cm"),
    plot.margin = ggplot2::margin(6, 6, 6, 6)
  )

energy_distance_panel <- function(rd, title, flagged = NULL) {
  panel_data <- rd$data
  p <- ggplot2::ggplot(panel_data, ggplot2::aes(x = Index, y = SignedDistance)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5,
                        linetype = "dashed", colour = "gray50")
  if (is.null(flagged)) {
    p <- p + ggplot2::geom_point(shape = 21, fill = NA, colour = "black",
                                 size = 2, stroke = 0.8, alpha = 1)
  } else {
    if (length(flagged) != nrow(panel_data)) stop("flagged must match the residual rows.")
    panel_data$Status <- factor(ifelse(flagged, "Flagged", "Not flagged"),
                                levels = c("Not flagged", "Flagged"))
    p <- p + ggplot2::geom_point(
      data = panel_data, ggplot2::aes(colour = Status),
      shape = 21, fill = NA, size = 2, stroke = 0.8, alpha = 1
    ) +
      ggplot2::scale_colour_manual(
        values = c("Not flagged" = "black", "Flagged" = "#D7191C"),
        drop = FALSE, name = NULL
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, order = 1))
  }
  p + ggplot2::labs(title = title, x = "Index", y = "Signed distance") + energy_theme
}

p1 <- energy_distance_panel(
  energy_original, "(a) Distance plot for original data"
)
energy_distance_plot <- energy_distance_panel(
  energy_signed, "(b) Distance plot for fitted residuals", energy_boxplot$flagged
)
energy_profile_plot <- plot_distance_profile_boxplot(
  energy_boxplot, title = "(c) Boxplot for fitted residuals", show_legend = TRUE
) + energy_theme +
  ggplot2::guides(
    fill = ggplot2::guide_legend(nrow = 1, order = 2),
    colour = ggplot2::guide_legend(nrow = 1, order = 3),
    linetype = ggplot2::guide_legend(nrow = 1, order = 3)
  )
p2 <- energy_distance_plot
energy_plot <- patchwork::wrap_plots(
  p1, energy_distance_plot, energy_profile_plot,
  nrow = 1, widths = c(1, 1, 1), guides = "collect"
) & ggplot2::theme(legend.position = "bottom")


cat("Number of residuals:", nrow(final_residual), "\n")
cat("Quartile fence:", format(energy_boxplot$fence, digits = 8), "\n")
cat("Number of outliers:", length(energy_outlier_indices), "\n")
cat("Outlier indices:", if (length(energy_outlier_indices)) {
  paste(energy_outlier_indices, collapse = ", ")
} else "none", "\n")
cat("Outlier year-months:", if (nrow(energy_outliers)) {
  paste(energy_outliers$YearMonth, collapse = ", ")
} else "none", "\n")
print(energy_outliers, row.names = FALSE)
print(energy_plot)


## prediction interval
PI_info_list <- get(load("energy_PI_info.RData"))
plot_list <- list()
for (k in 1:3) {
  current_inf <- PI_info_list[[k]]
  fit_residual <- current_inf[[1]]
  e_pred <- current_inf[[2]]
  e_true <- current_inf[[3]]
  PI_input <- objpi_proposed_calibrate(E=fit_residual,space='sphere')
  if (k==1){
    y_axis_title = "Profile quantile"
  } else{ 
    y_axis_title = ""
  }
  out <- objpi_prediction_boxplot(
    object = PI_input,
    E_pred = e_pred, 
    E_future = e_true,
    title = "",
    y_axis_title = y_axis_title
  )
  plot_list[[k]] <- out$plot
}



combined_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]]   +
  plot_layout(nrow = 1, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot


# ================================================================

# ================================================================
# ================================================================
## NYC bike sharing
# ================================================================
# ================================================================

NYC_holiday_calendar <- function(include_actual = TRUE,
                                 include_lincoln = FALSE) {
  holiday_names <- c("New Year's Day", "Martin Luther King Jr. Day",
                     "Presidents' Day", "Memorial Day", "Independence Day",
                     "Labor Day", "Columbus Day", "Election Day", "Veterans Day",
                     "Thanksgiving Day", "Christmas Day")
  calendar <- data.frame(
    Date = as.Date(c(
      "2018-01-01", "2018-01-15", "2018-02-19", "2018-05-28", "2018-07-04",
      "2018-09-03", "2018-10-08", "2018-11-06", "2018-11-12", "2018-11-22",
      "2018-12-25",
      "2019-01-01", "2019-01-21", "2019-02-18", "2019-05-27", "2019-07-04",
      "2019-09-02", "2019-10-14", "2019-11-05", "2019-11-11", "2019-11-28",
      "2019-12-25",
      "2020-01-01", "2020-01-20", "2020-02-17", "2020-05-25", "2020-07-03",
      "2020-09-07", "2020-10-12", "2020-11-03", "2020-11-11", "2020-11-26",
      "2020-12-25"
    )),
    Holiday = rep(holiday_names, 3L),
    Type = "Observed",
    stringsAsFactors = FALSE
  )
  if (include_actual) {
    calendar <- rbind(calendar, data.frame(
      Date = as.Date(c("2018-11-11", "2020-07-04")),
      Holiday = c("Veterans Day", "Independence Day"), Type = "Actual",
      stringsAsFactors = FALSE
    ))
  }
  if (include_lincoln) {
    calendar <- rbind(calendar, data.frame(
      Date = as.Date(c("2018-02-12", "2019-02-12", "2020-02-12")),
      Holiday = "Lincoln's Birthday", Type = "Floating",
      stringsAsFactors = FALSE
    ))
  }
  calendar <- calendar[order(calendar$Date), , drop = FALSE]
  rownames(calendar) <- NULL
  calendar
}

NY_trip_quantiles <- function(y_matrix, x = seq_len(ncol(y_matrix)),
                              p = seq(0, 1, length.out = 101L)) {
  y_matrix <- as.matrix(y_matrix)
  if (!is.numeric(y_matrix) || !nrow(y_matrix) || !ncol(y_matrix) ||
      any(!is.finite(y_matrix)) || any(y_matrix < 0) || any(rowSums(y_matrix) <= 0)) {
    stop("y_matrix must have finite nonnegative rows with positive totals. ",
         "For unavailable dates, subset the matrix and matching dates together first.")
  }
  if (length(x) != ncol(y_matrix) || any(!is.finite(x)) || any(diff(x) <= 0)) {
    stop("x must contain one finite, increasing support value per bin.")
  }
  if (length(p) < 2L || any(!is.finite(p)) || any(diff(p) <= 0) ||
      p[1L] != 0 || tail(p, 1L) != 1) {
    stop("p must be increasing and span [0, 1].")
  }
  row_to_quantile <- function(w) {
    keep <- w > 0
    support <- x[keep]
    cdf <- cumsum(w[keep])
    cdf <- cdf / tail(cdf, 1L)
    support[vapply(p, function(u) which(cdf >= u)[1L], integer(1))]
  }
  Q <- t(vapply(seq_len(nrow(y_matrix)), function(i) {
    row_to_quantile(y_matrix[i, ])
  }, numeric(length(p))))
  dimnames(Q) <- list(rownames(y_matrix), as.character(p))
  Q
}

NY_season <- function(dates) {
  dates <- as.Date(dates)
  if (!length(dates) || anyNA(dates)) {
    stop("Supply valid calendar dates to determine the season.")
  }
  month <- as.integer(format(dates, "%m"))
  season_names <- c("Winter", "Spring", "Summer", "Autumn")
  # New York seasons: Dec-Feb, Mar-May, Jun-Aug, and Sep-Nov.
  factor(season_names[(month %% 12L) %/% 3L + 1L], levels = season_names)
}

NY_seasonal_residual_test <- function(residuals, dates, p,
                                      reference = NULL, alpha = 0.05,
                                      n_randomizations = 100L, seed = NULL) {
  residuals <- as.matrix(residuals)
  n <- nrow(residuals)
  dates <- as.Date(dates)
  if (!is.numeric(residuals) || n < 2L || ncol(residuals) < 2L ||
      any(!is.finite(residuals)) || length(dates) != n || anyNA(dates)) {
    stop("Supply finite residual quantiles and matching calendar dates.")
  }
  season <- NY_season(dates)
  diagnostic_X <- droplevels(season)
  counts <- table(Season = season)
  cell_counts <- as.data.frame(counts, responseName = "N")
  sparse_cells <- cell_counts[cell_counts$N > 0L & cell_counts$N < 5L, , drop = FALSE]
  rownames(sparse_cells) <- NULL
  if (is.null(reference)) reference <- wasserstein_mean(residuals)
  if (!is.null(seed)) set.seed(seed)
  test <- residual_test(
    E = residuals, X = diagnostic_X, muY = reference, space = "wasserstein",
    method = "categorical", alpha = alpha, wasserstein_grid = p,
    return_fit = TRUE, n_randomizations = n_randomizations
  )
  test$data.name <- "Fitted bike residuals grouped by season"
  test$diagnostic_predictor <- "Season"
  test$null_hypothesis <- paste(
    "The conditional Frechet mean of the residuals equals the common",
    "residual reference in every observed season."
  )
  list(test = test, X = diagnostic_X, season = season,
       reference = reference, counts = counts, cell_counts = cell_counts,
       sparse_cells = sparse_cells,
       min_cell_n = min(cell_counts$N[cell_counts$N > 0L]))
}

NYC_categorical_residuals <- function(Q, dates = NULL,
                                      p = seq(0, 1, length.out = ncol(Q)),
                                      holiday_dates = NULL,
                                      alpha = 0.05, seed = 20261012L,
                                      n_randomizations = 100L) {
  needed <- c("residual_test", "metric_sub", "wasserstein_mean")
  if (!all(vapply(needed, exists, logical(1), mode = "function"))) {
    stop("Source the supplied basic_function.R before this function.")
  }
  Q <- as.matrix(Q)
  if (!is.numeric(Q) || nrow(Q) < 2L || ncol(Q) < 2L || any(!is.finite(Q)) ||
      any(Q[, -1L, drop = FALSE] < Q[, -ncol(Q), drop = FALSE])) {
    stop("Q must contain at least two finite, nondecreasing quantile rows.")
  }
  if (length(p) != ncol(Q) || any(!is.finite(p)) || any(diff(p) <= 0) ||
      p[1L] != 0 || tail(p, 1L) != 1) {
    stop("p must match the columns of Q and span [0, 1].")
  }
  if (is.null(dates)) {
    labels <- rownames(Q)
    if (!is.null(labels) && all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", labels))) {
      dates <- as.Date(labels)
    } else if (nrow(Q) == 731L) {
      dates <- seq(as.Date("2019-01-01"), as.Date("2020-12-31"), by = "day")
    } else {
      stop("Supply dates: one actual calendar date per row of Q.")
    }
  }
  dates <- as.Date(dates)
  if (length(dates) != nrow(Q) || anyNA(dates) || anyDuplicated(dates)) {
    stop("dates must contain one valid, distinct date per quantile row.")
  }
  if (is.null(holiday_dates)) {
    if (any(!format(dates, "%Y") %in% c("2018", "2019", "2020"))) {
      stop("The default holiday calendar covers 2018-2020; supply holiday_dates for other years.")
    }
    holiday_dates <- NYC_holiday_calendar()$Date
  }
  holiday_dates <- sort(unique(as.Date(holiday_dates)))
  if (anyNA(holiday_dates)) stop("holiday_dates must contain valid dates.")
  
  weekday_names <- c("Monday", "Tuesday", "Wednesday", "Thursday",
                     "Friday", "Saturday", "Sunday")
  weekday_number <- as.integer(format(dates, "%u"))
  weekday <- weekday_names[weekday_number]
  holiday <- dates %in% holiday_dates
  X <- factor(weekday, levels = weekday_names)
  season <- NY_season(dates)
  dimnames(Q) <- list(format(dates, "%Y-%m-%d"), as.character(p))
  reference <- wasserstein_mean(Q)
  
  # Reproducible randomizations for the original-data and residual tests.
  if (!is.null(seed)) set.seed(seed)
  fit_one <- function(X) {
    test_original <- residual_test(
      E = Q, X = X, muY = reference, space = "wasserstein",
      method = "categorical", alpha = alpha,
      wasserstein_grid = p, return_fit = TRUE,
      n_randomizations = n_randomizations
    )
    # This categorical fit returns means/group, rather than qout or yout.
    category_means <- do.call(rbind, test_original$fit$means)
    rownames(category_means) <- as.character(test_original$fit$categories)
    colnames(category_means) <- as.character(p)
    fitted <- category_means[test_original$fit$group, , drop = FALSE]
    dimnames(fitted) <- dimnames(Q)
    residuals <- t(vapply(seq_len(nrow(Q)), function(i) {
      metric_sub(nu1 = fitted[i, ], nu2 = Q[i, ], nu3 = reference,
                 space = "wasserstein", wasserstein_grid = p)
    }, numeric(ncol(Q))))
    dimnames(residuals) <- dimnames(Q)
    residual_reference <- wasserstein_mean(residuals)
    # The mean model uses day type only; its residual diagnostic uses season alone.
    seasonal_diagnostic <- NY_seasonal_residual_test(
      residuals, dates, p, reference = residual_reference, alpha = alpha,
      n_randomizations = n_randomizations
    )
    test_residuals <- seasonal_diagnostic$test
    # Keep test as an alias for the original-data test for existing callers.
    list(X = X, test = test_original, test_original = test_original,
         test_residuals = test_residuals, residual_reference = residual_reference,
         residual_test_predictor = seasonal_diagnostic$X,
         seasonal_diagnostic = seasonal_diagnostic,
         category_means = category_means, fitted = fitted,
         residuals = residuals, counts = table(X))
  }
  models <- list(day_type = fit_one(X))
  test_summary <- do.call(rbind, lapply(names(models), function(setting) {
    model <- models[[setting]]
    tests <- list(model$test_original, model$test_residuals)
    data.frame(
      Setting = setting, Data = c("Original data", "Fitted residuals"),
      TestPredictor = c("Day type", "Season"),
      N = nrow(Q),
      Groups = c(sum(model$counts > 0L), nlevels(model$residual_test_predictor)),
      MinGroupN = c(min(model$counts[model$counts > 0L]),
                    model$seasonal_diagnostic$min_cell_n),
      Dn = vapply(tests, function(test) test$Dn, numeric(1)),
      Statistic = vapply(tests, function(test) unname(test$statistic), numeric(1)),
      PValue = vapply(tests, function(test) test$p.value, numeric(1)),
      Reject = vapply(tests, function(test) test$reject, logical(1)),
      row.names = NULL
    )
  }))
  list(Q = Q, dates = dates, p = p, reference = reference, season = season,
       holiday_dates = holiday_dates,
       calendar = data.frame(Date = dates, Weekday = weekday, Holiday = holiday,
                             Season = season, DayType = X),
       day_type = models$day_type, test_summary = test_summary)
}

NY_outlier_data <- function(Q, final_residual,
                            p = seq(0, 1, length.out = ncol(Q)),
                            dates = NULL, kappa = 1.5,
                            indices = seq_len(nrow(Q))) {
  Q <- as.matrix(Q)
  final_residual <- as.matrix(final_residual)
  if (!is.numeric(Q) || !is.numeric(final_residual) ||
      !identical(dim(Q), dim(final_residual)) || nrow(Q) < 4L || ncol(Q) < 2L ||
      any(!is.finite(Q)) || any(!is.finite(final_residual))) {
    stop("Q and final_residual must be finite quantile matrices of the same size, ",
         "with at least four matching dates.")
  }
  if (any(Q[, -1L, drop = FALSE] < Q[, -ncol(Q), drop = FALSE]) ||
      any(final_residual[, -1L, drop = FALSE] <
          final_residual[, -ncol(final_residual), drop = FALSE])) {
    stop("Each original and residual quantile row must be nondecreasing.")
  }
  if (length(p) != ncol(Q) || any(!is.finite(p)) || any(diff(p) <= 0) ||
      p[1L] != 0 || tail(p, 1L) != 1) {
    stop("p must match the quantile columns and span [0, 1].")
  }
  if (length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("kappa must be one finite nonnegative value.")
  }
  if (is.null(dates)) {
    labels <- rownames(Q)
    if (!is.null(labels) && all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", labels))) {
      dates <- as.Date(labels)
    } else if (nrow(Q) == 731L) {
      dates <- seq(as.Date("2019-01-01"), as.Date("2020-12-31"), by = "day")
    } else {
      stop("Supply dates: one date per quantile row.")
    }
  }
  dates <- as.Date(dates)
  if (length(dates) != nrow(Q) || anyNA(dates) || anyDuplicated(dates)) {
    stop("dates must contain one valid, distinct date per quantile row.")
  }
  if (!is.numeric(indices) || length(indices) != nrow(Q) ||
      anyNA(indices) || any(!is.finite(indices)) || any(indices < 1L) ||
      any(indices != floor(indices)) || anyDuplicated(indices)) {
    stop("indices must contain one distinct positive original index per row.")
  }
  
  ny_boxplot <- distance_profile_boxplot(
    final_residual, space = "wasserstein", wasserstein_grid = p,
    kappa = kappa, fence = "quartile"
  )
  ny_original <- residual_distance_data(
    Q, space = "wasserstein", wasserstein_grid = p
  )
  ny_signed <- residual_distance_data(
    final_residual, space = "wasserstein", wasserstein_grid = p
  )
  ny_original$data$Index <- indices
  ny_signed$data$Index <- indices
  ny_calendar <- NYC_holiday_calendar()
  ny_results <- data.frame(
    Index = ny_signed$data$Index, Date = dates,
    Holiday = ny_calendar$Holiday[match(dates, ny_calendar$Date)],
    Distance = ny_signed$data$Distance,
    SignedDistance = ny_signed$data$SignedDistance,
    DistanceProfileScore = ny_boxplot$scores,
    Fence = ny_boxplot$fence, Outlier = ny_boxplot$flagged
  )
  list(boxplot = ny_boxplot, original = ny_original, signed = ny_signed,
       outlier_indices = indices[ny_boxplot$flagged], results = ny_results,
       outliers = ny_results[ny_results$Outlier, , drop = FALSE])
}

NY_distance_by_X_data <- function(rd, X) {
  df <- rd$data
  if (!is.data.frame(df) || !all(c("Index", "SignedDistance") %in% names(df))) {
    stop("rd must be returned by residual_distance_data().")
  }
  if (is.null(X) || !is.atomic(X) || !is.null(dim(X)) ||
      length(X) != nrow(df) || anyNA(X)) {
    stop("X must be a category vector/factor with one value per observation.")
  }
  if (!is.factor(X)) {
    labels <- as.character(X)
    weekday_levels <- c("Monday", "Tuesday", "Wednesday", "Thursday",
                        "Friday", "Saturday", "Sunday")
    category_levels <- if (all(labels %in% weekday_levels)) weekday_levels else unique(labels)
    X <- factor(labels, levels = category_levels)
  }
  df$X <- X
  df
}

NY_outlier_figure <- function(Q, final_residual,
                              p = seq(0, 1, length.out = ncol(Q)),
                              dates = NULL, kappa = 1.5, X = NULL,
                              x_axis = c("covariate", "index"),
                              indices = seq_len(nrow(Q)),
                              include_boxplot = TRUE) {
  x_axis <- match.arg(x_axis)
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages ggplot2 and patchwork are required for the figure.")
  }
  
  ny <- NY_outlier_data(
    Q, final_residual, p = p, dates = dates,
    kappa = kappa, indices = indices
  )
  
  if (x_axis == "covariate" || !is.null(X)) {
    ny$original$data <- NY_distance_by_X_data(ny$original, X)
    ny$signed$data <- NY_distance_by_X_data(ny$signed, X)
    ny$results$X <- ny$signed$data$X
  }
  ny$x_axis <- x_axis
  ny$outliers <- ny$results[ny$results$Outlier, , drop = FALSE]
  
  ny_theme <- ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 11, face = "plain", hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      ),
      axis.title = ggplot2::element_text(size = 10, face = "plain"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.width = grid::unit(1, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.spacing.x = grid::unit(0.15, "cm"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  
  ny_distance_panel <- function(rd, title, flagged = NULL) {
    panel_data <- rd$data
    mapping <- if (x_axis == "covariate") {
      ggplot2::aes(x = X, y = SignedDistance)
    } else {
      ggplot2::aes(x = Index, y = SignedDistance)
    }
    point_position <- if (x_axis == "covariate") {
      ggplot2::position_jitter(width = 0.15, height = 0, seed = 123)
    } else {
      "identity"
    }
    
    panel <- ggplot2::ggplot(panel_data, mapping) +
      ggplot2::geom_hline(
        yintercept = 0, linewidth = 0.5,
        linetype = "dashed", colour = "gray50"
      )
    if (is.null(flagged)) {
      panel <- panel + ggplot2::geom_point(
        shape = 21, fill = NA, colour = "black",
        size = 2, stroke = 0.8, alpha = 1,
        position = point_position
      )
    } else {
      if (length(flagged) != nrow(panel_data)) stop("flagged must match the residual rows.")
      panel_data$Status <- factor(ifelse(flagged, "Flagged", "Not flagged"),
                                  levels = c("Not flagged", "Flagged"))
      panel <- panel + ggplot2::geom_point(
        data = panel_data, ggplot2::aes(colour = Status),
        shape = 21, fill = NA, size = 2, stroke = 0.8, alpha = 1,
        position = point_position
      ) +
        ggplot2::scale_colour_manual(
          values = c("Not flagged" = "black", "Flagged" = "#D7191C"),
          drop = FALSE, name = NULL
        ) +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, order = 1))
    }
    panel <- panel +
      ggplot2::labs(
        title = title,
        x = if (x_axis == "covariate") "Day type" else "Index",
        y = "Signed distance"
      ) +
      ny_theme
    
    if (x_axis == "covariate") {
      panel <- panel + ggplot2::scale_x_discrete(drop = FALSE)
      if (nlevels(panel_data$X) > 3L) {
        panel <- panel + ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45, hjust = 1, vjust = 1
          )
        )
      }
    }
    panel
  }
  
  p1 <- ny_distance_panel(
    ny$original, "(a) Distance plot for original data"
  )
  p2 <- ny_distance_panel(
    ny$signed, "(b) Distance plot for fitted residuals",
    flagged = if (include_boxplot) ny$boxplot$flagged else NULL
  )
  if (!include_boxplot) p2 <- p2 + ggplot2::labs(y = NULL)
  ny$panels <- list(original = p1, distance = p2)
  if (include_boxplot) {
    p3 <- plot_distance_profile_boxplot(
      ny$boxplot, title = "(c) Boxplot for fitted residuals", show_legend = TRUE
    ) + ny_theme +
      ggplot2::guides(
        fill = ggplot2::guide_legend(nrow = 1, order = 2),
        colour = ggplot2::guide_legend(nrow = 1, order = 3),
        linetype = ggplot2::guide_legend(nrow = 1, order = 3)
      )
    ny$panels$profile <- p3
  }
  ny$include_boxplot <- include_boxplot
  ny$plot <- patchwork::wrap_plots(
    ny$panels, nrow = 1, widths = rep(1, length(ny$panels)), guides = "collect"
  ) + patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(6, 12, 6, 12))
  ) & ggplot2::theme(legend.position = if (include_boxplot) "bottom" else "none")
  
  ny
}


NY_categorical_fit <- function(Q, X, p, xout = X) {
  Q <- as.matrix(Q)
  if (!is.numeric(Q) || nrow(Q) < 2L || ncol(Q) < 2L ||
      any(!is.finite(Q)) ||
      any(Q[, -1L, drop = FALSE] < Q[, -ncol(Q), drop = FALSE])) {
    stop("Q must contain at least two finite, nondecreasing quantile rows.")
  }
  if (length(p) != ncol(Q) || any(!is.finite(p)) || any(diff(p) <= 0) ||
      p[1L] != 0 || tail(p, 1L) != 1) {
    stop("p must match the columns of Q and span [0, 1].")
  }
  if (length(X) != nrow(Q) || anyNA(X) || !length(xout) || anyNA(xout)) {
    stop("Supply one category per response and valid prediction categories.")
  }
  if (!is.factor(X)) X <- factor(X, levels = unique(X))
  prediction_categories <- match(as.character(xout), levels(X))
  if (anyNA(prediction_categories)) stop("Unknown prediction category.")
  n <- nrow(Q)
  fit <- objpi_regression_fit(
    xin = as.integer(X), Y = Q,
    xout = c(as.integer(X), prediction_categories),
    space = "wasserstein", design = "categorical", wasserstein_grid = p
  )
  fitted <- fit$fitted[seq_len(n), , drop = FALSE]
  predicted <- fit$fitted[n + seq_along(xout), , drop = FALSE]
  dimnames(fitted) <- dimnames(Q)
  colnames(predicted) <- colnames(Q)
  reference <- wasserstein_mean(Q)
  residuals <- objpi_residuals(Q, fitted, center = reference,
                               space = "wasserstein", wasserstein_grid = p)
  dimnames(residuals) <- dimnames(Q)
  list(X = X, reference = reference, fitted = fitted, predicted = predicted,
       residuals = residuals, counts = table(X))
}

NY_prediction_figure <- function(Q, X, p, dates, alpha = 0.05,
                                 seed = 20261013L) {
  day_types <- c("Weekday", "Weekend")
  categories <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
                   "Saturday", "Sunday")
  X <- factor(as.character(X), levels = categories)
  day_group <- ifelse(X %in% categories[1:5], "Weekday", "Weekend")
  eligible <- table(X)[as.character(X)] >= 2L
  if (length(X) != nrow(Q) || length(dates) != nrow(Q) || anyNA(X) ||
      !all(day_types %in% day_group[eligible])) {
    stop("Each illustrated group needs a day-of-week category with at least two observations for leave-one-out prediction.")
  }
  set.seed(seed)
  selected <- vapply(day_types, function(day_type) {
    candidates <- which(day_group == day_type & eligible)
    candidates[sample.int(length(candidates), 1L)]
  }, integer(1))
  plot_list <- fit_list <- calibration_list <- summary_list <-
    vector("list", length(day_types))
  for (k in seq_along(selected)) {
    i <- selected[k]
    training_indices <- setdiff(seq_len(nrow(Q)), i)
    fit <- NY_categorical_fit(Q[training_indices, , drop = FALSE],
                              X[training_indices], p, xout = X[i])
    fit$training_indices <- training_indices
    prediction_residual <- objpi_residuals(
      fit$predicted, fit$predicted, center = fit$reference,
      space = "wasserstein", wasserstein_grid = p
    )
    future_residual <- objpi_residuals(
      Q[i, , drop = FALSE], fit$predicted, center = fit$reference,
      space = "wasserstein", wasserstein_grid = p
    )
    calibrated <- objpi_proposed_calibrate(
      E = fit$residuals, center = fit$reference, space = "wasserstein",
      wasserstein_grid = p, alpha = alpha
    )
    panel <- objpi_prediction_boxplot(
      object = calibrated, E_pred = prediction_residual,
      E_future = future_residual,
      title = day_types[k],
      y_axis_title = if (k == 1L) "Profile quantile" else ""
    )
    error <- objpi_row_distances(Q[i, , drop = FALSE], fit$predicted,
                                 space = "wasserstein", wasserstein_grid = p)
    summary_list[[k]] <- data.frame(
      Index = unname(i), Date = as.Date(dates[i]), DayType = day_types[k],
      Category = as.character(X[i]),
      TrainingN = length(training_indices),
      CategoryTrainingN = sum(X[training_indices] == X[i]),
      PredictionError = error, SquaredPredictionError = error^2,
      MeanDistanceScore = panel$future_score,
      Threshold = calibrated$threshold, InsidePredictionSet = panel$future_inside
    )
    plot_list[[k]] <- panel$plot
    fit_list[[k]] <- fit
    calibration_list[[k]] <- calibrated
  }
  names(plot_list) <- names(fit_list) <- names(calibration_list) <- day_types
  combined <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = grid::unit(3, "line"),
                   legend.spacing.x = grid::unit(0.5, "cm"))
  list(indices = selected, summary = do.call(rbind, summary_list),
       fits = fit_list, calibration = calibration_list,
       panels = plot_list, plot = combined)
}

NY_kfold_prediction_error <- function(Q, X, p, dates, folds,
                                      indices = seq_len(nrow(Q)),
                                      dataset = "Original data") {
  # Fold zero marks observations retained in every training sample only.
  if (!is.factor(X)) X <- factor(X, levels = unique(X))
  n <- nrow(Q)
  if (length(X) != n || length(dates) != n || length(indices) != n ||
      length(folds) != n || anyNA(X) || anyNA(folds) ||
      !is.numeric(folds) || any(!is.finite(folds)) ||
      any(folds < 0L) || any(folds != floor(folds))) {
    stop("Q, X, dates, indices and nonnegative integer fold labels must match.")
  }
  evaluated <- which(folds > 0L)
  K <- as.integer(max(folds))
  if (K < 2L || !identical(sort(unique(as.integer(folds[evaluated]))), seq_len(K))) {
    stop("Validation fold labels must include 1 to K; zero is training only.")
  }
  predicted <- matrix(NA_real_, nrow = n, ncol = ncol(Q), dimnames = dimnames(Q))
  errors <- rep(NA_real_, n)
  for (k in seq_len(K)) {
    validation <- which(folds == k)
    training <- which(folds != k)
    fit <- objpi_regression_fit(
      xin = as.integer(X[training]), Y = Q[training, , drop = FALSE],
      xout = as.integer(X[validation]), space = "wasserstein",
      design = "categorical", wasserstein_grid = p
    )
    predicted[validation, ] <- fit$fitted
    errors[validation] <- objpi_row_distances(
      Q[validation, , drop = FALSE], fit$fitted,
      space = "wasserstein", wasserstein_grid = p
    )
  }
  if (any(!is.finite(errors[evaluated]))) stop("Some validation observations have invalid errors.")
  details <- data.frame(
    Dataset = dataset, Index = indices, Date = as.Date(dates), Category = X,
    Fold = folds, Evaluated = folds > 0L,
    PredictionError = errors, SquaredPredictionError = errors^2
  )
  summarize <- function(rows) {
    used <- rows[folds[rows] > 0L]
    data.frame(
      Dataset = dataset, N = length(used), K = K,
      MeanPredictionError = if (length(used)) mean(errors[used]) else NA_real_,
      MeanSquaredPredictionError = if (length(used)) mean(errors[used]^2) else NA_real_
    )
  }
  by_category <- do.call(rbind, lapply(levels(X), function(category) {
    cbind(Category = category, summarize(which(X == category)))
  }))
  fold_summary <- do.call(rbind, lapply(seq_len(K), function(k) {
    cbind(Fold = k, summarize(which(folds == k)))
  }))
  list(summary = summarize(seq_len(n)), details = details,
       by_category = by_category, fold_summary = fold_summary,
       fitted = predicted, folds = folds)
}

NY_season_stratified_folds <- function(X, dates, K = 10L, seed = 20261013L) {
  # Balance validation counts within season, using predictors only.
  if (length(X) != length(dates) || !length(X) || anyNA(X) ||
      length(K) != 1L || !is.numeric(K) || !is.finite(K) ||
      K < 2L || K != floor(K)) {
    stop("Supply matching categories/dates and an integer K >= 2.")
  }
  X <- if (is.factor(X)) droplevels(X) else factor(X, levels = unique(X))
  season <- NY_season(dates)
  joint <- interaction(X, season, drop = TRUE, sep = " / ")
  cell_n <- as.integer(table(joint)[as.integer(joint)])
  training_only <- cell_n == 1L
  evaluated <- which(!training_only)
  if (!length(evaluated)) stop("No weekday/season combination has two observations for validation.")
  season_n <- table(droplevels(season[evaluated]))
  K <- as.integer(K)
  if (any(season_n < K)) {
    stop("Each evaluated season needs at least K days; K is not reduced automatically.")
  }
  set.seed(seed)
  folds <- integer(length(X))
  for (s in names(season_n)) {
    rows <- which(season == s & !training_only)
    rows <- rows[sample.int(length(rows))]
    folds[rows] <- rep(sample.int(K), length.out = length(rows))
  }
  repairs <- data.frame(IndexA = integer(), IndexB = integer(),
                        Season = character(), FoldA = integer(), FoldB = integer())
  for (g in levels(droplevels(joint[evaluated]))) {
    rows <- which(joint == g)
    if (length(unique(folds[rows])) > 1L) next
    source_fold <- folds[rows[1L]]
    s <- season[rows[1L]]
    candidates <- which(season == s & folds > 0L & folds != source_fold)
    safe <- vapply(candidates, function(j) {
      other_rows <- setdiff(which(joint == joint[j]), j)
      any(folds[other_rows] != source_fold)
    }, logical(1L))
    candidates <- candidates[safe]
    if (!length(candidates)) {
      stop("Cannot preserve season balance and joint-category training coverage.")
    }
    i <- rows[sample.int(length(rows), 1L)]
    j <- candidates[sample.int(length(candidates), 1L)]
    target_fold <- folds[j]
    folds[c(i, j)] <- c(target_fold, source_fold)
    repairs <- rbind(repairs, data.frame(
      IndexA = i, IndexB = j, Season = as.character(s),
      FoldA = source_fold, FoldB = target_fold
    ))
  }
  list(folds = folds, K = K, training_only = which(training_only),
       season_counts = table(Season = season[evaluated], Fold = folds[evaluated]),
       joint_counts = table(Category = droplevels(joint[evaluated]), Fold = folds[evaluated]),
       repairs = repairs)
}

NY_compare_season_prediction <- function(Q, X, p, dates, folds = NULL,
                                         K = 10L, seed = 20261013L,
                                         indices = seq_len(nrow(Q)),
                                         dataset = "Original data") {
  Q <- as.matrix(Q)
  n <- nrow(Q)
  dates <- as.Date(dates)
  if (!is.numeric(Q) || n < 2L || ncol(Q) < 2L || any(!is.finite(Q)) ||
      any(Q[, -1L, drop = FALSE] < Q[, -ncol(Q), drop = FALSE]) ||
      length(X) != n || length(dates) != n || length(indices) != n ||
      anyNA(X) || anyNA(dates) || anyNA(indices) || anyDuplicated(indices)) {
    stop("Q, day types, dates, and distinct original indices must match.")
  }
  if (length(p) != ncol(Q) || any(!is.finite(p)) || any(diff(p) <= 0) ||
      p[1L] != 0 || tail(p, 1L) != 1) {
    stop("p must match the quantile columns and span [0,1].")
  }
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 2L || K != floor(K)) {
    stop("K must be an integer >= 2; it is never reduced automatically.")
  }
  X <- if (is.factor(X)) droplevels(X) else factor(X, levels = unique(X))
  season <- NY_season(dates)
  joint <- interaction(X, season, drop = TRUE, sep = " / ")
  cell_n <- as.integer(table(joint)[as.integer(joint)])
  holiday <- dates %in% NYC_holiday_calendar()$Date
  training_only_flag <- cell_n == 1L
  evaluated <- which(!training_only_flag)
  excluded <- which(training_only_flag)
  if (!length(evaluated)) {
    stop("No weekday/season category has two observations for validation.")
  }
  if (is.null(folds)) {
    plan <- NY_season_stratified_folds(X, dates, K = K, seed = seed)
    folds <- plan$folds
    fold_repairs <- plan$repairs
  } else {
    fold_repairs <- NULL
  }
  if (!is.numeric(folds) || length(folds) != n ||
      anyNA(folds[evaluated]) || any(!is.finite(folds[evaluated])) ||
      any(folds[evaluated] < 1L) ||
      any(folds[evaluated] != floor(folds[evaluated]))) {
    stop("Supply a positive integer fold label for every eligible observation.")
  }
  folds[excluded] <- 0L
  K <- as.integer(K)
  if (!identical(sort(unique(as.integer(folds[evaluated]))), seq_len(K))) {
    stop("Validation folds must include every label from 1 to the requested K.")
  }
  joint_folds <- table(droplevels(joint[evaluated]), folds[evaluated])
  if (any(apply(joint_folds, 1L, max) == rowSums(joint_folds))) {
    stop("A validation category is absent from training; stratify folds by day type and season.")
  }

  base <- NY_kfold_prediction_error(
    Q, X, p, dates, folds = folds, indices = indices, dataset = dataset
  )
  extended <- NY_kfold_prediction_error(
    Q, joint, p, dates, folds = folds, indices = indices, dataset = dataset
  )
  prediction_day <- base$fitted
  prediction_season <- extended$fitted
  error_day <- base$details$PredictionError
  error_season <- extended$details$PredictionError
  details <- data.frame(
    Dataset = dataset, Index = indices, Date = dates, DayType = X,
    Season = season, CellN = cell_n, Holiday = holiday, Evaluated = !training_only_flag,
    TrainingOnly = training_only_flag, Fold = folds,
    ExclusionReason = ifelse(cell_n == 1L, "Singleton combination: training only", NA_character_),
    ErrorDayType = error_day, ErrorDayTypeSeason = error_season,
    SquaredErrorDayType = error_day^2, SquaredErrorDayTypeSeason = error_season^2,
    SquaredErrorReduction = error_day^2 - error_season^2
  )
  summarize <- function(rows) {
    used <- rows[!training_only_flag[rows]]
    base_mse <- if (length(used)) mean(error_day[used]^2) else NA_real_
    season_mse <- if (length(used)) mean(error_season[used]^2) else NA_real_
    data.frame(
      Dataset = dataset, NTotal = length(rows), N = length(used),
      NExcluded = length(rows) - length(used), K = K,
      DayTypeMPE = if (length(used)) mean(error_day[used]) else NA_real_,
      DayTypeSeasonMPE = if (length(used)) mean(error_season[used]) else NA_real_,
      DayTypeMSPE = base_mse, DayTypeSeasonMSPE = season_mse,
      MSPEImprovement = base_mse - season_mse,
      MSPEImprovementPercent = if (is.finite(base_mse) && base_mse > 0)
        100 * (base_mse - season_mse) / base_mse else NA_real_,
      SeasonImprovesMSPE = season_mse < base_mse
    )
  }
  comparison <- summarize(seq_len(n))
  summary <- data.frame(
    Dataset = dataset, Model = c("Day type", "Day type + season"),
    NTotal = n, N = length(evaluated), NExcluded = length(excluded), K = K,
    MeanPredictionError = c(comparison$DayTypeMPE, comparison$DayTypeSeasonMPE),
    MeanSquaredPredictionError = c(comparison$DayTypeMSPE, comparison$DayTypeSeasonMSPE)
  )
  by_season <- do.call(rbind, lapply(levels(droplevels(season)), function(s) {
    cbind(Season = s, summarize(which(season == s)))
  }))
  by_category <- do.call(rbind, lapply(levels(X), function(x) {
    cbind(DayType = x, summarize(which(X == x)))
  }))
  by_fold <- do.call(rbind, lapply(seq_len(K), function(k) {
    cbind(Fold = k, summarize(which(folds == k)))
  }))
  fold_models <- lapply(seq_len(K), function(k) list(
    training_indices = indices[folds != k],
    validation_indices = indices[evaluated[folds[evaluated] == k]],
    engine = "within-category Frechet mean"
  ))
  cell_counts <- as.data.frame(table(DayType = X, Season = season), responseName = "N")
  cell_counts$KeptForTraining <- cell_counts$N > 0L
  cell_counts$UsedForPrediction <- cell_counts$N >= 2L
  training_only <- details[excluded, c("Dataset", "Index", "Date", "DayType",
                                       "Season", "CellN", "Holiday", "TrainingOnly", "ExclusionReason"), drop = FALSE]
  list(summary = summary, comparison = comparison, by_season = by_season,
       by_category = by_category, by_fold = by_fold, details = details,
       training_only = training_only, excluded = training_only,
       cell_counts = cell_counts, evaluated_indices = indices[evaluated],
       predictions = list(day_type = prediction_day, day_type_season = prediction_season),
       day_type_cv = base, day_type_season_cv = extended,
       fold_models = fold_models, folds = folds,
       fold_season_counts = table(Season = season[evaluated], Fold = folds[evaluated]),
       fold_joint_counts = joint_folds, fold_repairs = fold_repairs)
}

NY_season_prediction_plot <- function(comparison, title = NULL, y_title = "MSPE") {
  if (is.null(title)) {
    title <- if (comparison$comparison$Dataset == "Original data")
      "(a) Original data" else "(b) Data after outlier removal"
  }
  values <- rbind(cbind(Season = "Overall", comparison$comparison), comparison$by_season)
  plotting <- rbind(
    data.frame(Season = values$Season, Model = "Day type", MSE = values$DayTypeMSPE),
    data.frame(Season = values$Season, Model = "Day type + season", MSE = values$DayTypeSeasonMSPE)
  )
  plotting$Season <- factor(plotting$Season,
    levels = c("Overall", "Spring", "Summer", "Autumn", "Winter"))
  plotting$Model <- factor(plotting$Model, levels = c("Day type", "Day type + season"))
  ggplot2::ggplot(plotting, ggplot2::aes(x = Season, y = MSE, fill = Model)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8), width = 0.7,
      colour = "black", linewidth = 0.25, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c("Day type" = "#3B6FB6", "Day type + season" = "#D95F02")) +
    ggplot2::scale_y_continuous(limits = c(0, NA),
                              expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = y_title, fill = NULL, title = title) +
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 11, face = "plain", hjust = 0.5,
        margin = ggplot2::margin(b = 8)
      ),
      axis.title = ggplot2::element_text(size = 10, face = "plain"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

y_matrix <- get(load("NY_2019_trip_dist.RData"))
dates <- rownames(y_matrix)
if (is.null(dates)) {
  dates <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day")
  if (length(dates) != nrow(y_matrix)) {
    stop("The 2019 input has no date labels; supply the actual date for each row.")
  }
}
p <- seq(0, 1, length.out = 101L)
Q <- NY_trip_quantiles(y_matrix, x = seq_len(ncol(y_matrix)), p = p)

result <- NYC_categorical_residuals(Q, dates = dates, p = p)
ny_fitted <- result$day_type$fitted
ny_residuals <- result$day_type$residuals
ny_X <- result$day_type$X
dates <- result$dates

# Day-type axes 
distance_axis <- "covariate"

## Original responses: day-type effects. Fitted residuals: remaining seasonal effects.
print(result$test_summary, row.names = FALSE)
diagnostic <- result$day_type$seasonal_diagnostic
cat("\nDay-type fit: season counts for the residual test:\n")
print(diagnostic$counts)
if (nrow(diagnostic$sparse_cells)) {
  cat("Observed groups with fewer than five observations (retained in the test):\n")
  print(diagnostic$sparse_cells, row.names = FALSE)
}

ny_alpha <- 0.05
ny_seed <- 124L
ny_remove_mode <- "all" #first_n"
ny_n_remove <- 30L
ny_K <- 10L


ny_prediction <- NY_prediction_figure(Q, ny_X, p, dates,
                                      alpha = ny_alpha, seed = ny_seed)
cat("\nLeave-one-out prediction examples:\n")
print(ny_prediction$summary, row.names = FALSE)
print(ny_prediction$plot)

## Rank only detected outliers by descending mean distance to other residuals.
ny_eight <- NY_outlier_figure(Q, ny_residuals, p = p, dates = dates, X = ny_X,
                              x_axis = distance_axis)
ny_holiday <- ny_eight$results$Date %in% result$holiday_dates
ny_weekend <- as.integer(format(ny_eight$results$Date, "%u")) >= 6L
ny_holiday_outlier <- ny_eight$results$Outlier & ny_holiday
ny_holiday_nonoutlier <- !ny_eight$results$Outlier & ny_holiday
ny_holiday_outlier_count <- sum(ny_holiday_outlier)
ny_holiday_summary <- data.frame(
  Classification = c("Outlier", "Not outlier"),
  Total = c(ny_holiday_outlier_count, sum(ny_holiday_nonoutlier)),
  Weekdays = c(sum(ny_holiday_outlier & !ny_weekend),
               sum(ny_holiday_nonoutlier & !ny_weekend)),
  Weekends = c(sum(ny_holiday_outlier & ny_weekend),
               sum(ny_holiday_nonoutlier & ny_weekend))
)
ny_ranked_outliers <- ny_eight$outliers[
  order(-ny_eight$outliers$DistanceProfileScore, ny_eight$outliers$Index), , drop = FALSE
]
ny_ranked_outliers$Rank <- seq_len(nrow(ny_ranked_outliers))
ny_ranked_outliers <- ny_ranked_outliers[, c(
  "Rank", "Index", "Date", "Holiday", "X", "DistanceProfileScore", "Fence",
  "Outlier", "Distance", "SignedDistance"
)]
rownames(ny_ranked_outliers) <- NULL
ny_holiday_outliers <- ny_ranked_outliers[
  ny_ranked_outliers$Date %in% result$holiday_dates,
  c("Index", "Date", "Holiday"), drop = FALSE
]
ny_remove_mode <- match.arg(ny_remove_mode, c("first_n", "all"))
ny_removed_indices <- if (ny_remove_mode == "all") {
  ny_ranked_outliers$Index
} else {
  head(ny_ranked_outliers$Index, ny_n_remove)
}
ny_ranked_outliers$Removed <- ny_ranked_outliers$Index %in% ny_removed_indices
ny_removed_outliers <- ny_ranked_outliers[ny_ranked_outliers$Removed, , drop = FALSE]
if (ny_remove_mode == "first_n" && length(ny_removed_indices) < ny_n_remove) {
  warning("Only ", length(ny_removed_indices), " outliers were detected; removing all of them.")
}
cat("\nMonday-Sunday outlier ranking:\n")
cat("Mean distance score: sum of distances to all other residuals / (n - 1).\n")
cat("Quartile fence:", format(ny_eight$boxplot$fence, digits = 8), "\n")
cat("Number of detected outliers:", nrow(ny_ranked_outliers), "\n")
cat("Number of holiday observations flagged by the boxplot:", ny_holiday_outlier_count, "\n")
cat("Holiday observations by boxplot classification:\n")
print(ny_holiday_summary, row.names = FALSE)
cat("Holiday names among detected outliers:\n")
print(ny_holiday_outliers, row.names = FALSE)
cat("Removal option:", ny_remove_mode, "; number removed:", length(ny_removed_indices), "\n")
print(ny_ranked_outliers, row.names = FALSE)
cat("Removed original indices, in decreasing score order:",
    if (length(ny_removed_indices)) paste(ny_removed_indices, collapse = ", ") else "none",
    "\n")
print(ny_removed_outliers, row.names = FALSE)
print(ny_eight$plot)

ny_retained_indices <- setdiff(seq_len(nrow(Q)), ny_removed_indices)
Q_clean <- Q[ny_retained_indices, , drop = FALSE]
ny_X_clean <- ny_X[ny_retained_indices]
dates_clean <- dates[ny_retained_indices]
ny_clean_fit <- NY_categorical_fit(Q_clean, ny_X_clean, p)
ny_fitted_clean <- ny_clean_fit$fitted
ny_residuals_clean <- ny_clean_fit$residuals
ny_clean_figure <- NY_outlier_figure(
  Q_clean, ny_residuals_clean, p = p, dates = dates_clean, X = ny_X_clean,
  x_axis = distance_axis, indices = ny_retained_indices, include_boxplot = FALSE
)

cat("\nSeven-category counts before and after removal:\n")
print(rbind(Original = table(ny_X), AfterRemoval = table(ny_X_clean)))
print(ny_clean_figure$plot)

## Repeat the season-only residual diagnostic after the outlier refit.
ny_clean_seasonal_test <- NY_seasonal_residual_test(
  ny_residuals_clean, dates_clean, p,
  reference = wasserstein_mean(ny_residuals_clean),
  alpha = ny_alpha, seed = ny_seed, n_randomizations = 100L
)
ny_seasonal_tests <- list(
  "Original data" = result$day_type$seasonal_diagnostic,
  "After outlier removal" = ny_clean_seasonal_test
)
ny_seasonal_test_comparison <- do.call(rbind, lapply(names(ny_seasonal_tests), function(dataset) {
  diagnostic <- ny_seasonal_tests[[dataset]]
  test <- diagnostic$test
  data.frame(
    Dataset = dataset, Setting = "day_type", Data = "Fitted residuals",
    TestPredictor = "Season", N = length(diagnostic$X),
    Groups = nlevels(diagnostic$X), MinGroupN = diagnostic$min_cell_n,
    Dn = test$Dn, Statistic = unname(test$statistic),
    PValue = test$p.value, Reject = test$reject, row.names = NULL
  )
}))
cat("\nRemaining seasonal mean effects in the seven-category fitted residuals:\n")
print(ny_seasonal_test_comparison, row.names = FALSE)
cat("\nSeason counts after outlier removal:\n")
print(ny_clean_seasonal_test$counts)
if (nrow(ny_clean_seasonal_test$sparse_cells)) {
  cat("Observed groups with fewer than five observations (retained in the test):\n")
  print(ny_clean_seasonal_test$sparse_cells, row.names = FALSE)
}

## Additional season predictor
ny_season_prediction <- NY_compare_season_prediction(
  Q, ny_X, p, dates, K = ny_K, seed = ny_seed, dataset = "Original data"
)
ny_season_prediction_clean <- NY_compare_season_prediction(
  Q_clean, ny_X_clean, p, dates_clean, K = ny_K, seed = ny_seed,
  indices = ny_retained_indices, dataset = "After outlier removal"
)
stopifnot(ny_season_prediction$comparison$K == 10L,
          ny_season_prediction_clean$comparison$K == 10L)

# Preserve the baseline result objects using the same ten-fold comparison.
ny_cv_original <- ny_season_prediction$day_type_cv
ny_cv_clean <- ny_season_prediction_clean$day_type_cv
ny_mean_prediction_error <- rbind(ny_cv_original$summary, ny_cv_clean$summary)
ny_cv_by_category <- rbind(ny_cv_original$by_category, ny_cv_clean$by_category)
ny_cv_predictions <- rbind(ny_cv_original$details, ny_cv_clean$details)
ny_season_prediction_summary <- rbind(
  ny_season_prediction$summary, ny_season_prediction_clean$summary)
ny_season_prediction_improvement <- rbind(
  ny_season_prediction$comparison, ny_season_prediction_clean$comparison)
ny_season_prediction_by_season <- rbind(
  ny_season_prediction$by_season, ny_season_prediction_clean$by_season)
ny_season_prediction_by_daytype <- rbind(
  ny_season_prediction$by_category, ny_season_prediction_clean$by_category)
ny_season_prediction_details <- rbind(
  ny_season_prediction$details, ny_season_prediction_clean$details)
ny_season_prediction_training_only <- rbind(
  ny_season_prediction$training_only, ny_season_prediction_clean$training_only)
# Backward-compatible name: excluded from evaluation, retained in training.
ny_season_prediction_excluded <- ny_season_prediction_training_only
ny_season_prediction_fold_counts <- list(
  "Original data" = ny_season_prediction$fold_season_counts,
  "After outlier removal" = ny_season_prediction_clean$fold_season_counts)
cat("\nDoes adding season improve day-type point prediction?\n")
print(ny_season_prediction_summary, row.names = FALSE)
print(ny_season_prediction_improvement, row.names = FALSE)
cat("\nSingleton observations retained for training only:\n")
print(ny_season_prediction_training_only, row.names = FALSE)
cat("\nSeason counts in each validation fold:\n")
print(ny_season_prediction_fold_counts)
cat("\nPrediction improvement by season:\n")
print(ny_season_prediction_by_season, row.names = FALSE)
cat("\nPrediction improvement by day type:\n")
print(ny_season_prediction_by_daytype, row.names = FALSE)
ny_season_prediction_figure <- patchwork::wrap_plots(
  NY_season_prediction_plot(ny_season_prediction, title = "(a) Original data"),
  NY_season_prediction_plot(ny_season_prediction_clean,
                            title = "(b) Data after outlier removal", y_title = ""),
  nrow = 1, guides = "collect"
) & ggplot2::theme(legend.position = "bottom")
print(ny_season_prediction_figure)



