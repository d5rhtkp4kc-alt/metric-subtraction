metric_sub <- function(nu1, nu2, nu3,
                 space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                 tol = 1e-10,
                 check = TRUE,
                 wasserstein_grid = seq(0, 1, length.out = 1001)) {
  # Euclidean space
  if (space == "euclidean") {
    if (length(nu1) != length(nu2) ||
        length(nu1) != length(nu3)) {
      stop("nu1, nu2 and nu3 must have the same dimension.")
    }
    return(nu3 + nu2 - nu1)
  }
  

  # Sphere and simplex
  
  if (space %in% c("sphere", "simplex")) {
    nu1 <- as.numeric(nu1)
    nu2 <- as.numeric(nu2)
    nu3 <- as.numeric(nu3)
    
    if (length(nu1) != length(nu2) ||
        length(nu1) != length(nu3)) {
      stop("nu1, nu2 and nu3 must have the same dimension.")
    }
    
    norm_E <- function(x) sqrt(sum(x^2))
    
    if (space == "simplex") {
      nu1 <- sqrt(nu1)
      nu1 <- sqrt(nu1)
      nu2 <- sqrt(nu2)
    } 
    
    d <- length(nu1)
    I <- diag(d)
    
    cos_theta <- sum(nu1 * nu3)
    cos_theta <- max(-1, min(1, cos_theta))
    
    # nu1 = nu3: the transport is the identity
    if (abs(1 - cos_theta) < tol) {
      return(nu2)
    }
    
    # Antipodal case
    if (abs(1 + cos_theta) < tol) {
      stop("nu1 and nu3 are antipodal.")
    }
    
    theta <- acos(cos_theta)
    
    u1 <- nu1
    tmp <- nu3 - cos_theta * nu1
    u2 <- tmp / norm_E(tmp)
    Q <- outer(u2, u1) - outer(u1, u2)
    
    R <- I + sin(theta) * Q + (1 - cos(theta)) * (Q %*% Q)
    
    ans <- drop(R %*% nu2)
    ans <- ans / norm_E(ans)
    if (space == "simplex") {
      ans <- ans^2
    } 
    return(ans)
  }
  
  # SPD matrices 
  if (space == "spd") {
    
    if (!is.matrix(nu1) ||
        !is.matrix(nu2) ||
        !is.matrix(nu3)) {
      stop("For space = 'spd', nu1, nu2 and nu3 must be matrices.")
    }
    
    if (!all(dim(nu1) == dim(nu2)) ||
        !all(dim(nu1) == dim(nu3))) {
      stop("nu1, nu2 and nu3 must have the same dimension.")
    }
    
    if (nrow(nu1) != ncol(nu1)) {
      stop("SPD matrices must be square.")
    }
    
    if (check) {
      if (!isSymmetric(nu1, tol = 1e-8) ||
          !isSymmetric(nu2, tol = 1e-8) ||
          !isSymmetric(nu3, tol = 1e-8)) {
        stop("nu1, nu2 and nu3 must be symmetric.")
      }
    }

    L1 <- tryCatch(
      t(chol(nu1)),
      error = function(e) stop("nu1 is not positive definite.")
    )
    
    L2 <- tryCatch(
      t(chol(nu2)),
      error = function(e) stop("nu2 is not positive definite.")
    )
    
    L3 <- tryCatch(
      t(chol(nu3)),
      error = function(e) stop("nu3 is not positive definite.")
    )
    
    L4 <- matrix(
      0,
      nrow = nrow(L1),
      ncol = ncol(L1)
    )
    
    ind_lower <- lower.tri(L1)

    L4[ind_lower] <- L3[ind_lower] + L2[ind_lower] - L1[ind_lower]

    diag(L4) <- exp(log(diag(L3)) + log(diag(L2)) -  log(diag(L1)))
    S4 <- tcrossprod(L4)
    S4 <- (S4 + t(S4)) / 2
    return(S4)
  }
  

  # Wasserstein space
  if (space == "wasserstein") {
    
    is_fun <- c(
      is.function(nu1),
      is.function(nu2),
      is.function(nu3)
    )
    
    # Case 1: quantile-function inputs
    
    if (all(is_fun)) {
      
      p_grid <- wasserstein_grid
      
      if (any(p_grid < 0) || any(p_grid > 1)) {
        stop("wasserstein_grid must lie in [0,1].")
      }
      
      q1_grid <- nu1(p_grid)
      q2_grid <- nu2(p_grid)
      q3_grid <- nu3(p_grid)
      
      if (check) {
        
        if (any(!is.finite(q1_grid)) ||
            any(!is.finite(q2_grid)) ||
            any(!is.finite(q3_grid))) {
          stop(
            "Quantile functions must be finite on wasserstein_grid."
          )
        }
        
        if (any(diff(q1_grid) < -tol) ||
            any(diff(q2_grid) < -tol) ||
            any(diff(q3_grid) < -tol)) {
          stop(
            "Each input quantile function must be nondecreasing."
          )
        }
      }
      
      keep <- !duplicated(q1_grid, fromLast = TRUE)
      
      q1_unique <- q1_grid[keep]
      p_unique  <- p_grid[keep]
      
      if (length(q1_unique) < 2) {
        stop("nu1 does not provide an invertible quantile function.")
      }
      
      F1 <- function(x) {
        
        # F1(x) = 0 below the support and 1 above the support
        out <- approx(
          x = q1_unique,
          y = p_unique,
          xout = x,
          method = "linear",
          rule = 2,
          ties = "ordered"
        )$y
        
        pmin(1, pmax(0, out))
      }
      
      # Transported quantile
      q4 <- function(p) {
        
        if (any(p < 0) || any(p > 1)) {
          stop("p must lie in [0,1].")
        }
        
        z <- nu3(p)
        prob <- F1(z)
        
        nu2(prob)
      }
      
      if (check) {
        
        q4_grid <- q4(p_grid)
        
        if (any(diff(q4_grid) < -tol)) {
          stop(
            paste(
              "The transported quantile function is not",
              "nondecreasing on wasserstein_grid."
            )
          )
        }
      }
      
      return(q4)
    }
    
    
    
    # Case 2: discretized quantile vectors
    
    if (!any(is_fun)) {
      
      q1 <- as.numeric(nu1)
      q2 <- as.numeric(nu2)
      q3 <- as.numeric(nu3)
      
      if (length(q1) != length(q2) ||
          length(q1) != length(q3)) {
        stop(
          paste(
            "The three discretized quantile functions",
            "must have the same length."
          )
        )
      }
      

      if (length(wasserstein_grid) == length(q1)) {
        p_grid <- wasserstein_grid
      } else {
        p_grid <- seq(0, 1, length.out = length(q1))
      }
      
      if (check) {
        
        if (any(!is.finite(q1)) ||
            any(!is.finite(q2)) ||
            any(!is.finite(q3))) {
          stop("Quantile vectors must contain finite values.")
        }
        
        if (any(diff(q1) < -tol) ||
            any(diff(q2) < -tol) ||
            any(diff(q3) < -tol)) {
          stop(
            "Each input quantile vector must be nondecreasing."
          )
        }
      }
      
      keep <- !duplicated(q1, fromLast = TRUE)
      
      q1_unique <- q1[keep]
      p_unique  <- p_grid[keep]
      
      if (length(q1_unique) < 2) {
        stop("nu1 does not provide an invertible quantile function.")
      }
      
      F1_q3 <- approx(
        x = q1_unique,
        y = p_unique,
        xout = q3,
        method = "linear",
        rule = 2,
        ties = "ordered"
      )$y
      
      F1_q3 <- pmin(1, pmax(0, F1_q3))
      
      q4 <- approx(
        x = p_grid,
        y = q2,
        xout = F1_q3,
        method = "linear",
        rule = 2,
        ties = "ordered"
      )$y
      
      if (check && any(diff(q4) < -tol)) {
        stop(
          paste(
            "The transported quantile function",
            "Q2(F1(Q3(u))) is not nondecreasing."
          )
        )
      }
      
      return(q4)
    }
    
    
    stop(
      paste(
        "For Wasserstein space, either all three inputs must be",
        "quantile functions or all three must be numeric quantile vectors."
      )
    )
  }
  
  stop(
    "Unknown space. Use 'euclidean', 'sphere', 'spd', 'wasserstein', or 'simplex."
  )
}


norm_vec <- function(x) sqrt(sum(x^2))
normalize <- function(x) x / norm_vec(x)
clamp <- function(x, lower = -1, upper = 1) pmin(pmax(x, lower), upper)

## sphere

sphere_dist <- function(x, y) acos(clamp(sum(x * y)))

sphere_log <- function(x, mu, tol = 1e-10) {
  inner <- clamp(sum(x * mu))
  theta <- acos(inner)
  if (theta < tol) return(rep(0, length(mu)))
  tangent <- x - inner * mu
  theta * tangent / norm_vec(tangent)
}

sphere_exp <- function(v, mu, tol = 1e-10) {
  r <- norm_vec(v)
  if (r < tol) return(mu)
  normalize(cos(r) * mu + sin(r) * v / r)
}

sphere_intrinsic_mean <- function(Y, tol = 1e-10, max_iter = 500) {
  mu <- normalize(colMeans(Y))
  for (iter in seq_len(max_iter)) {
    log_mat <- t(apply(Y, 1, function(y) sphere_log(y, mu)))
    mean_log <- colMeans(log_mat)
    if (norm_vec(mean_log) < tol) break
    mu <- sphere_exp(mean_log, mu)
  }
  normalize(mu)
}

## Wasserstein

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

wasserstein_dist <- function(q1, q2, p_grid) {
  sqrt(trapz(p_grid, (q1 - q2)^2))
}

wasserstein_mean <- function(Q) colMeans(Q)




.validate_space <- function(space) {
  match.arg(space, c("euclidean", "sphere", "spd", "wasserstein", "simplex"))
}


# Euclidean/sphere/simplex/Wasserstein:
#   n x p matrix (one object per row) or a list of vectors.
# SPD:
#   list of p x p matrices or p x p x n array.
as_object_list <- function(X, space) {
  space <- .validate_space(space)
  
  if (space == "spd") {
    if (is.list(X)) return(X)
    if (length(dim(X)) == 3L && dim(X)[1] == dim(X)[2]) {
      return(lapply(seq_len(dim(X)[3]), function(i) X[, , i]))
    }
    stop("For space = 'spd', X must be a list of SPD matrices or a p x p x n array.")
  }
  
  if (is.list(X)) return(lapply(X, as.numeric))
  
  if (is.vector(X)) X <- matrix(X, nrow = 1)
  X <- as.matrix(X)
  lapply(seq_len(nrow(X)), function(i) as.numeric(X[i, ]))
}

# Cholesky distance for the SPD 
spd_cholesky_dist <- function(S1, S2) {
  L1 <- t(chol(S1))
  L2 <- t(chol(S2))
  ind_lower <- lower.tri(L1)
  
  sqrt(
    sum((L1[ind_lower] - L2[ind_lower])^2) +
      sum((log(diag(L1)) - log(diag(L2)))^2)
  )
}

# Fisher--Rao distance on the simplex 
simplex_dist <- function(x, y) {
  x <- pmax(as.numeric(x), 0)
  y <- pmax(as.numeric(y), 0)
  
  if (sum(x) <= 0 || sum(y) <= 0) {
    stop("Simplex observations must have positive total mass.")
  }
  
  x <- x / sum(x)
  y <- y / sum(y)
  
  acos(clamp(sum(sqrt(x * y))))
}

# Distance for the five spaces.
object_dist <- function(x, y,
                        space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                        wasserstein_grid = NULL) {
  space <- .validate_space(space)
  
  switch(
    space,
    euclidean = sqrt(sum((as.numeric(x) - as.numeric(y))^2)),
    sphere = sphere_dist(as.numeric(x), as.numeric(y)),
    spd = spd_cholesky_dist(x, y),
    wasserstein = {
      if (is.null(wasserstein_grid)) {
        wasserstein_grid <- seq(0, 1, length.out = length(x))
      }
      wasserstein_dist(as.numeric(x), as.numeric(y), wasserstein_grid)
    },
    simplex = simplex_dist(x, y)
  )
}

# Fréchet mean under the Cholesky metric.
spd_cholesky_mean <- function(S_list) {
  if (length(S_list) == 0L) stop("S_list must contain at least one SPD matrix.")
  
  L_list <- lapply(S_list, function(S) t(chol(S)))
  p <- nrow(L_list[[1]])
  ind_lower <- lower.tri(L_list[[1]])
  
  L_bar <- matrix(0, p, p)
  L_bar[ind_lower] <- Reduce(
    `+`, lapply(L_list, function(L) L[ind_lower])
  ) / length(L_list)
  
  diag(L_bar) <- exp(
    Reduce(`+`, lapply(L_list, function(L) log(diag(L)))) /
      length(L_list)
  )
  
  S_bar <- tcrossprod(L_bar)
  (S_bar + t(S_bar)) / 2
}

# Sample Fréchet mean for the five supported spaces.
object_frechet_mean <- function(X,
                                space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                                wasserstein_grid = NULL) {
  space <- .validate_space(space)
  X_list <- as_object_list(X, space)
  
  if (space == "spd") return(spd_cholesky_mean(X_list))
  
  M <- do.call(rbind, X_list)
  
  switch(
    space,
    euclidean = colMeans(M),
    sphere = sphere_intrinsic_mean(M),
    wasserstein = wasserstein_mean(M),
    simplex = {
      M <- pmax(M, 0)
      rs <- rowSums(M)
      if (any(rs <= 0)) stop("Simplex observations must have positive row sums.")
      M <- M / rs
      
      z_bar <- sphere_intrinsic_mean(sqrt(M))
      ans <- pmax(z_bar^2, 0)
      ans / sum(ans)
    }
  )
}

# fitted object-valued residuals
make_metric_residuals <- function(Y, fitted,
                                  space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                                  muY = NULL,
                                  wasserstein_grid = NULL,
                                  ...) {
  space <- .validate_space(space)
  Y_list <- as_object_list(Y, space)
  fitted_list <- as_object_list(fitted, space)
  
  if (length(Y_list) != length(fitted_list)) {
    stop("Y and fitted must contain the same number of objects.")
  }
  
  if (is.null(muY)) {
    muY <- object_frechet_mean(Y, space, wasserstein_grid)
  }
  
  E_hat <- lapply(seq_along(Y_list), function(i) {
    grid_i <- wasserstein_grid
    if (space == "wasserstein" && is.null(grid_i)) {
      grid_i <- seq(0, 1, length.out = length(Y_list[[i]]))
    }
    
    metric_sub(
      fitted_list[[i]], Y_list[[i]], muY,
      space = space,
      wasserstein_grid = if (is.null(grid_i)) seq(0, 1, length.out = 1001) else grid_i,
      ...
    )
  })
  
  list(residuals = E_hat, center = muY)
}

# signed residual distances
residual_distance_data <- function(E,
                                   space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                                   center = NULL,
                                   wasserstein_grid = NULL,
                                   a_minus = NULL,
                                   zero_tol = 1e-12,
                                   ...) {
  space <- .validate_space(space)
  E_list <- as_object_list(E, space)
  n <- length(E_list)
  
  if (n == 0L) stop("E must contain at least one residual object.")
  
  if (is.null(center)) {
    center <- object_frechet_mean(E, space, wasserstein_grid)
  }
  
  d_center <- vapply(
    E_list,
    function(e) object_dist(e, center, space, wasserstein_grid),
    numeric(1)
  )
  
  if (is.null(a_minus)) {
    j_minus <- which.max(d_center^2)
    a_minus <- E_list[[j_minus]]
  } else {
    j_minus <- NA_integer_
  }
  
  grid_use <- wasserstein_grid
  if (space == "wasserstein" && is.null(grid_use)) {
    grid_use <- seq(0, 1, length.out = length(E_list[[1]]))
  }
  
  a_plus <- metric_sub(
    a_minus, center, center,
    space = space,
    wasserstein_grid = if (is.null(grid_use)) seq(0, 1, length.out = 1001) else grid_use,
    ...
  )
  
  d_minus <- vapply(
    E_list,
    function(e) object_dist(e, a_minus, space, wasserstein_grid),
    numeric(1)
  )
  
  d_plus <- vapply(
    E_list,
    function(e) object_dist(e, a_plus, space, wasserstein_grid),
    numeric(1)
  )
  
  score <- d_minus^2 - d_plus^2
  sign_value <- ifelse(
    score < -zero_tol, -1,
    ifelse(score > zero_tol, 1, 0)
  )
  
  signed_distance <- sign_value * d_center
  
  group <- factor(
    ifelse(
      sign_value < 0, "Negative",
      ifelse(sign_value > 0, "Positive", "Boundary")
    ),
    levels = c("Negative", "Boundary", "Positive")
  )
  
  radius_minus <- object_dist(a_minus, center, space, wasserstein_grid)
  radius_plus <- object_dist(a_plus, center, space, wasserstein_grid)
  
  if (space == "wasserstein" && abs(radius_plus - radius_minus) > 1e-8) {
    warning(
      "For Wasserstein space, transport-defined a_plus need not be equidistant ",
      "from the center. Inspect radius_difference before interpreting the sign."
    )
  }
  
  out <- data.frame(
    Index = seq_len(n),
    Distance = d_center,
    Score = score,
    Sign = sign_value,
    SignedDistance = signed_distance,
    Group = group
  )
  
  structure(
    list(
      data = out,
      center = center,
      a_minus = a_minus,
      a_plus = a_plus,
      j_minus = j_minus,
      radius_minus = radius_minus,
      radius_plus = radius_plus,
      radius_difference = radius_plus - radius_minus,
      space = space
    ),
    class = "residual_distance_result"
  )
}

# Plot 
plot_residual_distance <- function(x,
                                   ordering = NULL,
                                   xlab = "Index",
                                   ylab = "Signed distance",
                                   title = NULL,
                                   mark_a_minus = TRUE,
                                   color_groups = TRUE,
                                   highlight_indices = NULL) {
  if (!inherits(x, "residual_distance_result")) {
    stop("x must be an object returned by residual_distance_data().")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  if (!is.logical(color_groups) || length(color_groups) != 1L ||
      is.na(color_groups)) {
    stop("color_groups must be TRUE or FALSE.")
  }
  
  df <- x$data
  
  if (!is.null(highlight_indices)) {
    if (!is.numeric(highlight_indices) || is.complex(highlight_indices) ||
        !is.null(dim(highlight_indices)) ||
        anyNA(highlight_indices) || any(!is.finite(highlight_indices)) ||
        any(highlight_indices != floor(highlight_indices)) ||
        any(highlight_indices < 1 | highlight_indices > nrow(df))) {
      stop("highlight_indices must contain integer indices between 1 and nrow(x$data).")
    }
    highlight_indices <- unique(as.integer(highlight_indices))
  }
  
  if (!is.null(ordering)) {
    if (length(ordering) != nrow(df)) {
      stop("ordering must have the same length as the number of residuals.")
    }
    df$Ordering <- ordering
    xvar <- "Ordering"
  } else {
    xvar <- "Index"
  }
  
  if (is.null(title)) {
    title <- paste0("Residual distance plot: ", x$space)
  }
  
  group_colors <- if (color_groups) {
    c(
      "Negative" = "#2C7BB6",
      "Boundary" = "gray55",
      "Positive" = "#D7191C"
    )
  } else {
    c("Negative" = "black", "Boundary" = "black", "Positive" = "black")
  }
  highlight_color <- "#7B3294"
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[xvar]], y = SignedDistance, color = Group
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.6,
      linetype = "dashed",
      color = "gray50"
    ) +
    ggplot2::geom_point(
      shape = 21,
      fill = NA,
      size = 2.6,
      stroke = 0.8,
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  # Draw highlights over the regular points with a fixed outline and fill.
  # This layer works independently of color_groups and the group color scale.
  if (length(highlight_indices) > 0L) {
    p <- p + ggplot2::geom_point(
      data = df[highlight_indices, , drop = FALSE],
      ggplot2::aes(x = .data[[xvar]], y = SignedDistance),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.2,
      stroke = 0.8,
      color = highlight_color,
      fill = highlight_color,
      alpha = 1,
      show.legend = FALSE
    )
  }
  
  if (mark_a_minus && !is.na(x$j_minus)) {
    marker_color <- if (x$j_minus %in% highlight_indices) {
      highlight_color
    } else {
      "#008837"
    }
    p <- p + ggplot2::geom_point(
      data = df[x$j_minus, , drop = FALSE],
      ggplot2::aes(x = .data[[xvar]], y = SignedDistance),
      inherit.aes = FALSE,
      shape = 15,
      size = 3.8,
      color = marker_color
    )
  }
  
  p
}

residual_distance_plot <- function(Y,
                                   fitted = NULL,
                                   space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                                   muY = NULL,
                                   center = NULL,
                                   wasserstein_grid = NULL,
                                   ordering = NULL,
                                   xlab = "Index",
                                   ylab = "Signed distance",
                                   title = NULL,
                                   mark_a_minus = TRUE,
                                   color_groups = TRUE,
                                   highlight_indices = NULL,
                                   ...) {
  space <- .validate_space(space)
  
  if (is.null(fitted)) {
    rd <- residual_distance_data(
      E = Y,
      space = space,
      center = center,
      wasserstein_grid = wasserstein_grid,
      ...
    )
  } else {
    residual_fit <- make_metric_residuals(
      Y = Y,
      fitted = fitted,
      space = space,
      muY = muY,
      wasserstein_grid = wasserstein_grid,
      ...
    )
    
    rd <- residual_distance_data(
      E = residual_fit$residuals,
      space = space,
      center = residual_fit$center,
      wasserstein_grid = wasserstein_grid,
      ...
    )
  }
  
  list(
    result = rd,
    plot = plot_residual_distance(
      rd,
      ordering = ordering,
      xlab = xlab,
      ylab = ylab,
      title = title,
      mark_a_minus = mark_a_minus,
      color_groups = color_groups,
      highlight_indices = highlight_indices
    )
  )
}

## Nadaraya-Watson Frechet regression

.nw_design <- function(x, name) {
  if (is.numeric(x) && is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
      any(!is.finite(x))) {
    stop(name, " must be a nonempty finite numeric vector or matrix.")
  }
  x
}

.nw_bandwidth <- function(bw, p) {
  if (!is.numeric(bw) || !length(bw) || any(!is.finite(bw)) ||
      any(bw <= 0) || !(length(bw) %in% c(1L, p))) {
    stop("bw must be positive and finite, with length 1 or p.")
  }
  rep(bw, length.out = p)
}

.nw_kernel <- function(kernel) {
  if (identical(kernel, "gausvar")) {
    stop("gausvar can be negative; use a nonnegative kernel for standard NW regression.")
  }
  match.arg(kernel, c("gauss", "epan", "rect", "quar"))
}

.nw_option_kernel <- function(optns) {
  # Accept the original spherical name 'kernel' and Wasserstein name 'ker'.
  kernel <- optns$kernel
  if (is.null(kernel)) kernel <- optns$ker
  if (is.null(kernel)) kernel <- "gauss"
  .nw_kernel(kernel)
}

.nw_weights <- function(X, Xout, bw, kernel, exclude = NULL) {
  bw <- .nw_bandwidth(bw, ncol(X))
  W <- matrix(0, nrow = nrow(Xout), ncol = nrow(X))
  for (j in seq_len(nrow(Xout))) {
    Z <- sweep(sweep(X, 2L, Xout[j, ], "-"), 2L, bw, "/")
    if (kernel == "gauss") {
      logw <- -0.5 * rowSums(Z^2)
    } else {
      if (kernel == "rect") {
        K <- 1 * (abs(Z) <= 1)
      } else {
        K <- pmax(1 - Z^2, 0)
        if (kernel == "quar") K <- K^2
      }
      logw <- rowSums(log(K))
    }
    # Exclude the held-out observation before rescaling Gaussian weights.
    if (!is.null(exclude)) logw[exclude[j]] <- -Inf
    a <- max(logw)
    if (!is.finite(a)) {
      stop("No positive, numerically representable kernel weights at prediction row ",
           j, ". Increase bw or change the prediction point/kernel.")
    }
    w <- exp(logw - a)
    W[j, ] <- w / sum(w)
  }
  W
}

NWFrechetWeights <- function(xin, xout = xin, bw, kernel = "gauss") {
  X <- .nw_design(xin, "xin")
  Xout <- .nw_design(xout, "xout")
  if (ncol(X) != ncol(Xout)) stop("xin and xout must have the same dimension.")
  .nw_weights(X, Xout, bw, .nw_kernel(kernel))
}

.nw_controls <- function(optns) {
  controls <- list(tol = 1e-7, maxit = 500L, nStarts = 1L)
  for (nm in names(controls)) {
    if (!is.null(optns[[nm]])) controls[[nm]] <- optns[[nm]]
    v <- controls[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a positive finite scalar.")
    }
  }
  for (nm in c("maxit", "nStarts")) {
    if (controls[[nm]] != floor(controls[[nm]])) stop(nm, " must be an integer.")
  }
  controls
}

.nw_sphere_stats <- function(y, Y, w, gradient = TRUE) {
  cprod <- pmax(-1, pmin(1, drop(Y %*% y)))
  V <- Y - tcrossprod(cprod, y)
  s <- sqrt(rowSums(V^2))
  theta <- atan2(s, cprod)
  value <- sum(w * theta^2)
  if (!gradient) return(value)
  L <- V
  regular <- s > 1e-12
  L[regular, ] <- V[regular, , drop = FALSE] * (theta[regular] / s[regular])
  antipodal <- !regular & cprod < 0
  if (any(antipodal)) {
    u <- numeric(length(y))
    u[which.min(abs(y))] <- 1
    u <- u - sum(u * y) * y
    u <- u / sqrt(sum(u^2))
    L[antipodal, ] <- tcrossprod(rep(pi, sum(antipodal)), u)
  }
  v <- colSums(L * w)
  v <- v - sum(v * y) * y
  list(value = value, direction = v)  # Riemannian gradient is -2 * v.
}

.nw_sphere_descent <- function(start, Y, w, controls) {
  y <- start / sqrt(sum(start^2))
  converged <- FALSE
  iterations <- 0L
  for (it in seq_len(controls$maxit)) {
    cur <- .nw_sphere_stats(y, Y, w)
    v <- cur$direction
    nv <- sqrt(sum(v^2))
    if (nv <= controls$tol) {
      converged <- TRUE
      break
    }
    step <- min(1, (pi / 2) / nv)
    accepted <- FALSE
    for (backtrack in seq_len(50L)) {
      candidate <- cos(step * nv) * y + sin(step * nv) * v / nv
      candidate <- candidate / sqrt(sum(candidate^2))
      fnew <- .nw_sphere_stats(candidate, Y, w, gradient = FALSE)
      if (fnew <= cur$value - 2e-4 * step * nv^2) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    y <- candidate
    iterations <- it
  }
  final <- .nw_sphere_stats(y, Y, w)
  converged <- converged || sqrt(sum(final$direction^2)) <= controls$tol
  list(y = y, value = final$value, converged = converged, iterations = iterations)
}

.nw_spherical_mean <- function(Y, w, controls) {
  active <- w > 0
  Y <- Y[active, , drop = FALSE]
  w <- w[active]
  z <- colSums(Y * w)
  nz <- sqrt(sum(z^2))
  start <- if (nz > 1e-12) z / nz else Y[which.max(w), ]
  starts <- matrix(start, nrow = 1L)
  if (controls$nStarts > 1L) {
    idx <- head(order(w, decreasing = TRUE), controls$nStarts - 1L)
    starts <- rbind(starts, Y[idx, , drop = FALSE])
  }
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    fit <- .nw_sphere_descent(starts[i, ], Y, w, controls)
    if (is.null(best) || fit$value < best$value ||
        (fit$value == best$value && fit$converged && !best$converged)) best <- fit
  }
  best
}

.nw_bw_grid <- function(bwGrid, p) {
  if (is.null(bwGrid)) stop("For bw = 'CV', supply optns$bwGrid.")
  if (is.numeric(bwGrid) && is.null(dim(bwGrid)) && p == 1L) {
    bwGrid <- matrix(bwGrid, ncol = 1L)
  }
  if (!is.matrix(bwGrid) || !is.numeric(bwGrid) || !nrow(bwGrid) ||
      ncol(bwGrid) != p || any(!is.finite(bwGrid)) || any(bwGrid <= 0)) {
    stop("bwGrid must contain positive bandwidths, with p columns for p predictors.")
  }
  bwGrid
}

.nw_cv <- function(X, bwGrid, kernel, loss) {
  if (nrow(X) < 2L) stop("Leave-one-out CV requires at least two observations.")
  grid <- .nw_bw_grid(bwGrid, ncol(X))
  scores <- rep(Inf, nrow(grid))
  for (g in seq_len(nrow(grid))) {
    scores[g] <- tryCatch({
      W <- .nw_weights(X, X, grid[g, ], kernel, exclude = seq_len(nrow(X)))
      loss(W)
    }, error = function(e) Inf)
  }
  if (!any(is.finite(scores))) {
    stop("All CV candidates failed. Check bwGrid (compact kernels need neighbors) and optimizer controls.")
  }
  idx <- which.min(scores)
  list(bw = as.numeric(grid[idx, ]), bandwidths = grid, scores = scores)
}

NWSpheReg <- function(xin, yin, xout = xin, optns = list()) {
  X <- .nw_design(xin, "xin")
  Xout <- .nw_design(xout, "xout")
  if (ncol(X) != ncol(Xout)) stop("xin and xout must have the same dimension.")
  if (!is.matrix(yin) || !is.numeric(yin) || nrow(yin) != nrow(X) ||
      ncol(yin) < 2L || any(!is.finite(yin))) {
    stop("yin must be a finite n by d numeric matrix, with d >= 2.")
  }
  norms <- sqrt(rowSums(yin^2))
  if (any(!is.finite(norms)) || any(norms == 0)) stop("yin contains an invalid or zero row.")
  if (any(abs(norms - 1) > 1e-6)) warning("Rows of yin have been normalized to unit length.")
  Y <- yin / norms
  kernel <- .nw_option_kernel(optns)
  controls <- .nw_controls(optns)
  bw <- optns$bw
  cv <- NULL
  if (identical(bw, "CV")) {
    cv <- .nw_cv(X, optns$bwGrid, kernel, function(W) {
      loss <- numeric(nrow(Y))
      for (i in seq_len(nrow(Y))) {
        fit <- .nw_spherical_mean(Y, W[i, ], controls)
        if (!fit$converged) return(Inf)
        loss[i] <- .nw_sphere_stats(fit$y, Y[i, , drop = FALSE], 1, gradient = FALSE)
      }
      mean(loss)
    })
    bw <- cv$bw
  }
  bw <- .nw_bandwidth(bw, ncol(X))
  W <- .nw_weights(X, Xout, bw, kernel)
  fits <- lapply(seq_len(nrow(W)), function(j) .nw_spherical_mean(Y, W[j, ], controls))
  yout <- do.call(rbind, lapply(fits, function(fit) fit$y))
  converged <- vapply(fits, function(fit) fit$converged, logical(1))
  if (any(!converged)) {
    warning("Some spherical fits did not meet the gradient tolerance; check $converged or increase maxit/nStarts.")
  }
  optns$bw <- bw
  optns$kernel <- kernel
  list(xout = xout, yout = yout, xin = xin, yin = Y, optns = optns,
       weights = W, effective_n = 1 / rowSums(W^2), converged = converged,
       objective = vapply(fits, function(fit) fit$value, numeric(1)),
       iterations = vapply(fits, function(fit) as.integer(fit$iterations), integer(1)),
       cv = cv)
}

NWWassReg <- function(xin, qin, xout = xin, optns = list()) {
  # Each row of qin is a quantile function on the same probability grid.
  X <- .nw_design(xin, "xin")
  Xout <- .nw_design(xout, "xout")
  if (ncol(X) != ncol(Xout)) stop("xin and xout must have the same dimension.")
  if (!is.matrix(qin) || !is.numeric(qin) || nrow(qin) != nrow(X) ||
      ncol(qin) < 2L || any(!is.finite(qin))) {
    stop("qin must be a finite n by m numeric matrix, with m >= 2.")
  }
  if (any(qin[, -1L, drop = FALSE] < qin[, -ncol(qin), drop = FALSE])) {
    stop("Every row of qin must be nondecreasing on the common probability grid.")
  }
  qSup <- optns$qSup
  if (is.null(qSup)) qSup <- seq(0, 1, length.out = ncol(qin))
  if (!is.numeric(qSup) || length(qSup) != ncol(qin) ||
      any(!is.finite(qSup)) || any(diff(qSup) <= 0) ||
      any(qSup < 0 | qSup > 1)) stop("qSup must be a strictly increasing probability grid matching qin.")
  lower <- if (is.null(optns$lower)) -Inf else optns$lower
  upper <- if (is.null(optns$upper)) Inf else optns$upper
  if (!is.numeric(lower) || !is.numeric(upper) ||
      length(lower) != 1L || length(upper) != 1L ||
      is.na(lower) || is.na(upper) || lower > upper ||
      lower == Inf || upper == -Inf) stop("Invalid lower/upper bounds.")
  bounded <- function(Q) pmax(pmin(Q, upper), lower)
  kernel <- .nw_option_kernel(optns)
  bw <- optns$bw
  cv <- NULL
  if (identical(bw, "CV")) {
    cv <- .nw_cv(X, optns$bwGrid, kernel, function(W) {
      E <- (bounded(W %*% qin) - qin)^2
      mean(drop((E[, -ncol(E), drop = FALSE] + E[, -1L, drop = FALSE]) %*%
                  (diff(qSup) / 2)))
    })
    bw <- cv$bw
  }
  bw <- .nw_bandwidth(bw, ncol(X))
  W <- .nw_weights(X, Xout, bw, kernel)
  optns$bw <- bw
  optns$ker <- kernel
  optns$qSup <- qSup
  list(xout = xout, qout = bounded(W %*% qin), qSup = qSup,
       xin = xin, qin = qin, optns = optns,
       weights = W, effective_n = 1 / rowSums(W^2), cv = cv)
}



## test

residual_test <- function(E, X, muY,
                          space = c("euclidean", "sphere", "spd",
                                    "wasserstein", "simplex"),
                          method = c("nw", "categorical"),
                          alpha = 0.05, optns = list(),
                          wasserstein_grid = NULL,
                          check = TRUE, return_fit = FALSE,
                          n_randomizations = 50L, cauchy_weights = NULL) {
  data_name <- paste(deparse(substitute(E)), collapse = "")
  space <- match.arg(space)
  method <- match.arg(method)
  if (!is.list(optns)) stop("optns must be a list.")
  required <- c("as_object_list", "object_dist", "object_frechet_mean")
  if (method == "nw") {
    required <- c(required, ".nw_design")
    if (space %in% c("sphere", "simplex")) {
      required <- c(required, "NWSpheReg")
    } else if (space == "wasserstein") {
      required <- c(required, "NWWassReg")
    } else {
      required <- c(required, ".nw_option_kernel", "NWFrechetWeights")
      if (identical(optns$bw, "CV")) required <- c(required, ".nw_cv")
    }
  }
  missing_helpers <- required[
    !vapply(required, exists, logical(1), mode = "function")
  ]
  if (length(missing_helpers)) {
    stop("Source your basic_function.R first. Missing: ",
         paste(missing_helpers, collapse = ", "))
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  flags <- list(check = check, return_fit = return_fit)
  if (any(!vapply(flags, function(z) {
    is.logical(z) && length(z) == 1L && !is.na(z)
  }, logical(1)))) {
    stop("check and return_fit must be TRUE or FALSE.")
  }
  if (method == "nw") {
    X <- .nw_design(X, "X")
    n <- nrow(X)
    if (is.null(optns$bw)) stop("Supply optns$bw for NW regression.")
  } else {
    if (!is.null(dim(X)) || !is.atomic(X) || anyNA(X)) {
      stop("For categorical regression, X must be a category vector/factor without NA.")
    }
    n <- length(X)
  }
  if (n < 2L) stop("At least two residual observations are required.")
  if (!is.numeric(n_randomizations) || length(n_randomizations) != 1L ||
      !is.finite(n_randomizations) || n_randomizations < 1 ||
      n_randomizations != floor(n_randomizations) ||
      n_randomizations > .Machine$integer.max) {
    stop("n_randomizations must be a positive integer.")
  }
  L <- as.integer(n_randomizations)
  if (is.null(cauchy_weights)) cauchy_weights <- rep(1, L)
  if (!is.numeric(cauchy_weights) || length(cauchy_weights) != L ||
      any(!is.finite(cauchy_weights)) || any(cauchy_weights < 0) ||
      !any(cauchy_weights > 0)) {
    stop("cauchy_weights must be a nonnegative finite vector of length n_randomizations with positive total weight.")
  }
  cauchy_weights <- cauchy_weights / max(cauchy_weights)
  cauchy_weights <- cauchy_weights / sum(cauchy_weights)
  if (space == "euclidean" && is.numeric(E) && is.null(dim(E))) {
    E <- matrix(E, ncol = 1L)
  }
  E_list <- as_object_list(E, space)
  if (length(E_list) != n) {
    stop("E and X must have the same number of observations.")
  }
  if (space == "spd") {
    d <- nrow(E_list[[1L]])
    if (check) {
      valid <- vapply(E_list, function(e) {
        is.matrix(e) && is.numeric(e) && all(dim(e) == c(d, d)) &&
          all(is.finite(e)) && isSymmetric(e, tol = 1e-8)
      }, logical(1))
      if (!all(valid) || !is.matrix(muY) || !is.numeric(muY) ||
          !all(dim(muY) == c(d, d)) || any(!is.finite(muY)) ||
          !isSymmetric(muY, tol = 1e-8)) {
        stop("Invalid SPD observations or null center.")
      }
      invisible(lapply(c(E_list, list(muY)), chol))
    }
  } else {
    sizes <- lengths(E_list)
    if (any(sizes != sizes[1L]) || length(muY) != sizes[1L]) {
      stop("All residual objects and muY must have the same dimension.")
    }
    M <- do.call(rbind, E_list)
    muY <- as.numeric(muY)
    if (any(!is.finite(M)) || any(!is.finite(muY))) {
      stop("Objects must be finite.")
    }
    if (space == "sphere") {
      if (ncol(M) < 2L) stop("Sphere observations need at least two coordinates.")
      norms <- sqrt(rowSums(M^2))
      center_norm <- sqrt(sum(muY^2))
      if (any(norms <= 0) || center_norm <= 0) {
        stop("Sphere objects cannot be zero.")
      }
      M <- M / norms
      muY <- muY / center_norm
    }
    if (space == "simplex") {
      if (any(M < 0) || any(muY < 0) || any(rowSums(M) <= 0) ||
          sum(muY) <= 0) {
        stop("Simplex objects must be nonnegative with positive total mass.")
      }
      M <- M / rowSums(M)
      muY <- muY / sum(muY)
    }
    if (space == "wasserstein") {
      if (is.null(wasserstein_grid)) wasserstein_grid <- optns$qSup
      if (is.null(wasserstein_grid)) {
        wasserstein_grid <- seq(0, 1, length.out = ncol(M))
      }
      if (!is.numeric(wasserstein_grid) || ncol(M) < 2L ||
          length(wasserstein_grid) != ncol(M) ||
          any(!is.finite(wasserstein_grid)) ||
          any(diff(wasserstein_grid) <= 0) ||
          abs(wasserstein_grid[1L]) > 1e-12 ||
          abs(tail(wasserstein_grid, 1L) - 1) > 1e-12) {
        stop("wasserstein_grid must match the quantiles and span [0,1].")
      }
      if (check && (any(M[, -1L, drop = FALSE] <
                        M[, -ncol(M), drop = FALSE]) || any(diff(muY) < 0))) {
        stop("Quantile vectors must be nondecreasing.")
      }
      optns$qSup <- wasserstein_grid
      lower <- if (is.null(optns$lower)) -Inf else optns$lower
      upper <- if (is.null(optns$upper)) Inf else optns$upper
      if (any(muY < lower | muY > upper)) {
        stop("The null center must satisfy the fitted quantile bounds.")
      }
    }
    E_list <- lapply(seq_len(n), function(i) as.numeric(M[i, ]))
  }
  distance2 <- function(a, b) {
    object_dist(a, b, space = space,
                wasserstein_grid = wasserstein_grid)^2
  }
  null_loss <- vapply(E_list, function(e) distance2(muY, e), numeric(1))
  W <- NULL
  bandwidth <- NULL
  if (method == "categorical") {
    categories <- unique(X)
    group <- match(X, categories)
    index <- lapply(seq_along(categories), function(k) which(group == k))
    means <- lapply(index, function(ii) {
      object_frechet_mean(E_list[ii], space = space,
                          wasserstein_grid = wasserstein_grid)
    })
    A_group <- vapply(index, function(ii) mean(null_loss[ii]), numeric(1))
    B_group <- vapply(seq_along(index), function(k) {
      mean(vapply(E_list[index[[k]]], function(e) {
        distance2(means[[k]], e)
      }, numeric(1)))
    }, numeric(1))
    A <- A_group[group]
    B <- B_group[group]
    regression_fit <- list(categories = categories, group = group, means = means)
  } else if (space %in% c("sphere", "simplex")) {
    fit <- NWSpheReg(xin = X, yin = if (space == "sphere") M else sqrt(M),
                     xout = X, optns = optns)
    if (is.null(fit$converged) || anyNA(fit$converged) ||
        any(!fit$converged)) {
      stop("NWSpheReg did not converge at all evaluation points.")
    }
    W <- fit$weights
    A <- drop(W %*% null_loss)
    fitted_objects <- if (space == "sphere") fit$yout else fit$yout^2
    fitted_list <- as_object_list(fitted_objects, space)
    if (length(fitted_list) != n) stop("Unexpected NWSpheReg output dimensions.")
    B <- vapply(seq_len(n), function(i) {
      sum(W[i, ] * vapply(E_list, function(e) {
        distance2(fitted_list[[i]], e)
      }, numeric(1)))
    }, numeric(1))
    regression_fit <- fit
    if (space == "simplex") regression_fit$composition_out <- fitted_objects
    bandwidth <- fit$optns$bw
  } else if (space == "wasserstein") {
    fit <- NWWassReg(xin = X, qin = M, xout = X, optns = optns)
    W <- fit$weights
    A <- drop(W %*% null_loss)
    fitted_centered <- sweep(fit$qout, 2L, muY, "-")
    bounded <- (!is.null(optns$lower) && any(is.finite(optns$lower))) ||
      (!is.null(optns$upper) && any(is.finite(optns$upper)))
    mean_centered <- if (bounded) {
      W %*% sweep(M, 2L, muY, "-")
    } else fitted_centered
    dq <- diff(wasserstein_grid)
    quadrature <- c(dq[1L] / 2,
                    (head(dq, -1L) + tail(dq, -1L)) / 2,
                    tail(dq, 1L) / 2)
    integrate_rows <- function(Z) drop(Z %*% quadrature)
    B <- A - 2 * integrate_rows(mean_centered * fitted_centered) +
      integrate_rows(fitted_centered^2)
    regression_fit <- fit
    bandwidth <- fit$optns$bw
  } else {
    if (space == "euclidean") {
      coordinates <- M
      reference <- muY
    } else {
      lower_index <- lower.tri(matrix(0, d, d))
      to_coordinates <- function(S) {
        Lchol <- t(chol(S))
        c(Lchol[lower_index], log(diag(Lchol)))
      }
      from_coordinates <- function(z) {
        Lchol <- matrix(0, d, d)
        n_lower <- sum(lower_index)
        Lchol[lower_index] <- z[seq_len(n_lower)]
        diag(Lchol) <- exp(z[n_lower + seq_len(d)])
        tcrossprod(Lchol)
      }
      coordinates <- do.call(rbind, lapply(E_list, to_coordinates))
      reference <- to_coordinates(muY)
    }
    centered <- sweep(coordinates, 2L, reference, "-")
    kernel <- .nw_option_kernel(optns)
    bandwidth <- optns$bw
    cv <- NULL
    if (identical(bandwidth, "CV")) {
      cv <- .nw_cv(X, optns$bwGrid, kernel, function(weights) {
        mean(rowSums((weights %*% centered - centered)^2))
      })
      bandwidth <- cv$bw
    }
    W <- NWFrechetWeights(xin = X, xout = X, bw = bandwidth, kernel = kernel)
    fitted_centered <- W %*% centered
    A <- drop(W %*% null_loss)
    B <- A - rowSums(fitted_centered^2)
    fitted_coordinates <- sweep(fitted_centered, 2L, reference, "+")
    fitted <- if (space == "euclidean") fitted_coordinates else {
      lapply(seq_len(n), function(i) from_coordinates(fitted_coordinates[i, ]))
    }
    regression_fit <- list(yout = fitted, weights = W, bw = bandwidth,
                           kernel = kernel, cv = cv)
  }
  if (length(A) != n || length(B) != n ||
      any(!is.finite(A)) || any(!is.finite(B))) {
    stop("Nonfinite or invalid objective values.")
  }
  tolerance <- 100 * .Machine$double.eps * pmax(1, abs(A), abs(B))
  if (any(A < -tolerance) || any(B < -tolerance)) {
    stop("Negative objective values.")
  }
  A <- pmax(A, 0)
  B <- pmax(B, 0)
  if (any(B > A + sqrt(.Machine$double.eps) * pmax(1, A, B))) {
    stop("The fitted objective exceeds the null objective. Check optimizer convergence and distance scaling.")
  }
  B <- pmin(B, A)
  objective_scale <- max(A, B)
  if (!is.finite(objective_scale) || objective_scale <= 0) {
    stop("Degenerate test variance.")
  }
  A_scaled <- A / objective_scale
  B_scaled <- B / objective_scale
  sigma_scaled <- sqrt(mean(A_scaled^2 + B_scaled^2) / 2)
  sigmahat <- objective_scale * sigma_scaled
  if (!is.finite(sigmahat) || sigmahat <= 0) stop("Degenerate test variance.")
  phi <- matrix(rgamma(as.double(n) * L, shape = 2, rate = 2),
                nrow = n, ncol = L)
  psi <- matrix(rgamma(as.double(n) * L, shape = 2, rate = 2),
                nrow = n, ncol = L)
  Dn_scaled <- mean(A_scaled - B_scaled)
  Tn_scaled <- Dn_scaled + drop(crossprod(A_scaled, phi - 1) -
                                  crossprod(B_scaled, psi - 1)) / n
  Zn <- sqrt(n) * Tn_scaled / sigma_scaled
  if (any(!is.finite(Zn))) stop("Nonfinite component test statistics.")
  component_log_p <- pnorm(Zn, lower.tail = FALSE, log.p = TRUE)
  component_log_q <- pnorm(Zn, lower.tail = TRUE, log.p = TRUE)
  active <- cauchy_weights > 0
  log_tail <- pmin(component_log_p[active], component_log_q[active])
  term_sign <- sign(Zn[active])
  log_abs_term <- rep(-Inf, sum(active))
  small_tail <- log_tail < log(1e-8)
  log_abs_term[small_tail] <- -log(pi) - log_tail[small_tail]
  intermediate <- !small_tail & log_tail < log(0.25)
  tail_probability <- exp(log_tail[intermediate])
  log_abs_term[intermediate] <- log(cospi(tail_probability)) -
    log(sinpi(tail_probability))
  central <- !small_tail & !intermediate & term_sign != 0
  log_abs_term[central] <- log(tanpi(
    -0.5 * expm1(log_tail[central] + log(2))
  ))
  log_weighted_term <- log(cauchy_weights[active]) + log_abs_term
  log_sum_exp <- function(z) {
    if (!length(z)) return(-Inf)
    zmax <- max(z)
    if (!is.finite(zmax)) return(zmax)
    zmax + log(sum(exp(z - zmax)))
  }
  log_positive <- log_sum_exp(log_weighted_term[term_sign > 0])
  log_negative <- log_sum_exp(log_weighted_term[term_sign < 0])
  if (is.infinite(log_positive) && log_positive > 0 &&
      is.infinite(log_negative) && log_negative > 0) {
    stop("Opposing infinite Cauchy terms; component statistics exceed numerical resolution.")
  }
  if (log_positive == log_negative) {
    cauchy_sign <- 0
    log_abs_cauchy <- -Inf
    Cn <- 0
    log_p_value <- -log(2)
  } else {
    cauchy_sign <- if (log_positive > log_negative) 1 else -1
    larger <- max(log_positive, log_negative)
    smaller <- min(log_positive, log_negative)
    log_abs_cauchy <- larger + log(-expm1(smaller - larger))
    Cn <- cauchy_sign * exp(log_abs_cauchy)
    log_positive_tail <- if (log_abs_cauchy > log(1e8)) {
      -log(pi) - log_abs_cauchy
    } else {
      pcauchy(exp(log_abs_cauchy), lower.tail = FALSE, log.p = TRUE)
    }
    log_p_value <- if (cauchy_sign > 0) {
      log_positive_tail
    } else {
      log1p(-exp(log_positive_tail))
    }
  }
  critical <- qcauchy(alpha, lower.tail = FALSE)
  answer <- list(
    statistic = c(Cauchy = Cn), parameter = c(n = n, L = L),
    p.value = exp(log_p_value), log.p.value = log_p_value,
    alternative = "the conditional Frechet mean differs from muY",
    method = paste0("One-sided Cauchy-combined residual test: ", space, " / ", method),
    data.name = data_name, reject = log_p_value < log(alpha),
    alpha = alpha, critical_value = critical,
    Tn = objective_scale * Tn_scaled,
    Dn = objective_scale * Dn_scaled, sigmahat = sigmahat,
    component_statistics = Zn,
    component_p.values = exp(component_log_p),
    component_log.p.values = component_log_p,
    n_randomizations = L, cauchy_weights = cauchy_weights,
    cauchy_sign = cauchy_sign, log_abs_cauchy = log_abs_cauchy,
    max_objective_weight = max(A_scaled^2 + B_scaled^2) /
      sum(A_scaled^2 + B_scaled^2),
    space = space, regression = method, muY = muY, bandwidth = bandwidth
  )
  if (return_fit) {
    answer$A <- A
    answer$B <- B
    answer$phi <- phi
    answer$psi <- psi
    answer$fit <- regression_fit
    answer$weights <- W
  }
  class(answer) <- "htest"
  answer
}



# Distance-profile boxplot and outlier detection

distance_profile_boxplot <- function(Y,
                                     space = c("euclidean", "sphere", "spd", "wasserstein", "simplex"),
                                     kappa = 1.5,
                                     fence = c("quartile", "adjusted"),
                                     central_level = 0.5,
                                     wasserstein_grid = NULL,
                                     u_grid = (seq_len(201L) - 0.5) / 201) {
  space <- .validate_space(space)
  fence <- match.arg(fence)
  
  D <- .boxplot_dist_matrix(Y, space, wasserstein_grid)
  n <- nrow(D)
  if (n < 4L) stop("At least four residuals are needed for the quartiles.")
  
  scores <- rowSums(D) / (n - 1)
  ordering <- order(scores)
  central <- seq_len(n) %in% ordering[seq_len(ceiling(central_level * n))]
  
  q1 <- .boxplot_quantile(scores, 0.25)
  q3 <- .boxplot_quantile(scores, 0.75)
  skewness <- if (fence == "adjusted") .boxplot_medcouple(scores) else 0
  fence_value <- q3 + kappa * exp(3 * skewness) * (q3 - q1)
  
  # Row i of the profile quantile matrix is hat G_i^{-1}(u) on u_grid: the ceil(u (n - 1))-th smallest of the n - 1 distances from residual i.
  sorted_distances <- t(apply(D, 1L, function(row) sort(row)[-1L]))
  rank_index <- pmin(pmax(ceiling(u_grid * (n - 1)), 1L), n - 1L)
  
  structure(
    list(
      distance = D,
      scores = scores,
      ordering = ordering,
      medoid = ordering[1L],
      central = central,
      quartiles = c(q1 = q1, q3 = q3),
      fence = fence_value,
      fence_type = fence,
      flagged = scores > fence_value,
      profile_quantiles = sorted_distances[, rank_index, drop = FALSE],
      u_grid = u_grid,
      kappa = kappa,
      space = space
    ),
    class = "distance_profile_boxplot"
  )
}


boxplot_figure <- function(Y,
                           space = c("euclidean", "sphere", "wasserstein"),
                           title = NULL,
                           highlight_indices = NULL,
                           highlight_label = "Highlighted",
                           kappa = 1.5,
                           fence = c("quartile", "adjusted"),
                           wasserstein_grid = NULL,
                           labels = NULL,
                           base_size = 8,
                           show_legend = FALSE,
                           show_ylab = TRUE) {
  space <- match.arg(space)
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required to assemble the figure.")
  }
  bp <- distance_profile_boxplot(Y, space, kappa = kappa, fence = match.arg(fence),
                                 wasserstein_grid = wasserstein_grid)
  panels <- list(
    profile = plot_distance_profile_boxplot(
      bp, title = title, show_legend = show_legend,
      ylab = if (show_ylab) "Profile quantile"
    ),
    distance = plot_boxplot_distances(
      Y, bp, space = space, highlight_indices = highlight_indices, highlight_label = highlight_label,
      show_legend = show_legend, wasserstein_grid = wasserstein_grid,
      ylab = if (show_ylab) "Signed distance"
    )
  )
  panels <- lapply(panels, .boxplot_panel_theme, base_size = base_size)
  column <- cowplot::plot_grid(plotlist = panels, ncol = 1L, align = "v",
                               labels = labels, label_size = base_size + 1)
  list(boxplot = bp, panels = panels, plot = column)
}


plot_distance_profile_boxplot <- function(bp,
                                          title = NULL,
                                          show_legend = FALSE,
                                          ylab = "Profile quantile") {
  Q <- bp$profile_quantiles
  u <- bp$u_grid
  envelope <- function(rows, name) {
    data.frame(u = u, lower = apply(Q[rows, , drop = FALSE], 2L, min),
               upper = apply(Q[rows, , drop = FALSE], 2L, max), Band = name)
  }
  bands <- rbind(envelope(!bp$flagged, "Non-flagged envelope"),
                 envelope(bp$central, "Central 50% envelope"))
  bands$Band <- factor(bands$Band, levels = c("Non-flagged envelope", "Central 50% envelope"))
  profile_lines <- function(rows, name) {
    data.frame(u = rep(u, length(rows)), value = as.vector(t(Q[rows, , drop = FALSE])),
               id = rep(rows, each = length(u)), Line = rep(name, length(u) * length(rows)),
               stringsAsFactors = FALSE)
  }
  lines <- rbind(profile_lines(which(bp$flagged), "Flagged profiles"),
                 profile_lines(bp$medoid, "Median profile"))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = bands,
                         ggplot2::aes(x = u, ymin = lower, ymax = upper, fill = Band)) +
    ggplot2::geom_line(data = lines,
                       ggplot2::aes(x = u, y = value, group = id, color = Line,
                                    linetype = Line, linewidth = Line)) +
    ggplot2::scale_fill_manual(values = c("Non-flagged envelope" = "gray85",
                                          "Central 50% envelope" = "gray55"), name = NULL) +
    ggplot2::scale_color_manual(values = c("Median profile" = "black",
                                           "Flagged profiles" = .boxplot_colors[["flagged"]]), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("Median profile" = "solid",
                                              "Flagged profiles" = "42"), name = NULL) +
    ggplot2::scale_linewidth_manual(values = c("Median profile" = 0.9,
                                               "Flagged profiles" = 0.4), guide = "none") +
    ggplot2::labs(title = title, x = "Probability level u", y = ylab)
  .boxplot_theme(p, legend = if (show_legend) c(0.02, 0.98))
}

plot_boxplot_distances <- function(Y, bp,
                                   space = c("euclidean", "sphere", "wasserstein"),
                                   highlight_indices = NULL,
                                   highlight_label = "Highlighted",
                                   show_legend = FALSE,
                                   wasserstein_grid = NULL,
                                   ylab = "Signed distance") {
  space <- match.arg(space)
  df <- residual_distance_data(Y, space = space, wasserstein_grid = wasserstein_grid)$data
  if (nrow(df) != length(bp$scores)) stop("Y and bp refer to different samples.")
  marks <- c("Not flagged", "Flagged", highlight_label)
  df$Mark <- factor(ifelse(bp$flagged, marks[2], marks[1]), levels = marks)
  highlighted <- df[seq_len(nrow(df)) %in% highlight_indices, , drop = FALSE]
  highlighted$Mark <- factor(rep(highlight_label, nrow(highlighted)), levels = marks)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Index, y = SignedDistance,
                                        color = Mark, shape = Mark)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(size = 1.2, stroke = 0.9) +
    ggplot2::geom_point(data = highlighted, size = 2.6, stroke = 0.6) +
    ggplot2::scale_color_manual(values = setNames(c("black", .boxplot_colors[["flagged"]],
                                                    .boxplot_colors[["highlight"]]), marks),
                                drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = setNames(c(21, 21, 0), marks), drop = FALSE, name = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = c(1.6, 1.6, 2.4), stroke = 0.8)
    )) +
    ggplot2::labs(x = "Index", y = ylab)
  .boxplot_theme(p, legend = if (show_legend) c(0.98, 0.02))
}


boxplot_legend <- function(panel, text_size = 7) {
  panel <- panel + ggplot2::theme(
    legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center",
    legend.text = ggplot2::element_text(size = text_size),
    legend.key.size = grid::unit(0.4, "cm"), legend.key.width = grid::unit(1, "cm"),
    legend.background = ggplot2::element_blank(), legend.box.background = ggplot2::element_blank()
  )
  suppressWarnings(cowplot::get_legend(panel))
}



.boxplot_colors <- c(flagged = "#D7191C", highlight = "#7B3294")


.boxplot_panel_theme <- function(p, base_size = 8) {
  p + ggplot2::theme(axis.title = ggplot2::element_text(size = base_size),
                     axis.text = ggplot2::element_text(size = base_size - 1),
                     plot.title = ggplot2::element_text(size = base_size + 1),
                     plot.margin = ggplot2::margin(2, 12, 1, 6))
}

.boxplot_theme <- function(p, legend = NULL) {
  p <- p + ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (is.null(legend)) return(p + ggplot2::theme(legend.position = "none"))
  p + ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = legend,
    legend.justification = legend,
    legend.background = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.text = ggplot2::element_text(size = 6),
    legend.spacing.y = grid::unit(0, "cm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0)
  )
}

# Pairwise distance matrix between the rows of Y.
.boxplot_dist_matrix <- function(Y, space, wasserstein_grid = NULL) {
  objects <- as_object_list(Y, space)
  n <- length(objects)
  if (space %in% c("spd", "simplex")) {
    D <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      D[i, j] <- D[j, i] <- object_dist(objects[[i]], objects[[j]], space)
    }
    return(D)
  }
  A <- do.call(rbind, objects)
  if (space == "sphere") return(acos(clamp(tcrossprod(A))))
  if (space == "wasserstein") {
    # The L2 distance of quantile functions by the trapezoidal rule is a Euclidean distance after scaling the columns by the trapezoid weights.
    if (is.null(wasserstein_grid)) wasserstein_grid <- seq(0, 1, length.out = ncol(A))
    h <- diff(wasserstein_grid)
    A <- sweep(A, 2L, sqrt(c(h, 0) / 2 + c(0, h) / 2), "*")
  }
  D <- sqrt(pmax(outer(rowSums(A^2), rowSums(A^2), "+") - 2 * tcrossprod(A), 0))
  diag(D) <- 0
  D
}


.boxplot_quantile <- function(x, alpha) sort(x)[ceiling(alpha * length(x))]


.boxplot_medcouple <- function(x) {
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("Package 'robustbase' is required for the skewness-adjusted fence.")
  }
  suppressMessages(robustbase::mc(x))
}




## prediction interval

.objpi_matrix <- function(Y, space) {
  M <- do.call(rbind, as_object_list(Y, space))
  if (!is.numeric(M) || !nrow(M) || !ncol(M) || any(!is.finite(M))) {
    stop("Responses must be finite numeric objects.")
  }
  M
}

.objpi_grid_weights <- function(grid, dimension) {
  if (is.null(grid)) grid <- seq(0, 1, length.out=dimension)
  if (length(grid)!=dimension || dimension<2L || any(!is.finite(grid)) ||
      any(diff(grid)<=0) || grid[1L]!=0 || grid[dimension]!=1) {
    stop("Use a common increasing Wasserstein probability grid from 0 to 1.")
  }
  delta <- diff(grid)
  c(delta[1L]/2, (head(delta, -1L)+tail(delta, -1L))/2,
    tail(delta, 1L)/2)
}


.objpi_local_linear_weights <- function(xin, xout, bandwidth,
                                        leave_one_out=FALSE) {
  if (!is.numeric(xin) || !is.numeric(xout) || !length(xout) ||
      length(xin)<2L || any(!is.finite(c(xin, xout))) ||
      length(unique(xin))<2L) stop("Use finite, nonconstant scalar predictors.")
  if (!is.numeric(bandwidth) || length(bandwidth)!=1L ||
      !is.finite(bandwidth) || bandwidth<=0) stop("Use a positive bandwidth.")
  if (leave_one_out && !identical(as.numeric(xin), as.numeric(xout))) {
    stop("Leave-one-out weights require xout=xin in the original order.")
  }
  z <- outer(xout, xin, function(x, u) (u-x)/bandwidth)
  log_kernel <- -.5*z^2
  if (leave_one_out) diag(log_kernel) <- -Inf
  log_kernel <- sweep(log_kernel, 1L, apply(log_kernel, 1L, max), "-")
  kernel <- exp(log_kernel)
  kernel <- kernel/rowSums(kernel)
  mean_z <- rowSums(kernel*z)
  centered_z <- sweep(z, 1L, mean_z, "-")
  variance_z <- rowSums(kernel*centered_z^2)
  if (any(!is.finite(variance_z)) ||
      any(variance_z<=.Machine$double.eps*pmax(1, mean_z^2))) {
    stop("Local-linear design is numerically singular; increase the bandwidth.")
  }
  kernel*(1-sweep(centered_z, 1L, mean_z/variance_z, "*"))
}


objpi_kfold_plan <- function(xin, n_folds=10L, trim=.2,
                             small_sample_loo=TRUE) {
  if (!is.numeric(xin) || length(xin)<3L || any(!is.finite(xin)) ||
      length(unique(xin))<3L) stop("Use at least three distinct finite predictors.")
  if (!is.numeric(n_folds) || length(n_folds)!=1L || !is.finite(n_folds) ||
      n_folds<2L || n_folds!=floor(n_folds)) stop("n_folds must be an integer >= 2.")
  if (!is.numeric(trim) || length(trim)!=1L || !is.finite(trim) ||
      trim<0 || trim>=.5) stop("trim must lie in [0, 0.5).")
  if (!is.logical(small_sample_loo) || length(small_sample_loo)!=1L ||
      is.na(small_sample_loo)) stop("small_sample_loo must be TRUE or FALSE.")
  limits <- range(xin)
  validation <- if (trim==0) seq_along(xin) else {
    which(xin>limits[1L]+trim*diff(limits) &
            xin<limits[2L]-trim*diff(limits))
  }
  m <- length(validation)
  if (m<2L) stop("Too few interior observations for K-fold validation.")
  k <- if (small_sample_loo && m<=30L) m else min(as.integer(n_folds), m)
  ordered <- validation[order(xin[validation])]
  labels <- unlist(lapply(seq_len(ceiling(m/k)), function(i) sample.int(k)),
                   use.names=FALSE)
  folds <- integer(length(xin))
  folds[ordered] <- labels[seq_len(m)]
  list(xin=as.numeric(xin), folds=folds, validation=validation,
       n_folds=as.integer(k), n_validation=m, trim=trim,
       small_sample_loo=small_sample_loo)
}


.objpi_gaussian_bw_range <- function(xin, xout) {
  x <- sort(unique(xin))
  if (length(x)<3L) stop("A fitting fold needs three distinct predictors.")
  lower <- 1.1*max(diff(x), x[2L]-min(xout),
                   max(xout)-x[length(x)-1L])/3
  upper <- diff(range(x))/3
  if (upper<lower) upper <- if (lower>upper*3/2) lower*1.01 else upper*3/2
  c(lower=lower, upper=upper)
}

.objpi_shared_bw_grid <- function(xin, bandwidth_grid, space, plan,
                                  predictor_domain=NULL) {
  if (!is.numeric(bandwidth_grid) || !length(bandwidth_grid) ||
      any(!is.finite(bandwidth_grid)) || any(bandwidth_grid<=0)) {
    stop("Use a positive numeric bandwidth grid.")
  }
  grid <- sort(unique(bandwidth_grid))
  limits <- NULL
  if (space=="sphere") {
    domain <- if (is.null(predictor_domain)) range(xin) else predictor_domain
    ranges <- .objpi_gaussian_bw_range(xin, domain)
    for (k in seq_len(plan$n_folds)) {
      held_out <- which(plan$folds==k)
      one <- .objpi_gaussian_bw_range(xin[-held_out], xin[held_out])
      ranges <- c(lower=max(ranges[1L], one[1L]),
                  upper=min(ranges[2L], one[2L]))
    }
    limits <- ranges
    grid <- grid[grid>ranges[1L]*(1+1e-10) & grid<ranges[2L]*(1-1e-10)]
    if (!length(grid)) {
      stop("No shared spherical bandwidth is admissible; extend fit_bw_constants.")
    }
  }
  list(grid=grid, requested_grid=sort(unique(bandwidth_grid)),
       admissible_range=limits)
}


.objpi_fixed_fit <- function(xin, Y, xout, space, regression, bandwidth,
                             wasserstein_grid=NULL, controls=NULL) {
  if (regression=="nw") {
    if (space=="euclidean") {
      fitted <- NWFrechetWeights(xin, xout, bw=bandwidth, kernel="gauss") %*% Y
      engine <- "NWFrechetWeights"
    } else if (space=="sphere") {
      optns <- .nw_controls(if (is.null(controls)) list() else controls)
      optns$bw <- bandwidth
      optns$kernel <- "gauss"
      prediction <- NWSpheReg(xin=xin, yin=Y, xout=xout, optns=optns)
      if (any(!prediction$converged)) stop("NW spherical regression failed to converge.")
      fitted <- prediction$yout
      engine <- "NWSpheReg"
    } else {
      prediction <- NWWassReg(xin=xin, qin=Y, xout=xout,
                              optns=list(bw=bandwidth, kernel="gauss",
                                         qSup=wasserstein_grid))
      fitted <- prediction$qout
      engine <- "NWWassReg"
    }
  } else if (space=="euclidean") {
    fitted <- .objpi_local_linear_weights(xin, xout, bandwidth) %*% Y
    engine <- "multivariate local-linear least squares"
  } else {
    if (!requireNamespace("frechet", quietly=TRUE)) {
      stop("Install the R package 'frechet' for local spherical/Wasserstein regression.")
    }
    if (space=="sphere") {
      bounds <- .objpi_gaussian_bw_range(xin, xout)
      if (bandwidth<bounds[1L] || bandwidth>bounds[2L]) {
        stop("LocSpheReg cannot honor this fixed bandwidth; revise the bandwidth grid.")
      }
      prediction <- frechet::LocSpheReg(xin=xin, yin=Y, xout=xout,
                                        optns=list(bw=bandwidth, kernel="gauss"))
      selected_bandwidth <- prediction$optns$bw
      fitted <- prediction$yout
      engine <- "frechet::LocSpheReg"
    } else {
      prediction <- frechet::LocDenReg(
        xin=xin, qin=Y, xout=xout,
        optns=list(bwReg=bandwidth, kernelReg="gauss",
                   qSup=wasserstein_grid, denLowerThreshold=FALSE))
      selected_bandwidth <- prediction$optns$bwReg
      fitted <- prediction$qout
      if (!isTRUE(all.equal(as.numeric(prediction$qSup),
                            as.numeric(wasserstein_grid), tolerance=1e-12))) {
        stop("LocDenReg returned a different quantile grid.")
      }
      engine <- "frechet::LocDenReg"
    }
    if (!isTRUE(all.equal(as.numeric(selected_bandwidth), as.numeric(bandwidth),
                          tolerance=1e-12))) {
      stop("frechet changed the fixed bandwidth during prediction.")
    }
  }
  fitted <- as.matrix(fitted)
  if (nrow(fitted)!=length(xout) || ncol(fitted)!=ncol(Y) ||
      any(!is.finite(fitted))) stop("Regression returned invalid predictions.")
  if (space=="sphere") {
    norms <- sqrt(rowSums(fitted^2))
    if (any(abs(norms-1)>1e-6)) stop("Spherical regression returned nonunit predictions.")
    fitted <- fitted/norms
  }
  if (space=="wasserstein" &&
      any(fitted[, -1L, drop=FALSE]-fitted[, -ncol(fitted), drop=FALSE] < -1e-8)) {
    stop("Wasserstein regression returned nonmonotone quantile predictions.")
  }
  list(fitted=fitted, engine=engine)
}

.objpi_kfold_cv <- function(xin, Y, space, regression, bandwidth_grid, plan,
                            wasserstein_grid=NULL, controls=NULL) {
  scores <- rep(Inf, length(bandwidth_grid))
  errors <- rep(NA_character_, length(bandwidth_grid))
  for (j in seq_along(bandwidth_grid)) {
    value <- tryCatch({
      losses <- rep(NA_real_, length(xin))
      for (k in seq_len(plan$n_folds)) {
        held_out <- which(plan$folds==k)
        fitting <- which(plan$folds!=k)
        prediction <- .objpi_fixed_fit(
          xin=xin[fitting], Y=Y[fitting, , drop=FALSE], xout=xin[held_out],
          space=space, regression=regression, bandwidth=bandwidth_grid[j],
          wasserstein_grid=wasserstein_grid, controls=controls)
        losses[held_out] <- objpi_row_distances(
          Y[held_out, , drop=FALSE], prediction$fitted, space,
          wasserstein_grid)^2
      }
      if (any(!is.finite(losses[plan$validation]))) stop("Invalid validation losses.")
      list(score=mean(losses[plan$validation]), error=NA_character_)
    }, error=function(e) list(score=Inf, error=conditionMessage(e)))
    scores[j] <- value$score
    errors[j] <- value$error
  }
  if (!any(is.finite(scores))) {
    stop("All K-fold bandwidth candidates failed: ", errors[which(!is.na(errors))[1L]])
  }
  selected <- which.min(scores)
  list(bandwidth=bandwidth_grid[selected], scores=scores, errors=errors,
       bandwidth_grid=bandwidth_grid, selected=selected,
       edge=selected %in% c(1L, length(bandwidth_grid)), plan=plan)
}


objpi_regression_fit <- function(
    xin, Y, xout, space, design=c("nw", "categorical"),
    wasserstein_grid=NULL, bandwidth="CV", bandwidth_grid=NULL,
    regression=c("local", "nw"), cv_plan=NULL, cv_folds=10L, cv_trim=.2,
    cv_small_sample_loo=TRUE, controls=NULL, predictor_domain=NULL) {
  design <- match.arg(design)
  space <- match.arg(space, c("euclidean", "sphere", "wasserstein"))
  regression <- match.arg(regression)
  Y <- .objpi_matrix(Y, space)
  if (!is.numeric(xin) || !is.numeric(xout) || !length(xout) ||
      length(xin)!=nrow(Y) || any(!is.finite(c(xin, xout)))) {
    stop("Use one finite scalar predictor per response and prediction point.")
  }
  if (design=="categorical") {
    categories <- sort(unique(xin))
    index <- match(xout, categories)
    if (anyNA(index)) stop("A prediction category is absent from the fitting sample.")
    means <- lapply(categories, function(x) {
      object_frechet_mean(Y[xin==x, , drop=FALSE], space=space,
                          wasserstein_grid=wasserstein_grid)
    })
    return(list(fitted=do.call(rbind, means)[index, , drop=FALSE],
                bandwidth=NA_real_, cv_edge=NA, cv=NULL,
                engine="within-category Frechet mean",
                bandwidth_selection="not applicable", regression=regression))
  }
  if (length(xin)<3L || length(unique(xin))<3L) {
    stop("Regression requires at least three distinct predictors.")
  }
  if (!is.null(controls) && !is.list(controls)) stop("controls must be a list or NULL.")
  if (!is.null(predictor_domain)) {
    if (!is.numeric(predictor_domain) || length(predictor_domain)!=2L ||
        any(!is.finite(predictor_domain)) || diff(predictor_domain)<=0 ||
        any(c(xin, xout)<predictor_domain[1L] |
            c(xin, xout)>predictor_domain[2L])) {
      stop("The predictor domain must contain all fitting and prediction points.")
    }
  }
  if (space=="sphere") {
    if (any(abs(rowSums(Y^2)-1)>1e-8)) stop("Spherical responses must have unit norm.")
    Y <- Y/sqrt(rowSums(Y^2))
  } else if (space=="wasserstein") {
    if (is.null(wasserstein_grid)) {
      wasserstein_grid <- seq(0, 1, length.out=ncol(Y))
    }
    .objpi_grid_weights(wasserstein_grid, ncol(Y))
    if (any(Y[, -1L, drop=FALSE]-Y[, -ncol(Y), drop=FALSE]<0)) {
      stop("Wasserstein responses must be increasing quantile functions.")
    }
  }
  use_cv <- identical(bandwidth, "CV")
  if (!use_cv && (!is.numeric(bandwidth) || length(bandwidth)!=1L ||
                  !is.finite(bandwidth) || bandwidth<=0)) {
    stop("bandwidth must be 'CV' or a positive scalar.")
  }
  cv <- NULL
  cv_edge <- NA
  if (use_cv) {
    if (is.null(cv_plan)) {
      cv_plan <- objpi_kfold_plan(xin, cv_folds, cv_trim, cv_small_sample_loo)
    }
    if (!identical(cv_plan$xin, as.numeric(xin)) ||
        length(cv_plan$folds)!=length(xin) ||
        any(!is.finite(cv_plan$folds)) ||
        !identical(sort(unique(cv_plan$folds[cv_plan$folds>0])),
                   seq_len(cv_plan$n_folds)) ||
        !identical(as.integer(which(cv_plan$folds>0)),
                   as.integer(cv_plan$validation))) {
      stop("The supplied K-fold plan does not match the fitting predictors.")
    }
    if (is.null(bandwidth_grid)) {
      bandwidth_grid <- exp(seq(log(.1), log(1), length.out=15L))*length(xin)^(-1/5)
    }
    grid <- .objpi_shared_bw_grid(xin, bandwidth_grid, space, cv_plan,
                                  predictor_domain)
    cv <- .objpi_kfold_cv(xin, Y, space, regression, grid$grid, cv_plan,
                          wasserstein_grid, controls)
    cv$requested_grid <- grid$requested_grid
    cv$admissible_range <- grid$admissible_range
    bandwidth <- cv$bandwidth
    cv_edge <- cv$edge
  }
  prediction <- .objpi_fixed_fit(xin, Y, xout, space, regression, bandwidth,
                                 wasserstein_grid, controls)
  selection <- if (use_cv) "shared K-fold CV using fitting sample only" else "fixed"
  list(fitted=prediction$fitted, bandwidth=as.numeric(bandwidth),
       cv_edge=cv_edge, cv=cv, engine=prediction$engine,
       bandwidth_selection=selection, regression=regression)
}


objpi_local_fit <- function(
    xin, Y, xout, space, design=c("nw", "categorical"),
    wasserstein_grid=NULL, bandwidth="CV", bandwidth_grid=NULL,
    cv_plan=NULL, cv_folds=10L, cv_trim=.2, cv_small_sample_loo=TRUE,
    controls=NULL, predictor_domain=NULL) {
  objpi_regression_fit(
    xin=xin, Y=Y, xout=xout, space=space, design=match.arg(design),
    wasserstein_grid=wasserstein_grid, bandwidth=bandwidth,
    bandwidth_grid=bandwidth_grid, regression="local", cv_plan=cv_plan,
    cv_folds=cv_folds, cv_trim=cv_trim, cv_small_sample_loo=cv_small_sample_loo,
    controls=controls, predictor_domain=predictor_domain)
}

objpi_nw_fit <- function(
    xin, Y, xout, space, design=c("nw", "categorical"),
    wasserstein_grid=NULL, bandwidth="CV", bandwidth_grid=NULL,
    cv_plan=NULL, cv_folds=10L, cv_trim=.2, cv_small_sample_loo=TRUE,
    controls=NULL, predictor_domain=NULL) {
  objpi_regression_fit(
    xin=xin, Y=Y, xout=xout, space=space, design=match.arg(design),
    wasserstein_grid=wasserstein_grid, bandwidth=bandwidth,
    bandwidth_grid=bandwidth_grid, regression="nw", cv_plan=cv_plan,
    cv_folds=cv_folds, cv_trim=cv_trim, cv_small_sample_loo=cv_small_sample_loo,
    controls=controls, predictor_domain=predictor_domain)
}

objpi_distance_matrix <- function(A, B=NULL,
                                  space=c("euclidean", "sphere", "wasserstein"),
                                  wasserstein_grid=NULL) {
  space <- match.arg(space)
  A <- .objpi_matrix(A, space)
  same <- is.null(B)
  B <- if (same) A else .objpi_matrix(B, space)
  if (ncol(A)!=ncol(B)) stop("Object dimensions do not match.")
  if (space=="sphere") {
    dot_products <- tcrossprod(A, B)
    D <- acos(pmax(pmin(dot_products, 1), -1))
  } else {
    origin <- A[1L, ]
    AA <- sweep(A, 2L, origin, "-")
    BB <- sweep(B, 2L, origin, "-")
    if (space=="wasserstein") {
      weights <- sqrt(.objpi_grid_weights(wasserstein_grid, ncol(A)))
      AA <- sweep(AA, 2L, weights, "*")
      BB <- sweep(BB, 2L, weights, "*")
    }
    squared <- outer(rowSums(AA^2), rowSums(BB^2), "+") - 2*tcrossprod(AA, BB)
    D <- sqrt(pmax(squared, 0))
  }
  if (same) diag(D) <- 0
  D
}

objpi_row_distances <- function(A, B,
                                space=c("euclidean", "sphere", "wasserstein"),
                                wasserstein_grid=NULL) {
  space <- match.arg(space)
  A <- .objpi_matrix(A, space)
  B <- .objpi_matrix(B, space)
  if (!identical(dim(A), dim(B))) stop("Use one fitted object per response.")
  if (space=="sphere") return(acos(pmax(-1, pmin(1, rowSums(A*B)))))
  squared <- (A-B)^2
  if (space=="wasserstein") {
    return(sqrt(drop(squared %*% .objpi_grid_weights(wasserstein_grid, ncol(A)))))
  }
  sqrt(rowSums(squared))
}

objpi_residuals <- function(Y, fitted, center, space, wasserstein_grid=NULL) {
  do.call(rbind, make_metric_residuals(
    Y, fitted, space=space, muY=center,
    wasserstein_grid=wasserstein_grid)$residuals)
}

objpi_order_quantile <- function(scores, alpha=.05, conformal=TRUE) {
  if (!length(scores) || any(!is.finite(scores))) stop("Scores must be finite.")
  if (length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1) {
    stop("alpha must lie strictly between zero and one.")
  }
  m <- length(scores)
  rank <- ceiling((m+as.integer(conformal))*(1-alpha)-1e-12)
  if (rank>m) return(Inf)
  sort(scores, partial=rank)[rank]
}

objpi_proposed_calibrate <- function(Y=NULL, fitted=NULL, center=NULL, space,
                                     wasserstein_grid=NULL, alpha=.05, E=NULL) {
  space <- match.arg(space, c("euclidean", "sphere", "wasserstein"))
  if (is.null(E)) {
    if (is.null(Y) || is.null(fitted)) {
      stop("Supply calibration residuals E, or both Y and fitted.")
    }
    Y <- .objpi_matrix(Y, space)
    if (is.null(center)) {
      center <- object_frechet_mean(Y, space, wasserstein_grid)
    }
    E <- objpi_residuals(Y, fitted, center, space, wasserstein_grid)
  }
  E <- .objpi_matrix(E, space)
  if (nrow(E)<2L) stop("At least two calibration residuals are needed.")
  if (space=="sphere" && any(abs(rowSums(E^2)-1)>1e-6)) {
    stop("Spherical calibration residuals must be unit vectors.")
  }
  distances <- objpi_distance_matrix(E, space=space,
                                     wasserstein_grid=wasserstein_grid)
  scores <- rowSums(distances)/(nrow(E)-1L)
  list(residuals=E, center=center, space=space,
       wasserstein_grid=wasserstein_grid, calibration_scores=scores,
       threshold=objpi_order_quantile(scores, alpha, conformal=FALSE))
}

objpi_proposed_scores <- function(object, Y=NULL, fitted=NULL, E=NULL) {
  if (is.null(E)) {
    if (is.null(Y) || is.null(fitted) || is.null(object$center)) {
      stop("Supply residuals E, or Y and fitted with a stored calibration center.")
    }
    E <- objpi_residuals(Y, fitted, object$center, object$space,
                         object$wasserstein_grid)
  }
  E <- .objpi_matrix(E, object$space)
  calibration <- .objpi_matrix(object$residuals, object$space)
  if (ncol(E)==1L && nrow(E)==ncol(calibration) && ncol(calibration)>1L) {
    E <- t(E)
  }
  if (object$space=="sphere" && any(abs(rowSums(E^2)-1)>1e-6)) {
    stop("Spherical prediction residuals must be unit vectors.")
  }
  rowMeans(objpi_distance_matrix(E, calibration, object$space,
                                 object$wasserstein_grid))
}


objpi_lm_calibrate <- function(Y, fitted, space,
                               wasserstein_grid=NULL, alpha=.05) {
  scores <- objpi_row_distances(Y, fitted, space, wasserstein_grid)
  list(space=space, wasserstein_grid=wasserstein_grid,
       calibration_scores=scores,
       threshold=objpi_order_quantile(scores, alpha, conformal=TRUE))
}

objpi_lm_scores <- function(object, Y, fitted) {
  objpi_row_distances(Y, fitted, object$space, object$wasserstein_grid)
}

objpi_zm_weights <- function(xin, xout, bandwidth, design=c("nw", "categorical")) {
  design <- match.arg(design)
  if (!length(xin) || !length(xout) || any(!is.finite(c(xin, xout)))) {
    stop("Predictors must be finite scalar values.")
  }
  if (design=="categorical") {
    W <- 1*outer(xout, xin, "==")
    if (any(rowSums(W)==0)) stop("A prediction category is absent from training.")
  } else {
    if (length(bandwidth)!=1L || !is.finite(bandwidth) || bandwidth<=0) {
      stop("Use a positive profile bandwidth.")
    }
    W <- pmax(1-(outer(xout, xin, "-")/bandwidth)^2, 0)
    empty <- which(rowSums(W)==0)
    for (i in empty) W[i, which.min(abs(xin-xout[i]))] <- 1
  }
  W/rowSums(W)
}

objpi_zm_inverse_profiles <- function(distances, weights, probabilities) {
  if (!identical(dim(distances), dim(weights))) {
    stop("Use one conditional weight vector per distance profile.")
  }
  Q <- matrix(0, nrow(distances), length(probabilities))
  for (i in seq_len(nrow(distances))) {
    positive <- which(weights[i, ]>0)
    ordered <- positive[order(distances[i, positive])]
    cumulative <- cumsum(weights[i, ordered])
    cumulative[length(cumulative)] <- 1
    index <- findInterval(probabilities, cumulative, left.open=TRUE)+1L
    Q[i, ] <- distances[i, ordered[pmin(index, length(ordered))]]
  }
  Q
}

objpi_zm_transport_costs <- function(Q, Q_train, weights) {
  if (nrow(Q)!=nrow(weights) || nrow(Q_train)!=ncol(weights) ||
      ncol(Q)!=ncol(Q_train)) stop("Profile dimensions do not match.")
  vapply(seq_len(nrow(Q)), function(i) {
    positive <- which(weights[i, ]>0)
    difference <- sweep(Q_train[positive, , drop=FALSE], 2L, Q[i, ], "-")
    sum(weights[i, positive]*rowMeans(abs(difference)))
  }, numeric(1))
}

objpi_zm_conditional_cdf <- function(costs, training_costs, weights) {
  if (length(costs)!=nrow(weights) || length(training_costs)!=ncol(weights)) {
    stop("Cost and weight dimensions do not match.")
  }
  scores <- vapply(seq_along(costs), function(i) {
    sum(weights[i, ]*(training_costs<=costs[i]))
  }, numeric(1))
  scores[abs(scores-1)<1e-12] <- 1
  pmax(0, pmin(1, scores))
}

objpi_zm_fit <- function(X, Y, space,
                         design=c("nw", "categorical"), wasserstein_grid=NULL,
                         bandwidth=.5*(length(X)/log(length(X)))^(-1/5),
                         profile_points=201L) {
  design <- match.arg(design)
  Y <- .objpi_matrix(Y, space)
  if (length(X)!=nrow(Y) || length(X)<2L) stop("Training dimensions do not match.")
  if (length(profile_points)!=1L || !is.finite(profile_points) ||
      profile_points<2 || profile_points!=floor(profile_points)) {
    stop("Use at least two profile integration points.")
  }
  probabilities <- (seq_len(profile_points)-.5)/profile_points
  W <- objpi_zm_weights(X, X, bandwidth, design)
  D <- objpi_distance_matrix(Y, space=space, wasserstein_grid=wasserstein_grid)
  Q <- objpi_zm_inverse_profiles(D, W, probabilities)
  costs <- objpi_zm_transport_costs(Q, Q, W)
  list(X=X, Y=Y, space=space, design=design,
       wasserstein_grid=wasserstein_grid, bandwidth=bandwidth,
       probabilities=probabilities, training_profiles=Q, training_costs=costs)
}

objpi_zm_scores <- function(object, X, Y) {
  Y <- .objpi_matrix(Y, object$space)
  if (length(X)!=nrow(Y)) stop("Prediction dimensions do not match.")
  W <- objpi_zm_weights(object$X, X, object$bandwidth, object$design)
  D <- objpi_distance_matrix(Y, object$Y, object$space, object$wasserstein_grid)
  Q <- objpi_zm_inverse_profiles(D, W, object$probabilities)
  costs <- objpi_zm_transport_costs(Q, object$training_profiles, W)
  objpi_zm_conditional_cdf(costs, object$training_costs, W)
}

objpi_zm_calibrate <- function(object, X, Y, alpha=.05) {
  object$calibration_scores <- objpi_zm_scores(object, X, Y)
  object$threshold <- objpi_order_quantile(object$calibration_scores, alpha,
                                           conformal=TRUE)
  object
}

objpi_prediction_boxplot_data <- function(object, fitted = NULL, future = NULL,
                                          candidates = NULL,
                                          u_grid = (seq_len(201L) - 0.5) / 201,
                                          E_pred = NULL, E_future = NULL,
                                          E_candidates = NULL) {
  required <- c("residuals", "center", "space", "threshold")
  if (!all(required %in% names(object))) {
    stop("Use an object returned by objpi_proposed_calibrate().")
  }
  if (length(object$threshold) != 1L || is.na(object$threshold) || object$threshold < 0) {
    stop("The prediction threshold must be nonnegative, possibly Inf.")
  }
  space <- match.arg(object$space, c("euclidean", "sphere", "wasserstein"))
  calibration <- .objpi_matrix(object$residuals, space)
  if (nrow(calibration) < 2L) stop("At least two calibration residuals are required.")
  if (!length(u_grid) || any(!is.finite(u_grid)) ||
      any(u_grid < 0 | u_grid > 1) || is.unsorted(u_grid)) {
    stop("u_grid must be increasing probability levels in [0, 1].")
  }
  validate_objects <- function(Y, name, single = FALSE) {
    Y <- .objpi_matrix(Y, space)
    if (single && ncol(Y) == 1L && nrow(Y) == ncol(calibration) && ncol(calibration) > 1L) {
      Y <- t(Y)
    }
    if (ncol(Y) != ncol(calibration)) stop(name, " has an incompatible dimension.")
    if (single && nrow(Y) != 1L) stop(name, " must contain one object.")
    if (space == "sphere" && any(abs(rowSums(Y^2) - 1) > 1e-6)) {
      stop(name, " must contain unit vectors for spherical data.")
    }
    Y
  }
  validate_objects(calibration, "Calibration residuals")
  if (!is.null(fitted) && !is.null(E_pred)) {
    stop("Supply fitted or E_pred for the prediction, not both.")
  }
  if (!is.null(future) && !is.null(E_future)) {
    stop("Supply future or E_future for the known future point, not both.")
  }
  if (!is.null(candidates) && !is.null(E_candidates)) {
    stop("Supply candidates or E_candidates, not both.")
  }
  if (is.null(E_pred)) {
    if (is.null(fitted) || is.null(object$center)) {
      stop("Supply E_pred directly, or fitted with a stored calibration center.")
    }
    fitted <- validate_objects(fitted, "fitted", single = TRUE)
    prediction_residual <- objpi_residuals(fitted, fitted, object$center, space,
                                           object$wasserstein_grid)
    input_mode <- "responses"
  } else {
    prediction_residual <- validate_objects(E_pred, "E_pred", single = TRUE)
    input_mode <- "residuals"
  }
  fitted_distances <- objpi_distance_matrix(prediction_residual, calibration, space,
                                            object$wasserstein_grid)
  if (!is.null(E_candidates)) {
    candidate_residuals <- validate_objects(E_candidates, "E_candidates")
    candidate_source <- "Supplied candidate residuals"
  } else if (!is.null(candidates)) {
    if (is.null(fitted) || is.null(object$center)) {
      stop("Use E_candidates when supplying prediction residuals directly.")
    }
    candidates <- validate_objects(candidates, "candidates")
    fitted_rows <- fitted[rep(1L, nrow(candidates)), , drop = FALSE]
    candidate_residuals <- objpi_residuals(candidates, fitted_rows, object$center,
                                           space, object$wasserstein_grid)
    candidate_source <- "Supplied candidate objects"
  } else {
    candidate_residuals <- rbind(calibration, prediction_residual)
    candidate_source <- "Calibration residuals and the prediction residual"
  }
  distances <- objpi_distance_matrix(candidate_residuals, calibration, space,
                                     object$wasserstein_grid)
  scores <- rowMeans(distances)
  inside <- scores <= object$threshold
  m <- ncol(distances)
  ranks <- pmax(1L, pmin(m, ceiling(u_grid * m)))
  profiles <- t(apply(distances, 1L, sort))[, ranks, drop = FALSE]
  accepted <- which(inside)
  envelope <- data.frame(u = numeric(0), lower = numeric(0), upper = numeric(0))
  if (length(accepted)) {
    envelope <- data.frame(
      u = u_grid,
      lower = apply(profiles[accepted, , drop = FALSE], 2L, min),
      upper = apply(profiles[accepted, , drop = FALSE], 2L, max)
    )
  }
  future_profile <- future_score <- future_inside <- NULL
  future_residual <- NULL
  if (!is.null(E_future)) {
    future_residual <- validate_objects(E_future, "E_future", single = TRUE)
  } else if (!is.null(future)) {
    if (is.null(fitted) || is.null(object$center)) {
      stop("Use E_future when supplying prediction residuals directly.")
    }
    future <- validate_objects(future, "future", single = TRUE)
    future_residual <- objpi_residuals(future, fitted, object$center, space,
                                       object$wasserstein_grid)
  }
  if (!is.null(future_residual)) {
    future_distances <- objpi_distance_matrix(future_residual, calibration, space,
                                              object$wasserstein_grid)
    future_profile <- sort(as.numeric(future_distances))[ranks]
    future_score <- mean(future_distances)
    future_inside <- future_score <= object$threshold
  }
  list(
    u_grid = u_grid, profile_quantiles = profiles, scores = scores, inside = inside,
    threshold = object$threshold, envelope = envelope,
    fitted_profile = sort(as.numeric(fitted_distances))[ranks],
    fitted_score = mean(fitted_distances),
    future_profile = future_profile, future_score = future_score,
    future_inside = future_inside, candidate_source = candidate_source,
    prediction_residual = prediction_residual, future_residual = future_residual,
    input_mode = input_mode
  )
}

objpi_prediction_boxplot <- function(object, fitted = NULL, future = NULL,
                                     candidates = NULL,
                                     title = "Prediction-set boxplot",
                                     base_size = 10,
                                     u_grid = (seq_len(201L) - 0.5) / 201,
                                     E_pred = NULL, E_future = NULL,
                                     E_candidates = NULL,
                                     y_axis_title="Profile quantile") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  d <- objpi_prediction_boxplot_data(
    object = object, fitted = fitted, future = future, candidates = candidates,
    u_grid = u_grid, E_pred = E_pred, E_future = E_future, E_candidates = E_candidates
  )
  outside <- which(!d$inside)
  levels <- c("Outside prediction set", "Prediction", "Known future point")
  lines <- data.frame(
    u = rep(d$u_grid, length(outside)),
    value = as.vector(t(d$profile_quantiles[outside, , drop = FALSE])),
    id = rep(seq_along(outside), each = length(d$u_grid)),
    Profile = rep(levels[1L], length(outside) * length(d$u_grid))
  )
  lines <- rbind(lines, data.frame(
    u = d$u_grid, value = d$fitted_profile,
    id = length(outside) + 1L, Profile = levels[2L]
  ))
  if (!is.null(d$future_profile)) {
    lines <- rbind(lines, data.frame(
      u = d$u_grid, value = d$future_profile,
      id = length(outside) + 2L, Profile = levels[3L]
    ))
  }
  lines$Profile <- factor(lines$Profile, levels = levels)
  p <- ggplot2::ggplot()
  if (nrow(d$envelope)) {
    band <- d$envelope
    band$Band <- "Prediction set"
    p <- p + ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = u, ymin = lower, ymax = upper, fill = Band)
    )
  }
  p <- p +
    ggplot2::geom_line(
      data = lines,
      ggplot2::aes(x = u, y = value, group = id, colour = Profile,
                   linetype = Profile, linewidth = Profile)
    ) +
    ggplot2::scale_fill_manual(
      values = c("Prediction set" = "grey80"), name = NULL
    ) +
    ggplot2::scale_colour_manual(
      values = c("Outside prediction set" = "grey60", "Prediction" = "#2166AC",
                 "Known future point" = "black"), name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Outside prediction set" = "dotted", "Prediction" = "dashed",
                 "Known future point" = "solid"), name = NULL
    ) +
    ggplot2::scale_linewidth_manual(
      values = c("Outside prediction set" = 0.4, "Prediction" = 0.9,
                 "Known future point" = 0.9), guide = "none"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1),
      colour = ggplot2::guide_legend(order = 2),
      linetype = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(title = title, x = "Probability level", y = y_axis_title) +
    ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain", hjust = 0.5),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      legend.position = "bottom", legend.direction = "horizontal",
      legend.box = "horizontal", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  list(plot = p, data = d, future_score = d$future_score,
       future_inside = d$future_inside)
}
