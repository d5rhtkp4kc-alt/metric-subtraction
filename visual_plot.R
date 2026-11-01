library(ggplot2)
library(ggforce)
library(patchwork)
library(showtext)
library(sysfonts)
library(grid)


font_add_google("Noto Sans Math", "noto_math")
showtext_auto()

# ============================================================
# Figure 1. subtraction
# ============================================================

title_style <- element_text(
  hjust = 0.5,
  size = 14,
  family = "noto_math",
  margin = margin(t = 5, b = 10)
)

common_theme <- theme_void() +
  theme(
    plot.title = title_style,
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )


# 1.1 Panel (a): Euclidean space


df_euclid <- data.frame(
  x = c(0, 3, 1, 2),
  y = c(0, 2, 3, -1),
  label = c("italic(c)", "italic(a)", "italic(b)", "italic(d)")
)

p1 <- ggplot(df_euclid) +
  geom_vline(xintercept = 0, color = "gray75", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "gray75", linewidth = 0.6) +
  
  # red transported vectors
  geom_segment(
    aes(x = 1, y = 3, xend = 3, yend = 2),
    color = "firebrick",
    linewidth = 1.4,
    arrow = arrow(length = unit(0.28, "cm"), type = "closed")
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = 2, yend = -1),
    color = "firebrick",
    linewidth = 1.4,
    arrow = arrow(length = unit(0.28, "cm"), type = "closed")
  ) +
  
  # blue dashed comparison segments
  geom_segment(
    aes(x = 0, y = 0, xend = 1, yend = 3),
    color = "#1E40FF",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = 2, y = -1, xend = 3, yend = 2),
    color = "#1E40FF",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  geom_point(aes(x, y), size = 3.5, color = "firebrick") +
  geom_text(
    aes(x, y, label = label),
    parse = TRUE,
    size = 6,
    vjust = -1.1
  ) +
  
  labs(title = "(a) Euclidean space") +
  coord_fixed(
    xlim = c(-1.25, 4.25),
    ylim = c(-1.50, 4.00)
  ) +
  common_theme



norm_vec <- function(x) sqrt(sum(x^2))

geo_dist <- function(x, y) {
  acos(max(min(sum(x * y), 1), -1))
}

get_rotation_matrix <- function(g_from, g_to) {
  p <- length(g_from)
  theta <- geo_dist(g_from, g_to)
  
  if (theta < 1e-10) return(diag(p))
  
  u1 <- g_from
  proj <- g_to - sum(g_to * u1) * u1
  u2 <- proj / norm_vec(proj)
  Q <- (u2 %*% t(u1)) - (u1 %*% t(u2))
  R <- diag(p) + sin(theta) * Q + (1 - cos(theta)) * (Q %*% Q)
  return(R)
}

get_smooth_geodesic <- function(p1, p2, n = 500) {
  p1 <- as.numeric(p1 / norm_vec(p1))
  p2 <- as.numeric(p2 / norm_vec(p2))
  theta <- acos(sum(p1 * p2))
  
  if (theta < 1e-10) {
    return(data.frame(
      V1 = rep(p1[1], n),
      V2 = rep(p1[2], n),
      V3 = rep(p1[3], n)
    ))
  }
  
  t <- seq(0, 1, length.out = n)
  path <- sapply(t, function(ti) {
    (sin((1 - ti) * theta) * p1 + sin(ti * theta) * p2) / sin(theta)
  })
  as.data.frame(t(path))
}


# 1.2 Panel (b): Sphere


g1_3d <- c(cos(2.2), 0, sin(2.2))   # b
g3_3d <- c(cos(3.6), 0, sin(3.6))   # c
q2    <- c(-0.1, -0.6, 0.5)
g2_3d <- as.numeric(q2 / norm_vec(q2))  # a

R_mat <- get_rotation_matrix(g1_3d, g3_3d)
g4_3d <- drop(R_mat %*% g2_3d)           # d

path_red_top <- get_smooth_geodesic(g1_3d, g2_3d)  # b -> a
path_red_bot <- get_smooth_geodesic(g3_3d, g4_3d)  # c -> d
path_blue1   <- get_smooth_geodesic(g1_3d, g3_3d)  # b -> c
path_blue2   <- get_smooth_geodesic(g2_3d, g4_3d)  # a -> d

t_eq <- seq(0, pi, length.out = 1000)
eq_front <- data.frame(x = cos(t_eq + pi), y = 0.3 * sin(t_eq + pi))
eq_back  <- data.frame(x = cos(t_eq),      y = 0.3 * sin(t_eq))

p2 <- ggplot() +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "black", linewidth = 0.7) +
  geom_path(data = eq_back, aes(x = x, y = y),
            color = "gray70", linewidth = 0.6, linetype = "dashed") +
  geom_path(data = eq_front, aes(x = x, y = y),
            color = "gray35", linewidth = 0.6) +
  
  geom_path(data = path_blue1, aes(x = V1, y = V3),
            color = "#1E40FF", linewidth = 1.2, linetype = "dashed") +
  geom_path(data = path_blue2, aes(x = V1, y = V3),
            color = "#1E40FF", linewidth = 1.2, linetype = "dashed") +
  
  geom_path(data = path_red_top, aes(x = V1, y = V3),
            color = "firebrick", linewidth = 1.4,
            arrow = arrow(length = unit(0.24, "cm"), type = "closed")) +
  geom_path(data = path_red_bot, aes(x = V1, y = V3),
            color = "firebrick", linewidth = 1.4,
            arrow = arrow(length = unit(0.24, "cm"), type = "closed")) +
  
  geom_point(
    aes(
      x = c(g1_3d[1], g2_3d[1], g3_3d[1], g4_3d[1]),
      y = c(g1_3d[3], g2_3d[3], g3_3d[3], g4_3d[3])
    ),
    color = "firebrick",
    size = 3
  ) +
  
  annotate("text", x = g1_3d[1] - 0.12, y = g1_3d[3] + 0.03, label = "italic(b)", parse = TRUE, size = 6) +
  annotate("text", x = g2_3d[1] + 0.12, y = g2_3d[3],        label = "italic(a)", parse = TRUE, size = 6) +
  annotate("text", x = g3_3d[1] - 0.12, y = g3_3d[3] - 0.05, label = "italic(c)", parse = TRUE, size = 6) +
  annotate("text", x = g4_3d[1] + 0.12, y = g4_3d[3],        label = "italic(d)", parse = TRUE, size = 6) +
  
  labs(title = "(b) Sphere") +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  common_theme


# 1.3 Panel (c): General metric space


theta <- seq(0, 2 * pi, length.out = 500)

r <- 1 + 0.10 * sin(3 * theta) + 0.06 * cos(5 * theta)

metric_boundary <- data.frame(
  x = 1.40 * r * cos(theta),
  y = 1.00 * r * sin(theta)
)

metric_pts <- data.frame(
  name = c("b", "a", "c", "d"),
  x    = c(-0.90,  0.07, -0.65,  0.70),
  y    = c( 0.00,  0.69, -0.62, -0.39)
)

b <- metric_pts[metric_pts$name == "b", ]
a <- metric_pts[metric_pts$name == "a", ]
c <- metric_pts[metric_pts$name == "c", ]
d <- metric_pts[metric_pts$name == "d", ]

p3 <- ggplot() +
  geom_polygon(
    data = metric_boundary,
    aes(x = x, y = y),
    fill = "#D4D8DD",
    color = NA
  ) +
  
  # blue dashed reference curves
  geom_curve(
    aes(x = b$x, y = b$y, xend = c$x, yend = c$y),
    curvature = 0.25,
    color = "#1E40FF",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_curve(
    aes(x = a$x, y = a$y, xend = d$x, yend = d$y),
    curvature = -0.22,
    color = "#1E40FF",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  # red geodesic transport curves
  geom_curve(
    aes(x = b$x, y = b$y, xend = a$x, yend = a$y),
    curvature = -0.25,
    color = "firebrick",
    linewidth = 1.5,
    arrow = arrow(length = unit(0.26, "cm"), type = "closed")
  ) +
  geom_curve(
    aes(x = c$x, y = c$y, xend = d$x, yend = d$y),
    curvature = -0.20,
    color = "firebrick",
    linewidth = 1.5,
    arrow = arrow(length = unit(0.26, "cm"), type = "closed")
  ) +
  
  geom_point(data = metric_pts, aes(x = x, y = y), size = 3.5, color = "firebrick") +
  
  annotate("text", x = b$x - 0.12, y = b$y - 0.18, label = "italic(b)", parse = TRUE, size = 6) +
  annotate("text", x = a$x + 0.16, y = a$y + 0.1, label = "italic(a)", parse = TRUE, size = 6) +
  annotate("text", x = c$x - 0.10, y = c$y - 0.18, label = "italic(c)", parse = TRUE, size = 6) +
  annotate("text", x = d$x + 0.16, y = d$y - 0.1, label = "italic(d)", parse = TRUE, size = 6) +
  
  labs(title = "(c) General metric space") +
  coord_fixed(
    xlim = c(-1.75, 1.75),
    ylim = c(-1.75, 1.75),
    clip = "off"
  ) +
  common_theme


# Combine 


(p1 + p2 + p3) + plot_layout(ncol = 3, widths = c(1, 1, 1))



# ============================================================
# Figure 2. signed distance plot
# ============================================================


title_style <- element_text(
  hjust = 0.5, size = 14, family = "noto_math",
  margin = margin(t = 5, b = 10)
)

common_theme <- theme_void() +
  theme(
    plot.title = title_style,
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

scatter_theme <- theme_classic() +
  theme(
    plot.title = title_style,
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    axis.title = element_text(size = 13, family = "noto_math"),
    axis.text = element_text(size = 11, family = "noto_math"),
    legend.position = "none"
  )



norm_vec <- function(x) {sqrt(sum(x^2))}
normalize <- function(x) {x / norm_vec(x)}
clamp <- function(x, lower = -1, upper = 1) {pmin(pmax(x, lower), upper)}

cross_product <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

sphere_dist <- function(x, y) {acos(clamp(sum(x * y)))}

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

slerp <- function(x, y, t) {
  theta <- sphere_dist(x, y)
  if (theta < 1e-10) {
    out <- matrix(rep(x, times = length(t)), nrow = length(t), byrow = TRUE)
    colnames(out) <- c("V1", "V2", "V3")
    return(out)
  }
  out <- t(sapply(t, function(tt) {
    z <- sin((1 - tt) * theta) / sin(theta) * x +
      sin(tt * theta) / sin(theta) * y
    normalize(z)
  }))
  colnames(out) <- c("V1", "V2", "V3")
  out
}

get_smooth_geodesic <- function(x, y, n_grid = 300) {
  temp <- slerp(x, y, seq(0, 1, length.out = n_grid))
  data.frame(V1 = temp[, 1], V2 = temp[, 2], V3 = temp[, 3])
}


# Generate spherical observations


n <- 120

##dgp.setting <- 'two group'
dgp.setting <- 'periodic'
##dgp.setting <- 'other'

if (dgp.setting=='two group'){
  mu_true <- normalize(c(0.10, -0.05, 1.00))
  
  e1 <- c(1, 0, 0)
  e1 <- e1 - sum(e1 * mu_true) * mu_true
  e1 <- normalize(e1)
  
  e2 <- c(0, 1, 0)
  e2 <- e2 - sum(e2 * mu_true) * mu_true - sum(e2 * e1) * e1
  e2 <- normalize(e2)
  
  latent_group <- sample(c(-1, 1), n, replace = TRUE)
  z1 <- 0.45 * latent_group + rnorm(n, sd = 0.20)
  z2 <- rnorm(n, sd = 0.30)
  
  R_hat <- matrix(NA_real_, nrow = n, ncol = 3)
  
  for (i in seq_len(n)) {
    tangent_i <- z1[i] * e1 + z2[i] * e2
    R_hat[i, ] <- sphere_exp(tangent_i, mu_true)
  }
  colnames(R_hat) <- c("V1", "V2", "V3")  
} else if (dgp.setting=='periodic'){
  period <- 60
  mu0 <- normalize(c(0.10, -0.05, 1.00))
  e1 <- c(1, 0, 0)
  e1 <- e1 - sum(e1 * mu0) * mu0
  e1 <- normalize(e1)
  
  e2 <- c(0, 1, 0)
  e2 <- e2 - sum(e2 * mu0) * mu0 - sum(e2 * e1) * e1
  e2 <- normalize(e2)

  mu_phase <- matrix(NA_real_, nrow = period, ncol = 3)
  amp <- 0.35
  for (j in seq_len(period)) {
    
    theta_j <- 2 * pi * (j - 1) / period
    
    ## tangent direction determining periodic mean
    tangent_mu <- amp * (cos(theta_j) * e1 +sin(theta_j) * e2)
    
    mu_phase[j, ] <- sphere_exp(tangent_mu, mu0)
  }
  
  colnames(mu_phase) <- c("V1", "V2", "V3")

  phase <- ((seq_len(n) - 1) %% period) + 1
  
  mu_true <- mu_phase[phase, ]
  
  R_hat <- matrix(NA_real_, nrow = n, ncol = 3)
  
  for (i in seq_len(n)) {
    mu_i <- mu_true[i, ]
    a1 <- c(1, 0, 0)
    if (abs(sum(a1 * mu_i)) > 0.95) {
      a1 <- c(0, 1, 0)
    }
    
    a1 <- a1 - sum(a1 * mu_i) * mu_i
    a1 <- normalize(a1)
    
    a2 <- cross_product(mu_i, a1)
    a2 <- normalize(a2)
    latent_group <- sample(c(-1, 1), 1)
    
    z1 <- 0.45 * latent_group + rnorm(1, sd = 0.20)
    z2 <- rnorm(1, sd = 0.30)
    
    tangent_i <- z1 * a1 + z2 * a2
    
    R_hat[i, ] <- sphere_exp(tangent_i, mu_i)
  }
  R_hat=mu_true
  colnames(R_hat) <- c("V1", "V2", "V3")
} else {
  mu_true <- normalize(c(0.10, -0.05, 1.00))
  R_hat <- movMF::rmovMF(n, 10 * mu_true)
  colnames(R_hat) <- c("V1", "V2", "V3")
}


mu_hat <- sphere_intrinsic_mean(R_hat)


## Signed distance based on sample mean


distance_sample <- apply(R_hat, 1, function(x) sphere_dist(x, mu_hat))
j_minus_sample <- which.max(distance_sample^2)
a_minus_sample <- R_hat[j_minus_sample, ]
a_plus_sample <- sphere_exp(-sphere_log(a_minus_sample, mu_hat), mu_hat)

d_minus_sample <- apply(R_hat, 1, function(x) sphere_dist(x, a_minus_sample))
d_plus_sample <- apply(R_hat, 1, function(x) sphere_dist(x, a_plus_sample))

score_sample <- d_minus_sample^2 - d_plus_sample^2
sign_sample <- ifelse(score_sample < 0, -1, 1)
signed_distance_sample <- sign_sample * distance_sample

group_sample <- factor(
  ifelse(sign_sample < 0, "Negative", "Positive"),
  levels = c("Negative", "Positive")
)


## Misspecified reference mean


mu_ref <- sphere_intrinsic_mean(rbind(a_minus_sample, mu_hat))

distance_ref <- apply(R_hat, 1, function(x) sphere_dist(x, mu_ref))
j_minus_ref <- which.max(distance_ref^2)
a_minus_ref <- R_hat[j_minus_ref, ]
a_plus_ref <- sphere_exp(-sphere_log(a_minus_ref, mu_ref), mu_ref)

d_minus_ref <- apply(R_hat, 1, function(x) sphere_dist(x, a_minus_ref))
d_plus_ref <- apply(R_hat, 1, function(x) sphere_dist(x, a_plus_ref))

score_ref <- d_minus_ref^2 - d_plus_ref^2
sign_ref <- ifelse(score_ref < 0, -1, 1)
signed_distance_ref <- sign_ref * distance_ref

group_ref <- factor(
  ifelse(sign_ref < 0, "Negative", "Positive"),
  levels = c("Negative", "Positive")
)


## plot 


log_mat <- t(apply(R_hat, 1, function(x) sphere_log(x, mu_hat)))
Sigma_log <- crossprod(log_mat) / n
eig <- eigen(Sigma_log, symmetric = TRUE)

pc1 <- eig$vectors[, 1]
pc2 <- eig$vectors[, 2]

pc1 <- pc1 - sum(pc1 * mu_hat) * mu_hat
pc1 <- normalize(pc1)

pc2 <- pc2 - sum(pc2 * mu_hat) * mu_hat - sum(pc2 * pc1) * pc1
pc2 <- normalize(pc2)

view_angle <- 55 * pi / 180
view_direction <- normalize(cos(view_angle) * mu_hat + sin(view_angle) * pc2)

vertical_axis <- mu_hat - sum(mu_hat * view_direction) * view_direction
vertical_axis <- normalize(vertical_axis)

horizontal_axis <- cross_product(vertical_axis, view_direction)
horizontal_axis <- normalize(horizontal_axis)

if (sum(horizontal_axis * pc1) < 0) horizontal_axis <- -horizontal_axis

project_to_2d <- function(X) {
  if (is.null(dim(X))) X <- matrix(X, nrow = 1)
  data.frame(
    x = as.vector(X %*% horizontal_axis),
    y = as.vector(X %*% vertical_axis)
  )
}


R_proj <- project_to_2d(R_hat)

plot_sample <- data.frame(
  Index = seq_len(n),
  x = R_proj$x,
  y = R_proj$y,
  SignedDistance = signed_distance_sample,
  group = group_sample
)

plot_ref <- data.frame(
  Index = seq_len(n),
  SignedDistance = signed_distance_ref,
  group = group_ref
)

sample_mean_proj <- project_to_2d(mu_hat)
minus_sample_proj <- project_to_2d(a_minus_sample)
plus_sample_proj <- project_to_2d(a_plus_sample)

reference_dat <- data.frame(
  x = c(sample_mean_proj$x, minus_sample_proj$x, plus_sample_proj$x),
  y = c(sample_mean_proj$y, minus_sample_proj$y, plus_sample_proj$y),
  reference = c("Sample mean", "a-", "a+")
)



path_minus <- get_smooth_geodesic(a_minus_sample, mu_hat)
path_plus <- get_smooth_geodesic(mu_hat, a_plus_sample)
path_reference <- rbind(path_minus, path_plus[-1, ])
path_matrix <- as.matrix(path_reference[, c("V1", "V2", "V3")])
path_proj <- project_to_2d(path_matrix)


t_eq <- seq(0, pi, length.out = 1000)
eq_back <- data.frame(x = cos(t_eq), y = 0.30 * sin(t_eq))
eq_front <- data.frame(x = cos(t_eq + pi), y = 0.30 * sin(t_eq + pi))


group_cols <- c(
  "Negative" = "#2C7BB6",
  "Positive" = "#D7191C"
)

reference_cols <- c(
  "Sample mean" = "goldenrod4",
  "a-" = "#008837",
  "a+" = "lightyellow4"
)


## Figure 2.1: Sphere


p_sphere <- ggplot() +
  geom_circle(
    aes(x0 = 0, y0 = 0, r = 1),
    color = "black", linewidth = 0.7
  ) +
  geom_path(
    data = eq_back, aes(x = x, y = y),
    color = "gray72", linewidth = 0.6, linetype = "dashed"
  ) +
  geom_path(
    data = eq_front, aes(x = x, y = y),
    color = "gray38", linewidth = 0.6
  ) +
  geom_path(
    data = path_proj, aes(x = x, y = y),
    color = "gray35", linewidth = 0.75, linetype = "dashed"
  ) +
  geom_point(
    data = subset(plot_sample, group == "Negative"),
    aes(x = x, y = y),
    shape = 21, color = group_cols["Negative"], fill = NA,
    size = 2.6, stroke = 0.8, alpha = 0.9
  ) +
  geom_point(
    data = subset(plot_sample, group == "Positive"),
    aes(x = x, y = y),
    shape = 21, color = group_cols["Positive"], fill = NA,
    size = 2.6, stroke = 0.8, alpha = 0.9
  ) +
  
  geom_point(
    data = subset(reference_dat, reference == "Sample mean"),
    aes(x = x, y = y),
    shape = 16, size = 3.7, color = reference_cols["Sample mean"]
  ) +
  geom_point(
    data = subset(reference_dat, reference == "a-"),
    aes(x = x, y = y),
    shape = 15, size = 3.8, color = reference_cols["a-"]
  ) +
  geom_point(
    data = subset(reference_dat, reference == "a+"),
    aes(x = x, y = y),
    shape = 17, size = 4.0, color = reference_cols["a+"]
  ) +
  labs(title = "(a) Observations on sphere") +
  coord_fixed(
    xlim = c(-1.08, 1.08),
    ylim = c(-1.08, 1.08),
    clip = "off"
  ) +
  common_theme


## Figure 2.2: Signed distance using sample mean


p_sample <- ggplot(
  plot_sample,
  aes(x = Index, y = SignedDistance, color = group)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    shape = 21, fill = NA,
    size = 2.6, stroke = 0.8, alpha = 0.9
  ) +
  geom_point(
    data = plot_sample[j_minus_sample, ],
    aes(x = Index, y = SignedDistance),
    inherit.aes = FALSE,
    shape = 15, size = 3.8,
    color = reference_cols["a-"]
  ) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "(b) Distance plot using sample mean",
    x = "Index",
    y = "Signed distance"
  ) +
  scatter_theme


## Figure 2.3: Signed distance using reference mean


p_ref <- ggplot(
  plot_ref,
  aes(x = Index, y = SignedDistance, color = group)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    shape = 21, fill = NA,
    size = 2.6, stroke = 0.8, alpha = 0.9
  ) +
  geom_point(
    data = plot_ref[j_minus_ref, ],
    aes(x = Index, y = SignedDistance),
    inherit.aes = FALSE,
    shape = 15, size = 3.8,
    color = reference_cols["a-"]
  ) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "(c) Distance plot using misspecified mean",
    x = "Index",
    y = "Signed distance"
  ) +
  scatter_theme


legend_dat <- data.frame(
  x = c(2.2, 2.95, 3.75, 4.90, 5.55),
  y = rep(1, 5),
  label = c("Negative", "Positive", "Sample mean", "a-", "a+"),
  color = c("#2C7BB6", "#D7191C", "goldenrod4", "#008837", "lightyellow4"),
  shape = c(21, 21, 16, 15, 17)
)

p_legend <- ggplot() +
  geom_point(
    data = legend_dat[1:2, ],
    aes(x = x, y = y),
    shape = 21,
    color = legend_dat$color[1:2],
    fill = NA,
    size = 3.3,
    stroke = 0.9
  ) +
  geom_point(
    data = legend_dat[3:5, ],
    aes(x = x, y = y),
    color = legend_dat$color[3:5],
    shape = legend_dat$shape[3:5],
    size = 3.3
  ) +
  geom_text(
    data = legend_dat,
    aes(x = x + 0.15, y = y, label = label),
    hjust = 0,
    size = 4.2,
    family = "noto_math"
  ) +
  coord_cartesian(
    xlim = c(0.5, 7.5),
    ylim = c(0.5, 1.5),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

## Combine


top_row <- p_sphere + p_sample + p_ref +
  plot_layout(widths = c(0.95, 1.25, 1.25))

combined_plot <- top_row / p_legend +
  plot_layout(heights = c(1, 0.12))

combined_plot


## Diagnostics


cat("\nSample mean construction\n")
cat("Index of a_-:", j_minus_sample, "\n")
cat("d(a_-, sample mean) =",
    round(sphere_dist(a_minus_sample, mu_hat), 4), "\n")

cat("\nReference mean construction\n")
cat("Index of a_-:", j_minus_ref, "\n")
cat("d(a_-, reference mean) =",
    round(sphere_dist(a_minus_ref, mu_ref), 4), "\n")
cat("d(sample mean, reference mean) =",
    round(sphere_dist(mu_hat, mu_ref), 4), "\n")




# ================================================================
# Figure 3: distance-profile boxplots of residuals with planted outliers
# ================================================================


source("basic_function.R")

library(ggplot2)
library(cowplot)

n <- 200L
k <- 10L
delta <- 3
kappa <- 1.5
seed <- 1L
p_grid <- seq(0, 1, length.out = 501L)
vmf_kappa <- 10
outlier_label <- "Outlier"
out_dir <- "../outputs/boxplot_figure_R"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

.sim_parameters <- function(space, p_grid, kappa) {
  if (space == "euclidean") {
    Sigma <- matrix(0.3, 5L, 5L)
    diag(Sigma) <- 1
    
    return(list(
      mu_true = c(1, -1, 0.5, 2, -0.5),
      reference = c(-1, 0.5, 1, -0.5, 2),
      Sigma = Sigma,
      direction = rep(1 / sqrt(5), 5),
      s0 = sqrt(sum(diag(Sigma)))
    ))
  }
  
  if (space == "sphere") {
    V <- integrate(
      function(z) {
        acos(pmax(-1, pmin(1, z)))^2 *
          kappa * exp(kappa * (z - 1)) / (-expm1(-2 * kappa))
      },
      -1, 1,
      rel.tol = 1e-10,
      subdivisions = 1000L
    )$value
    
    return(list(
      mu_true = normalize(c(1, 0, 0)),
      reference = normalize(c(1, 1, 0)),
      direction = c(0, 0, 1),
      kappa = kappa,
      s0 = sqrt(V)
    ))
  }
  
  q_ref <- 0.15 + 0.7 * p_grid
  basis <- approx(
    p_grid, sin(2 * pi * p_grid),
    xout = q_ref, rule = 2
  )$y
  
  list(
    mu_true = p_grid,
    reference = q_ref,
    p_grid = p_grid,
    s0 = sqrt(0.05^2 / 3 + (0.04^2 / 3) * trapz(p_grid, basis^2))
  )
}

.sim_generate_null <- function(n, space, model) {
  if (space == "euclidean") {
    Y <- MASS::mvrnorm(n, mu = model$mu_true, Sigma = model$Sigma)
  } else if (space == "sphere") {
    Y <- movMF::rmovMF(n, model$kappa * model$mu_true)
  } else {
    A <- runif(n, -0.05, 0.05)
    B <- runif(n, -0.04, 0.04)
    
    Y <- outer(A, rep(1, length(model$p_grid))) +
      outer(B, sin(2 * pi * model$p_grid))
    Y <- sweep(Y, 2L, model$p_grid, "+")
  }
  
  residuals <- lapply(seq_len(n), function(i) {
    metric_sub(
      model$mu_true, Y[i, ], model$reference,
      space = space,
      wasserstein_grid = if (space == "wasserstein") {
        model$p_grid
      } else {
        seq(0, 1, length.out = 1001L)
      }
    )
  })
  
  do.call(rbind, residuals)
}

.sim_add_signal <- function(E0, g, delta, space, model) {
  shift <- delta * model$s0 * g
  
  if (space == "wasserstein") {
    return(sweep(E0, 1L, shift, "+"))
  }
  
  if (space == "euclidean") {
    return(E0 + outer(shift, model$direction))
  }
  
  shifted <- lapply(seq_len(nrow(E0)), function(i) {
    conditional_center <- sphere_exp(
      shift[i] * model$direction,
      model$reference
    )
    
    metric_sub(
      model$reference, E0[i, ], conditional_center,
      space = "sphere"
    )
  })
  
  do.call(rbind, shifted)
}

spaces <- c("euclidean", "sphere", "wasserstein")
space_titles <- c(
  euclidean = "Euclidean",
  sphere = "Sphere",
  wasserstein = "Wasserstein"
)

outliers <- seq_len(k)
signal <- as.numeric(seq_len(n) <= k)

set.seed(seed)

samples <- lapply(spaces, function(space) {
  model <- .sim_parameters(space, p_grid, vmf_kappa)
  E0 <- .sim_generate_null(n, space, model)
  .sim_add_signal(E0, signal, delta, space, model)
})

names(samples) <- spaces

column_widths <- c(1.11, 1, 1.065)
row_heights <- c(1, 0.08, 0.07, 1, 0.08)
fig_width_mm <- 192
fig_height_mm <- 113

columns <- lapply(seq_along(spaces), function(j) {
  space <- spaces[j]
  
  boxplot_figure(
    samples[[space]],
    space = space,
    title = space_titles[[space]],
    highlight_indices = outliers,
    highlight_label = outlier_label,
    kappa = kappa,
    wasserstein_grid = if (space == "wasserstein") p_grid else NULL,
    show_ylab = (j == 1L),
    show_legend = FALSE
  )
})

names(columns) <- spaces

aligned <- lapply(columns, function(col) {
  align_plots(
    col$panels$profile + theme(legend.position = "none"),
    col$panels$distance + theme(legend.position = "none"),
    align = "v",
    axis = "lr"
  )
})

panel_row <- function(row) {
  plot_grid(
    plotlist = lapply(aligned, `[[`, row),
    nrow = 1L,
    align = "h",
    axis = "tb",
    rel_widths = column_widths
  )
}

with_legend <- boxplot_figure(
  samples$euclidean,
  space = "euclidean",
  highlight_indices = outliers,
  highlight_label = outlier_label,
  kappa = kappa,
  show_legend = TRUE
)

legend_row <- function(row) {
  p <- with_legend$panels[[row]] +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
  g <- ggplotGrob(p)
  idx <- which(grepl("^guide-box", g$layout$name))
  idx <- idx[vapply(g$grobs[idx], inherits, logical(1), "gtable")]
  
  if (!length(idx)) {
    stop(
      "No legend was generated for panel ", row,
      ". Check aes(), show.legend, and guides() inside boxplot_figure()."
    )
  }
  
  legend <- g$grobs[[idx[1L]]]
  legend$vp <- grid::viewport(
    x = 0.5,
    y = 0.5,
    width = sum(legend$widths),
    height = sum(legend$heights),
    just = c(0.5, 0.5)
  )
  
  legend
}

figure <- plot_grid(
  panel_row(1L),
  legend_row(1L),
  NULL,
  panel_row(2L),
  legend_row(2L),
  ncol = 1L,
  rel_heights = row_heights
) +
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave(
  file.path(out_dir, "fig_boxplot.png"),
  plot = figure,
  width = fig_width_mm,
  height = fig_height_mm,
  units = "mm",
  dpi = 300,
  bg = "white"
)

## Summary 

summary_table <- do.call(rbind, lapply(spaces, function(space) {
  bp <- columns[[space]]$boxplot
  data.frame(space = space, n = n, k = k, delta = delta, kappa = kappa,
             flagged = sum(bp$flagged),
             outliers_flagged = sum(bp$flagged[outliers]),
             clean_flagged = sum(bp$flagged[-outliers]),
             fence = unname(bp$fence))
}))
write.csv(summary_table, file.path(out_dir, "results_figure.csv"), row.names = FALSE)
print(summary_table, row.names = FALSE, digits = 4)



