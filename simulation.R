source('basic_function.R')

library(movMF)
library(ggplot2)
library(dplyr)

# ================================================================
# ================================================================
# Part 1
# ================================================================
# ================================================================

rep_num <- 1000
sample_sizes <- c(1:10) * 100
result <- data.frame()

## ================================================================
## 1. Real
## ================================================================

d <- 5

mu_true <- c(1, -1, 0.5, 2, -0.5)
ref_point <- c(-1, 0.5, 1, -0.5, 2)
Sigma <- matrix(0.3, d, d)
diag(Sigma) <- 1

for (n_sample in sample_sizes) {
  dist_vec <- numeric(rep_num)
  for (rep in 1:rep_num) {
    set.seed(rep)
    y <- MASS::mvrnorm(n_sample, mu = mu_true, Sigma = Sigma)
    err_matrix <- matrix(NA, nrow = n_sample, ncol = d)
    
    for (i in 1:n_sample) {
      err_matrix[i, ] <- metric_sub(mu_true, y[i, ], ref_point,
                                    space = "euclidean")
    }
    
    err_mean <- colMeans(err_matrix)
    dist_vec[rep] <- sum((err_mean - ref_point)^2)
  }
  
  result <- rbind(result,data.frame(space = "Euclidean", n = n_sample, rep = 1:rep_num, dist2 = dist_vec))
}

## ================================================================
## 2. Sphere
## ================================================================

mu_true <- normalize(c(1, 0, 0))
ref_point <- normalize(c(1, 1, 0))
kappa <- 15

for (n_sample in sample_sizes) {
  dist_vec <- numeric(rep_num)
  for (rep in 1:rep_num) {
    set.seed(rep)
    y <- movMF::rmovMF(n_sample, kappa * mu_true)
    err_matrix <- matrix(NA, nrow = n_sample, ncol = 3)
    
    for (i in 1:n_sample) {
      err_matrix[i, ] <- metric_sub(mu_true, y[i, ], ref_point, space = "sphere")
    }
    err_mean <- sphere_intrinsic_mean(err_matrix)
    dist_vec[rep] <- sphere_dist(err_mean, ref_point)^2
  }
  result <- rbind(result, data.frame(space = "Sphere", n = n_sample, rep = 1:rep_num, dist2 = dist_vec))
}


## ================================================================
## 3. Wasserstein
## ================================================================

p_grid <- seq(0, 1, length.out = 501)

## true Frechet mean: Uniform(0, 1)
q_mu <- p_grid

## reference point: Uniform(0.15, 0.85)
q_ref <- 0.15 + 0.7 * p_grid

for (n_sample in sample_sizes) {
  dist_vec <- numeric(rep_num)
  
  for (rep in 1:rep_num) {
    set.seed(rep)
    
    A <- runif(n_sample, -0.05, 0.05)
    B <- runif(n_sample, -0.04, 0.04)
    
    q_y <- matrix(NA, nrow = n_sample, ncol = length(p_grid))
    
    for (i in 1:n_sample) {
      q_y[i, ] <- p_grid + A[i] + B[i] * sin(2 * pi * p_grid)
    }
    
    err_matrix <- matrix(NA, nrow = n_sample, ncol = length(p_grid))
    
    for (i in 1:n_sample) {
      err_matrix[i, ] <- metric_sub(q_mu, q_y[i, ], q_ref,
                                    space = "wasserstein",
                                    wasserstein_grid = p_grid)
    }
    
    err_mean <- wasserstein_mean(err_matrix)
    dist_vec[rep] <- wasserstein_dist(err_mean, q_ref, p_grid)^2
  }
  
  result <- rbind(result, data.frame(space = "Wasserstein", n = n_sample, rep = 1:rep_num, dist2 = dist_vec))
}


## ================================================================
## Summary
## ================================================================

result_mean <- aggregate(dist2 ~ space + n, result, mean)
result_sd <- aggregate(dist2 ~ space + n, result, sd)
result_summary <- merge(result_mean, result_sd,by = c("space", "n"),suffixes = c("_mean", "_sd"))
print(result_summary)

result_mean <- result %>%
  group_by(space, n) %>%
  summarise(MSE = mean(dist2), .groups = "drop") %>%
  mutate(space = recode(space, "Real" = "Euclidean"))

result_mean$space <- factor(
  result_mean$space,
  levels = c("Euclidean", "Sphere", "Wasserstein"),
  labels = c("(a) Euclidean space", "(b) Sphere", "(c) Wasserstein space")
)

ggplot(result_mean, aes(x = n, y = MSE)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(~ space, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = sample_sizes) +
  labs(x = "Sample size", y = "MSE") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 15),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    panel.background = element_blank(),
    plot.margin = margin(5, 15, 5, 5),
    legend.title = element_blank(),
    legend.position = "bottom"
  )


# ================================================================
# ================================================================
# Part 2.1 known residuals
# ================================================================
# ================================================================
basic_file <- "basic_function.R"  
sample_sizes <- c(200L, 500L, 1000L)
signal_strengths <- seq(0,1,by=0.1)
rep_num <- 1000L
alpha <- 0.05
seed <- 123L
metric_spaces <- c("euclidean","sphere","wasserstein")
regression_methods <- c("nw","categorical")
p_grid <- seq(0,1,length.out=501L)
kappa <- 10  
bandwidth_grid_constants <- exp(seq(log(0.1),log(1),length.out=15L))
bandwidth_exponent <- 1/6
available_cores <- parallel::detectCores(logical=FALSE)
if (is.na(available_cores)) available_cores <- 2L
n_workers <- max(1L,available_cores-1L)

.sim_parameters <- function(space,p_grid,kappa) {
  if (space=="euclidean") {
    Sigma <- matrix(0.3,5L,5L)
    diag(Sigma) <- 1
    return(list(mu_true=c(1,-1,0.5,2,-0.5),
                reference=c(-1,0.5,1,-0.5,2),Sigma=Sigma,
                direction=rep(1/sqrt(5),5),s0=sqrt(sum(diag(Sigma)))))
  }
  if (space=="sphere") {
    V <- integrate(function(z) {
      acos(pmax(-1,pmin(1,z)))^2 * kappa*exp(kappa*(z-1))/(-expm1(-2*kappa))
    },-1,1,rel.tol=1e-10,subdivisions=1000L)$value
    return(list(mu_true=normalize(c(1,0,0)),reference=normalize(c(1,1,0)),
                direction=c(0,0,1),kappa=kappa,s0=sqrt(V)))
  }
  q_ref <- 0.15+0.7*p_grid
  basis <- approx(p_grid,sin(2*pi*p_grid),xout=q_ref,rule=2)$y
  list(mu_true=p_grid,reference=q_ref,p_grid=p_grid,
       s0=sqrt(0.05^2/3+(0.04^2/3)*trapz(p_grid,basis^2)))
}

.sim_generate_null <- function(n,space,model) {
  if (space=="euclidean") {
    Y <- MASS::mvrnorm(n,mu=model$mu_true,Sigma=model$Sigma)
  } else if (space=="sphere") {
    Y <- movMF::rmovMF(n,model$kappa*model$mu_true)
  } else {
    A <- runif(n,-0.05,0.05)
    B <- runif(n,-0.04,0.04)
    Y <- outer(A,rep(1,length(model$p_grid))) +
      outer(B,sin(2*pi*model$p_grid))
    Y <- sweep(Y,2L,model$p_grid,"+")
  }
  residuals <- lapply(seq_len(n),function(i) {
    metric_sub(model$mu_true,Y[i,],model$reference,space=space,
               wasserstein_grid=if(space=="wasserstein") model$p_grid else
                 seq(0,1,length.out=1001L))
  })
  do.call(rbind,residuals)
}

.sim_add_signal <- function(E0,g,delta,space,model) {
  shift <- delta*model$s0*g
  if (space=="wasserstein") return(sweep(E0,1L,shift,"+"))
  if (space=="euclidean") return(E0+outer(shift,model$direction))
  shifted <- lapply(seq_len(nrow(E0)),function(i) {
    conditional_center <- sphere_exp(shift[i]*model$direction,model$reference)
    metric_sub(model$reference,E0[i,],conditional_center,space="sphere")
  })
  do.call(rbind,shifted)
}


.sim_one_replication <- function(task, model, signals, methods, alpha,
                                 bandwidth_grid_constants,
                                 bandwidth_exponent, p_grid,
                                 n_randomizations = 50L,
                                 cauchy_weights = NULL) {
  if (!is.numeric(n_randomizations) || length(n_randomizations) != 1L ||
      !is.finite(n_randomizations) || n_randomizations < 1 ||
      n_randomizations != floor(n_randomizations) ||
      n_randomizations > .Machine$integer.max) {
    stop("n_randomizations must be a positive integer.")
  }
  L <- as.integer(n_randomizations)
  n <- task$n
  space <- task$space
  U <- runif(n)
  E0 <- .sim_generate_null(n, space, model)
  h_grid <- bandwidth_grid_constants * n^(-bandwidth_exponent)
  rows <- vector("list", length(signals) * length(methods))
  cursor <- 0L
  for (method in methods) {
    X <- if (method == "nw") U else ifelse(U <= 0.5, -1, 1)
    g <- if (method == "nw") sqrt(3) * (2 * X - 1) else X
    for (delta in signals) {
      E <- .sim_add_signal(E0, g, delta, space, model)
      fit_options <- if (method == "nw") {
        list(bw = "CV", bwGrid = h_grid, kernel = "gauss",
             tol = 1e-6, maxit = 500L)
      } else list()
      out <- residual_test(
        E = E, X = X, muY = model$reference,
        space = space, method = method,
        alpha = alpha, optns = fit_options,
        wasserstein_grid = if (space == "wasserstein") p_grid else NULL,
        check = TRUE, return_fit = FALSE,
        n_randomizations = L, cauchy_weights = cauchy_weights
      )
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        space = space, method = method, n = n,
        signal = delta, replication = task$replication,
        reject = out$reject,
        Z = unname(out$statistic),
        p_value = out$p.value, Dn = out$Dn,
        sigmahat = out$sigmahat,
        bandwidth = if (method == "nw") unname(out$bandwidth) else NA_real_,
        n_randomizations = out$n_randomizations
      )
    }
  }
  do.call(rbind, rows)
}


run_residual_test_simulation <- function(
    basic_file="basic_function.R",
    sample_sizes=c(200L,500L,1000L),signals=seq(0,1,by=0.1),repetitions=200L,
    spaces=c("euclidean","sphere","wasserstein"),methods=c("nw","categorical"),
    alpha=0.05,seed=123L,n_workers=2L,
    p_grid=seq(0,1,length.out=501L),kappa=10,
    bandwidth_grid_constants=exp(seq(log(0.1),log(1),length.out=15L)),
    bandwidth_exponent=1/6) {
  packages <- c("MASS","movMF","foreach","doSNOW","ggplot2")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("Install required packages: ",paste(missing,collapse=", "))
  if (!file.exists(basic_file) && identical(basic_file,"basic_function.R") &&
      file.exists("basic_function(1).R")) basic_file <- "basic_function(1).R"
  basic_path <- normalizePath(basic_file,mustWork=TRUE)
  source(basic_path,local=.GlobalEnv)
  if (!length(sample_sizes) || any(!is.finite(sample_sizes)) ||
      any(sample_sizes<2 | sample_sizes!=floor(sample_sizes))) stop("Invalid sample sizes.")
  if (!length(signals) || any(!is.finite(signals)) || any(signals<0) ||
      !any(signals==0)) stop("Signals must be nonnegative and include zero.")
  if (length(repetitions)!=1L || !is.finite(repetitions) || repetitions<2 ||
      repetitions!=floor(repetitions)) stop("Invalid number of replications.")
  if (length(n_workers)!=1L || !is.finite(n_workers) || n_workers<1 ||
      n_workers!=floor(n_workers)) stop("n_workers must be a positive integer.")
  if (!length(spaces) || any(!spaces %in% c("euclidean","sphere","wasserstein"))) {
    stop("Cover euclidean, sphere, and wasserstein.")
  }
  if (!length(methods) || any(!methods %in% c("nw","categorical"))) stop("Invalid methods.")
  if (length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1) stop("Invalid alpha.")
  if (length(kappa)!=1L || !is.finite(kappa) || kappa<=0) stop("kappa must be positive.")
  if (length(p_grid)<2L || any(!is.finite(p_grid)) || any(diff(p_grid)<=0) ||
      abs(p_grid[1L])>1e-12 || abs(tail(p_grid,1L)-1)>1e-12) stop("Invalid quantile grid.")
  if (!is.numeric(bandwidth_grid_constants) || !length(bandwidth_grid_constants) ||
      any(!is.finite(bandwidth_grid_constants)) || any(bandwidth_grid_constants<=0) ||
      length(unique(bandwidth_grid_constants))<2L || length(bandwidth_exponent)!=1L ||
      !is.finite(bandwidth_exponent) || bandwidth_exponent<=0 ||
      bandwidth_exponent>=1/4) {
    stop("Use at least two positive CV grid constants and an exponent in (0,1/4).")
  }
  bandwidth_grid_constants <- sort(unique(bandwidth_grid_constants))
  sample_sizes <- sort(unique(sample_sizes))
  shape_values <- c(16,17,15,3,7,8,0,1,2,4,5,6,9,10,11,12,13,14,18:25)
  if (length(sample_sizes)>length(shape_values)) {
    stop("At most 26 sample sizes are supported by the distinct point shapes.")
  }
  shape_values <- setNames(shape_values[seq_along(sample_sizes)],as.character(sample_sizes))
  signals <- sort(unique(signals))
  spaces <- unique(spaces)
  methods <- unique(methods)
  models <- setNames(lapply(spaces,.sim_parameters,p_grid=p_grid,kappa=kappa),spaces)
  max_g <- if("nw" %in% methods) sqrt(3) else 1
  if ("sphere" %in% spaces && max(signals)*models$sphere$s0*max_g>=pi/2) {
    stop("Reduce the spherical signal grid to keep center shifts below pi/2.")
  }
  tasks <- expand.grid(space=spaces,n=sample_sizes,
                       replication=seq_len(repetitions),stringsAsFactors=FALSE)
  n_tasks <- nrow(tasks)
  n_workers <- min(as.integer(n_workers),n_tasks)
  
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
  cl <- NULL
  pb <- NULL
  on.exit({
    if (!is.null(pb)) close(pb)
    if (!is.null(cl)) snow::stopCluster(cl)
    foreach::registerDoSEQ()
    do.call(RNGkind,as.list(old_kind))
    if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else {
      if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) {
        rm(".Random.seed",envir=.GlobalEnv)
      }
    }
  },add=TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list",n_tasks)
  streams[[1L]] <- .Random.seed
  if (n_tasks>1L) for (i in 2L:n_tasks) {
    streams[[i]] <- parallel::nextRNGStream(streams[[i-1L]])
  }
  
  cl <- snow::makeCluster(n_workers,type="SOCK")
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl,function(basic_path) {
    source(basic_path,local=.GlobalEnv)
    RNGkind("L'Ecuyer-CMRG")
    NULL
  },basic_path))
  simulation_helpers <- list(.sim_one_replication=.sim_one_replication,
                             .sim_generate_null=.sim_generate_null,
                             .sim_add_signal=.sim_add_signal)
  invisible(snow::clusterCall(cl,function(helpers) {
    list2env(helpers,envir=.GlobalEnv)
    NULL
  },simulation_helpers))
  pb <- utils::txtProgressBar(min=0,max=n_tasks,style=3)
  progress <- function(n) utils::setTxtProgressBar(pb,n)
  snow_options <- list(progress=progress,preschedule=FALSE)
  `%dopar%` <- foreach::`%dopar%`
  result <- foreach::foreach(task_id=seq_len(n_tasks),.combine=rbind,
                             .multicombine=TRUE,.inorder=FALSE,.errorhandling="stop",
                             .packages=c("MASS","movMF"),.options.snow=snow_options,
                             .noexport=c("residual_test",".sim_one_replication",
                                         ".sim_generate_null",".sim_add_signal")) %dopar% {
                                           assign(".Random.seed",streams[[task_id]],envir=.GlobalEnv)
                                           task <- tasks[task_id,,drop=FALSE]
                                           .sim_one_replication(task,models[[task$space]],signals,methods,alpha,
                                                                bandwidth_grid_constants,bandwidth_exponent,p_grid)
                                         }
  close(pb)
  pb <- NULL
  snow::stopCluster(cl)
  cl <- NULL
  foreach::registerDoSEQ()
  rownames(result) <- NULL
  result <- result[order(match(result$space,spaces),match(result$method,methods),
                         result$n,result$signal,result$replication),]
  rownames(result) <- NULL
  groups <- split(result,interaction(result$space,result$method,result$n,
                                     result$signal,drop=TRUE))
  result_summary <- do.call(rbind,lapply(groups,function(z) {
    R <- nrow(z)
    rejected <- sum(z$reject)
    p <- rejected/R
    selected_h <- z$bandwidth[is.finite(z$bandwidth)]
    data.frame(space=z$space[1L],method=z$method[1L],n=z$n[1L],signal=z$signal[1L],
               repetitions=R,rejections=rejected,percent_rejected=100*p,
               mc_se_percent=100*sqrt(p*(1-p)/R),
               bandwidth_mean=if(length(selected_h)) mean(selected_h) else NA_real_,
               bandwidth_sd=if(length(selected_h)>1L) sd(selected_h) else NA_real_,
               bandwidth_min=if(length(selected_h)) min(selected_h) else NA_real_,
               bandwidth_max=if(length(selected_h)) max(selected_h) else NA_real_)
  }))
  result_summary <- result_summary[order(match(result_summary$method,methods),
                                         match(result_summary$space,spaces),result_summary$n,result_summary$signal),]
  rownames(result_summary) <- NULL
  if (any(result_summary$repetitions!=repetitions)) {
    stop("Unexpected replication counts; no failed runs are silently discarded.")
  }
  plot_data <- result_summary
  plot_data$space <- factor(plot_data$space,
                            levels=c("euclidean","sphere","wasserstein"),
                            labels=c("Euclidean space","Sphere","Wasserstein space"))
  plot_data$method <- factor(plot_data$method,levels=c("nw","categorical"),
                             labels=c("Continuous X","Discrete X"))
  plot_data$sample_size <- factor(plot_data$n,levels=sample_sizes)
  p <- ggplot2::ggplot(plot_data,ggplot2::aes(x=signal,y=percent_rejected,
                                              colour=sample_size,linetype=sample_size,
                                              shape=sample_size,group=sample_size)) +
    ggplot2::geom_hline(yintercept=100*alpha,linetype="dashed",colour="grey45") +
    ggplot2::geom_line(linewidth=0.8) + ggplot2::geom_point(size=2) +
    ggplot2::scale_shape_manual(values=shape_values) +
    ggplot2::facet_grid(
      method~space,
      labeller=ggplot2::labeller(
        method=ggplot2::as_labeller(
          c("Continuous X"="plain(Continuous)~italic(X)",
            "Discrete X"="plain(Discrete)~italic(X)"),
          default=ggplot2::label_parsed))) +
    ggplot2::scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20)) +
    ggplot2::labs(x="Signal strength",y="Percentage of rejections (%)",colour="Sample size",
                  linetype="Sample size",shape="Sample size") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="bottom",legend.key.width=grid::unit(2,"cm"),
                   strip.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_rect(fill=NA,colour="black"))
  list(replications=result,summary=result_summary,plot=p,
       settings=list(sample_sizes=sample_sizes,signals=signals,repetitions=repetitions,
                     spaces=spaces,methods=methods,alpha=alpha,seed=seed,n_workers=n_workers,
                     kappa=kappa,bandwidth_selection="leave-one-out CV",
                     bandwidth_grid_constants=bandwidth_grid_constants,
                     bandwidth_exponent=bandwidth_exponent))
}

simulation <- run_residual_test_simulation(
  basic_file=basic_file,sample_sizes=sample_sizes,
  signals=signal_strengths,repetitions=rep_num,spaces=metric_spaces,
  methods=regression_methods,alpha=alpha,seed=seed,n_workers=n_workers,
  p_grid=p_grid,kappa=kappa,bandwidth_grid_constants=bandwidth_grid_constants,
  bandwidth_exponent=bandwidth_exponent)
result <- simulation$replications
result_summary <- simulation$summary
print(result_summary,row.names=FALSE,digits=4)
print(simulation$plot)



# ================================================================
# ================================================================
# Part 2.2 fitted residuals
# ================================================================
# ================================================================

.frsim_load_basic <- function(basic_file) {
  path <- normalizePath(basic_file, mustWork=TRUE)
  source(path, local=.GlobalEnv)
  required <- c("metric_sub", "object_frechet_mean", "object_dist",
                "make_metric_residuals", "NWFrechetWeights", "NWSpheReg",
                "NWWassReg", "residual_test", ".nw_cv")
  if (!all(vapply(required, exists, logical(1), mode="function"))) {
    stop("The supplied basic function file is missing required functions.")
  }
  path
}

.frsim_g <- function(x, design) {
  if (design == "nw") sqrt(3)*(2*x-1) else x
}

.frsim_model <- function(space, p_grid, sphere_radius) {
  if (space == "euclidean") {
    Sigma <- matrix(.3, 5L, 5L)
    diag(Sigma) <- 1
    return(list(space=space, reference=c(-1,.5,1,-.5,2), Sigma=Sigma,
                v=rep(1/sqrt(5),5), s0=sqrt(sum(diag(Sigma)))))
  }
  if (space == "sphere") {
    return(list(space=space, reference=normalize(c(1,1,0)),
                v=c(0,0,1), v2=normalize(c(-1,1,0)),
                radius=sphere_radius, s0=sphere_radius/sqrt(2)))
  }
  b1 <- sin(pi*p_grid)/pi
  b2 <- sin(2*pi*p_grid)/(2*pi)
  b1[c(1L,length(b1))] <- b2[c(1L,length(b2))] <- 0
  a_bound <- .05
  b_bound <- .04
  s0 <- sqrt(a_bound^2/(6*pi^2) + b_bound^2/(24*pi^2))
  list(space=space, reference=.15+.7*p_grid, p_grid=p_grid,
       b1=b1, b2=b2, a_bound=a_bound, b_bound=b_bound, s0=s0,
       # v1 = sqrt(2)*sin(pi*u), with L2 norm exactly 1.
       tau=s0*pi*sqrt(2))
}

.frsim_mean <- function(x, design, model) {
  g <- .frsim_g(x, design)
  if (model$space == "euclidean") {
    return(sweep(outer(model$s0*g,model$v),2L,model$reference,"+"))
  }
  if (model$space == "sphere") {
    return(do.call(rbind,lapply(g,function(z) {
      sphere_exp(model$s0*z*model$v,model$reference)
    })))
  }
  sweep(outer(model$tau*g,model$b1),2L,model$reference,"+")
}

.frsim_response <- function(n, design, model, X=NULL, noise=NULL) {
  if (is.null(X)) X <- if (design == "nw") runif(n) else 2*rbinom(n,1,.5)-1
  if (length(X)!=n) stop("The response-generating predictor must have length n.")
  M <- .frsim_mean(X,design,model)
  if (model$space == "euclidean") {
    if (is.null(noise)) noise <- matrix(rnorm(n*5L),n,5L) %*% chol(model$Sigma)
    Y <- M+noise
  } else if (model$space == "sphere") {
    if (is.null(noise)) {
      radius <- model$radius*sqrt(runif(n))
      angle <- runif(n,0,2*pi)
      noise <- outer(radius*cos(angle),model$v) +
        outer(radius*sin(angle),model$v2)
    }
    Y <- do.call(rbind,lapply(seq_len(n),function(i) {
      Z <- sphere_exp(noise[i,],model$reference)
      metric_sub(model$reference,Z,M[i,],space="sphere")
    }))
  } else {
    if (is.null(noise)) {
      A <- runif(n,-model$a_bound,model$a_bound)
      B <- runif(n,-model$b_bound,model$b_bound)
      noise <- outer(A,model$b1)+outer(B,model$b2)
    }
    Y <- M+noise
  }
  list(X=X,Y=Y,true_mean=M,noise=noise)
}

.frsim_fit <- function(xin,Y,xout,space,design,p_grid,optns) {
  if (design == "categorical") {
    categories <- sort(unique(xin))
    idx <- match(xout,categories)
    if (anyNA(idx)) stop("A training category is absent; increase n_train.")
    means <- lapply(categories,function(x) {
      object_frechet_mean(Y[xin==x,,drop=FALSE],space=space,
                          wasserstein_grid=p_grid)
    })
    return(list(fitted=do.call(rbind,means)[idx,,drop=FALSE],
                bandwidth=NA_real_,cv_edge=NA))
  }
  if (space == "sphere") {
    fit <- NWSpheReg(xin,Y,xout,optns=optns)
    if (any(!fit$converged)) stop("First-stage spherical regression failed to converge.")
    fitted <- fit$yout
    bandwidth <- fit$optns$bw
    cv <- fit$cv
  } else if (space == "wasserstein") {
    optns$qSup <- p_grid
    fit <- NWWassReg(xin,Y,xout,optns=optns)
    fitted <- fit$qout
    bandwidth <- fit$optns$bw
    cv <- fit$cv
  } else {
    bandwidth <- optns$bw
    cv <- NULL
    if (identical(bandwidth,"CV")) {
      cv <- .nw_cv(matrix(xin,ncol=1L),optns$bwGrid,"gauss",function(W) {
        mean(rowSums((W %*% Y-Y)^2))
      })
      bandwidth <- cv$bw
    }
    W <- NWFrechetWeights(xin,xout,bw=bandwidth,kernel="gauss")
    fitted <- W %*% Y
  }
  edge <- if (is.null(cv)) NA else {
    selected <- which.min(cv$scores)
    selected %in% c(1L,length(cv$scores))
  }
  list(fitted=fitted,bandwidth=as.numeric(bandwidth),cv_edge=edge)
}

.frsim_mse <- function(A,B,space,p_grid) {
  mean(vapply(seq_len(nrow(A)),function(i) {
    object_dist(A[i,],B[i,],space=space,wasserstein_grid=p_grid)^2
  },numeric(1)))
}

.frsim_options <- function(n,constants,exponent,bw_mode,controls) {
  h <- sort(unique(constants))*n^(-exponent)
  c(list(bw=if (bw_mode=="CV") "CV" else h[ceiling(length(h)/2)],
         bwGrid=h,kernel="gauss"),controls)
}

.frsim_one_task <- function(task, settings, stream) {
  assign(".Random.seed", stream, envir = .GlobalEnv)
  n <- as.integer(task$n)
  n_train <- as.integer(ceiling(n * settings$training_multiplier))
  space <- task$space
  design <- task$design
  model <- settings$models[[space]]
  n_randomizations <- settings$n_randomizations
  if (is.null(n_randomizations)) n_randomizations <- 1L
  train <- .frsim_response(n_train, design, model)
  evaluation <- .frsim_response(n, design, model)
  reference <- object_frechet_mean(
    train$Y, space = space, wasserstein_grid = settings$p_grid
  )
  oracle <- do.call(rbind, make_metric_residuals(
    evaluation$Y, evaluation$true_mean, space = space, muY = reference,
    wasserstein_grid = settings$p_grid
  )$residuals)
  test_options <- .frsim_options(
    n, settings$test_bw_constants, settings$test_bw_exponent,
    settings$bw_mode, settings$controls
  )
  fit_options <- .frsim_options(
    n_train, settings$fit_bw_constants, settings$fit_bw_exponent,
    settings$bw_mode, settings$controls
  )
  fit <- .frsim_fit(
    train$X, train$Y, evaluation$X, space, design,
    settings$p_grid, fit_options
  )
  baseline_mse <- .frsim_mse(
    fit$fitted, evaluation$true_mean, space, settings$p_grid
  )
  corruption_order <- sample.int(n)
  altered_predictors <- if (design == "nw") {
    (evaluation$X + settings$continuous_shift) %% 1
  } else -evaluation$X
  row <- function(test, rate, type, mse = 0, residual_mse = 0,
                  clean_mse = 0, bw_fit = NA_real_, cv_edge = NA,
                  realized = 0) {
    data.frame(
      space = space, design = design, n = n, n_train = n_train,
      corruption = rate, realized_corruption = realized,
      continuous_shift = if (design == "nw") settings$continuous_shift else NA_real_,
      replication = task$replication, residual_type = type,
      reject = test$reject, Z = unname(test$statistic), p_value = test$p.value,
      Dn = test$Dn, sqrt_n_Dn = sqrt(n) * test$Dn, sigmahat = test$sigmahat,
      prediction_mse = mse, clean_prediction_mse = clean_mse,
      residual_mse = residual_mse,
      bandwidth_fit = bw_fit,
      bandwidth_test = if (is.null(test$bandwidth)) NA_real_ else unname(test$bandwidth),
      fit_cv_edge = cv_edge
    )
  }
  oracle_test <- residual_test(
    oracle, evaluation$X, muY = reference,
    space = space, method = design, alpha = settings$alpha,
    optns = test_options, wasserstein_grid = settings$p_grid,
    n_randomizations = n_randomizations,
    cauchy_weights = settings$cauchy_weights
  )
  rows <- list(row(oracle_test, rate = 0, type = "oracle"))
  for (p in settings$corruption_rates) {
    count <- floor(p * n + 1e-9)
    chosen <- head(corruption_order, count)
    X_response <- evaluation$X
    X_response[chosen] <- altered_predictors[chosen]
    if (sum(X_response != evaluation$X) != count) {
      stop("Corruption count does not match its definition.")
    }
    evaluation_p <- .frsim_response(
      n, design, model, X = X_response, noise = evaluation$noise
    )
    E <- do.call(rbind, make_metric_residuals(
      evaluation_p$Y, fit$fitted, space = space, muY = reference,
      wasserstein_grid = settings$p_grid
    )$residuals)
    oracle_p <- if (p == 0) oracle else do.call(rbind, make_metric_residuals(
      evaluation_p$Y, evaluation_p$true_mean, space = space, muY = reference,
      wasserstein_grid = settings$p_grid
    )$residuals)
    test <- residual_test(
      E, evaluation$X, muY = reference, space = space,
      method = design, alpha = settings$alpha, optns = test_options,
      wasserstein_grid = settings$p_grid,
      n_randomizations = n_randomizations,
      cauchy_weights = settings$cauchy_weights
    )
    rows[[length(rows) + 1L]] <- row(
      test, p, "fitted",
      mse = .frsim_mse(fit$fitted, evaluation_p$true_mean, space, settings$p_grid),
      residual_mse = .frsim_mse(E, oracle_p, space, settings$p_grid),
      clean_mse = baseline_mse, bw_fit = fit$bandwidth,
      cv_edge = fit$cv_edge, realized = count / n
    )
  }
  do.call(rbind, rows)
}

.frsim_summary <- function(result) {
  grouping <- interaction(result$space,result$design,result$n,
                          result$corruption,result$residual_type,drop=TRUE)
  groups <- split(result,grouping)
  do.call(rbind,lapply(groups,function(z) {
    R <- nrow(z)
    q <- mean(z$reject)
    # Wilson binomial interval, including boundary probabilities 0 and 1.
    t <- qnorm(.975)
    denominator <- 1+t^2/R
    center <- (q+t^2/(2*R))/denominator
    half <- t*sqrt(q*(1-q)/R+t^2/(4*R^2))/denominator
    mean_na <- function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
    data.frame(space=z$space[1],design=z$design[1],n=z$n[1],
               n_train=z$n_train[1],corruption=z$corruption[1],
               realized_corruption=z$realized_corruption[1],
               continuous_shift=z$continuous_shift[1],
               residual_type=z$residual_type[1],repetitions=R,rejections=sum(z$reject),
               percent_rejected=100*q,mc_se_percent=100*sqrt(q*(1-q)/R),
               lower_95=100*(center-half),upper_95=100*(center+half),
               prediction_mse=mean(z$prediction_mse),
               clean_prediction_mse=mean(z$clean_prediction_mse),
               residual_mse=mean(z$residual_mse),
               sqrt_n_Dn_mean=mean(z$sqrt_n_Dn),
               bandwidth_fit_mean=mean_na(z$bandwidth_fit),
               bandwidth_test_mean=mean_na(z$bandwidth_test),
               fit_cv_edge_fraction=mean_na(z$fit_cv_edge))
  }))
}

run_fitted_residual_simulation <- function(
    basic_file="basic_function.R",sample_sizes=c(200L,500L,1000L),
    corruption_rates=seq(0,1,.1),repetitions=200L,
    spaces=c("euclidean","sphere","wasserstein"),
    designs=c("nw","categorical"),training_multiplier=1,
    alpha=.05,seed=123L,n_workers=1L,
    p_grid=seq(0,1,length.out=501L),sphere_radius=.35,
    fit_bw_constants=exp(seq(log(.1),log(1),length.out=15L)),
    test_bw_constants=exp(seq(log(.1),log(1),length.out=15L)),
    fit_bw_exponent=1/5,test_bw_exponent=1/6,bw_mode=c("CV","fixed"),
    controls=list(tol=1e-6,maxit=500L,nStarts=1L),progress=TRUE,
    continuous_shift=.5) {
  packages <- c("foreach","doSNOW","snow","ggplot2")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) {
    stop("Install required packages: ",paste(missing,collapse=", "))
  }
  basic_path <- .frsim_load_basic(basic_file)
  bw_mode <- match.arg(bw_mode)
  integer_ok <- function(x,minimum=1L) length(x)>0 && is.numeric(x) &&
    all(is.finite(x)) && all(x>=minimum & x==floor(x))
  if (!integer_ok(sample_sizes,20L)) stop("Use integer sample sizes >=20.")
  if (!integer_ok(repetitions,2L) || length(repetitions)!=1L ||
      !integer_ok(n_workers) || length(n_workers)!=1L) stop("Invalid repetitions/workers.")
  if (!is.logical(progress) || length(progress)!=1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.")
  }
  if (length(unique(sample_sizes))>20L) stop("At most 20 sample sizes are supported.")
  if (!length(corruption_rates) || any(!is.finite(corruption_rates)) ||
      any(corruption_rates<0 | corruption_rates>1) ||
      !any(corruption_rates==0)) stop("Corruption rates must be in [0,1] and include zero.")
  if (!length(spaces) || any(!spaces %in% c("euclidean","sphere","wasserstein")) ||
      !length(designs) || any(!designs %in% c("nw","categorical"))) {
    stop("Invalid metric spaces/predictor designs.")
  }
  scalar_between <- function(x,a,b) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x>a && x<b
  if (!scalar_between(alpha,0,1) ||
      !scalar_between(training_multiplier,0,Inf)) stop("Invalid alpha/training multiplier.")
  if (any(ceiling(sample_sizes*training_multiplier)<20)) stop("n_train must be >=20.")
  if (!integer_ok(seed,0) || length(seed)!=1L) stop("Use a nonnegative integer seed.")
  if (!scalar_between(continuous_shift,0,1)) {
    stop("continuous_shift must be in (0,1); a shift of 0 or 1 would leave X unchanged after wrapping.")
  }
  if (!scalar_between(sphere_radius,0,pi/4)) {
    stop("sphere_radius must be in (0, pi/4).")
  }
  max_abs_g <- if ("nw" %in% designs) sqrt(3) else 1
  if ("sphere" %in% spaces &&
      sphere_radius*(1+max_abs_g/sqrt(2))>=pi/4) {
    stop("Decrease sphere_radius: all spherical responses must lie in a ball of radius < pi/4.")
  }
  if (length(p_grid)<3L || any(!is.finite(p_grid)) || any(diff(p_grid)<=0) ||
      abs(p_grid[1])>1e-12 || abs(tail(p_grid,1)-1)>1e-12) stop("Invalid probability grid.")
  for (values in list(fit_bw_constants,test_bw_constants)) {
    if (!is.numeric(values) || length(unique(values))<2L ||
        any(!is.finite(values)) || any(values<=0)) stop("Supply >=2 positive CV constants.")
  }
  if (!scalar_between(fit_bw_exponent,1/6,1/4) ||
      !scalar_between(test_bw_exponent,0,1/4)) {
    stop("Use fit exponent in (1/6,1/4), and test exponent in (0,1/4).")
  }
  .nw_controls(controls)
  spaces <- unique(spaces)
  designs <- unique(designs)
  sample_sizes <- sort(unique(sample_sizes))
  models <- setNames(lapply(spaces,.frsim_model,p_grid=p_grid,
                            sphere_radius=sphere_radius),spaces)
  if ("wasserstein" %in% spaces) {
    model <- models[["wasserstein"]]
    derivative_lower_bound <- .7-model$tau*max_abs_g-model$a_bound-model$b_bound
    if (derivative_lower_bound<=0) {
      stop("Wasserstein quantiles must be strictly increasing on the predictor domain [0,1].")
    }
  }
  settings <- list(sample_sizes=sample_sizes,
                   corruption_rates=sort(unique(corruption_rates)),repetitions=repetitions,
                   spaces=spaces,designs=designs,training_multiplier=training_multiplier,
                   response_strength=1,corruption_scope="evaluation response-generating predictors",
                   continuous_shift=continuous_shift,
                   continuous_corruption="fixed shift modulo 1",
                   fitting_predictors="clean",
                   evaluation_predictors="clean",
                   alpha=alpha,seed=seed,p_grid=p_grid,sphere_radius=sphere_radius,
                   fit_bw_constants=sort(unique(fit_bw_constants)),
                   test_bw_constants=sort(unique(test_bw_constants)),
                   fit_bw_exponent=fit_bw_exponent,test_bw_exponent=test_bw_exponent,
                   bw_mode=bw_mode,controls=controls,models=models)
  tasks <- expand.grid(space=spaces,design=designs,n=sample_sizes,
                       replication=seq_len(repetitions),stringsAsFactors=FALSE)
  n_tasks <- nrow(tasks)
  n_workers <- min(as.integer(n_workers),n_tasks)
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed",envir=.GlobalEnv)
  cl <- NULL
  pb <- NULL
  on.exit({
    if (!is.null(pb)) close(pb)
    if (!is.null(cl)) snow::stopCluster(cl)
    foreach::registerDoSEQ()
    do.call(RNGkind,as.list(old_kind))
    if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else {
      if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
        rm(".Random.seed",envir=.GlobalEnv)
    }
  },add=TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list",n_tasks)
  streams[[1L]] <- .Random.seed
  if (n_tasks>1L) for (i in 2:n_tasks) {
    streams[[i]] <- parallel::nextRNGStream(streams[[i-1L]])
  }
  cl <- snow::makeCluster(n_workers,type="SOCK")
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl,function(path) {
    source(path,local=.GlobalEnv)
    RNGkind("L'Ecuyer-CMRG")
    NULL
  },basic_path))
  helper_env <- environment(.frsim_one_task)
  helpers <- ls(helper_env,pattern="^\\.frsim_",all.names=TRUE)
  snow::clusterExport(cl,helpers,envir=helper_env)
  snow_options <- list(preschedule=FALSE)
  if (progress) {
    pb <- utils::txtProgressBar(min=0,max=n_tasks,style=3)
    progress_callback <- function(n) utils::setTxtProgressBar(pb,n)
    snow_options$progress <- progress_callback
  }
  `%dopar%` <- foreach::`%dopar%`
  result <- foreach::foreach(task_id=seq_len(n_tasks),.combine=rbind,
                             .multicombine=TRUE,.inorder=FALSE,.errorhandling="stop",
                             .options.snow=snow_options,.noexport=helpers) %dopar% {
                               .frsim_one_task(tasks[task_id,,drop=FALSE],settings,streams[[task_id]])
                             }
  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }
  snow::stopCluster(cl)
  cl <- NULL
  foreach::registerDoSEQ()
  result <- result[order(result$space,result$design,result$n,
                         result$corruption,result$residual_type,result$replication),]
  rownames(result) <- NULL
  summary <- .frsim_summary(result)
  summary <- summary[order(summary$space,summary$design,summary$n,
                           summary$corruption,summary$residual_type),]
  rownames(summary) <- NULL
  if (any(summary$repetitions!=repetitions)) stop("Unexpected replication counts.")
  settings$n_workers <- n_workers
  settings$parallel_backend <- "foreach / doSNOW"
  settings$progress_bar <- progress
  settings$models <- NULL
  p <- .frsim_plot(summary,settings)
  structure(list(replications=result,summary=summary,plot=p,settings=settings),
            class="fitted_residual_simulation")
}

.frsim_plot <- function(result_summary,settings) {
  plot_data <- result_summary[result_summary$residual_type=="fitted",,drop=FALSE]
  plot_data$signal <- plot_data$corruption
  plot_data$space <- factor(plot_data$space,
                            levels=c("euclidean","sphere","wasserstein"),
                            labels=c("Euclidean space","Sphere","Wasserstein space"))
  plot_data$method <- factor(plot_data$design,levels=c("nw","categorical"),
                             labels=c("Continuous X","Discrete X"))
  plot_data$sample_size <- factor(plot_data$n,levels=settings$sample_sizes)
  shape_values <- c(16,17,15,3,7,8,0,1,2,4,5,6,9,10,11,12,13,14,18:25)
  shape_values <- setNames(shape_values[seq_along(settings$sample_sizes)],
                           as.character(settings$sample_sizes))
  ggplot2::ggplot(plot_data,ggplot2::aes(x=signal,y=percent_rejected,
                                         colour=sample_size,linetype=sample_size,shape=sample_size,group=sample_size)) +
    ggplot2::geom_hline(yintercept=100*settings$alpha,linetype="dashed",colour="grey45") +
    ggplot2::geom_line(linewidth=.8) + ggplot2::geom_point(size=2) +
    ggplot2::scale_shape_manual(values=shape_values) +
    ggplot2::facet_grid(method~space,
                        labeller=ggplot2::labeller(method=ggplot2::as_labeller(
                          c("Continuous X"="plain(Continuous)~italic(X)",
                            "Discrete X"="plain(Discrete)~italic(X)"),
                          default=ggplot2::label_parsed))) +
    ggplot2::scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20)) +
    ggplot2::labs(x="Signal strength",y="Percentage of rejections (%)",
                  colour="Sample size",linetype="Sample size",shape="Sample size") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="bottom",legend.key.width=grid::unit(2,"cm"),
                   strip.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_rect(fill=NA,colour="black"))
}

plot_fitted_residual_simulation <- function(simulation) {
  p <- simulation$plot
  if (is.null(p)) p <- .frsim_plot(simulation$summary,simulation$settings)
  print(p)
  invisible(p)
}

validate_fitted_residual_design <- function(basic_file="basic_function(2).R",
                                            p_grid=seq(0,1,length.out=501L),
                                            continuous_shift=.5,sphere_radius=.35) {
  if (!is.numeric(continuous_shift) || length(continuous_shift)!=1L ||
      !is.finite(continuous_shift) || continuous_shift<=0 || continuous_shift>=1) {
    stop("continuous_shift must be in (0,1).")
  }
  .frsim_load_basic(basic_file)
  report <- list()
  add <- function(check,error,tolerance) {
    report[[length(report)+1L]] <<- data.frame(check=check,error=error,
                                               tolerance=tolerance,passed=error<=tolerance)
  }
  for (space in c("euclidean","sphere","wasserstein")) {
    model <- .frsim_model(space,p_grid,sphere_radius)
    m <- drop(.frsim_mean((.8+continuous_shift) %% 1,"nw",model))
    if (space=="euclidean") {
      v <- c(.1,-.2,.15,.05,-.1)
      Y <- rbind(m+v,m-v)
    } else if (space=="sphere") {
      v <- .12*model$v+.08*model$v2
      Z <- rbind(sphere_exp(v,model$reference),sphere_exp(-v,model$reference))
      Y <- do.call(rbind,lapply(1:2,function(i) {
        metric_sub(model$reference,Z[i,],m,space="sphere")
      }))
    } else {
      v <- .04*model$b1-.03*model$b2
      Y <- rbind(m+v,m-v)
      add("Wasserstein endpoints",max(abs(Y[,c(1,length(p_grid))]-
                                            matrix(c(.15,.85),2,2,byrow=TRUE))),1e-12)
      add("Wasserstein nondecreasing quantiles",max(0,-min(t(apply(Y,1,diff)))),0)
    }
    fitted <- matrix(m,2,length(m),byrow=TRUE)
    E <- do.call(rbind,make_metric_residuals(Y,fitted,space=space,
                                             muY=model$reference,wasserstein_grid=p_grid)$residuals)
    mean_Y <- object_frechet_mean(Y,space,wasserstein_grid=p_grid)
    mean_E <- object_frechet_mean(E,space,wasserstein_grid=p_grid)
    add(paste(space,"antithetic conditional mean"),
        object_dist(mean_Y,m,space,wasserstein_grid=p_grid),1e-6)
    add(paste(space,"antithetic residual center"),
        object_dist(mean_E,model$reference,space,wasserstein_grid=p_grid),1e-6)
  }
  q0 <- .15+.7*p_grid
  bad <- metric_sub(q0+.04,q0+.04,q0,space="wasserstein",wasserstein_grid=p_grid)
  bad_error <- max(abs(bad-q0))
  if (bad_error<.039) stop("The expected shifted-support counterexample was not detected.")
  checks <- do.call(rbind,report)
  if (!all(checks$passed)) stop("A deterministic generator/centering check failed.")
  list(checks=checks,shifted_support_counterexample_max_error=bad_error,
       sphere_response_ball_radius=sphere_radius*(1+sqrt(3/2)),
       wasserstein_derivative_lower_bound=.7-sqrt(.05^2+.04^2/4)-.05-.04)
}


basic_file <- "basic_function.R"
sample_sizes <- c(200L,500L,1000L)
signal_strengths <- seq(0,1,by=.1)  
continuous_shift <- .5  # Selected continuous X_i -> (X_i + .5) %% 1.
sphere_radius <- .35   
rep_num <- 1000L
alpha <- .05
seed <- 123L
available_cores <- parallel::detectCores(logical=FALSE)
if (is.na(available_cores)) available_cores <- 2L
n_workers <- max(1L,available_cores-1L)

simulation <- run_fitted_residual_simulation(
  basic_file=basic_file,
  sample_sizes=sample_sizes,
  corruption_rates=signal_strengths,
  repetitions=rep_num,
  spaces=c("euclidean","sphere","wasserstein"),
  designs=c("nw","categorical"),
  training_multiplier=1,
  continuous_shift=continuous_shift,
  sphere_radius=sphere_radius,
  alpha=alpha,
  seed=seed,
  n_workers=n_workers,
  progress=TRUE)

result <- simulation$replications
result_summary <- simulation$summary
print(result_summary,row.names=FALSE,digits=4)
print(simulation$plot)


# ================================================================
# ================================================================
# Part 4 Prediction interval
# ================================================================
# ================================================================
.pisi_subset <- function(data, index) {
  list(X=data$X[index], Y=data$Y[index, , drop=FALSE])
}

.pisi_in_prediction_range <- function(X, settings) {
  X>=settings$prediction_x_range[1L] & X<=settings$prediction_x_range[2L]
}

.pisi_future <- function(design, model, settings) {
  if (design=="categorical") {
    return(.frsim_response(settings$n_future, design, model))
  }
  X <- runif(settings$n_future, min=settings$prediction_x_range[1L],
             max=settings$prediction_x_range[2L])
  .frsim_response(settings$n_future, design, model, X=X)
}

.pisi_regression <- function(X, Y, Xout, space, design, settings, cv_plan=NULL) {
  bandwidth_grid <- settings$fit_bw_constants*length(X)^(-settings$fit_bw_exponent)
  bandwidth <- if (settings$bw_mode=="CV") "CV" else {
    bandwidth_grid[ceiling(length(bandwidth_grid)/2)]
  }
  objpi_regression_fit(
    xin=X, Y=Y, xout=Xout, space=space, design=design,
    wasserstein_grid=settings$p_grid, bandwidth=bandwidth,
    bandwidth_grid=bandwidth_grid, regression=settings$regression,
    cv_plan=cv_plan, cv_folds=settings$cv_folds, cv_trim=settings$cv_trim,
    cv_small_sample_loo=settings$cv_small_sample_loo, controls=settings$controls,
    predictor_domain=if (design=="nw") settings$predictor_domain else NULL)
}

.pisi_cv_plan <- function(X, design, settings) {
  if (design=="categorical" || settings$bw_mode=="fixed") return(NULL)
  objpi_kfold_plan(X, settings$cv_folds, settings$cv_trim,
                   settings$cv_small_sample_loo)
}

.pisi_candidates <- function(x, design, model, count, settings) {
  true_mean <- drop(.frsim_mean(x, design, model))
  if (model$space=="euclidean") {
    direction <- matrix(rnorm(count*5L), count, 5L)
    direction <- direction/sqrt(rowSums(direction^2))
    displacement <- direction*(settings$euclidean_search_radius*runif(count)^(1/5))
    return(sweep(displacement, 2L, true_mean, "+"))
  }
  if (model$space=="sphere") {
    z <- runif(count, -1, 1)
    theta <- runif(count, 0, 2*pi)
    return(cbind(sqrt(1-z^2)*cos(theta), sqrt(1-z^2)*sin(theta), z))
  }
  a <- runif(count, -settings$wasserstein_search_bounds[1L],
             settings$wasserstein_search_bounds[1L])
  b <- runif(count, -settings$wasserstein_search_bounds[2L],
             settings$wasserstein_search_bounds[2L])
  sweep(outer(a, model$b1)+outer(b, model$b2), 2L, true_mean, "+")
}

.pisi_one <- function(task, settings, stream) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", stream, envir=.GlobalEnv)
  space <- task$space
  design <- task$design
  n <- task$n
  n_train <- floor(settings$training_fraction*n)
  n_cal_available <- n-n_train
  model <- .frsim_model(space, settings$p_grid, settings$sphere_radius)
  sample <- .frsim_response(n, design, model)
  training_index <- sample.int(n, n_train)
  calibration_index <- setdiff(seq_len(n), training_index)
  training <- .pisi_subset(sample, training_index)
  calibration <- .pisi_subset(sample, calibration_index)
  if (design=="nw") {
    retained <- which(.pisi_in_prediction_range(calibration$X, settings))
    calibration <- .pisi_subset(calibration, retained)
  }
  n_cal <- length(calibration$X)
  future <- .pisi_future(design, model, settings)
  x_size <- if (design=="categorical") c(-1, 1) else settings$size_predictors
  # Construct BOTH plans before either regression engine can consume RNG.
  training_plan <- .pisi_cv_plan(training$X, design, settings)
  proposed_plan <- .pisi_cv_plan(sample$X, design, settings)
  # Isolate search-candidate draws from any stochastic regression optimization.
  candidate_stream <- parallel::nextRNGSubStream(stream)
  fitting_stream <- parallel::nextRNGSubStream(candidate_stream)
  assign(".Random.seed", fitting_stream, envir=.GlobalEnv)
  fit <- .pisi_regression(training$X, training$Y,
                          c(calibration$X, future$X, x_size), space, design,
                          settings, cv_plan=training_plan)
  i_cal <- seq_len(n_cal)
  i_future <- n_cal+seq_len(settings$n_future)
  i_size <- n_cal+settings$n_future+seq_along(x_size)
  fitted_cal <- fit$fitted[i_cal, , drop=FALSE]
  fitted_future <- fit$fitted[i_future, , drop=FALSE]
  lm <- if (n_cal>0L) {
    objpi_lm_calibrate(calibration$Y, fitted_cal, space,
                       settings$p_grid, settings$alpha)
  } else {
    list(space=space, wasserstein_grid=settings$p_grid,
         calibration_scores=numeric(0), threshold=Inf)
  }
  profile_bandwidth <- if (is.null(settings$profile_bandwidth)) {
    settings$profile_bw_constant*(n_train/log(n_train))^(-1/5)
  } else settings$profile_bandwidth
  zm <- objpi_zm_fit(training$X, training$Y, space, design,
                     settings$p_grid, profile_bandwidth, settings$profile_points)
  if (n_cal>0L) {
    zm <- objpi_zm_calibrate(zm, calibration$X, calibration$Y, settings$alpha)
  } else {
    zm$calibration_scores <- numeric(0)
    zm$threshold <- Inf
  }
  proposed_fit <- .pisi_regression(sample$X, sample$Y,
                                   c(sample$X, future$X, x_size), space, design,
                                   settings, cv_plan=proposed_plan)
  i_proposed_sample <- seq_len(n)
  i_proposed_future <- n+seq_len(settings$n_future)
  i_proposed_size <- n+settings$n_future+seq_along(x_size)
  fitted_proposed_sample <- proposed_fit$fitted[i_proposed_sample, , drop=FALSE]
  fitted_proposed_future <- proposed_fit$fitted[i_proposed_future, , drop=FALSE]
  center <- object_frechet_mean(sample$Y, space,
                                wasserstein_grid=settings$p_grid)
  proposed <- objpi_proposed_calibrate(sample$Y, fitted_proposed_sample, center, space,
                                       settings$p_grid, settings$alpha)
  hit <- cbind(
    Proposed=objpi_proposed_scores(proposed, future$Y, fitted_proposed_future)<=proposed$threshold,
    LM=objpi_lm_scores(lm, future$Y, fitted_future)<=lm$threshold,
    ZM=objpi_zm_scores(zm, future$X, future$Y)<=zm$threshold)
  sizes <- matrix(0, length(x_size), ncol(hit), dimnames=list(NULL, colnames(hit)))
  assign(".Random.seed", candidate_stream, envir=.GlobalEnv)
  for (i in seq_along(x_size)) {
    candidates <- .pisi_candidates(x_size[i], design, model,
                                   settings$n_candidates, settings)
    fitted_candidates <- matrix(fit$fitted[i_size[i], ], nrow(candidates),
                                ncol(candidates), byrow=TRUE)
    fitted_proposed_candidates <- matrix(proposed_fit$fitted[i_proposed_size[i], ],
                                         nrow(candidates), ncol(candidates), byrow=TRUE)
    sizes[i, "Proposed"] <- mean(objpi_proposed_scores(
      proposed, candidates, fitted_proposed_candidates)<=proposed$threshold)
    sizes[i, "LM"] <- mean(objpi_lm_scores(lm, candidates, fitted_candidates)<=lm$threshold)
    sizes[i, "ZM"] <- mean(objpi_zm_scores(
      zm, rep(x_size[i], nrow(candidates)), candidates)<=zm$threshold)
  }
  bins <- if (design=="categorical") ifelse(future$X<0, 1L, 2L) else {
    findInterval(future$X, settings$prediction_bins,
                 rightmost.closed=TRUE, all.inside=TRUE)
  }
  conditional <- do.call(rbind, lapply(sort(unique(bins)), function(bin) {
    data.frame(space=space, design=design, n=n, replication=task$replication,
               bin=bin, method=colnames(hit),
               bin_lower=if (design=="categorical") c(-1,1)[bin] else {
                 settings$prediction_bins[bin]
               },
               bin_upper=if (design=="categorical") c(-1,1)[bin] else {
                 settings$prediction_bins[bin+1L]
               },
               coverage=colMeans(hit[bins==bin, , drop=FALSE]),
               n_future=sum(bins==bin), row.names=NULL)
  }))
  replication <- data.frame(
    space=space, design=design, regression=settings$regression, n=n,
    n_training=c(n, n_train, n_train),
    n_calibration=c(0L, n_cal, n_cal), n_scoring=c(n, n_cal, n_cal),
    n_calibration_available=c(0L, n_cal_available, n_cal_available),
    replication=task$replication, method=colnames(hit),
    coverage=colMeans(hit), relative_size=colMeans(sizes),
    threshold=c(proposed$threshold, lm$threshold, zm$threshold),
    zm_saturated=c(NA_real_, NA_real_, as.numeric(zm$threshold>=1)),
    fit_bandwidth=c(proposed_fit$bandwidth, fit$bandwidth, NA_real_),
    fit_cv_edge=c(proposed_fit$cv_edge, fit$cv_edge, NA),
    fit_engine=c(proposed_fit$engine, fit$engine, NA_character_),
    fit_bandwidth_selection=c(proposed_fit$bandwidth_selection,
                              fit$bandwidth_selection, NA_character_),
    fit_cv_folds=c(if (is.null(proposed_plan)) NA_integer_ else proposed_plan$n_folds,
                   if (is.null(training_plan)) NA_integer_ else training_plan$n_folds,
                   NA_integer_),
    profile_bandwidth=c(NA_real_, NA_real_,
                        if (design=="categorical") NA_real_ else profile_bandwidth),
    row.names=NULL)
  list(replication=replication, conditional=conditional,
       cv_plans=list(proposed=proposed_plan, LM=training_plan))
}

.pisi_summary <- function(result, keys=c("space", "design", "n", "method"),
                          variables=c("coverage", "relative_size", "zm_saturated")) {
  groups <- split(result, interaction(result[keys], drop=TRUE, lex.order=TRUE))
  summary <- do.call(rbind, lapply(groups, function(group) {
    row <- group[1L, keys, drop=FALSE]
    row$repetitions <- nrow(group)
    for (variable in variables) {
      values <- group[[variable]]
      values <- values[is.finite(values)]
      row[[variable]] <- if (length(values)) mean(values) else NA_real_
    }
    row
  }))
  rownames(summary) <- NULL
  summary
}

prediction_interval_latex <- function(summary, settings,
                                      file="prediction_interval_table.tex", print_table=TRUE) {
  if (!identical(settings$proposed_sample_usage, "full sample")) {
    stop("This table requires results from the updated full-sample proposed-method simulation.")
  }
  format_entry <- function(row, variable) {
    base::sprintf("%.2f", 100*row[[variable]])
  }
  fitting_caption <- if (is.null(settings$regression)) "" else {
    base::paste0(" with ", if (settings$regression=="nw") {
      "Nadaraya--Watson"
    } else "local Fr\\'echet", " regression")
  }
  table_label <- if (is.null(settings$regression)) {
    "\\label{tab:prediction-interval-simulation}"
  } else base::paste0("\\label{tab:prediction-interval-", settings$regression, "}")
  evaluation_caption <- if ("nw" %in% settings$designs &&
                            !is.null(settings$prediction_x_range)) {
    base::sprintf(" Continuous coverage and RSV are evaluated for $X\\in[%.2f,%.2f]$.",
                  settings$prediction_x_range[1L], settings$prediction_x_range[2L])
  } else ""
  table_rows <- base::c(
    "\\begin{table}[tbp]", "\\centering", "\\small", "\\setlength{\\tabcolsep}{4pt}",
    base::sprintf(base::paste0("\\caption{Prediction-set simulation%s at the nominal %.0f\\%% level, ",
                               "based on %d replications.%s}"),
                  fitting_caption, 100*(1-settings$alpha), settings$repetitions,
                  evaluation_caption),
    table_label,
    "\\begin{tabular}{lllrrrrrr}", "\\toprule",
    " & & & \\multicolumn{2}{c}{Proposed, fitted} & \\multicolumn{2}{c}{Lugosi--Matabuena} & \\multicolumn{2}{c}{Zhou--M\\\"uller} \\\\",
    "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
    "Space & $X$ & $n$ & EC & RSV & EC & RSV & EC & RSV \\\\", "\\midrule")
  labels <- base::c(euclidean="$\\mathbb{R}^{5}$", sphere="$\\mathbb{S}^{2}$",
                    wasserstein="$\\mathcal{W}_{2}$")
  space_rows <- length(settings$designs)*length(settings$sample_sizes)
  for (space_index in seq_along(settings$spaces)) {
    space <- settings$spaces[space_index]
    if (space_index>1L) table_rows <- base::c(table_rows, "\\midrule")
    for (design_index in seq_along(settings$designs)) {
      design <- settings$designs[design_index]
      if (design_index>1L) table_rows <- base::c(table_rows, "\\cmidrule(lr){2-9}")
      for (n_index in seq_along(settings$sample_sizes)) {
        n <- settings$sample_sizes[n_index]
        block <- summary[summary$space==space & summary$design==design & summary$n==n, , drop=FALSE]
        values <- base::unlist(lapply(base::c("Proposed", "LM", "ZM"), function(method) {
          row <- block[block$method==method, , drop=FALSE]
          if (nrow(row)!=1L) stop("The table requires one complete summary per method/setting.")
          base::c(format_entry(row, "coverage"), format_entry(row, "relative_size"))
        }), use.names=FALSE)
        space_cell <- if (design_index==1L && n_index==1L) {
          base::sprintf("\\multirow{%d}{*}{%s}", space_rows, labels[[space]])
        } else ""
        design_cell <- if (n_index==1L) {
          base::sprintf("\\multirow{%d}{*}{%s}", length(settings$sample_sizes),
                        if (design=="nw") "Continuous" else "Categorical")
        } else ""
        table_rows <- base::c(table_rows,
                              base::paste0(base::paste(base::c(space_cell, design_cell, n, values),
                                                       collapse=" & "), " \\\\"))
      }
    }
  }
  table_rows <- base::c(table_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  latex <- base::paste(base::as.character(table_rows), collapse="\n")
  if (!is.null(file)) {
    dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
    base::writeLines(latex, con=file, useBytes=TRUE)
  }
  if (print_table) base::cat(latex, "\n", sep="")
  invisible(latex)
}

run_prediction_interval_simulation <- function(
    basic_file="basic_function.R", sample_sizes=c(200L, 500L, 1000L),
    repetitions=1000L, spaces=c("euclidean", "sphere", "wasserstein"),
    designs=c("nw", "categorical"), training_fraction=.5, alpha=.05,
    n_future=1000L, n_candidates=1024L, seed=123L, n_workers=1L,
    p_grid=seq(0, 1, length.out=501L), sphere_radius=.35,
    fit_bw_constants=exp(seq(log(.1), log(1), length.out=15L)),
    fit_bw_exponent=1/5, bw_mode=c("CV", "fixed"),
    controls=list(tol=1e-6, maxit=500L, nStarts=1L),
    profile_bandwidth=NULL, profile_bw_constant=.5, profile_points=201L,
    size_predictors=NULL, euclidean_search_radius=6,
    wasserstein_search_bounds=c(.15, .12), progress=TRUE,
    latex_file="prediction_interval_table.tex", print_latex=TRUE,
    regression=c("local", "nw"), cv_folds=10L, cv_trim=.2,
    cv_small_sample_loo=TRUE, prediction_x_range=c(0.02,0.98)) {
  regression <- match.arg(regression)
  if (missing(latex_file)) {
    latex_file <- paste0("prediction_interval_", regression, "_table.tex")
  }
  packages <- c("foreach", "doSNOW", "snow")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse=", "))
  basic_path <- normalizePath(basic_file, mustWork=TRUE)
  source(basic_path, local=.GlobalEnv)
  required <- c(".frsim_model", ".frsim_g", ".frsim_mean", ".frsim_response",
                "objpi_regression_fit", "objpi_kfold_plan", "objpi_proposed_calibrate",
                "objpi_lm_calibrate", "objpi_zm_fit", "objpi_zm_calibrate")
  absent <- required[!vapply(required, exists, logical(1), mode="function")]
  if (length(absent)) stop("Load Part 2.2 and the new basic functions first. Missing: ",
                           paste(absent, collapse=", "))
  if (regression=="nw") {
    required_nw <- c("NWFrechetWeights", "NWSpheReg", "NWWassReg", ".nw_controls")
    absent <- required_nw[!vapply(required_nw, exists, logical(1), mode="function")]
    if (length(absent)) stop("Original NW functions are missing: ",
                             paste(absent, collapse=", "))
  }
  integer_ok <- function(x, minimum=1L) {
    is.numeric(x) && length(x)>0L && all(is.finite(x)) &&
      all(x>=minimum & x==floor(x))
  }
  if (!integer_ok(sample_sizes, 40L) || !integer_ok(repetitions, 2L) ||
      length(repetitions)!=1L || !integer_ok(n_workers) || length(n_workers)!=1L ||
      !integer_ok(n_future, 10L) || length(n_future)!=1L ||
      !integer_ok(n_candidates, 2L) || length(n_candidates)!=1L ||
      !integer_ok(profile_points, 2L) || length(profile_points)!=1L ||
      !integer_ok(cv_folds, 2L) || length(cv_folds)!=1L) {
    stop("Invalid sample sizes, replication counts, or numerical grid sizes.")
  }
  scalar_positive <- function(x) length(x)==1L && is.finite(x) && x>0
  if (!is.numeric(prediction_x_range) || length(prediction_x_range)!=2L ||
      any(!is.finite(prediction_x_range)) || prediction_x_range[1L]<=0 ||
      prediction_x_range[2L]>=1 || diff(prediction_x_range)<=0) {
    stop("Use an increasing prediction_x_range strictly inside (0,1).")
  }
  if (is.null(size_predictors)) {
    size_predictors <- seq(prediction_x_range[1L], prediction_x_range[2L],
                           length.out=5L)
  }
  if (!is.numeric(cv_trim) || length(cv_trim)!=1L || !is.finite(cv_trim) ||
      cv_trim<0 || cv_trim>=.5 ||
      !is.logical(cv_small_sample_loo) || length(cv_small_sample_loo)!=1L ||
      is.na(cv_small_sample_loo)) stop("Invalid K-fold validation settings.")
  if (!scalar_positive(training_fraction) || training_fraction>=1 ||
      any(floor(training_fraction*sample_sizes)<10) ||
      any(sample_sizes-floor(training_fraction*sample_sizes)<10)) {
    stop("Each split must contain at least ten observations.")
  }
  if (!scalar_positive(alpha) || alpha>=1 || !integer_ok(seed, 0L) ||
      length(seed)!=1L || seed>.Machine$integer.max ||
      !scalar_positive(profile_bw_constant) ||
      (!is.null(profile_bandwidth) && !scalar_positive(profile_bandwidth)) ||
      !scalar_positive(euclidean_search_radius) || !scalar_positive(sphere_radius) ||
      !scalar_positive(fit_bw_exponent)) stop("Invalid scalar setting.")
  if (!is.numeric(fit_bw_constants) || !length(fit_bw_constants) ||
      any(!is.finite(fit_bw_constants)) || any(fit_bw_constants<=0)) {
    stop("Use positive bandwidth-grid constants.")
  }
  if (!is.numeric(size_predictors) || !length(size_predictors) ||
      any(!is.finite(size_predictors)) ||
      any(size_predictors<prediction_x_range[1L] |
          size_predictors>prediction_x_range[2L])) {
    stop("Continuous size-evaluation predictors must lie inside prediction_x_range.")
  }
  if (length(wasserstein_search_bounds)!=2L ||
      any(!is.finite(wasserstein_search_bounds)) || any(wasserstein_search_bounds<=0)) {
    stop("Use two positive Wasserstein search bounds.")
  }
  if (!is.logical(progress) || length(progress)!=1L || is.na(progress) ||
      !is.logical(print_latex) || length(print_latex)!=1L || is.na(print_latex)) {
    stop("progress and print_latex must be TRUE or FALSE.")
  }
  if (!is.null(latex_file) &&
      (!is.character(latex_file) || length(latex_file)!=1L || is.na(latex_file) ||
       !nzchar(latex_file))) stop("Use one LaTeX filename or NULL.")
  .objpi_grid_weights(p_grid, length(p_grid))
  if (!is.null(controls) && !is.list(controls)) stop("controls must be a list or NULL.")
  if (regression=="nw") .nw_controls(if (is.null(controls)) list() else controls)
  bw_mode <- match.arg(bw_mode)
  spaces <- match.arg(spaces, c("euclidean", "sphere", "wasserstein"), several.ok=TRUE)
  designs <- match.arg(designs, c("nw", "categorical"), several.ok=TRUE)
  spaces <- unique(spaces)
  designs <- unique(designs)
  sample_sizes <- sort(unique(sample_sizes))
  needs_frechet <- regression=="local" && "nw" %in% designs &&
    any(c("sphere", "wasserstein") %in% spaces)
  if (needs_frechet && !requireNamespace("frechet", quietly=TRUE)) {
    stop("Install required package: frechet (install.packages('frechet')).")
  }
  if ("wasserstein" %in% spaces) {
    model <- .frsim_model("wasserstein", p_grid, sphere_radius)
    if (.7-model$tau*sqrt(3)-sum(wasserstein_search_bounds)<=0) {
      stop("Wasserstein search quantiles must remain strictly increasing.")
    }
  }
  settings <- list(sample_sizes=sample_sizes, repetitions=repetitions,
                   spaces=spaces, designs=designs, training_fraction=training_fraction,
                   sample_size_definition="n = total observed sample for each method",
                   proposed_sample_usage="full sample",
                   conformal_sample_usage="training/calibration split of the same observed sample",
                   prediction_x_range=prediction_x_range,
                   prediction_bins=seq(prediction_x_range[1L], prediction_x_range[2L], length.out=6L),
                   continuous_future_predictor="Uniform on prediction_x_range",
                   continuous_coverage_target="coverage averaged over X in prediction_x_range",
                   continuous_calibration="retain assigned split-calibration points inside prediction_x_range",
                   observed_predictor="Part 2.2: Unif(0,1) for continuous predictors",
                   regression=regression,
                   continuous_mean_regression=if (regression=="nw") {
                     "Nadaraya-Watson Frechet regression for Proposed and LM"
                   } else "local-linear Frechet regression for Proposed and LM",
                   continuous_fit_engines=if (regression=="nw") {
                     c(euclidean="NWFrechetWeights", sphere="NWSpheReg", wasserstein="NWWassReg")
                   } else c(euclidean="multivariate local-linear least squares",
                            sphere="frechet::LocSpheReg", wasserstein="frechet::LocDenReg"),
                   continuous_fit_cv="shared K-fold CV using fitting sample only",
                   cv_folds=cv_folds, cv_trim=cv_trim, cv_small_sample_loo=cv_small_sample_loo,
                   cv_loss="mean squared object distance on the held-out interior observations",
                   cv_fold_policy="fixed across bandwidths and regression choices",
                   predictor_domain=c(0, 1),
                   fit_bw_constants_usage="common K-fold grid; same admissible spherical grid for NW/local",
                   controls_usage="NW spherical optimization; local package optimizers use defaults",
                   paired_evaluation="same data, conformal split, CV plans and search candidates for equal seeds/settings",
                   categorical_mean_regression="within-category Frechet mean",
                   frechet_version=if (needs_frechet) as.character(utils::packageVersion("frechet")) else NA_character_,
                   alpha=alpha, n_future=n_future, n_candidates=n_candidates,
                   seed=seed, p_grid=p_grid, sphere_radius=sphere_radius,
                   fit_bw_constants=sort(unique(fit_bw_constants)), fit_bw_exponent=fit_bw_exponent,
                   bw_mode=bw_mode, controls=controls, profile_bandwidth=profile_bandwidth,
                   profile_bw_constant=profile_bw_constant, profile_points=profile_points,
                   size_predictors=size_predictors, euclidean_search_radius=euclidean_search_radius,
                   wasserstein_search_bounds=wasserstein_search_bounds,
                   generator="Part 2.2, no predictor corruption",
                   zm_reference="https://github.com/hg-zh/Object-Conformal/tree/main/Simulation/Proposed",
                   zm_commit="1893f14361bf4645720f7220a776467cc893c3e9")
  tasks <- expand.grid(space=spaces, design=designs, n=sample_sizes,
                       replication=seq_len(repetitions), stringsAsFactors=FALSE)
  n_tasks <- nrow(tasks)
  n_workers <- min(as.integer(n_workers), n_tasks)
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir=.GlobalEnv)
  cluster <- NULL
  bar <- NULL
  on.exit({
    if (!is.null(bar)) close(bar)
    if (!is.null(cluster)) snow::stopCluster(cluster)
    foreach::registerDoSEQ()
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir=.GlobalEnv) else {
      if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) {
        rm(".Random.seed", envir=.GlobalEnv)
      }
    }
  }, add=TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  streams <- vector("list", n_tasks)
  streams[[1L]] <- .Random.seed
  if (n_tasks>1L) for (i in 2:n_tasks) {
    streams[[i]] <- parallel::nextRNGStream(streams[[i-1L]])
  }
  basic_env <- environment(objpi_proposed_calibrate)
  basic_names <- ls(basic_env, pattern="^\\.?objpi_", all.names=TRUE)
  simulation_env <- environment(.pisi_one)
  simulation_names <- ls(simulation_env, pattern="^\\.pisi_", all.names=TRUE)
  generator_names <- c(".frsim_model", ".frsim_g", ".frsim_mean", ".frsim_response")
  generator_env <- environment(.frsim_response)
  helpers <- c(mget(basic_names, envir=basic_env),
               mget(simulation_names, envir=simulation_env),
               setNames(lapply(generator_names, function(name) {
                 get(name, envir=generator_env, mode="function", inherits=TRUE)
               }), generator_names))
  cluster <- snow::makeCluster(n_workers, type="SOCK")
  doSNOW::registerDoSNOW(cluster)
  invisible(snow::clusterCall(cluster, function(path) {
    source(path, local=.GlobalEnv)
    NULL
  }, basic_path))
  invisible(snow::clusterCall(cluster, function(functions) {
    list2env(functions, envir=.GlobalEnv)
    NULL
  }, helpers))
  snow_options <- list(preschedule=FALSE)
  if (progress) {
    bar <- utils::txtProgressBar(min=0, max=n_tasks, style=3)
    snow_options$progress <- function(completed) utils::setTxtProgressBar(bar, completed)
  }
  `%dopar%` <- foreach::`%dopar%`
  output <- foreach::foreach(task_id=seq_len(n_tasks), .inorder=TRUE,
                             .errorhandling="stop", .options.snow=snow_options,
                             .packages=if (needs_frechet) "frechet" else character(),
                             .noexport=names(helpers)) %dopar% {
                               .pisi_one(tasks[task_id, , drop=FALSE], settings, streams[[task_id]])
                             }
  if (!is.null(bar)) {
    close(bar)
    bar <- NULL
  }
  snow::stopCluster(cluster)
  cluster <- NULL
  foreach::registerDoSEQ()
  replications <- do.call(rbind, lapply(output, `[[`, "replication"))
  conditional <- do.call(rbind, lapply(output, `[[`, "conditional"))
  cv_plans <- lapply(output, function(value) value$cv_plans)
  replications <- replications[order(match(replications$space, spaces),
                                     match(replications$design, designs), replications$n,
                                     replications$replication, match(replications$method, c("Proposed", "LM", "ZM"))), ]
  rownames(replications) <- NULL
  summary <- .pisi_summary(replications)
  conditional_summary <- .pisi_summary(conditional,
                                       keys=c("space", "design", "n", "method", "bin", "bin_lower", "bin_upper"),
                                       variables="coverage")
  if (any(summary$repetitions!=repetitions)) stop("Unexpected replication counts.")
  settings$n_workers <- n_workers
  settings$parallel_backend <- "foreach / doSNOW"
  settings$basic_file_md5 <- unname(tools::md5sum(basic_path))
  latex <- prediction_interval_latex(summary, settings, latex_file, print_latex)
  invisible(list(latex=latex, summary=summary, replications=replications,
                 conditional_summary=conditional_summary,
                 conditional_replications=conditional, cv_plans=cv_plans,
                 cv_tasks=tasks, settings=settings,
                 session_info=utils::sessionInfo()))
}


run_prediction_interval_nw <- function(
    ..., latex_file="prediction_interval_nw_table.tex") {
  run_prediction_interval_simulation(..., regression="nw", latex_file=latex_file)
}

run_prediction_interval_local <- function(
    ..., latex_file="prediction_interval_local_table.tex") {
  run_prediction_interval_simulation(..., regression="local", latex_file=latex_file)
}


run_prediction_interval_nw(basic_file="basic_function.R", n_workers=n_workers)

