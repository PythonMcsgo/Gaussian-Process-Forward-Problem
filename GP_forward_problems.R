library(mvtnorm)
options(width = 100, digits = 6)

script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

output_dir <- script_directory()
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Color set
paper_palette <- function(n = 160L, palette = c("viridis", "inferno")) {
  palette <- match.arg(palette)
  anchors <- if (palette == "viridis") {
    c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  } else {
    c("#000004", "#420A68", "#932667", "#DD513A", "#FCA50A", "#FCFFA4")
  }
  grDevices::colorRampPalette(anchors)(n)
}

# ----------------------------------------------
# Kernel derivatives and stable conditioning
# -------------------------------------------------------

pair_diff <- function(x1, x2) outer(x1, x2, "-")

se_kernel_1d <- function(x1, x2, sigma_f = 1, ell = 1) {
  stopifnot(is.numeric(x1), is.numeric(x2), sigma_f > 0, ell > 0)
  d <- pair_diff(x1, x2)
  sigma_f^2 * exp(-0.5 * (d / ell)^2)
}

se_derivatives_1d <- function(x1, x2, sigma_f = 1, ell = 1) {
  d <- pair_diff(x1, x2)
  k <- se_kernel_1d(x1, x2, sigma_f, ell)
  list(
    k = k,
    dx = -d * k / ell^2,
    dy = d * k / ell^2,
    dxx = (d^2 - ell^2) * k / ell^4,
    dyy = (d^2 - ell^2) * k / ell^4,
    dxy = (ell^2 - d^2) * k / ell^4,
    dxxy = d * (d^2 - 3 * ell^2) * k / ell^6,
    dxyy = d * (3 * ell^2 - d^2) * k / ell^6,
    dxxyy = ((d / ell)^4 - 6 * (d / ell)^2 + 3) * k / ell^4
  )
}

stable_chol <- function(K, relative_jitter = 1e-10, max_tries = 8L) {
  stopifnot(is.matrix(K), nrow(K) == ncol(K), relative_jitter > 0)
  K <- 0.5 * (K + t(K))
  scale <- max(mean(diag(K)), 1)
  for (attempt in 0:(max_tries - 1L)) {
    jitter <- relative_jitter * scale * 10^attempt
    R <- tryCatch(chol(K + jitter * diag(nrow(K))), error = identity)
    if (!inherits(R, "error")) return(list(R = R, jitter = jitter))
  }
  stop("Cholesky factorisation failed after adaptive jitter escalation.")
}

condition_gp <- function(K_obs, y, K_star_obs, prior_variance,
                         observation_variance = rep(0, length(y)),
                         relative_jitter = 1e-10) {
  stopifnot(length(y) == nrow(K_obs), ncol(K_star_obs) == length(y),
            length(prior_variance) == nrow(K_star_obs),
            length(observation_variance) == length(y),
            all(observation_variance >= 0))
  K_likelihood <- K_obs + diag(observation_variance, nrow(K_obs))
  fac <- stable_chol(K_likelihood, relative_jitter)
  alpha <- backsolve(fac$R, forwardsolve(t(fac$R), y))
  projected <- forwardsolve(t(fac$R), t(K_star_obs))
  variance <- pmax(prior_variance - colSums(projected^2), 0)
  mean <- drop(K_star_obs %*% alpha)
  sd <- sqrt(variance)
  list(
    mean = mean,
    sd = sd,
    lower95 = mean - 1.96 * sd,
    upper95 = mean + 1.96 * sd,
    jitter = fac$jitter
  )
}

# -----------------------------------------------------
# One-dimensional forward solvers
# -----------------------------------------------------

solve_poisson_gp <- function(x_bc, g_bc, x_col, f_col, x_star,
                              sigma_f = 1, ell = 0.45,
                              noise_sd_forcing = 0,
                              relative_jitter = 1e-10) {
  stopifnot(length(x_bc) == length(g_bc), length(x_col) == length(f_col))
  d_bc <- se_derivatives_1d(x_bc, x_bc, sigma_f, ell)
  d_bf <- se_derivatives_1d(x_bc, x_col, sigma_f, ell)
  d_ff <- se_derivatives_1d(x_col, x_col, sigma_f, ell)
  K_obs <- rbind(
    cbind(d_bc$k, -d_bf$dyy),
    cbind(-t(d_bf$dyy), d_ff$dxxyy)
  )
  d_sb <- se_derivatives_1d(x_star, x_bc, sigma_f, ell)
  d_sf <- se_derivatives_1d(x_star, x_col, sigma_f, ell)
  fit <- condition_gp(
    K_obs, c(g_bc, f_col), cbind(d_sb$k, -d_sf$dyy),
    rep(sigma_f^2, length(x_star)),
    c(rep(0, length(x_bc)), rep(noise_sd_forcing^2, length(x_col))),
    relative_jitter
  )
  fit$n_boundary_initial <- length(x_bc)
  fit$n_collocation <- length(x_col)
  fit
}

solve_first_order_gp <- function(t_ic, g_ic, t_col, f_col, t_star, alpha,
                                  sigma_f = 1.5, ell = 1.3,
                                  noise_sd_forcing = 0,
                                  relative_jitter = 1e-10) {
  stopifnot(length(t_ic) == length(g_ic), length(t_col) == length(f_col))
  d_ii <- se_derivatives_1d(t_ic, t_ic, sigma_f, ell)
  d_if <- se_derivatives_1d(t_ic, t_col, sigma_f, ell)
  d_ff <- se_derivatives_1d(t_col, t_col, sigma_f, ell)
  K_uL <- d_if$dy + alpha * d_if$k
  K_LL <- d_ff$dxy + alpha * (d_ff$dx + d_ff$dy) + alpha^2 * d_ff$k
  K_obs <- rbind(cbind(d_ii$k, K_uL), cbind(t(K_uL), K_LL))
  d_si <- se_derivatives_1d(t_star, t_ic, sigma_f, ell)
  d_sf <- se_derivatives_1d(t_star, t_col, sigma_f, ell)
  fit <- condition_gp(
    K_obs, c(g_ic, f_col), cbind(d_si$k, d_sf$dy + alpha * d_sf$k),
    rep(sigma_f^2, length(t_star)),
    c(rep(0, length(t_ic)), rep(noise_sd_forcing^2, length(t_col))),
    relative_jitter
  )
  fit$n_boundary_initial <- length(t_ic)
  fit$n_collocation <- length(t_col)
  fit
}

cov_euler_blocks <- function(x1, x2, p, q, sigma_f, ell) {
  d <- se_derivatives_1d(x1, x2, sigma_f, ell)
  K_Mu <- x1^2 * d$dxx + (p * x1) * d$dx + q * d$k
  K_uM <- sweep(d$dyy, 2, x2^2, "*") +
    sweep(d$dy, 2, p * x2, "*") + q * d$k
  K_MM <- outer(x1^2, x2^2) * d$dxxyy +
    p * outer(x1^2, x2) * d$dxxy +
    p * outer(x1, x2^2) * d$dxyy +
    p^2 * outer(x1, x2) * d$dxy +
    q * x1^2 * d$dxx +
    q * sweep(d$dyy, 2, x2^2, "*") +
    p * q * x1 * d$dx +
    p * q * sweep(d$dy, 2, x2, "*") +
    q^2 * d$k
  list(K_Mu = K_Mu, K_uM = K_uM, K_MM = K_MM)
}

solve_euler_gp <- function(t_bc, g_bc, t_col, f_col, t_star, p, q,
                            sigma_f = 1.5, ell = 0.9,
                            noise_sd_forcing = 0,
                            relative_jitter = 1e-10) {
  stopifnot(length(t_bc) == length(g_bc), length(t_col) == length(f_col))
  K_uu <- se_kernel_1d(t_bc, t_bc, sigma_f, ell)
  K_uM <- cov_euler_blocks(t_bc, t_col, p, q, sigma_f, ell)$K_uM
  K_MM <- cov_euler_blocks(t_col, t_col, p, q, sigma_f, ell)$K_MM
  K_obs <- rbind(cbind(K_uu, K_uM), cbind(t(K_uM), K_MM))
  K_su <- se_kernel_1d(t_star, t_bc, sigma_f, ell)
  K_sM <- cov_euler_blocks(t_star, t_col, p, q, sigma_f, ell)$K_uM
  fit <- condition_gp(
    K_obs, c(g_bc, f_col), cbind(K_su, K_sM),
    rep(sigma_f^2, length(t_star)),
    c(rep(0, length(t_bc)), rep(noise_sd_forcing^2, length(t_col))),
    relative_jitter
  )
  fit$n_boundary_initial <- length(t_bc)
  fit$n_collocation <- length(t_col)
  fit
}

# ---------------------------------
# Heat-equation forward solver
# ---------------------------------
se_kernel_2d <- function(z1, z2, sigma_f = 1, ell_x = 0.22, ell_t = 0.35) {
  stopifnot(is.matrix(z1), is.matrix(z2), ncol(z1) == 2, ncol(z2) == 2)
  dx <- pair_diff(z1[, 1], z2[, 1])
  dt <- pair_diff(z1[, 2], z2[, 2])
  sigma_f^2 * exp(-0.5 * (dx / ell_x)^2 - 0.5 * (dt / ell_t)^2)
}

cov_uL_heat <- function(z1, z2, kappa, sigma_f, ell_x, ell_t) {
  dx <- pair_diff(z1[, 1], z2[, 1])
  dt <- pair_diff(z1[, 2], z2[, 2])
  k <- se_kernel_2d(z1, z2, sigma_f, ell_x, ell_t)
  k * (dt / ell_t^2 - kappa * (dx^2 / ell_x^4 - 1 / ell_x^2))
}

cov_LL_heat <- function(z1, z2, kappa, sigma_f, ell_x, ell_t) {
  dx <- pair_diff(z1[, 1], z2[, 1])
  dt <- pair_diff(z1[, 2], z2[, 2])
  k <- se_kernel_2d(z1, z2, sigma_f, ell_x, ell_t)
  k * (1 / ell_t^2 - dt^2 / ell_t^4 +
       kappa^2 * (dx^4 / ell_x^8 - 6 * dx^2 / ell_x^6 + 3 / ell_x^4))
}

solve_heat_gp <- function(z_boundary_initial, g_boundary_initial,
                           z_col, f_col, z_star, kappa = 1,
                           sigma_f = 1, ell_x = 0.22, ell_t = 0.35,
                           noise_sd_forcing = 0,
                           relative_jitter = 1e-9) {
  stopifnot(nrow(z_boundary_initial) == length(g_boundary_initial),
            nrow(z_col) == length(f_col))
  K_uu <- se_kernel_2d(
    z_boundary_initial, z_boundary_initial, sigma_f, ell_x, ell_t
  )
  K_uL <- cov_uL_heat(
    z_boundary_initial, z_col, kappa, sigma_f, ell_x, ell_t
  )
  K_LL <- cov_LL_heat(z_col, z_col, kappa, sigma_f, ell_x, ell_t)
  K_obs <- rbind(cbind(K_uu, K_uL), cbind(t(K_uL), K_LL))
  K_star <- cbind(
    se_kernel_2d(z_star, z_boundary_initial, sigma_f, ell_x, ell_t),
    cov_uL_heat(z_star, z_col, kappa, sigma_f, ell_x, ell_t)
  )
  fit <- condition_gp(
    K_obs, c(g_boundary_initial, f_col), K_star,
    rep(sigma_f^2, nrow(z_star)),
    c(rep(0, nrow(z_boundary_initial)),
      rep(noise_sd_forcing^2, nrow(z_col))),
    relative_jitter
  )
  fit$n_boundary_initial <- nrow(z_boundary_initial)
  fit$n_collocation <- nrow(z_col)
  fit
}

# ------------------------------------------------
# Empirical-Bayes hyperparameter estimation
# ------------------------------------------------

gp_negative_log_marginal <- function(K_obs, y,
                                     observation_variance = rep(0, length(y)),
                                     relative_jitter = 1e-10) {
  stopifnot(length(observation_variance) == length(y),
            all(observation_variance >= 0))
  fac <- stable_chol(
    K_obs + diag(observation_variance, nrow(K_obs)), relative_jitter
  )
  alpha <- backsolve(fac$R, forwardsolve(t(fac$R), y))
  0.5 * sum(y * alpha) + sum(log(diag(fac$R))) +
    0.5 * length(y) * log(2 * pi)
}

optimise_empirical_bayes <- function(objective, start, lower, upper,
                                     maxit = 200L) {
  stopifnot(length(start) == length(lower), length(lower) == length(upper),
            all(start > 0), all(lower > 0), all(upper > lower))
  fit <- stats::optim(
    par = log(start), fn = objective, method = "L-BFGS-B",
    lower = log(lower), upper = log(upper),
    control = list(maxit = maxit, factr = 1e7)
  )
  algorithm <- "L-BFGS-B"
  if (fit$convergence != 0L) {
    retry <- stats::nlminb(
      start = fit$par, objective = objective,
      lower = log(lower), upper = log(upper),
      control = list(iter.max = maxit, eval.max = 3L * maxit,
                     rel.tol = 1e-7, x.tol = 1e-7)
    )
    if (is.finite(retry$objective) &&
        (retry$convergence == 0L || retry$objective <= fit$value)) {
      fit <- list(
        par = retry$par, value = retry$objective,
        convergence = retry$convergence, counts = retry$evaluations
      )
      algorithm <- "nlminb retry"
    }
  }
  list(
    par = exp(fit$par),
    nll = fit$value,
    value = fit$value,
    convergence = fit$convergence,
    counts = fit$counts,
    algorithm = algorithm
  )
}

deterministic_multistarts <- function(start, lower, upper) {
  p <- length(start)
  fractions <- if (p == 2L) {
    rbind(
      c(0.15, 0.15), c(0.15, 0.85), c(0.85, 0.15), c(0.85, 0.85),
      c(0.50, 0.25), c(0.25, 0.50), c(0.75, 0.50), c(0.50, 0.75)
    )
  } else if (p == 3L) {
    as.matrix(expand.grid(rep(list(c(0.2, 0.8)), 3)))
  } else {
    stop("Only two- and three-parameter empirical-Bayes problems are supported.")
  }
  log_lower <- log(lower)
  span <- log(upper) - log_lower
  generated <- t(apply(fractions, 1, function(frac) exp(log_lower + frac * span)))
  unique(rbind(start, generated))
}

optimise_multistart_empirical_bayes <- function(objective, start, lower, upper,
                                                maxit = 200L) {
  starts <- deterministic_multistarts(start, lower, upper)
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optimise_empirical_bayes(objective, starts[i, ], lower, upper, maxit)
  })
  values <- vapply(fits, function(x) x$nll, numeric(1))
  convergence <- vapply(fits, function(x) x$convergence, integer(1))
  eligible <- which(is.finite(values) & convergence == 0L)
  if (!length(eligible)) eligible <- which(is.finite(values))
  best <- eligible[which.min(values[eligible])]
  result <- fits[[best]]
  result$n_starts <- nrow(starts)
  result$best_start <- starts[best, ]
  result$all_nll <- values
  result$all_convergence <- convergence
  result
}

poisson_system <- function(x_bc, g_bc, x_col, f_col, sigma_f, ell,
                           noise_sd_forcing = 0) {
  d_bc <- se_derivatives_1d(x_bc, x_bc, sigma_f, ell)
  d_bf <- se_derivatives_1d(x_bc, x_col, sigma_f, ell)
  d_ff <- se_derivatives_1d(x_col, x_col, sigma_f, ell)
  list(
    K = rbind(
      cbind(d_bc$k, -d_bf$dyy),
      cbind(-t(d_bf$dyy), d_ff$dxxyy)
    ),
    y = c(g_bc, f_col),
    observation_variance = c(
      rep(0, length(x_bc)), rep(noise_sd_forcing^2, length(x_col))
    )
  )
}

first_order_system <- function(t_ic, g_ic, t_col, f_col, alpha,
                               sigma_f, ell, noise_sd_forcing = 0) {
  d_ii <- se_derivatives_1d(t_ic, t_ic, sigma_f, ell)
  d_if <- se_derivatives_1d(t_ic, t_col, sigma_f, ell)
  d_ff <- se_derivatives_1d(t_col, t_col, sigma_f, ell)
  K_uL <- d_if$dy + alpha * d_if$k
  K_LL <- d_ff$dxy + alpha * (d_ff$dx + d_ff$dy) + alpha^2 * d_ff$k
  list(
    K = rbind(cbind(d_ii$k, K_uL), cbind(t(K_uL), K_LL)),
    y = c(g_ic, f_col),
    observation_variance = c(
      rep(0, length(t_ic)), rep(noise_sd_forcing^2, length(t_col))
    )
  )
}

euler_system <- function(t_bc, g_bc, t_col, f_col, p, q, sigma_f, ell,
                         noise_sd_forcing = 0) {
  K_uu <- se_kernel_1d(t_bc, t_bc, sigma_f, ell)
  K_uM <- cov_euler_blocks(t_bc, t_col, p, q, sigma_f, ell)$K_uM
  K_MM <- cov_euler_blocks(t_col, t_col, p, q, sigma_f, ell)$K_MM
  list(
    K = rbind(cbind(K_uu, K_uM), cbind(t(K_uM), K_MM)),
    y = c(g_bc, f_col),
    observation_variance = c(
      rep(0, length(t_bc)), rep(noise_sd_forcing^2, length(t_col))
    )
  )
}

heat_system <- function(z_boundary_initial, g_boundary_initial,
                        z_col, f_col, kappa, sigma_f, ell_x, ell_t,
                        noise_sd_forcing = 0) {
  K_uu <- se_kernel_2d(
    z_boundary_initial, z_boundary_initial, sigma_f, ell_x, ell_t
  )
  K_uL <- cov_uL_heat(
    z_boundary_initial, z_col, kappa, sigma_f, ell_x, ell_t
  )
  K_LL <- cov_LL_heat(z_col, z_col, kappa, sigma_f, ell_x, ell_t)
  list(
    K = rbind(cbind(K_uu, K_uL), cbind(t(K_uL), K_LL)),
    y = c(g_boundary_initial, f_col),
    observation_variance = c(
      rep(0, nrow(z_boundary_initial)),
      rep(noise_sd_forcing^2, nrow(z_col))
    )
  )
}

safe_log_objective <- function(system_builder, relative_jitter = 1e-10) {
  force(system_builder)
  function(log_parameters) {
    value <- tryCatch({
      sys <- system_builder(exp(log_parameters))
      gp_negative_log_marginal(
        sys$K, sys$y, sys$observation_variance, relative_jitter
      )
    }, error = function(e) Inf)
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }
}

# -----------------------------------------------
# Numerical designs and forward solutions
# ------------------------------------------------


manual_p <- c(sigma_f = 1.0, ell = 0.45)
lower_p <- c(sigma_f = 0.05, ell = 0.03)
upper_p <- c(sigma_f = 20.0, ell = 3.00)
manual_1 <- c(sigma_f = 1.5, ell = 1.30)
lower_1 <- c(sigma_f = 0.05, ell = 0.10)
upper_1 <- c(sigma_f = 20.0, ell = 8.00)
manual_e <- c(sigma_f = 1.5, ell = 0.90)
lower_e <- c(sigma_f = 0.05, ell = 0.05)
upper_e <- c(sigma_f = 20.0, ell = 5.00)
manual_h <- c(sigma_f = 1.0, ell_x = 0.18, ell_t = 0.50)
start_h <- c(sigma_f = 1.0, ell_x = 0.20, ell_t = 1.00)
lower_h <- c(sigma_f = 0.05, ell_x = 0.03, ell_t = 0.05)
upper_h <- c(sigma_f = 20.0, ell_x = 1.00, ell_t = 3.00)


noise_ratio <- 0.01

set.seed(20161387)

# Poisson equation: two boundary values and four interior midpoint constraints.
x_star_p <- seq(0, 1, length.out = 300)
x_bc_p <- c(0, 1)
x_col_p <- (seq_len(4) - 0.5) / 4
truth_p <- x_star_p - x_star_p^2
g_bc_p <- c(0, 0)
f_p_exact <- rep(2, length(x_col_p))
noise_sd_f_p <- noise_ratio * max(abs(f_p_exact))
f_p_noisy <- f_p_exact + rnorm(length(f_p_exact), 0,noise_sd_f_p)

run_poisson_case <- function(f_data, noise_sd_forcing) {
  builder <- function(par) poisson_system(
    x_bc_p, g_bc_p, x_col_p, f_data, par[1], par[2], noise_sd_forcing
  )
  objective <- safe_log_objective(builder)
  eb <- optimise_empirical_bayes(objective, manual_p, lower_p, upper_p)
  ms <- optimise_multistart_empirical_bayes(objective, manual_p, lower_p, upper_p)
  solve_with <- function(parameters) solve_poisson_gp(
    x_bc_p, g_bc_p, x_col_p, f_data, x_star_p,
    sigma_f = parameters[1], ell = parameters[2],
    noise_sd_forcing = noise_sd_forcing
  )
  list(
    g_boundary_initial = g_bc_p, f_data = f_data,
    noise_sd_forcing = noise_sd_forcing,
    objective = objective, eb = eb, ms = ms,
    fixed_fit = solve_with(manual_p), eb_fit = solve_with(eb$par),
    ms_fit = solve_with(ms$par)
  )
}

result_p_nf <- run_poisson_case(f_p_exact, 0)
result_p_noisy <- run_poisson_case(f_p_noisy, noise_sd_f_p)

# First-order ODE: one initial value and four midpoint constraints.
t_star_1 <- seq(0, 4, length.out = 300)
t_ic_1 <- 0
t_col_1 <- (seq_len(4) - 0.5) * 4 / 4
truth_first <- function(t) (2 * sin(t) - cos(t)) / 5 + 1.2 * exp(-2 * t)
truth_1 <- truth_first(t_star_1)
g_ic_1 <- 1
f_1_exact <- sin(t_col_1)
noise_sd_f_1 <- noise_ratio * max(abs(f_1_exact))
f_1_noisy <- f_1_exact + rnorm(length(f_1_exact),0, noise_sd_f_1)


run_first_order_case <- function(f_data, noise_sd_forcing) {
  builder <- function(par) first_order_system(
    t_ic_1, g_ic_1, t_col_1, f_data, alpha = 2,
    sigma_f = par[1], ell = par[2],
    noise_sd_forcing = noise_sd_forcing
  )
  objective <- safe_log_objective(builder)
  eb <- optimise_empirical_bayes(objective, manual_1, lower_1, upper_1)
  ms <- optimise_multistart_empirical_bayes(objective, manual_1, lower_1, upper_1)
  solve_with <- function(parameters) solve_first_order_gp(
    t_ic_1, g_ic_1, t_col_1, f_data, t_star_1, alpha = 2,
    sigma_f = parameters[1], ell = parameters[2],
    noise_sd_forcing = noise_sd_forcing
  )
  list(
    g_boundary_initial = g_ic_1, f_data = f_data,
    noise_sd_forcing = noise_sd_forcing,
    objective = objective, eb = eb, ms = ms,
    fixed_fit = solve_with(manual_1), eb_fit = solve_with(eb$par),
    ms_fit = solve_with(ms$par)
  )
}

result_1_nf <- run_first_order_case(f_1_exact, 0)
result_1_noisy <- run_first_order_case(f_1_noisy, noise_sd_f_1)

# Euler--Cauchy ODE: two boundary values and four midpoint constraints.
t_star_e <- seq(0.5, 3, length.out = 300)
t_bc_e <- c(1, 2)
t_col_e <- 0.5 + (seq_len(4) - 0.5) * (3 - 0.5) / 4
truth_euler <- function(t) t^2 / 3 + 2 * t / 9 + 4 / (9 * t)
truth_e <- truth_euler(t_star_e)
g_bc_e <- c(1, 2)
f_e_exact <- t_col_e^2
noise_sd_f_e <- noise_ratio * max(abs(f_e_exact))
f_e_noisy <- f_p_exact + rnorm(length(f_e_exact),0,noise_sd_f_e)

run_euler_case <- function(f_data, noise_sd_forcing) {
  builder <- function(par) euler_system(
    t_bc_e, g_bc_e, t_col_e, f_data, p = 1, q = -1,
    sigma_f = par[1], ell = par[2],
    noise_sd_forcing = noise_sd_forcing
  )
  objective <- safe_log_objective(builder)
  eb <- optimise_empirical_bayes(objective, manual_e, lower_e, upper_e)
  ms <- optimise_multistart_empirical_bayes(objective, manual_e, lower_e, upper_e)
  solve_with <- function(parameters) solve_euler_gp(
    t_bc_e, g_bc_e, t_col_e, f_data, t_star_e, p = 1, q = -1,
    sigma_f = parameters[1], ell = parameters[2],
    noise_sd_forcing = noise_sd_forcing
  )
  list(
    g_boundary_initial = g_bc_e, f_data = f_data,
    noise_sd_forcing = noise_sd_forcing,
    objective = objective, eb = eb, ms = ms,
    fixed_fit = solve_with(manual_e), eb_fit = solve_with(eb$par),
    ms_fit = solve_with(ms$par)
  )
}

result_e_nf <- run_euler_case(f_e_exact, 0)
result_e_noisy <- run_euler_case(f_e_noisy, noise_sd_f_e)

# Heat equation: 10 initial values, 16 boundary values.
kappa <- 0.5
x_ic_h <- seq(0, 1, length.out = 12)[2:11]
t_bc_h <- seq(0, 1, length.out = 8)
z_boundary_initial_h <- rbind(
  cbind(x_ic_h, 0), cbind(0, t_bc_h), cbind(1, t_bc_h)
)
heat_truth <- function(x, t) exp(-t) * sin(2 * pi * x)
heat_forcing <- function(x, t) (2 * pi^2 - 1) * exp(-t) * sin(2 * pi * x)
g_boundary_initial_h <- c(
  sin(2 * pi * x_ic_h), rep(0, 2 * length(t_bc_h))
)
z_col_h <- as.matrix(expand.grid(
  x = seq(0.04, 0.96, length.out = 8),
  t = seq(0.02, 1, length.out = 8)
))
f_h_exact <- heat_forcing(z_col_h[, 1], z_col_h[, 2])
grid_h <- expand.grid(x = seq(0, 1, length.out = 60),
                      t = seq(0, 1, length.out = 60))
z_star_h <- as.matrix(grid_h)
truth_h <- heat_truth(grid_h$x, grid_h$t)
noise_sd_f_h <- noise_ratio * max(abs(f_h_exact))
f_h_noisy <- f_h_exact + rnorm(length(f_h_exact),0,noise_sd_f_h)




run_heat_case <- function(f_data, noise_sd_forcing) {
  builder <- function(par) heat_system(
    z_boundary_initial_h, g_boundary_initial_h, z_col_h, f_data, kappa,
    sigma_f = par[1], ell_x = par[2], ell_t = par[3],
    noise_sd_forcing = noise_sd_forcing
  )
  objective <- safe_log_objective(builder)
  eb <- optimise_empirical_bayes(
    objective, start_h, lower_h, upper_h, maxit = 600L
  )
  ms <- optimise_multistart_empirical_bayes(
    objective, start_h, lower_h, upper_h, maxit = 600L
  )
  solve_with <- function(parameters) solve_heat_gp(
    z_boundary_initial_h, g_boundary_initial_h, z_col_h, f_data,
    z_star_h, kappa,
    sigma_f = parameters[1], ell_x = parameters[2], ell_t = parameters[3],
    noise_sd_forcing = noise_sd_forcing
  )
  list(
    g_boundary_initial = g_boundary_initial_h, f_data = f_data,
    noise_sd_forcing = noise_sd_forcing,
    objective = objective, eb = eb, ms = ms,
    fixed_fit = solve_with(manual_h), eb_fit = solve_with(eb$par),
    ms_fit = solve_with(ms$par)
  )
}

result_h_nf <- run_heat_case(f_h_exact, 0)
result_h_noisy <- run_heat_case(f_h_noisy, noise_sd_f_h)

# -----------------------------------------------------------------------------
# Empirical-Bayes collocation studies and numerical summaries
# -----------------------------------------------------------------------------

summarise_fit <- function(problem, fit, truth) {
  covered <- truth >= fit$lower95 & truth <= fit$upper95
  data.frame(
    problem = problem,
    max_abs_error = max(abs(fit$mean - truth)),
    rmse = sqrt(mean((fit$mean - truth)^2)),
    mean_posterior_sd = mean(fit$sd),
    max_posterior_sd = max(fit$sd),
    pointwise_95_coverage = mean(covered),
    effective_jitter = fit$jitter,
    check.names = FALSE
  )
}

study_row <- function(problem, regime, n_col, optimisation, fit, truth,
                      noise_sd_forcing) {
  data.frame(
    problem = problem,
    data_regime = regime,
    collocation_points = n_col,
    boundary_initial_noise_sd = 0,
    forcing_noise_sd = noise_sd_forcing,
    sigma_f = unname(optimisation$par[1]),
    ell = unname(optimisation$par[2]),
    negative_log_marginal_likelihood = optimisation$value,
    rmse = sqrt(mean((fit$mean - truth)^2)),
    max_abs_error = max(abs(fit$mean - truth)),
    pointwise_95_coverage = mean(truth >= fit$lower95 & truth <= fit$upper95),
    optimizer_convergence = optimisation$convergence,
    check.names = FALSE
  )
}

run_poisson_study <- function(n_col, regime) {
  x_col <- (seq_len(n_col) - 0.5) / n_col
  g_exact <- c(0, 0)
  f_exact <- rep(2, n_col)
  if (regime == "Noisy forcing") {
    sd_f <- noise_ratio * max(abs(f_exact))
    
    # N = 4 equals the main experiment
    if (n_col == 4){
      f_data <- f_p_noisy
    }else{
      f_data <- f_exact + rnorm(length(f_exact),0,sd_f)
    }
  } else {
    f_data <- f_exact
    sd_f <- 0
  }
  builder <- function(par) poisson_system(
    x_bc_p, g_exact, x_col, f_data, par[1], par[2], sd_f
  )
  opt <- optimise_empirical_bayes(
    safe_log_objective(builder), manual_p, lower_p, upper_p
  )
  fit <- solve_poisson_gp(
    x_bc_p, g_exact, x_col, f_data, x_star_p,
    sigma_f = opt$par[1], ell = opt$par[2],
    noise_sd_forcing = sd_f
  )
  study_row("Poisson equation", regime, n_col, opt, fit, truth_p, sd_f)
}
study_p <- do.call(rbind, lapply(c("Noise-free", "Noisy forcing"), function(regime) {
  do.call(rbind, lapply(1:30, run_poisson_study, regime = regime))
}))

run_first_order_study <- function(n_col, regime) {
  t_col <- (seq_len(n_col) - 0.5) * 4 / n_col
  g_exact <- 1
  f_exact <- sin(t_col)
  if (regime == "Noisy forcing") {
    sd_f <- noise_ratio * max(abs(f_exact))
    
    # N = 4 equals the main experiment
    if (n_col == 4){
      f_data <- f_1_noisy
    }else{
      f_data <- f_exact + rnorm(length(f_exact),0,sd_f)
    }
  } else {
    f_data <- f_exact
    sd_f <- 0
  }
  builder <- function(par) first_order_system(
    t_ic_1, g_exact, t_col, f_data, alpha = 2,
    sigma_f = par[1], ell = par[2],
    noise_sd_forcing = sd_f
  )
  opt <- optimise_empirical_bayes(
    safe_log_objective(builder), manual_1, lower_1, upper_1
  )
  fit <- solve_first_order_gp(
    t_ic_1, g_exact, t_col, f_data, t_star_1, alpha = 2,
    sigma_f = opt$par[1], ell = opt$par[2],
    noise_sd_forcing = sd_f
  )
  study_row("First-order ODE", regime, n_col, opt, fit, truth_1, sd_f)
}
study_1 <- do.call(rbind, lapply(c("Noise-free", "Noisy forcing"), function(regime) {
  do.call(rbind, lapply(1:30, run_first_order_study, regime = regime))
}))

run_euler_study <- function(n_col, regime) {
  t_col <- 0.5 + (seq_len(n_col) - 0.5) * (3 - 0.5) / n_col
  g_exact <- c(1, 2)
  f_exact <- t_col^2
  if (regime == "Noisy forcing") {
    sd_f <- noise_ratio * max(abs(f_exact))
    
    # N = 4 equals the main experiment
    if (n_col == 4){
      f_data <- f_e_noisy
    }else{
      f_data <- f_exact + rnorm(length(f_exact),0,sd_f)
    }
  } else {
    f_data <- f_exact
    sd_f <- 0
  }
  builder <- function(par) euler_system(
    t_bc_e, g_exact, t_col, f_data, p = 1, q = -1,
    sigma_f = par[1], ell = par[2],
    noise_sd_forcing = sd_f
  )
  opt <- optimise_empirical_bayes(
    safe_log_objective(builder), manual_e, lower_e, upper_e
  )
  fit <- solve_euler_gp(
    t_bc_e, g_exact, t_col, f_data, t_star_e, p = 1, q = -1,
    sigma_f = opt$par[1], ell = opt$par[2],
    noise_sd_forcing = sd_f
  )
  study_row("Euler-Cauchy ODE", regime, n_col, opt, fit, truth_e, sd_f)
}
study_e <- do.call(rbind, lapply(c("Noise-free", "Noisy forcing"), function(regime) {
  do.call(rbind, lapply(1:30, run_euler_study, regime = regime))
}))
ode_collocation_study <- rbind(study_p, study_1, study_e)

comparison_row <- function(problem, regime, method, n_col, fit, truth,
                           parameters, nll, convergence, n_starts,
                           noise_sd_forcing) {
  metrics <- summarise_fit(problem, fit, truth)
  data.frame(
    problem = problem,
    data_regime = regime,
    method = method,
    collocation_points = n_col,
    boundary_initial_noise_sd = 0,
    forcing_noise_sd = noise_sd_forcing,
    sigma_f = unname(parameters[1]),
    ell_x = unname(parameters[2]),
    ell_t = if (length(parameters) >= 3) unname(parameters[3]) else NA_real_,
    negative_log_marginal_likelihood = nll,
    rmse = metrics$rmse,
    max_abs_error = metrics$max_abs_error,
    mean_posterior_sd = metrics$mean_posterior_sd,
    max_posterior_sd = metrics$max_posterior_sd,
    pointwise_95_coverage = metrics$pointwise_95_coverage,
    effective_jitter = metrics$effective_jitter,
    optimizer_convergence = convergence,
    n_starts = n_starts,
    check.names = FALSE
  )
}

case_comparison <- function(problem, regime, n_col, result, truth, manual) {
  rbind(
    comparison_row(
      problem, regime, "Fixed/manual", n_col, result$fixed_fit, truth,
      manual, result$objective(log(manual)), NA_integer_, 0,
      result$noise_sd_forcing
    ),
    comparison_row(
      problem, regime, "Empirical Bayes", n_col, result$eb_fit, truth,
      result$eb$par, result$eb$value, result$eb$convergence, 1,
      result$noise_sd_forcing
    ),
    comparison_row(
      problem, regime, "Multi-start empirical Bayes", n_col,
      result$ms_fit, truth, result$ms$par, result$ms$value,
      result$ms$convergence, result$ms$n_starts,
      result$noise_sd_forcing
    )
  )
}

hyperparameter_comparison <- rbind(
  case_comparison("Poisson equation", "Noise-free", 4, result_p_nf, truth_p, manual_p),
  case_comparison("Poisson equation", "Noisy forcing", 4, result_p_noisy, truth_p, manual_p),
  case_comparison("First-order ODE", "Noise-free", 4, result_1_nf, truth_1, manual_1),
  case_comparison("First-order ODE", "Noisy forcing", 4, result_1_noisy, truth_1, manual_1),
  case_comparison("Euler-Cauchy ODE", "Noise-free", 4, result_e_nf, truth_e, manual_e),
  case_comparison("Euler-Cauchy ODE", "Noisy forcing", 4, result_e_noisy, truth_e, manual_e),
  case_comparison("Heat equation", "Noise-free", 64, result_h_nf, truth_h, manual_h),
  case_comparison("Heat equation", "Noisy forcing", 64, result_h_noisy, truth_h, manual_h)
)

summary_case <- function(problem, regime, result, truth) {
  out <- summarise_fit(problem, result$eb_fit, truth)
  data.frame(
    problem = out$problem,
    data_regime = regime,
    boundary_initial_noise_sd = 0,
    forcing_noise_sd = result$noise_sd_forcing,
    out[, setdiff(names(out), "problem"), drop = FALSE],
    check.names = FALSE
  )
}

numerical_summary <- rbind(
  summary_case("Poisson equation", "Noise-free", result_p_nf, truth_p),
  summary_case("Poisson equation", "Noisy forcing", result_p_noisy, truth_p),
  summary_case("First-order ODE", "Noise-free", result_1_nf, truth_1),
  summary_case("First-order ODE", "Noisy forcing", result_1_noisy, truth_1),
  summary_case("Euler-Cauchy ODE", "Noise-free", result_e_nf, truth_e),
  summary_case("Euler-Cauchy ODE", "Noisy forcing", result_e_noisy, truth_e),
  summary_case("Heat equation", "Noise-free", result_h_nf, truth_h),
  summary_case("Heat equation", "Noisy forcing", result_h_noisy, truth_h)
)

design_row <- function(problem, regime, n_boundary_initial, n_col,
                       n_prediction, result) {
  data.frame(
    problem = problem,
    data_regime = regime,
    hyperparameter_method = "Empirical Bayes (single start)",
    boundary_initial_constraints = n_boundary_initial,
    interior_collocation = n_col,
    total_conditioning = n_boundary_initial + n_col,
    prediction_points = n_prediction,
    boundary_initial_noise_sd = 0,
    forcing_noise_sd = result$noise_sd_forcing,
    sigma_f = result$eb$par[1],
    ell_x = result$eb$par[2],
    ell_t = if (length(result$eb$par) == 3) result$eb$par[3] else NA_real_,
    check.names = FALSE
  )
}

design_summary <- rbind(
  design_row("Poisson equation", "Noise-free", 2, 4, 300, result_p_nf),
  design_row("Poisson equation", "Noisy forcing", 2, 4, 300, result_p_noisy),
  design_row("First-order ODE", "Noise-free", 1, 4, 300, result_1_nf),
  design_row("First-order ODE", "Noisy forcing", 1, 4, 300, result_1_noisy),
  design_row("Euler-Cauchy ODE", "Noise-free", 2, 4, 300, result_e_nf),
  design_row("Euler-Cauchy ODE", "Noisy forcing", 2, 4, 300, result_e_noisy),
  design_row("Heat equation", "Noise-free", 26, 64, 3600, result_h_nf),
  design_row("Heat equation", "Noisy forcing", 26, 64, 3600, result_h_noisy)
)

utils::write.csv(numerical_summary, file.path(output_dir, "numerical_summary.csv"), row.names = FALSE)
utils::write.csv(design_summary, file.path(output_dir, "design_summary.csv"), row.names = FALSE)
utils::write.csv(hyperparameter_comparison,
                 file.path(output_dir, "hyperparameter_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(ode_collocation_study,
                 file.path(output_dir, "ode_collocation_study.csv"),
                 row.names = FALSE)

# -----------------------------------------------------------------------------
# Figure 1: one-dimensional posterior means and pointwise bands
# -----------------------------------------------------------------------------

draw_band_panel <- function(x, truth, fit, value_x, value_y, title, xlab) {
  ylim <- range(truth, fit$lower95, fit$upper95)
  pad <- 0.04 * diff(ylim)
  if (!is.finite(pad) || pad == 0) pad <- 0.04
  ylim <- ylim + c(-pad, pad)
  graphics::plot(x, truth, type = "n", ylim = ylim, xlab = xlab, ylab = "u",
                 main = title, family = "serif", bty = "l")
  graphics::grid(col = "#D9D9D9", lty = 1, lwd = 0.8)
  graphics::polygon(
    c(x, rev(x)), c(fit$lower95, rev(fit$upper95)),
    col = grDevices::adjustcolor("#BDD7EE", 0.70), border = NA
  )
  graphics::lines(x, fit$mean, col = "#0055A4", lwd = 2.0, lty = 1)
  graphics::lines(x, truth, col = "#B2182B", lwd = 1.8, lty = 2)
  graphics::points(value_x, value_y, pch = 19, col = "black", cex = 0.82)
}

grDevices::png(
  file.path(figure_dir, "figure_1_1d_posterior.png"),
  width = 2700, height = 2250, res = 240, type = "cairo"
)
graphics::layout(
  matrix(c(1, 1, 2, 3, 4, 5, 6, 7), nrow = 4, byrow = TRUE),
  heights = c(0.16, 1, 1, 1)
)
old_par <- graphics::par(family = "serif", las = 1)
graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::legend(
  "center", horiz = TRUE, bty = "n", cex = 0.84,
  legend = c("95% posterior interval", "GP posterior mean",
             "Analytical solution", "Exact boundary/initial constraints"),
  fill = c("#BDD7EE", NA, NA, NA), border = c(NA, NA, NA, NA),
  lty = c(NA, 1, 2, NA), lwd = c(NA, 2.0, 1.8, NA),
  col = c(NA, "#0055A4", "#B2182B", "black"),
  pch = c(NA, NA, NA, 19), pt.cex = 0.82
)
graphics::par(mar = c(4.0, 4.2, 3.0, 1.2),
              cex.axis = 0.80, cex.lab = 0.90, cex.main = 0.92)
draw_band_panel(x_star_p, truth_p, result_p_nf$eb_fit, x_bc_p, g_bc_p,
                "Poisson: noise-free", "x")
draw_band_panel(x_star_p, truth_p, result_p_noisy$eb_fit, x_bc_p, g_bc_p,
                "Poisson: 1% noisy forcing", "x")
draw_band_panel(t_star_1, truth_1, result_1_nf$eb_fit, t_ic_1, g_ic_1,
                "First-order ODE: noise-free", "t")
draw_band_panel(t_star_1, truth_1, result_1_noisy$eb_fit, t_ic_1, g_ic_1,
                "First-order ODE: 1% noisy forcing", "t")
draw_band_panel(t_star_e, truth_e, result_e_nf$eb_fit, t_bc_e, g_bc_e,
                "Euler--Cauchy ODE: noise-free", "t")
draw_band_panel(t_star_e, truth_e, result_e_noisy$eb_fit, t_bc_e, g_bc_e,
                "Euler--Cauchy ODE: 1% noisy forcing", "t")
graphics::par(old_par)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# Figure 4: empirical-Bayes RMSE for N = 1, ..., 30 in all three 1D problems
# -----------------------------------------------------------------------------

grDevices::png(
  file.path(figure_dir, "figure_4_ode_rmse_empirical_bayes.png"),
  width = 2700, height = 1050, res = 240, type = "cairo"
)
old_par <- graphics::par(
  mfrow = c(1, 3), mar = c(4.5, 4.6, 3.2, 1.0), family = "serif", las = 1,
  cex.axis = 0.88, cex.lab = 0.94, cex.main = 0.98
)
study_names <- unique(ode_collocation_study$problem)
panel_titles <- c("Poisson equation", "First-order ODE", "Euler--Cauchy ODE")
for (j in seq_along(study_names)) {
  dat_nf <- ode_collocation_study[
    ode_collocation_study$problem == study_names[j] &
      ode_collocation_study$data_regime == "Noise-free", ]
  dat_noisy <- ode_collocation_study[
    ode_collocation_study$problem == study_names[j] &
      ode_collocation_study$data_regime == "Noisy forcing", ]
  all_rmse <- c(dat_nf$rmse, dat_noisy$rmse)
  graphics::plot(
    dat_nf$collocation_points, dat_nf$rmse,
    type = "n", log = "y", xaxt = "n",
    ylim = c(min(all_rmse) * 0.70, max(all_rmse) * 1.45),
    xlab = "Operator collocation points, N",
    ylab = "RMSE (log scale)", main = panel_titles[j], bty = "l"
  )
  graphics::axis(1, at = dat_nf$collocation_points)
  graphics::grid(col = "#D9D9D9", lty = 1, lwd = 0.8)
  graphics::abline(v = 4, col = "#777777", lty = 3, lwd = 1.0)
  graphics::lines(dat_nf$collocation_points, dat_nf$rmse,
                  col = "#0055A4", lwd = 2.2, lty = 1)
  graphics::points(dat_nf$collocation_points, dat_nf$rmse,
                   pch = 19, col = "#0055A4", cex = 0.82)
  graphics::lines(dat_noisy$collocation_points, dat_noisy$rmse,
                  col = "#D95F02", lwd = 2.2, lty = 2)
  graphics::points(dat_noisy$collocation_points, dat_noisy$rmse,
                   pch = 17, col = "#D95F02", cex = 0.86)
  graphics::points(4, dat_nf$rmse[dat_nf$collocation_points == 4],
                   pch = 21, bg = "#0055A4", col = "black", cex = 1.10)
  graphics::points(4, dat_noisy$rmse[dat_noisy$collocation_points == 4],
                   pch = 24, bg = "#D95F02", col = "black", cex = 1.10)
  if (j == 1) {
    graphics::legend(
      "bottomleft",
      legend = c("Noise-free EB", "Noisy-forcing EB", "Main design: N = 4"),
      col = c("#0055A4", "#D95F02", "#777777"),
      lty = c(1, 2, 3), lwd = c(2.2, 2.2, 1.0),
      pch = c(19, 17, NA), bty = "n", cex = 0.70
    )
  }
}
graphics::par(old_par)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# Figure 2: heat-equation exact field, posterior mean, and absolute error
# -----------------------------------------------------------------------------

nx_h <- length(unique(grid_h$x))
nt_h <- length(unique(grid_h$t))
x_grid_h <- sort(unique(grid_h$x))
t_grid_h <- sort(unique(grid_h$t))
as_field <- function(v) matrix(v, nrow = nx_h, ncol = nt_h)

plot_heat_field <- function(values, main, zlim, palette) {
  graphics::image(
    x_grid_h, t_grid_h, as_field(values), col = paper_palette(160, palette),
    zlim = zlim, xlab = "x", ylab = "t", main = main, useRaster = TRUE,
    xaxs = "i", yaxs = "i", asp = 1
  )
  graphics::box()
}

plot_colorbar <- function(zlim, palette, label) {
  yy <- seq(zlim[1], zlim[2], length.out = 160)
  graphics::image(1, yy, matrix(yy, nrow = 1), col = paper_palette(160, palette),
                  zlim = zlim, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  ticks <- pretty(zlim, n = 4)
  ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  graphics::axis(4, at = ticks, las = 2, cex.axis = 0.62, tck = -0.32)
  graphics::mtext(label, side = 3, line = 0.2, cex = 0.62, font = 2)
  graphics::box()
}

fit_h_nf <- result_h_nf$eb_fit
fit_h_noisy <- result_h_noisy$eb_fit
field_lim <- range(truth_h, fit_h_nf$mean, fit_h_noisy$mean)
error_h_nf <- abs(fit_h_nf$mean - truth_h)
error_h_noisy <- abs(fit_h_noisy$mean - truth_h)
error_lim <- c(0, max(error_h_nf, error_h_noisy))

grDevices::png(
  file.path(figure_dir, "figure_2_heat_posterior_fields.png"),
  width = 3000, height = 1900, res = 240, type = "cairo"
)
graphics::layout(matrix(1:12, nrow = 2, byrow = TRUE),
                 widths = rep(c(1, 0.11), 3))
old_par <- graphics::par(family = "serif", oma = c(1.4, 0.5, 0.7, 0.5), las = 1)

graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(truth_h, "Noise-free: exact solution", field_lim, "viridis")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(field_lim, "viridis", "u")
graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(fit_h_nf$mean, "Noise-free: posterior mean", field_lim, "viridis")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(field_lim, "viridis", "u")
graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(error_h_nf, "Noise-free: absolute error", error_lim, "inferno")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(error_lim, "inferno", "|error|")
graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(truth_h, "Noisy forcing: exact solution", field_lim, "viridis")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(field_lim, "viridis", "u")
graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(fit_h_noisy$mean, "Noisy forcing: posterior mean", field_lim, "viridis")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(field_lim, "viridis", "u")
graphics::par(mar = c(3.3, 3.4, 2.7, 0.4)); plot_heat_field(error_h_noisy, "Noisy forcing: absolute error", error_lim, "inferno")
graphics::par(mar = c(3.3, 0.1, 2.7, 1.9)); plot_colorbar(error_lim, "inferno", "|error|")
graphics::mtext(
  "Viridis: exact solution and posterior mean; Inferno: absolute error.",
  side = 1, outer = TRUE, line = 0.15, cex = 0.72
)
graphics::par(old_par)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# Figure 3: heat-equation posterior bands at three time slices
# -----------------------------------------------------------------------------

slice_times <- c(0.05, 0.25, 0.75)
slice_indices <- vapply(slice_times, function(tt) which.min(abs(t_grid_h - tt)), integer(1))

grDevices::png(
  file.path(figure_dir, "figure_3_heat_posterior_slices.png"),
  width = 2700, height = 1500, res = 240, type = "cairo"
)
graphics::layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), nrow = 3, byrow = TRUE),
                 heights = c(0.18, 1, 1))
old_par <- graphics::par(family = "serif", las = 1)
graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::legend(
  "center", horiz = TRUE, bty = "n", cex = 0.82,
  legend = c("95% posterior interval", "GP posterior mean",
             "Analytical solution"),
  fill = c("#BDD7EE", NA, NA), border = c(NA, NA, NA),
  lty = c(NA, 1, 2), lwd = c(NA, 2.0, 1.8),
  col = c(NA, "#0055A4", "#B2182B")
)
graphics::par(mar = c(3.8, 4.0, 3.0, 1.0),
              cex.axis = 0.84, cex.lab = 0.94, cex.main = 0.98)
slice_ylim <- lapply(slice_indices, function(idx) {
  rows <- seq_len(nx_h) + (idx - 1L) * nx_h
  lim <- range(
    truth_h[rows], fit_h_nf$lower95[rows], fit_h_nf$upper95[rows],
    fit_h_noisy$lower95[rows], fit_h_noisy$upper95[rows]
  )
  pad <- 0.04 * diff(lim)
  if (!is.finite(pad) || pad == 0) pad <- 0.01
  lim + c(-pad, pad)
})
slice_fits <- list("Noise-free" = fit_h_nf, "Noisy forcing" = fit_h_noisy)
for (regime in names(slice_fits)) {
  fit_slice <- slice_fits[[regime]]
  for (j in seq_along(slice_indices)) {
    idx <- slice_indices[j]
    rows <- seq_len(nx_h) + (idx - 1L) * nx_h
    graphics::plot(
      x_grid_h, truth_h[rows], type = "n", ylim = slice_ylim[[j]],
      xlab = "x", ylab = "u(x,t)",
      main = sprintf("%s, t = %.3f", regime, t_grid_h[idx]), bty = "l"
    )
    graphics::grid(col = "#D9D9D9", lty = 1, lwd = 0.8)
    graphics::polygon(
      c(x_grid_h, rev(x_grid_h)),
      c(fit_slice$lower95[rows], rev(fit_slice$upper95[rows])),
      col = grDevices::adjustcolor("#BDD7EE", 0.70), border = NA
    )
    graphics::lines(x_grid_h, fit_slice$mean[rows],
                    col = "#0055A4", lwd = 2.0, lty = 1)
    graphics::lines(x_grid_h, truth_h[rows],
                    col = "#B2182B", lwd = 1.8, lty = 2)
  }
}
graphics::par(old_par)
grDevices::dev.off()

cat("\nDesign summary\n")
print(design_summary, row.names = FALSE)
cat("\nNumerical summary\n")
print(numerical_summary, row.names = FALSE)
cat("\nHyperparameter comparison\n")
print(hyperparameter_comparison, row.names = FALSE)
cat("\nOne-dimensional empirical-Bayes collocation study\n")
print(ode_collocation_study, row.names = FALSE)
cat(sprintf("\nFigures and CSV summaries written to: %s\n", output_dir))


##########Indometh example ##############
# Four observations (indices 1, 4, 7, 10) are used 
# SSbiexp is fitted only to these four observations

data("Indometh", package = "datasets")
indometh <- as.data.frame(Indometh)
indometh$Subject <- as.character(indometh$Subject)
train_index <- c(1L, 4L, 7L, 10L)
t_grid <- seq(0, 8, length.out = 401L)

# ordinary GP
fit_sparse_gp <- function(t, y, t_new) {
  objective <- function(log_par) {
    par <- exp(log_par)
    gp_negative_log_marginal(
      se_kernel_1d(t, t, par[1], par[2]), y,
      observation_variance = rep(par[3]^2, length(y))
    )
  }
  eb <- optimise_multistart_empirical_bayes(
    objective,
    start = c(max(y), 0.8, max(0.02, 0.05 * stats::sd(y))),
    lower = c(0.03, 0.03, 0.005),
    upper = c(30, 20, 3)
  )
  par <- eb$par
  fit <- condition_gp(
    se_kernel_1d(t, t, par[1], par[2]), y,
    se_kernel_1d(t_new, t, par[1], par[2]),
    rep(par[1]^2, length(t_new)),
    rep(par[3]^2, length(y))
  )
  fit$sigma_f <- par[1]
  fit$ell <- par[2]
  fit$sigma_y <- par[3]
  fit
}

# SSbiexp fit parameters
fit_four_point_ssbiexp <- function(dat) {
  total_start <- 1.5 * max(dat$conc)
  start_grid <- expand.grid(
    fast = c(1, 3, 8), slow = c(0.05, 0.15, 0.4),
    fast_fraction = c(0.4, 0.7)
  )
  fits <- lapply(seq_len(nrow(start_grid)), function(i) {
    s <- start_grid[i, ]
    try(
      stats::nls(
        conc ~ stats::SSbiexp(time, A1, lrc1, A2, lrc2),
        data = dat,
        start = list(
          A1 = total_start * s$fast_fraction,
          lrc1 = log(s$fast),
          A2 = total_start * (1 - s$fast_fraction),
          lrc2 = log(s$slow)
        ),
        algorithm = "port",
        lower = c(A1 = 0, lrc1 = log(0.01),
                  A2 = 0, lrc2 = log(0.01)),
        upper = c(A1 = 30, lrc1 = log(30),
                  A2 = 30, lrc2 = log(30)),
        control = stats::nls.control(maxiter = 1000, warnOnly = TRUE)
      ),
      silent = TRUE
    )
  })
  valid <- vapply(
    fits,
    function(x) !inherits(x, "try-error") &&
      isTRUE(x$convInfo$isConv) && all(is.finite(stats::coef(x))),
    logical(1)
  )
  if (!any(valid)) stop("The four-point SSbiexp fit did not converge.")
  fits <- fits[valid]
  sse <- vapply(fits, function(x) sum(stats::residuals(x)^2), numeric(1))
  fits[[which.min(sse)]]
}

#ode-constrained GP
fit_sparse_ode_gp <- function(t, y, t_new, sigma_y, sigma_f, ell,
                              A_fast, A_slow,
                              lambda_fast, lambda_slow) {
  a <- lambda_fast + lambda_slow
  b <- lambda_fast * lambda_slow
  tc <- (seq_len(64L) - 0.5) * 8 / 64
  ode_mean <- function(x) {
    A_fast * exp(-lambda_fast * x) + A_slow * exp(-lambda_slow * x)
  }

  ddd <- se_derivatives_1d(t, t, sigma_f, ell)
  ddc <- se_derivatives_1d(t, tc, sigma_f, ell)
  dcc <- se_derivatives_1d(tc, tc, sigma_f, ell)

  K_data_L <- ddc$dyy + a * ddc$dy + b * ddc$k
  K_LL <- dcc$dxxyy + a * dcc$dxxy + b * dcc$dxx +
    a * dcc$dxyy + a^2 * dcc$dxy + a * b * dcc$dx +
    b * dcc$dyy + a * b * dcc$dy + b^2 * dcc$k
  K_obs <- rbind(
    cbind(ddd$k, K_data_L),
    cbind(t(K_data_L), K_LL)
  )

  dnd <- se_derivatives_1d(t_new, t, sigma_f, ell)
  dnc <- se_derivatives_1d(t_new, tc, sigma_f, ell)
  K_new_obs <- cbind(
    dnd$k, dnc$dyy + a * dnc$dy + b * dnc$k
  )
  operator_sd <- 1e-5 * sigma_f
  fit <- condition_gp(
    K_obs,
    c(y - ode_mean(t), rep(0, length(tc))),
    K_new_obs,
    rep(sigma_f^2, length(t_new)),
    c(rep(sigma_y^2, length(y)), rep(operator_sd^2, length(tc))),
    relative_jitter = 1e-8
  )
  fit$mean <- ode_mean(t_new) + fit$mean
  fit$lower95 <- fit$mean - 1.96 * fit$sd
  fit$upper95 <- fit$mean + 1.96 * fit$sd
  fit
}

subjects <- sort(unique(indometh$Subject))
indometh_result <- vector("list", length(subjects))
names(indometh_result) <- subjects

for (id in subjects) {
  dat <- indometh[indometh$Subject == id, ]
  dat <- dat[order(dat$time), ]

  train_dat <- dat[train_index, , drop = FALSE]
  ss_fit <- fit_four_point_ssbiexp(train_dat)
  cf <- stats::coef(ss_fit)
  component <- data.frame(
    A = c(unname(cf["A1"]), unname(cf["A2"])),
    lambda = exp(c(unname(cf["lrc1"]), unname(cf["lrc2"])))
  )
  component <- component[order(component$lambda, decreasing = TRUE), ]
  Af <- component$A[1]
  As <- component$A[2]
  lf <- component$lambda[1]
  ls <- component$lambda[2]
  all_times <- c(t_grid, dat$time)
  gp <- fit_sparse_gp(
    train_dat$time, train_dat$conc, all_times
  )
  ode_gp <- fit_sparse_ode_gp(
    train_dat$time, train_dat$conc, all_times,
    gp$sigma_y, gp$sigma_f, gp$ell,
    A_fast = Af, A_slow = As,
    lambda_fast = lf, lambda_slow = ls
  )
  obs_rows <- length(t_grid) + seq_len(nrow(dat))
  check <- setdiff(seq_len(nrow(dat)), train_index)

  indometh_result[[id]] <- list(
    data = dat,
    gp = gp,
    ode_gp = ode_gp,
    ss_curve = Af * exp(-lf * t_grid) + As * exp(-ls * t_grid),
    metrics = data.frame(
      Subject = id,
      model = c("Ordinary GP", "ODE-GP"),
      check_RMSE = c(
        sqrt(mean((dat$conc[check] - gp$mean[obs_rows][check])^2)),
        sqrt(mean((dat$conc[check] - ode_gp$mean[obs_rows][check])^2))
      )
    )
  )
}

indometh_metrics <- do.call(
  rbind, lapply(indometh_result, function(x) x$metrics)
)
utils::write.csv(
  indometh_metrics,
  file.path(output_dir, "Indometh_GP_vs_ODE_GP_metrics.csv"),
  row.names = FALSE
)

grDevices::png(
  file.path(output_dir, "Indometh_GP_vs_ODE_GP.png"),
  width = 2500, height = 1550, res = 220, type = "cairo"
)
old_par <- graphics::par(
  mfrow = c(2, 3), mar = c(4.1, 4.3, 3.0, 1.0), las = 1
)
for (id in subjects) {
  x <- indometh_result[[id]]
  dat <- x$data
  n_grid <- length(t_grid)
  ylim <- range(
    0, dat$conc, x$gp$upper95[seq_len(n_grid)],
    x$ode_gp$upper95[seq_len(n_grid)]
  )
  graphics::plot(
    dat$time, dat$conc, type = "n", ylim = ylim,
    xlab = "Time (hour)", ylab = "Concentration (mcg/ml)",
    main = paste("Subject", id), bty = "l"
  )
  graphics::grid(col = "#DDDDDD")
  graphics::polygon(
    c(t_grid, rev(t_grid)),
    c(x$gp$lower95[seq_len(n_grid)],
      rev(x$gp$upper95[seq_len(n_grid)])),
    border = NA, col = grDevices::adjustcolor("#9ECAE1", 0.55)
  )
  # Conditional ODE-GP uncertainty given the four-point plug-in SSbiexp
  # parameters.  It does not include uncertainty in those fitted parameters.
  graphics::polygon(
    c(t_grid, rev(t_grid)),
    c(x$ode_gp$lower95[seq_len(n_grid)],
      rev(x$ode_gp$upper95[seq_len(n_grid)])),
    border = NA, col = grDevices::adjustcolor("#A1D99B", 0.45)
  )
  graphics::lines(t_grid, x$gp$mean[seq_len(n_grid)],
                  col = "#2171B5", lwd = 2)
  graphics::lines(t_grid, x$ode_gp$mean[seq_len(n_grid)],
                  col = "#238B45", lwd = 2.2)
  graphics::lines(t_grid, x$ss_curve,
                  col = "#D95F0E", lwd = 1.4, lty = 2)
  graphics::points(dat$time[-train_index], dat$conc[-train_index],
                   pch = 1, col = "#CB181D", lwd = 1.3)
  graphics::points(dat$time[train_index], dat$conc[train_index], pch = 19)
  if (id == subjects[1]) {
    graphics::legend(
      "topright",
      c("Ordinary GP", "ODE-GP", "4-point SSbiexp",
        "4 training", "7 check", "95% bands"),
      col = c("#2171B5", "#238B45", "#D95F0E",
              "black", "#CB181D", "#74C476"),
      lty = c(1, 1, 2, NA, NA, 1),
      lwd = c(2, 2, 1.4, NA, NA, 7),
      pch = c(NA, NA, NA, 19, 1, NA),
      bty = "n", cex = 0.58
    )
  }
}
graphics::par(old_par)
grDevices::dev.off()

cat("\nIndometh ordinary GP versus ODE-GP\n")
print(indometh_metrics, row.names = FALSE)
