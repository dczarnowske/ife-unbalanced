######################################### Appendix #################################################
# Comparison of Matrix Completion Algorithms
#
####################################################################################################

### load packages and define functions

# load package
library(bench)
library(InteractiveEffects)
library(pbmcapply) # Note: Unix OS required
library(softImpute)

# cross validation for soft imputation
cv_nu <- function(w, observed, r_max, n_trials = 30L, n_rep = 5L) {
  n <- nrow(w)
  t <- ncol(w)
  observed <- which(observed)
  nt_obs <- length(observed)
  scaling <- sqrt((n * t) / nt_obs)
  nu_max <- lambda0(w)
  nu_grid <- exp(seq(
    log(nu_max * 0.99), log(nu_max * 0.001),
    length = n_trials
  )) * scaling
  cv_mse <- matrix(NA_real_, nrow = length(nu_grid), ncol = n_rep)
  for (k in seq.int(n_rep)) {
    observed_train <- sample(observed, floor(nt_obs^2 / (n * t)))
    observed_test <- setdiff(observed, observed_train)
    w_train <- matrix(NA_real_, n, t)
    w_train[observed_train] <- w[observed_train]
    fit <- NULL
    for (j in seq_along(nu_grid)) {
      fit <- suppressWarnings(softImpute(
        w_train,
        rank.max = r_max, lambda = nu_grid[[j]],
        type = "svd", warm.start = fit
      ))
      w_hat <- impute(fit, row(w)[observed_test], col(w)[observed_test])
      cv_mse[j, k] <- mean((w[observed_test] - w_hat)^2)
    }
  }
  nu_grid[[which.min(rowMeans(cv_mse))]]
}

# plug-in approach for soft imputation
plugin_nu <- function(w, observed, r_max, n_iter = 15L, n_sim = 1000L) {
  n <- nrow(w)
  t <- ncol(w)
  sigma_hat <- sd(as.vector(w), na.rm = TRUE)
  z_ast <- sapply(seq.int(n_sim), function(r, n, t, observed, sigma_hat) {
    z <- matrix(rnorm(n * t, 0.0, sigma_hat), n, t)
    z[!observed] <- 0.0
    norm(z, "2")
  }, n = n, t = t, observed = observed, sigma_hat = sigma_hat)
  nu <- 1.1 * quantile(z_ast, 0.95)
  fit <- NULL
  for (k in seq.int(n_iter)) {
    nu_old <- nu
    fit <- softImpute(w, rank.max = r_max, lambda = nu, type = "svd", warm.start = fit)
    w_hat <- tcrossprod(fit[["u"]], fit[["v"]] * fit[["d"]])
    u_hat <- w - w_hat
    sigma_hat <- sqrt(mean(u_hat^2, na.rm = TRUE))
    z_ast <- sapply(seq.int(n_sim), function(r, n, t, observed, sigma_hat) {
      z <- matrix(rnorm(n * t, 0.0, sigma_hat), n, t)
      z[!observed] <- 0.0
      norm(z, "2")
    }, n = n, t = t, observed = observed, sigma_hat = sigma_hat)
    nu <- 1.1 * quantile(z_ast, 0.95)
    if (abs(nu - nu_old) / nu_old < 0.01) break
  }
  names(nu) <- NULL
  nu
}

# data generation
simulate_data <- function(n, t, psi = 0.0, error_term = c("iid_norm", "xshet_chisq")) {
  # validity checks
  n <- as.integer(n)
  if (n < 1L) stop("Number of cross-sectional units should be at least one.", call. = FALSE)
  t <- as.integer(t)
  if (t < 1L) stop("Number of time periods should be at least one.", call. = FALSE)
  psi <- as.double(psi)
  if (psi < 0.0 || psi > 0.5) stop("'psi' has to be in [0, 0.5].", call. = FALSE)
  error_term <- match.arg(error_term)

  # set some fixed parameters
  burnin <- 1000L
  r <- 2L
  rho_f <- 0.5
  sigma_f <- 0.5

  # auxiliary variables
  id <- rep.int(seq.int(n), t)
  time <- rep(seq.int(t), each = n)

  # generate dependent variable
  if (error_term == "iid_norm") {
    E <- matrix(rnorm(n * (t + burnin + 1L), 0.0, sqrt(3.0)), n, t + burnin + 1L)
  } else {
    E <- t(sapply(seq.int(n), function(x, t) {
      if (x %% 2L) {
        u <- sqrt(2.0) * (rchisq(t, 5L) - 5.0) / sqrt(10.0)
      } else {
        u <- 2.0 * (rchisq(t, 5L) - 5.0) / sqrt(10.0)
      }
      u
    }, t = t + burnin + 1L))
  }
  lambda <- matrix(rnorm(n * r, 1.0, 1.0), n, r)
  ft <- matrix(NA_real_, r, t + burnin + 1L)
  ft[, 1L] <- rnorm(r, 0.0, sigma_f)
  Theta <- matrix(NA_real_, n, t + burnin + 1L)
  Theta[, 1L] <- as.vector(lambda %*% ft[, 1])
  Y <- matrix(NA_real_, n, t + burnin + 1L)
  Y[, 1L] <- Theta[, 1L] + E[, 1L]
  for (s in seq.int(2L, t + burnin + 1L)) {
    ft[, s] <- rho_f * ft[, s - 1L] + rnorm(r, 0.0, sqrt(1 - rho_f^2) * sigma_f)
    Theta[, s] <- as.vector(lambda %*% ft[, s])
    Y[, s] <- Theta[, s] + E[, s]
  }
  Theta <- Theta[, -seq.int(burnin)]
  Theta <- Theta[, -1L]
  Y <- Y[, -seq.int(burnin)]
  Y <- Y[, -1L]

  # transform matrices to vectors
  y <- as.vector(Y)
  theta <- as.vector(Theta)
  lambda_long <- rep(lambda[, 1L], t)

  # generate balanced panel data set
  df <- data.frame(id, time, y, theta)
  df <- df[order(time, lambda_long, id), ]

  # determine which data to drop if requested
  if (psi > 0.0) {
    # time series of type 1 individuals start randomly at (1, t2 + 1) and last for t2 periods
    alpha <- 2.0 * psi
    n1 <- floor(alpha * n)
    n2 <- n - n1
    t1 <- floor(t / 2)
    keep <- as.vector(rbind(
      t(sapply(sample.int(t1 + 1L, n1, replace = TRUE), function(x) {
        c(rep(FALSE, x - 1L), rep(TRUE, t1), rep(FALSE, t - t1 - x + 1L))
      })),
      matrix(TRUE, n2, t)
    ))
    df[["keep"]] <- keep
  } else {
    df[["keep"]] <- TRUE
  }

  # return panel data set
  df
}

# simulation function (returns one iteration for each configuration)
sim_func <- function(x, configs) {
  res_list <- mapply(
    function(psi, nbar) {
      # generate synthetic data
      seed <- sample(.Machine[["integer.max"]], 1L)
      n <- round(nbar / (1.0 - psi))
      df <- simulate_data(n, n, psi, "xshet_chisq")
      df <- panel_data(df, xs_id = "id", time_id = "time")[["pd"]]
      pd <- panel_data(subset(df, keep == TRUE), xs_id = "id", time_id = "time")
      
      # true gamma
      r <- 2L
      theta_true <- df[["theta"]]
      
      # prepare matrix to complete
      observed <- pd[["observed"]]
      gamma <- matrix(NA_real_, n, n)
      gamma[observed] <- pd[["pd"]][["y"]]
      
      # selection of \nu
      nu_cv <- cv_nu(gamma, observed, r)
      nu_pi <- plugin_nu(gamma, observed, r)
      
      # bias and rmse - em algorithm
      gamma_ast <- InteractiveEffects:::matrix_completion(gamma, r, observed, 0.0)
      theta_hat <- as.vector(InteractiveEffects:::principal_components(gamma_ast, r)[["theta"]])
      bias_em_all <- mean(abs((theta_true - theta_hat)))
      bias_em_obs <- mean(abs((theta_true - theta_hat)[observed]))
      bias_em_miss <- mean(abs((theta_true - theta_hat)[!observed]))
      rmse_em_all <- sqrt(mean(((theta_true - theta_hat))^2))
      rmse_em_obs <- sqrt(mean(((theta_true - theta_hat)[observed])^2))
      rmse_em_miss <- sqrt(mean(((theta_true - theta_hat)[!observed])^2))
      
      # bias and rmse - redebias algorithm - cv
      gamma_ast <- InteractiveEffects:::matrix_completion(gamma, r, observed, nu_cv)
      theta_hat <- as.vector(InteractiveEffects:::principal_components(gamma_ast, r)[["theta"]])
      bias_rdcv_all <- mean(abs((theta_true - theta_hat)))
      bias_rdcv_obs <- mean(abs((theta_true - theta_hat)[observed]))
      bias_rdcv_miss <- mean(abs((theta_true - theta_hat)[!observed]))
      rmse_rdcv_all <- sqrt(mean(((theta_true - theta_hat))^2))
      rmse_rdcv_obs <- sqrt(mean(((theta_true - theta_hat)[observed])^2))
      rmse_rdcv_miss <- sqrt(mean(((theta_true - theta_hat)[!observed])^2))
      
      # bias and rmse - redebias algorithm - pi
      gamma_ast <- InteractiveEffects:::matrix_completion(gamma, r, observed, nu_pi)
      theta_hat <- as.vector(InteractiveEffects:::principal_components(gamma_ast, r)[["theta"]])
      bias_rdpi_all <- mean(abs((theta_true - theta_hat)))
      bias_rdpi_obs <- mean(abs((theta_true - theta_hat)[observed]))
      bias_rdpi_miss <- mean(abs((theta_true - theta_hat)[!observed]))
      rmse_rdpi_all <- sqrt(mean(((theta_true - theta_hat))^2))
      rmse_rdpi_obs <- sqrt(mean(((theta_true - theta_hat)[observed])^2))
      rmse_rdpi_miss <- sqrt(mean(((theta_true - theta_hat)[!observed])^2))
      
      # measure time - mean out of 30 replications
      time <- summary(mark(
        EM    = InteractiveEffects:::matrix_completion(gamma, r, observed, 0.0),
        DRMC1 = InteractiveEffects:::matrix_completion(gamma, r, observed, nu_cv),
        DRMC2 = InteractiveEffects:::matrix_completion(gamma, r, observed, nu_pi),
        check = FALSE
      ))[["median"]]
      time_rdcv_em <- as.double(time[[2L]] / time[[1L]])
      time_rdpi_em <- as.double(time[[3L]] / time[[1L]])
      
      # return data.frame with results
      data.frame(
        psi = psi,
        n = n,
        t = n,
        nbar = nbar,
        tbar = nbar,
        seed = seed,
        r = 2L,
        nu_cv = nu_cv,
        nu_pi = nu_pi,
        bias_em_all = bias_em_all,
        bias_em_obs = bias_em_obs,
        bias_em_miss = bias_em_miss,
        bias_rdcv_all = bias_rdcv_all,
        bias_rdcv_obs = bias_rdcv_obs,
        bias_rdcv_miss = bias_rdcv_miss,
        bias_rdpi_all = bias_rdpi_all,
        bias_rdpi_obs = bias_rdpi_obs,
        bias_rdpi_miss = bias_rdpi_miss,
        rmse_em_all = rmse_em_all,
        rmse_em_obs = rmse_em_obs,
        rmse_em_miss = rmse_em_miss,
        rmse_rdcv_all = rmse_rdcv_all,
        rmse_rdcv_obs = rmse_rdcv_obs,
        rmse_rdcv_miss = rmse_rdcv_miss,
        rmse_rdpi_all = rmse_rdpi_all,
        rmse_rdpi_obs = rmse_rdpi_obs,
        rmse_rdpi_miss = rmse_rdpi_miss,
        time_rdcv_em = time_rdcv_em,
        time_rdpi_em = time_rdpi_em,
        stringsAsFactors = FALSE
      )
    },
    psi = configs[["psi"]],
    nbar = configs[["nbar"]],
    SIMPLIFY = FALSE
  )
  res <- Reduce(rbind, res_list)
  rownames(res) <- seq.int(length(res_list))
  res
}


### start simulation

# choose initial seed
set.seed(4321L)
text <- "app"

# different configurations
R <- 1000L
nbar <- c(60L, 120L, 180L)
psi <- 0.4

# generate data.frame with all possible combinations
configs <- expand.grid(
  psi              = psi,
  nbar             = nbar,
  stringsAsFactors = FALSE
)

# start simulation
res_list <- pbmclapply(seq.int(R), sim_func, configs = configs, mc.cores = 30L)
res <- Reduce(rbind, res_list)

# Save results
saveRDS(res, file = paste0("results/mc_comparison_", text, ".rds"), compress = "xz")