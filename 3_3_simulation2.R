######################################### Online Supplement ########################################
# Additional Monte Carlo Simulation
#
####################################################################################################

### load packages and define functions

# load package
library(pbmcapply) # Note: Unix OS required
library(InteractiveEffects)

# data generation
simulate_data <- function(n, t, seed, errorterm = c("homo", "homo_fattail", "xshet", "xshet_tshet"), missingdata = c("none", "pattern1", "pattern2", "pattern3"), psi = NULL) {
  # validity checks
  n <- as.integer(n)
  if (n < 1L) stop("Number of cross-sectional units should be at least one.", call. = FALSE)
  t <- as.integer(t)
  if (t < 1L) stop("Number of time periods should be at least one.", call. = FALSE)
  seed <- as.integer(seed)
  if (seed < 0L) stop("'seed' has to be a positive integer.", call. = FALSE)
  psi <- as.double(psi)
  if (psi < 0.0 || psi > 0.5) stop("'psi' has to be in [0, 0.5].", call. = FALSE)
  errorterm <- match.arg(errorterm)
  missingdata <- match.arg(missingdata)

  # fixed parameters
  set.seed(seed)
  beta <- 1.0
  r <- 2L
  rhoe <- 0.5
  theta <- 2.0

  # panel identifiers
  id <- rep.int(seq.int(n), t)
  time <- rep(seq.int(t), each = n)

  # generate factor structure and regressor
  lambda <- matrix(rnorm(n * r, 1.0), n, r)
  chi <- matrix(rnorm(n * r, 1.0), n, r)
  f <- matrix(rnorm((t + 1L) * r), t + 1L, r)
  X <- 1.0 + matrix(rnorm(n * t), n, t) + tcrossprod(lambda + chi, f[-1L, ] + f[-t - 1L, ])

  # generate idiosyncratic error term
  if (errorterm == "homo") {
    # a) homoskedasticity
    E <- matrix(rnorm(n * t, 0.0, 2.0), n, t)
  } else if (errorterm == "homo_fattail") {
    # b) fat tail
    E <- matrix(rt(n * t, 5L), n, t) / sqrt(5.0 / 12.0)
  } else if (errorterm == "xshet") {
    # c) heteroskedasticity in cross-sectional dimension
    E <- t(sapply(seq.int(n), function(x, t) {
      if (x %% 2L) {
        u <- rnorm(t, 0.0, sqrt(2.0)) # Odd
      } else {
        u <- rnorm(t, 0.0, sqrt(6.0)) # Even
      }
      u
    }, t = t))
  } else {
    # d) heteroskedasticity in both dimensions
    E1 <- t(sapply(seq.int(n), function(x, t) {
      if (x %% 2L) {
        u <- rnorm(t) # Odd
      } else {
        u <- rnorm(t, 0.0, sqrt(3.0)) # Even
      }
      u
    }, t = t))
    E2 <- sapply(seq.int(t), function(x, n) {
      if (x %% 2L) {
        u <- rnorm(n) # Odd
      } else {
        u <- rnorm(n, 0.0, sqrt(3.0)) # Even
      }
      u
    }, n = n)
    E <- E1 + E2
  }

  # generate dependent variable
  Y <- beta * X + tcrossprod(lambda, f[-1L, ]) + E

  # transform matrices to vectors
  y <- as.vector(Y)
  x <- as.vector(X)
  errorterm <- as.vector(E)

  # generate balanced panel data set
  df <- data.frame(id, time, y, x, errorterm)

  # randomly drop data if requested
  if (missingdata != "none") {
    # check validity of 'psi'
    if (is.null(psi)) {
      stop("'psi' has to be specified.", call. = FALSE)
    } else {
      psi <- as.double(psi)
      if (psi < 0.0 || psi > 0.5) {
        stop("'psi' has to be in [0, 0.5].", call. = FALSE)
      }
    }

    # declare missing observations based on requested missing data pattern
    if (missingdata == "pattern1") {
      # randomly drop observations
      keep <- sort(sample.int(n * t, round(n * t * (1.0 - psi))))
    } else {
      # set additional parameters
      alpha <- 2.0 * psi
      n1 <- round(alpha * n)
      n2 <- n - n1
      t1 <- round(t / 2)

      # different missing data patterns
      if (missingdata == "pattern2") {
        # all time series start at t = 1 but type 1 individuals drop out earlier at t = t1.
        keep <- as.vector(rbind(
          matrix(c(rep(TRUE, t1), rep(FALSE, t - t1)), n1, t, byrow = TRUE),
          matrix(TRUE, n2, t)
        ))
      } else {
        # time series of type 1 individuals start randomly at (1, t2 + 1) and last for t2 periods
        keep <- as.vector(rbind(
          t(sapply(sample.int(t1 + 1L, n1, replace = TRUE), function(x) {
            c(rep(FALSE, x - 1L), rep(TRUE, t1), rep(FALSE, t - t1 - x + 1L))
          })),
          matrix(TRUE, n2, t)
        ))
      }
    }

    # drop observations from balanced panel data set
    df <- df[keep, ]
  }

  # return panel data set
  df
}

# simulation function (returns one iteration for each configuration)
sim_func <- function(x, configs) {
  res_list <- mapply(
    function(missingdata, psi, errorterm, nbar, tbar) {
      # generate synthetic data
      seed <- sample(.Machine[["integer.max"]], 1L)
      n <- round(nbar / (1.0 - psi))
      t <- round(tbar / (1.0 - psi))
      df <- simulate_data(
        n, t, seed,
        errorterm = errorterm,
        missingdata = missingdata,
        psi = psi
      )
      df <- panel_data(df, xs_id = "id", time_id = "time")

      # debiased LS interactive effects estimator (em algorithm for imputation)
      fit <- lm_ie(y ~ x, df, 2L, nu = 0.0)
      fit <- bias_correction(fit, L = 0L, xs_het = TRUE, ts_het = TRUE)
      betals <- coef(fit)
      sebetals <- sqrt(diag(vcov(fit)))

      # estimate number of factors using different estimators
      rmax <- round(12.0 * (min(nbar, tbar) / 100.0)^0.25)
      fit <- lm_ie(y ~ x, df, rmax, nu = 0.0)
      nof <- number_of_factors(fit)

      # return data.frame with results
      data.frame(
        missingdata = missingdata,
        psi = psi,
        errorterm = errorterm,
        nbar = nbar,
        tbar = tbar,
        seed = seed,
        beta_ls = betals,
        se_beta_ls = sebetals,
        r = 2L,
        rmax = rmax,
        as.list(nof),
        stringsAsFactors = FALSE
      )
    },
    missingdata = configs[["missingdata"]],
    psi = configs[["psi"]],
    errorterm = configs[["errorterm"]],
    nbar = configs[["nbar"]],
    tbar = configs[["tbar"]],
    SIMPLIFY = FALSE
  )
  res <- Reduce(rbind, res_list)
  rownames(res) <- seq.int(length(res_list))
  res
}


### start simulation

# choose initial seed
set.seed(4321L)
text <- "supp"

# different configurations
R <- 500L
errorterm <- c("homo", "homo_fattail", "xshet", "xshet_tshet")
nbar <- c(120L, 240L)
missingdata <- c("none", paste0("pattern", seq.int(3L)))
psi <- c(0.0, 0.2, 0.4)
tbar <- c(24L, 48L, 96L)

# generate data.frame with all possible combinations
configs <- expand.grid(
  missingdata      = missingdata,
  psi              = psi,
  errorterm        = errorterm,
  nbar             = nbar,
  tbar             = tbar,
  stringsAsFactors = FALSE
)

# drop duplicated configurations
configs <- subset(configs, !(missingdata == "none" & psi != 0.0))
configs <- subset(configs, !(missingdata %in% paste0("pattern", seq.int(3L)) & psi == 0.0))
configs <- configs[order(
  configs[["missingdata"]],
  configs[["psi"]],
  configs[["errorterm"]],
  configs[["nbar"]],
  configs[["tbar"]]
), ]

# start simulation
res_list <- pbmclapply(seq.int(R), sim_func, configs = configs, mc.cores = 30L)
res <- Reduce(rbind, res_list)

# Save results
saveRDS(res, file = paste0("results/ife_simulation2_", text, ".rds"), compress = "xz")
