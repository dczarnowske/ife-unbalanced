######################################### Section 4 ################################################
# Monte Carlo Simulation
#
####################################################################################################

### load packages and define functions

# load package
library(pbmcapply) # Note: Unix OS required
library(InteractiveEffects)

# simulation function (returns one iteration for each configuration)
sim_func <- function(x, configs) {
  res_list <- mapply(
    function(psi, rho, nbar, tbar, L, rmax) {
      # generate data
      seed <- sample(.Machine[["integer.max"]], 1L)
      n <- round(nbar / (1.0 - psi))
      t <- round(tbar / (1.0 - psi))
      df <- synthetic_data(n, t, seed, rho = rho, psi = psi)
      df <- panel_data(df, xs_id = "id", time_id = "time")
      
      # debiased LS interactive effects estimator (em algorithm for imputation)
      fit <- lm_ie(y ~ y_lag, df, 1L, nu = 0.0)
      fit <- bias_correction(fit, L = L, xs_het = TRUE, ts_het = TRUE)
      beta_em <- coef(fit)
      se_beta_em <- sqrt(diag(vcov(fit)))
      
      # debiased LS interactive effects estimator (redebias algorithm for imputation);
      # plug-in approach for \nu
      fit <- lm_ie(y ~ y_lag, df, 1L)
      fit <- bias_correction(fit, L = L, xs_het = TRUE, ts_het = TRUE)
      beta_rd <- coef(fit)
      se_beta_rd <- sqrt(diag(vcov(fit)))

      # estimate number of factors using different estimators
      fit <- lm_ie(y ~ y_lag, df, rmax, nu = 0.0)
      nof_rescale_em <- number_of_factors(fit, rescale = TRUE)
      names(nof_rescale_em) <- paste0(names(nof_rescale_em), "_rescale_em")
      nof_em <- number_of_factors(fit, rescale = FALSE)
      names(nof_em) <- paste0(names(nof_em), "_em")
      
      # estimate number of factors using different estimators
      fit <- lm_ie(y ~ y_lag, df, rmax)
      nof_rescale_rd <- number_of_factors(fit, rescale = TRUE)
      names(nof_rescale_rd) <- paste0(names(nof_rescale_rd), "_rescale_rd")
      nof_rd <- number_of_factors(fit, rescale = FALSE)
      names(nof_rd) <- paste0(names(nof_rd), "_rd")

      # return results
      data.frame(
        psi = psi,
        rho = rho,
        nbar = nbar,
        tbar = tbar,
        seed = seed,
        L = L,
        beta_em = beta_em,
        se_beta_em = se_beta_em,
        beta_rd = beta_rd,
        se_beta_rd = se_beta_rd,
        r = 1,
        rmax = rmax,
        as.list(nof_rescale_em),
        as.list(nof_em),
        as.list(nof_rescale_rd),
        as.list(nof_rd),
        stringsAsFactors = FALSE
      )
    },
    psi = configs[["psi"]],
    rho = configs[["rho"]],
    nbar = configs[["nbar"]],
    tbar = configs[["tbar"]],
    L = configs[["L"]],
    rmax = configs[["rmax"]],
    SIMPLIFY = FALSE
  )
  res <- Reduce(rbind, res_list)
  rownames(res) <- seq.int(length(res_list))
  res
}


### start simulation

# choose initial seed
set.seed(1234L)
text <- "main"

# different configurations
R <- 1000L
nbar <- 100L
tbar <- c(5L, 10L, 20L, 40L, 80L)
rho <- c(0.3, 0.9)
psi <- c(0.0, 0.2, 0.4)

# generate config data with all possible combinations
configs <- expand.grid(
  psi              = psi,
  nbar             = nbar,
  tbar             = tbar,
  rho              = rho,
  stringsAsFactors = FALSE
)
configs <- configs[order(
  configs[["psi"]],
  configs[["rho"]],
  configs[["nbar"]],
  configs[["tbar"]]
), ]

# add bandwidth for bias correction
configs[configs[["tbar"]] == 5, "L"] <- 2L
configs[configs[["tbar"]] == 10, "L"] <- 3L
configs[configs[["tbar"]] == 20, "L"] <- 4L
configs[configs[["tbar"]] == 40, "L"] <- 5L
configs[configs[["tbar"]] == 80, "L"] <- 6L
configs[configs[["tbar"]] == 5, "rmax"] <- 2L
configs[configs[["tbar"]] == 10, "rmax"] <- 5L
configs[configs[["tbar"]] == 20, "rmax"] <- 10L
configs[configs[["tbar"]] == 40, "rmax"] <- 10L
configs[configs[["tbar"]] == 80, "rmax"] <- 10L

# start simulation
res_list <- pbmclapply(seq.int(R), sim_func, configs = configs, mc.cores = 30L)
res <- Reduce(rbind, res_list)

# save results
saveRDS(res, file = paste0("results/ife_simulation_", text, ".rds"), compress = "xz")
