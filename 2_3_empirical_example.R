######################################### Appendix #################################################
# Empirical Illustration - Acemoglu et al. (2019)
#
####################################################################################################

### Load packages and define functions

# required packages
library(data.table)
library(forcats)
library(ggplot2)
library(ggthemes)
library(haven)
library(InteractiveEffects)

# compute short-run, persistence, and long-run effect (and standard errors)
effects <- function(mod, se = FALSE) {
  beta <- coef(mod)
  eff <- c(
    beta[[1L]],
    sum(beta[-1L]),
    beta[[1L]] / (1.0 - sum(beta[-1L]))
  )
  if (se) {
    vcov <- vcov(mod)
    jacpers <- c(0.0, rep(1.0, length(beta) - 1L))
    jaclong <- c(
      1.0 / (1.0 - sum(beta[-1L])),
      rep(beta[[1L]] / (1.0 - sum(beta[-1L]))^2, length(beta) - 1L)
    )
    stderr <- sqrt(c(
      vcov(mod)[1L, 1L],
      as.double(jacpers %*% vcov %*% jacpers),
      as.double(jaclong %*% vcov %*% jaclong)
    ))
    eff <- cbind(eff, stderr)
  }
  eff
}

### import and prepare data

# import data
dt <- read_dta("DDCGdata_final.dta")
setDT(dt)

# generate lagged outcome variables
setkey(dt, wbcode, year)
dt[, y_lag1 := shift(y, 1L), by = wbcode]
dt[, y_lag2 := shift(y, 2L), by = wbcode]
dt[, y_lag3 := shift(y, 3L), by = wbcode]
dt[, y_lag4 := shift(y, 4L), by = wbcode]

# only keep complete observations
dt <- dt[is.finite(y) & is.finite(dem)]
dt <- dt[, .(y, y_lag1, y_lag2, y_lag3, y_lag4, dem, wbcode, year)]
dt[, gdp := exp(y / 100.0)]

# compute some summary statistics
nt <- nrow(dt)
n <- dt[, length(unique(wbcode))]
t <- dt[, length(unique(year))]
n_bar <- ceiling(nt / t)
t_bar <- ceiling(nt / n)

# split data sets and drop initial periods
dt1 <- dt[, .(y, y_lag1, dem, wbcode, year)]
dt1 <- na.omit(dt1)
pd1 <- panel_data(
  data = dt1,
  xs_id = "wbcode",
  time_id = "year"
)
dt2 <- dt[, .(y, y_lag1, y_lag2, dem, wbcode, year)]
dt2 <- na.omit(dt2)
pd2 <- panel_data(
  data = dt2,
  xs_id = "wbcode",
  time_id = "year"
)
dt3 <- dt[, .(y, y_lag1, y_lag2, y_lag3, y_lag4, dem, wbcode, year)]
dt3 <- na.omit(dt3)
pd3 <- panel_data(
  data = dt3,
  xs_id = "wbcode",
  time_id = "year"
)
rm(dt1, dt2, dt3)

### estimation

## additive fixed effects

# within estimator
fit_fe1 <- lm_fe(
  formula = y ~ dem + y_lag1,
  data = pd1,
  additive_effects = "xs+time"
)
fit_fe2 <- lm_fe(
  formula = y ~ dem + y_lag1 + y_lag2,
  data = pd2,
  additive_effects = "xs+time"
)
fit_fe3 <- lm_fe(
  formula = y ~ dem + y_lag1 + y_lag2 + y_lag3 + y_lag4,
  data = pd3,
  additive_effects = "xs+time"
)

## interactive fixed effects
load("results/ife_empirical_example.RData")

## generate figures

# upper bound bandwidth parameter
L_max <- 8L

# bias-corrected results
fit_fe1_bc_bw <- Reduce(rbind, lapply(lapply(seq.int(L_max), bias_correction, fit = fit_fe1), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 1",
    estimator = "FE",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_fe1_bc_bw[["L"]] <- seq.int(L_max)
fit_ife1_bc_bw1 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife1[[1L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 1",
    estimator = "R = 1",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife1_bc_bw1[["L"]] <- seq.int(L_max)
fit_ife1_bc_bw2 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife1[[2L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 1",
    estimator = "R = 2",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife1_bc_bw2[["L"]] <- seq.int(L_max)
fit_ife1_bc_bw3 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife1[[3L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 1",
    estimator = "R = 3",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife1_bc_bw3[["L"]] <- seq.int(L_max)
tab31 <- rbind(fit_fe1_bc_bw, fit_ife1_bc_bw1, fit_ife1_bc_bw2, fit_ife1_bc_bw3)
fit_fe2_bc_bw <- Reduce(rbind, lapply(lapply(seq.int(L_max), bias_correction, fit = fit_fe2), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 2",
    estimator = "FE",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_fe2_bc_bw[["L"]] <- seq.int(L_max)
fit_ife2_bc_bw1 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife2[[1L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 2",
    estimator = "R = 1",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife2_bc_bw1[["L"]] <- seq.int(L_max)
fit_ife2_bc_bw2 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife2[[2L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 2",
    estimator = "R = 2",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife2_bc_bw2[["L"]] <- seq.int(L_max)
fit_ife2_bc_bw3 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife2[[3L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 2",
    estimator = "R = 3",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife2_bc_bw3[["L"]] <- seq.int(L_max)
tab32 <- rbind(fit_fe2_bc_bw, fit_ife2_bc_bw1, fit_ife2_bc_bw2, fit_ife2_bc_bw3)
fit_fe3_bc_bw <- Reduce(rbind, lapply(lapply(seq.int(L_max), bias_correction, fit = fit_fe3), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 3",
    estimator = "FE",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_fe3_bc_bw[["L"]] <- seq.int(L_max)
fit_ife3_bc_bw1 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife2[[1L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 3",
    estimator = "R = 1",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife3_bc_bw1[["L"]] <- seq.int(L_max)
fit_ife3_bc_bw2 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife3[[2L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 3",
    estimator = "R = 2",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife3_bc_bw2[["L"]] <- seq.int(L_max)
fit_ife3_bc_bw3 <- Reduce(rbind, lapply(lapply(
  seq.int(L_max),
  bias_correction,
  fit = fit_list_ife3[[3L]],
  xs_het = TRUE,
  ts_het = TRUE
), function(x) {
  res <- effects(x)
  data.frame(
    specification = "Specification 3",
    estimator = "R = 3",
    beta = res[[1L]],
    gamma = res[[2L]],
    phi = res[[3L]]
  )
}))
fit_ife3_bc_bw3[["L"]] <- seq.int(L_max)
tab33 <- rbind(fit_fe3_bc_bw, fit_ife3_bc_bw1, fit_ife3_bc_bw2, fit_ife3_bc_bw3)
tab3 <- rbind(tab31, tab32, tab33)
setDT(tab3)
setkeyv(tab3, c("specification", "estimator"))
tab3[, specification := factor(specification)]
tab3[, estimator := factor(estimator, levels = c("FE", paste0("R = ", 1:3)))]
tab3 <- melt(tab3, id.vars = c("specification", "estimator", "L"), variable.name = "effect")
setattr(tab3[["effect"]], "levels", c("short-run", "persistence", "long-run"))

# generate figure 4
ggplot(tab3, aes(x = L, y = value, shape = estimator, linetype = estimator)) +
  facet_wrap(specification ~ effect, scales = "free") +
  geom_point() +
  geom_line(alpha = 0.75) +
  labs(
    x        = "Bandwidth Parameter - L",
    y        = "Estimate",
    shape    = "Estimator",
    linetype = "Estimator"
  ) +
  scale_x_continuous(n.breaks = L_max) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave("figures/app_figure4.pdf", width = 12, height = 12)


# bandwidth parameter (rule-of-thumbs)
L_rot <- round(0.75 * sqrt(t_bar))

# bias-corrected results
fit_list_ife1_bc <- lapply(
  fit_list_ife1, bias_correction,
  L = L_rot, xs_het = TRUE, ts_het = TRUE
)
tab41 <- as.data.frame(sapply(fit_list_ife1_bc, function(x) {
  as.vector(t(effects(x)))
}))
tab41[["specification"]] <- "Specification 1"
tab41[["effect"]] <- c("short-run", "persistence", "long-run")
fit_list_ife2_bc <- lapply(
  fit_list_ife2, bias_correction,
  L = L_rot, xs_het = TRUE, ts_het = TRUE
)
tab42 <- as.data.frame(sapply(fit_list_ife2_bc, function(x) {
  as.vector(t(effects(x)))
}))
tab42[["specification"]] <- "Specification 2"
tab42[["effect"]] <- c("short-run", "persistence", "long-run")
fit_list_ife3_bc <- lapply(
  fit_list_ife3, bias_correction,
  L = L_rot, xs_het = TRUE, ts_het = TRUE
)
tab43 <- as.data.frame(sapply(fit_list_ife3_bc, function(x) {
  as.vector(t(effects(x)))
}))
tab43[["specification"]] <- "Specification 3"
tab43[["effect"]] <- c("short-run", "persistence", "long-run")
tab4 <- rbind(tab41, tab42, tab43)
colnames(tab4) <- c(as.character(seq.int(length(fit_list_ife1))), "specification", "effect")
setDT(tab4)
setkeyv(tab4, c("specification", "effect"))
tab4[, specification := factor(specification)]
tab4[, effect := factor(effect, levels = c("short-run", "persistence", "long-run"))]
tab4 <- melt(tab4, id.vars = c("specification", "effect"), variable.name = "R", variable.factor = FALSE)
tab4[, R := as.integer(R)]

# generate figure 5
ggplot(tab4, aes(x = R, y = value, shape = specification, linetype = specification)) +
  facet_wrap(. ~ effect, scales = "free_y") +
  geom_point() +
  geom_line(alpha = 0.75) +
  labs(
    x        = "Number of Factors - R",
    y        = "Estimate",
    shape    = NULL,
    linetype = NULL
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave("figures/app_figure5.pdf", width = 12, height = 5)
