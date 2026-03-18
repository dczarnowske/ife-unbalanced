######################################### Section 5 ################################################
# Empirical Illustration - Acemoglu et al. (2019)
#
####################################################################################################

### load packages and define functions


# required packages
library(data.table)
library(forcats)
library(ggplot2)
library(ggthemes)
library(haven)
library(InteractiveEffects)
library(pbmcapply) # Note: Unix OS required
library(xtable)

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

# format numbers in table
table_formater <- function(x, digits) {
  format(round(x, digits), digits = digits, nsmall = digits, trim = TRUE)
}

# generate data for scree plot
scree_plot_data <- function(fit, r_max) {
  # validity check
  if (!inherits(fit, "lm_ie")) stop("'fit' has to be of class lm_ie.", call. = FALSE)

  # extract required quantities from 'fit'
  n <- fit[["mf"]][["pdim"]][[1L]]
  t <- fit[["mf"]][["pdim"]][[2L]]
  nt <- n * t
  scaling <- length(fit[["mf"]][["y"]])^2 / nt
  observed <- fit[["mf"]][["observed"]]

  # compute pure factor model (impute missing values with zeros)
  if (is.null(observed)) {
    gamma <- matrix(fit[["pfm"]], n, t)
  } else {
    gamma <- matrix(0.0, n, t)
    gamma[observed] <- fit[["pfm"]]
  }

  # compute eigenvalues of the sample covariance
  if (n > t) {
    mu <- eigen(crossprod(gamma) / scaling, symmetric = TRUE, only.values = TRUE)[["values"]]
  } else {
    mu <- eigen(tcrossprod(gamma) / scaling, symmetric = TRUE, only.values = TRUE)[["values"]]
  }

  # parallel analysis based on Buja and Eyuboglu
  # compute quantile of eigenvalues of randomly permuted data
  np <- 199L
  epsilon <- 0.05
  rmu <- apply(sapply(seq.int(np), function(b, w, n, t, scaling) {
    rgamma <- apply(gamma, 2L, sample)
    if (n > t) {
      eigen(crossprod(rgamma) / scaling, symmetric = TRUE, only.values = TRUE)[["values"]]
    } else {
      eigen(tcrossprod(rgamma) / scaling, symmetric = TRUE, only.values = TRUE)[["values"]]
    }
  }, w = w, n = n, t = t, scaling = scaling), 1L, max)
  rmu <- rmu * (1.0 + epsilon)

  # generate data.frame
  nof <- seq.int(r_max)
  data.frame(nof = nof, sigma = sqrt(mu[nof]), rsigma = sqrt(rmu[nof]))
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
dt[, sum(dem)]
dt[, .N - sum(dem)]
dt[dem == 0, mean(gdp)]
dt[dem == 1, mean(gdp)]
dt[, mean(dem), by = wbcode][0.0 < V1 & V1 < 1.0, .N]
dt[, switcher := fifelse(0 < mean(dem) & mean(dem) < 1, 1, 0), by = wbcode]
dt[, ti := .N, by = wbcode]
summary(dt[ti != t])
dt[ti != t, length(unique(wbcode))]

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

# set maximum number of factor
r_max <- 5L

# estimate different specifications and store results
fit_list_ife1 <- pbmclapply(
  seq.int(r_max), lm_ie,
  formula = y ~ dem + y_lag1,
  data = pd1,
  additive_effects = "xs+time",
  nu = 0.0,
  mc.cores = r_max
)
fit_list_ife2 <- pbmclapply(
  seq.int(r_max), lm_ie,
  formula = y ~ dem + y_lag1 + y_lag2,
  data = pd2,
  additive_effects = "xs+time",
  nu = 0.0,
  mc.cores = r_max
)
fit_list_ife3 <- pbmclapply(
  seq.int(r_max), lm_ie,
  formula = y ~ dem + y_lag1 + y_lag2 + y_lag3 + y_lag4,
  data = pd3,
  additive_effects = "xs+time",
  nu = 0.0,
  mc.cores = r_max
)
save(
  fit_list_ife1, fit_list_ife2, fit_list_ife3,
  file = "results/ife_empirical_example.RData",
  compress = "xz"
)

## generate tables and figures

# figure 2 (motiviation)
fig1 <- dt
fig1[, control := fifelse(switcher == 0 & dem == 0, 1, 0)]
fig1_1 <- fig1[control == 1] 
fig1_1 <- fig1_1[, .(y_control = mean(y)), by = year]
fig1_2 <- fig1[switcher == 1]
fig1_2[, dem_lag1 := shift(dem), by = wbcode]
fig1_2 <- fig1_2[, .(wbcode, year, y, dem, dem_lag1)]
fig1_2[, transition := fifelse(dem - dem_lag1 == 1, 1, 0)]
fig1_2[, transition := shift(transition), by = wbcode]
fig1_2 <- na.omit(fig1_2)
fig1_2[, treat_clean := fifelse(sum(transition) == 1, 1, 0) , by = wbcode]
fig1_2 <- fig1_2[treat_clean == 1]
fig1_2[, year_treat := max(transition * year), by = wbcode]
fig1_2[, time_to_treat := year - year_treat]
fig1 <- merge(fig1_2, fig1_1, by = "year", all.x = TRUE)
fig1[, y_diff := y - y_control]
fig1[, n := .N, by = time_to_treat]
fig1 <- fig1[, .(y_diff = mean(y_diff)), by = time_to_treat]
setkey(fig1, time_to_treat)
y_diff_0 <- fig1[time_to_treat == 0, y_diff]
fig1[, y_diff := y_diff - y_diff_0]
fig1 <- fig1[time_to_treat >= - 15 & time_to_treat <= 15]
ggplot(fig1, aes(x = time_to_treat, y = y_diff)) +
  geom_line(alpha = 0.75) +
  geom_point() +
  labs(
    x        = "Years around Democratization",
    y        = "Change in GDP per capita (in log points)"
  ) +
  scale_x_continuous(breaks = seq.int(fig1[, min(time_to_treat)], fig1[, max(time_to_treat)])) +
  scale_y_continuous(n.breaks = 8) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave("figures/main_figure2.pdf", width = 12, height = 5)

# bandwidth parameter (rule-of-thumbs)
L_rot <- round(0.75 * sqrt(t_bar))

# table 3 (number of factors)
set.seed(1234L)
fit_ife1 <- fit_list_ife1[[r_max]]
fit_ife2 <- fit_list_ife2[[r_max]]
fit_ife3 <- fit_list_ife3[[r_max]]
tab1 <- rbind(
  number_of_factors(fit_ife1),
  number_of_factors(fit_ife2),
  number_of_factors(fit_ife3)
)[, c(5L, seq.int(7L, 11L))]
tab1 <- as.data.frame(tab1)
tab1 <- cbind(spec = c("$p = 1$", "$p = 2$", "$p = 4$"), tab1)
colnames(tab1) <- c(
  "Specification",
  "$\\text{IC}_{2}$", "$\\text{BIC}_{3}$", "ER", "GR", "ED", "PA"
)
print(
  xtable(tab1),
  file = "tables/main_table3.txt",
  booktabs = TRUE,
  comment = FALSE,
  include.rownames = FALSE,
  only.contents = TRUE,
  sanitize.text.function = function(x) {
    x
  }
)

# figure 3 (scree plots)
limit <- 2 * r_max
fig2 <- rbind(
  scree_plot_data(fit_ife1, limit),
  scree_plot_data(fit_ife2, limit),
  scree_plot_data(fit_ife3, limit)
)
fig2[["specification"]] <- rep(paste("Specification", seq.int(3L)), each = limit)
setDT(fig2)
fig2 <- melt(fig2, id.vars = c("nof", "specification"))
fig2[, specification := factor(specification)]
nlvls <- c("observed" = "sigma", "randomized" = "rsigma")
fig2[, variable := fct_recode(variable, !!!nlvls)]
ggplot(fig2, aes(x = nof, y = value, shape = variable, linetype = variable)) +
  facet_grid(. ~ specification) +
  geom_line(alpha = 0.75) +
  geom_point() +
  labs(
    x        = "Factor Number",
    y        = "Singular Value",
    shape    = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(breaks = 1:limit) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave("figures/main_figure3.pdf", width = 12, height = 5)


# table 4 (bias-corrected results)
fit_fe1_bc <- bias_correction(fit_fe1, L_rot)
ab <- c(0.959, 0.477, 0.946, 0.009, 17.608, 10.609)
hhk <- c(0.781, 0.455, 0.938, 0.011, 12.644, 8.282)
fit_ife1_bc1 <- bias_correction(fit_list_ife1[[1L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife1_bc2 <- bias_correction(fit_list_ife1[[2L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife1_bc3 <- bias_correction(fit_list_ife1[[3L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
tab21 <- cbind(
  as.vector(t(effects(fit_fe1_bc, TRUE))),
  ab,
  hhk,
  as.vector(t(effects(fit_ife1_bc1, TRUE))),
  as.vector(t(effects(fit_ife1_bc2, TRUE))),
  as.vector(t(effects(fit_ife1_bc3, TRUE)))
)
fit_fe2_bc <- bias_correction(fit_fe2, L_rot)
ab <- c(0.797, 0.417, 0.946, 0.009, 14.882, 9.152)
hhk <- c(0.582, 0.387, 0.941, 0.010, 9.929, 7.258)
fit_ife2_bc1 <- bias_correction(fit_list_ife2[[1L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife2_bc2 <- bias_correction(fit_list_ife2[[2L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife2_bc3 <- bias_correction(fit_list_ife2[[3L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
tab22 <- cbind(
  as.vector(t(effects(fit_fe2_bc, TRUE))),
  ab,
  hhk,
  as.vector(t(effects(fit_ife2_bc1, TRUE))),
  as.vector(t(effects(fit_ife2_bc2, TRUE))),
  as.vector(t(effects(fit_ife2_bc3, TRUE)))
)
fit_fe3_bc <- bias_correction(fit_fe3, L_rot)
ab <- c(0.875, 0.374, 0.947, 0.009, 16.448, 8.436)
hhk <- c(1.178, 0.370, 0.953, 0.009, 25.032, 10.581)
fit_ife3_bc1 <- bias_correction(fit_list_ife3[[1L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife3_bc2 <- bias_correction(fit_list_ife3[[2L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
fit_ife3_bc3 <- bias_correction(fit_list_ife3[[3L]], L = L_rot, xs_het = TRUE, ts_het = TRUE)
tab23 <- cbind(
  as.vector(t(effects(fit_fe3_bc, TRUE))),
  ab,
  hhk,
  as.vector(t(effects(fit_ife3_bc1, TRUE))),
  as.vector(t(effects(fit_ife3_bc2, TRUE))),
  as.vector(t(effects(fit_ife3_bc3, TRUE)))
)

# bind panels and format results
effectlabels <- c(
  "Democracy", "",
  "Persistence of", "\\quad GDP process",
  "Long-run effect", "\\quad of democracy"
)
tab2 <- rbind(tab21, tab22, tab23)
tab2 <- table_formater(tab2, 3L)
tab2[seq.int(2L, nrow(tab2), by = 2L), ] <-
  apply(tab2[seq.int(2L, nrow(tab2), by = 2L), ], 1:2, function(x) {
    paste0("(", x, ")")
  })
tab2 <- cbind(rep(effectlabels, 3L), tab2)

# print table 4
estlabels <- "&FE&AB&HHK&\\multicolumn{3}{c}{IFE}\\\\\n"
speclabels <- paste0("Specification ", 1:3, " - $p = ", c(1:2, 4L), "$")
factorlabels <- paste0("$R = ", 1:3, "$")
pos <- as.list(c(numeric(4L), seq(from = 0L, by = 6L, length.out = 3L)))
hlines <- nrow(tab2)
command <- c(
  "\\toprule\n",
  estlabels,
  "\\cmidrule(lr){5-7}\n",
  paste0("&&&&", paste0(factorlabels, collapse = "&"), "\\\\\n"),
  paste0("\\midrule\n&\\multicolumn{6}{c}{", speclabels, "}\\\\\n\\cmidrule(lr){2-7}\n")
)
rows <- list(pos, command)
print(
  xtable(tab2),
  file = "tables/main_table4.txt",
  add.to.row = rows,
  booktabs = TRUE,
  comment = FALSE,
  hline.after = hlines,
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents = TRUE,
  sanitize.text.function = function(x) {
    x
  }
)
