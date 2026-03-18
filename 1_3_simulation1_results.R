######################################### Section 4 ################################################
# Results from Monte Carlo Simulations
#
####################################################################################################

### import and prepare data set with results

# load packages
library(data.table)
library(xtable)

# import results
dt <- readRDS("results/ife_simulation_main.rds")
setDT(dt)

# helper variables
criteria <- paste0(c("ic2", "bic3", "er", "gr", "ed", "pa"), "_rescale_em")
configs <- c("psi", "rho", "tbar")
rholabels <- c("$\\rho = 0.3$", "$\\rho = 0.9$")
properties <- c("bias", "ratio", "size")

# sort results
setkeyv(dt, configs)

### finite sample properties
table1 <- dt[, .(
  L       = first(L),
  bias    = mean((beta_em - rho) / rho) * 100.0,
  ratio   = mean(se_beta_em) / sd(beta_em),
  size    = mean(abs(beta_em - rho) / se_beta_em > qnorm(0.975))
), by = configs]

# generate table 1
tab_bal <- table1[psi == 0.0]
tab_unbal <- table1[psi != 0.0]
tab1 <- rbind(tab_bal, tab_unbal)
tab1 <- tab1[, lapply(.SD, function(x) {
  paste0(format(round(x, 3L), digits = 3L, nsmall = 3L), collapse = " / ")
}), .SDcols = properties, by = c("rho", "tbar", "L")]
tab1 <- tab1[, rho := NULL]

# print table 1
explabels <- "$\\psi = 0.0 \\; / \\; \\psi = 0.2 \\; / \\; \\psi = 0.4$"
proplabels <- c("Bias", "Ratio", "Size")
pos <- as.list(c(numeric(4L), seq(from = 0L, by = 5L, length.out = 2L)))
hlines <- nrow(tab1)
command <- c(
  "\\toprule\n",
  paste0("$\\overline{T}$&$L$&\\multicolumn{3}{c}{", explabels, "}\\\\\n"),
  "\\cmidrule(lr){3-5}\n",
  paste0("&&", paste0(proplabels, collapse = "&"), "\\\\\n"),
  paste0("\\midrule\n&&\\multicolumn{3}{c}{", rholabels, "}\\\\\n\\cmidrule(lr){3-5}\n")
)
rows <- list(pos, command)
print(
  xtable(tab1),
  file = "tables/main_table1.txt",
  add.to.row       = rows,
  booktabs         = TRUE,
  comment          = FALSE,
  hline.after      = hlines,
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents    = TRUE
)

### estimating the number of factors
table2 <- dt[, lapply(.SD, mean), .SDcols = c("rmax", criteria), by = configs]

# generate table 2
tab_bal <- table2[psi == 0.0]
tab_unbal <- table2[psi != 0.0]
tab2 <- rbind(tab_bal, tab_unbal)
tab2 <- tab2[, lapply(.SD, function(x) {
  paste0(format(round(x, 3L), digits = 3L, nsmall = 3L), collapse = " / ")
}), .SDcols = criteria, by = c("rho", "tbar", "rmax")]
tab2 <- tab2[, rho := NULL]
tab2[, rmax := as.integer(rmax)]
tab2_1 <- tab2[, c("tbar", "rmax", criteria[1:3]), with = FALSE]
setnames(tab2_1, criteria[1:3], paste0("c", 1:3))
tab2_2 <- tab2[, c("tbar", "rmax", criteria[4:6]), with = FALSE]
setnames(tab2_2, criteria[4:6], paste0("c", 1:3))
tab2 <- rbindlist(list(tab2_1, tab2_2))

# print table 2
explabels <- "$\\psi = 0.0 \\; / \\; \\psi = 0.2 \\; / \\; \\psi = 0.4$"
critlabels <- c("$\\text{IC}_{2}$", "$\\text{BIC}_{3}$", "ER", "GR", "ED", "PA")
pos <- as.list(c(numeric(5L), c(5L, 10L, 10L, 15L)))
hlines <- nrow(tab2)
command <- c(
  "\\toprule\n",
  paste0("$\\overline{T}$&$\\overline{R}$&\\multicolumn{3}{c}{", explabels, "}\\\\\n"),
  "\\cmidrule(lr){3-5}\n",
  paste0("&&", paste0(critlabels[1:3], collapse = "&"), "\\\\\n"),
  paste0("\\midrule\n&&\\multicolumn{3}{c}{", rholabels, "}\\\\\n\\cmidrule(lr){3-5}\n"),
  paste0("&&", paste0(critlabels[4:6], collapse = "&"), "\\\\\n"),
  paste0("\\midrule\n&&\\multicolumn{3}{c}{", rholabels, "}\\\\\n\\cmidrule(lr){3-5}\n")
)
rows <- list(pos, command)
print(
  xtable(tab2),
  file = "tables/main_table2.txt",
  add.to.row       = rows,
  booktabs         = TRUE,
  comment          = FALSE,
  hline.after      = hlines,
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents    = TRUE
)
