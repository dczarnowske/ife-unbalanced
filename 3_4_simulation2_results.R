######################################### Online Supplement ########################################
# Results from Additional Monte Carlo Simulations
#
####################################################################################################

### import and prepare data set with results

# load packages
library(data.table)
library(xtable)

# import results
dt <- readRDS("results/ife_simulation2_supp.rds")
setDT(dt)

# helper variables
criteria <- c("ic2", "bic3", "er", "gr", "ed", "pa")
criteriast <- paste0(criteria, "_st")
configs <- c("missingdata", "psi", "errorterm", "nbar", "tbar")
errorterms <- c("homo", "homo_fattail", "xshet", "xshet_tshet")
properties <- c("bias", "ratio", "size")

# sort results
setkeyv(dt, configs)

### finite sample properties
table1 <- dt[, .(
  bias = (mean(beta_ls) - 1.0) * 100.0,
  ratio = mean(se_beta_ls) / sd(beta_ls),
  size = mean(abs(beta_ls - 1.0) / se_beta_ls > qnorm(0.975))
),
by = configs
]

# generate table 1
table1list <- lapply(paste0("pattern", seq.int(3L)), function(x) {
  tabbal <- table1[missingdata == "none"]
  tabbal[, "missingdata" := NULL]
  tabunbal <- table1[missingdata == eval(x)]
  tabunbal[, "missingdata" := NULL]
  tab <- rbind(tabbal, tabunbal)
  tab <- tab[, lapply(.SD, function(x) {
    paste0(format(round(x, 2L), digits = 2L, nsmall = 2L), collapse = " / ")
  }), .SDcols = properties, by = c("errorterm", "nbar", "tbar")]
  tab <- tab[, errorterm := NULL]
  tab
})

# print table 9-11
sink("tables/supp_table9-11.txt")
invisible(lapply(table1list, function(x) {
  explabels <- "$\\psi = 0.0 \\; / \\; \\psi = 0.2 \\; / \\; \\psi = 0.4$"
  proplabels <- c("Bias", "Ratio", "Size")
  errorlabels <- c(
    "Homoskedastic",
    "Homoskedastic with Fat Tails",
    "Heteroskedastic across Units",
    "Heteroskedastic across Units and over Time"
  )
  pos <- as.list(c(numeric(4L), seq(from = 0L, by = 6L, length.out = 4L)))
  hlines <- nrow(x)
  command <- c(
    "\\toprule\n",
    paste0("$\\overline{N}$&$\\overline{T}$&\\multicolumn{3}{c}{", explabels, "}\\\\\n"),
    "\\cmidrule(lr){3-5}\n",
    paste0("&&", paste0(proplabels, collapse = "&"), "\\\\\n"),
    paste0("\\midrule\n&&\\multicolumn{3}{c}{", errorlabels, "}\\\\\n\\cmidrule(lr){3-5}\n")
  )
  rows <- list(pos, command)
  cat("\n\n\n--------------------------- Start New Table---------------------------\n\n\n")
  print(xtable(x),
    add.to.row       = rows,
    booktabs         = TRUE,
    comment          = FALSE,
    hline.after      = hlines,
    include.rownames = FALSE,
    include.colnames = FALSE,
    only.contents    = TRUE
  )
}))
sink()

### estimating the number of factors
table2 <- dt[, lapply(.SD, mean), .SDcols = criteria, by = configs]

# generate table 2
table2list <- lapply(paste0("pattern", seq.int(3L)), function(x) {
  tabbal <- table2[missingdata == "none"]
  tabbal[, "missingdata" := NULL]
  tabunbal <- table2[missingdata == eval(x)]
  tabunbal[, "missingdata" := NULL]
  tab <- rbind(tabbal, tabunbal)
  tab <- tab[, lapply(.SD, function(x) {
    paste0(format(round(x, 2L), digits = 2L, nsmall = 2L), collapse = " / ")
  }), .SDcols = criteria, by = c("errorterm", "nbar", "tbar")]
  tab <- tab[, errorterm := NULL]
  tab
})

# print table 12-14
sink("tables/supp_table12-14.txt")
invisible(lapply(table2list, function(x) {
  explabels <- "$\\psi = 0.0 \\; / \\; \\psi = 0.2 \\; / \\; \\psi = 0.4$"
  critlabels <- c("$\\text{IC}_{2}$", "$\\text{BIC}_{3}$", "ER", "GR", "ED", "PA")
  errorlabels <- c(
    "Homoskedastic",
    "Homoskedastic with Fat Tails",
    "Heteroskedastic across Units",
    "Heteroskedastic across Units and over Time"
  )
  pos <- as.list(c(numeric(4L), seq(from = 0L, by = 6L, length.out = 4L)))
  hlines <- nrow(x)
  command <- c(
    "\\toprule\n",
    paste0("$\\overline{N}$&$\\overline{T}$&\\multicolumn{6}{c}{", explabels, "}\\\\\n"),
    "\\cmidrule(lr){3-8}\n",
    paste0("&&", paste0(critlabels, collapse = "&"), "\\\\\n"),
    paste0("\\midrule\n&&\\multicolumn{6}{c}{", errorlabels, "}\\\\\n\\cmidrule(lr){3-8}\n")
  )
  rows <- list(pos, command)
  cat("\n\n\n--------------------------- Start New Table---------------------------\n\n\n")
  print(xtable(x),
    add.to.row       = rows,
    booktabs         = TRUE,
    comment          = FALSE,
    hline.after      = hlines,
    include.rownames = FALSE,
    include.colnames = FALSE,
    only.contents    = TRUE
  )
}))
sink()
