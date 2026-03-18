######################################### Appendix #################################################
# Results from Comparison of Matrix Completion Algorithms
#
####################################################################################################

### import and prepare data set with results

# load packages
library(data.table)
library(xtable)

# import results
dt <- readRDS("results/mc_comparison_app.rds")
setDT(dt)

# helper variables
configs <- "n"

# sort results
setkeyv(dt, configs)

# format numbers in table
table_formater <- function(x, digits) {
  format(round(x, digits), digits = digits, nsmall = digits, trim = TRUE)
}


### tuning parameter selection and computation time
tab1 <- dt[, .(
  nbar = first(nbar),
  nu_cv = mean(nu_cv),
  nu_pi = mean(nu_pi),
  time_rdcv_em = mean(time_rdcv_em),
  time_rdpi_em = mean(time_rdpi_em)
), by = configs][-1]
tab1[, n := as.integer(n)]
cols <- c("nu_cv", "nu_pi", "time_rdcv_em", "time_rdpi_em")
tab1[, (cols) := lapply(.SD, table_formater, 4L), .SDcols = cols]

# print table 5
pos <- as.list(numeric(3L))
hlines <- nrow(tab1)
command <- c(
  "\\toprule\n",
  "$N = T$&$\\overline{N}=\\overline{T}$&$\\nu_{\\text{CV}}$&$\\nu_{\\text{PI}}$&Time CV / Time EM&Time PI / Time EM\\\\\n",
  "\\midrule\n"
)
rows <- list(pos, command)
print(
  xtable(tab1),
  file = "tables/app_table5.txt",
  add.to.row = rows,
  booktabs = TRUE,
  comment = FALSE,
  hline.after = hlines,
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents = TRUE
)


### prediction performance
tab2_1 <- dt[, .(
  nbar      = first(nbar),
  bias_em   = mean(bias_em_all),
  bias_rdcv = mean(bias_rdcv_all),
  bias_rdpi = mean(bias_rdpi_all),
  rmse_em   = mean(rmse_em_all),
  rmse_rdcv = mean(rmse_rdcv_all),
  rmse_rdpi = mean(rmse_rdpi_all)
), by = configs]
tab2_2 <- dt[, .(
  nbar      = first(nbar),
  bias_em   = mean(bias_em_obs),
  bias_rdcv = mean(bias_rdcv_obs),
  bias_rdpi = mean(bias_rdpi_obs),
  rmse_em   = mean(rmse_em_obs),
  rmse_rdcv = mean(rmse_rdcv_obs),
  rmse_rdpi = mean(rmse_rdpi_obs)
), by = configs]
tab2_3 <- dt[, .(
  nbar      = first(nbar),
  bias_em   = mean(bias_em_miss),
  bias_rdcv = mean(bias_rdcv_miss),
  bias_rdpi = mean(bias_rdpi_miss),
  rmse_em   = mean(rmse_em_miss),
  rmse_rdcv = mean(rmse_rdcv_miss),
  rmse_rdpi = mean(rmse_rdpi_miss)
), by = configs]
tab2 <- rbindlist(list(tab2_1, tab2_2, tab2_3))
tab2[, n := as.integer(n)]
cols <- c("bias_em", "bias_rdcv", "bias_rdpi", "rmse_em", "rmse_rdcv", "rmse_rdpi")
tab2[, (cols) := lapply(.SD, table_formater, 4L), .SDcols = cols]

# print table 6
proplabels <- paste0("\\multicolumn{3}{c}{", c("Bias", "RMSE"), "}", collapse = "&")
estlabels <- rep(c("EM", "CV", "PI"), 2L)
entlabels <- c("All Entries", "Observed Entries", "Missing Entries")
pos <- as.list(c(numeric(4L), seq(from = 0L, by = 3L, length.out = 3L)))
hlines <- nrow(tab2)
command <- c(
  "\\toprule\n",
  paste0("$N = T$&$\\overline{N}=\\overline{T}$&", proplabels, "\\\\\n"),
  "\\cmidrule(lr){3-8}\n",
  paste0("&&", paste0(estlabels, collapse = "&"), "\\\\\n"),
  paste0("\\midrule\n&&\\multicolumn{6}{c}{", entlabels, "}\\\\\n\\cmidrule(lr){3-8}\n")
)
rows <- list(pos, command)
print(
  xtable(tab2),
  file = "tables/app_table6.txt",
  add.to.row = rows,
  booktabs = TRUE,
  comment = FALSE,
  hline.after = hlines,
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents = TRUE
)
