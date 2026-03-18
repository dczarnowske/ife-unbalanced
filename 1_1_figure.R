######################################### Section 4 ################################################
# Figure "Pattern of Missing Observations"
#
####################################################################################################

# required packages
library(data.table)
library(ggplot2)
library(ggthemes)

# general settings
set.seed(1805L)
nbar <- 100L
tbar <- 20L
psi <- c(0.0, 0.2, 0.4)
n <- round(nbar / (1.0 - psi))
t <- round(tbar / (1.0 - psi))
alpha <- 2.0 * psi
n1 <- floor(alpha * n)
n2 <- n - n1
t1 <- floor(t / 2)
t2 <- t
id2 <- sample(n[[2L]])
id3 <- sample(n[[3L]])

# psi = 0.0
dt1 <- data.frame(
  i = rep(seq.int(n[[1L]]), each = t[[1L]]),
  t = rep(seq.int(t[[1L]]), n[[1L]])
)
setDT(dt1)

# psi = 0.2
s <- rep(sample(seq.int(0L, t1[[2L]]), n1[[2L]], replace = TRUE), each = t1[[2L]])
df1 <- data.frame(
  i = rep(id2[seq.int(n1[[2L]])], each = t1[[2L]]),
  t = s + rep.int(seq.int(t1[[2L]]), n1[[2L]])
)
df2 <- data.frame(
  i = rep(id2[-seq.int(n1[[2L]])], each = t2[[2L]]),
  t = rep(seq.int(t2[[2L]]), n2[[2L]])
)
dt2 <- rbind(df1, df2)
setDT(dt2)

# psi = 0.4
s <- rep(sample(seq.int(0L, t1[[3L]]), n1[[3L]], replace = TRUE), each = t1[[3L]])
df1 <- data.frame(
  i = rep(id3[seq.int(n1[[3L]])], each = t1[[3L]]),
  t = s + rep.int(seq.int(t1[[3L]]), n1[[3L]])
)
df2 <- data.frame(
  i = rep(id3[-seq.int(n1[[3L]])], each = t2[[3L]]),
  t = rep(seq.int(t2[[3L]]), n2[[3L]])
)
dt3 <- rbind(df1, df2)
setDT(dt3)

# generate plot
setkey(dt1, i)
setkey(dt2, i)
setkey(dt3, i)
dt1[, type := "psi = 0.0"]
dt2[, type := "psi = 0.2"]
dt3[, type := "psi = 0.4"]
dt <- rbind(dt1, dt2, dt3)
ggplot(dt, aes(t, i)) +
  geom_point(shape = 20) +
  facet_grid(~type) +
  theme_classic() +
  xlab("Time Period") +
  ylab("Unit")
ggsave("figures/main_figure1.pdf", width = 12, height = 5)
