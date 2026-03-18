######################################### Online Supplement ########################################
# Figure "Patterns of Randomly Missing Observations"
#
####################################################################################################

# general settings
library(data.table)
library(ggplot2)
library(ggthemes)
set.seed(1805L)
phi <- 0.4
n <- 30L
n1 <- n * phi
n2 <- n * (1.0 - phi)
t <- 30L
t1 <- t / 2L
t2 <- t
id <- sample(n)

# pattern 1
df <- data.frame(
  i = rep(id[seq.int(n)], each = t),
  t = rep.int(seq.int(t), n)
)
dt1 <- df[sample(n * t, n * t * (1.0 - phi)), ]
setDT(dt1)

# pattern 2
df1 <- data.frame(
  i = rep(id[seq.int(n1)], each = t1),
  t = rep.int(seq.int(t1), n1)
)
df2 <- data.frame(
  i = rep(id[seq.int(to = n1 + n2, length.out = n2)], each = t2),
  t = rep.int(seq.int(t2), n2)
)
dt2 <- rbind(df1, df2)
setDT(dt2)

# pattern 3
s <- rep(sample(seq.int(0L, t1), n1, replace = TRUE), each = t1)
df1 <- data.frame(
  i = rep(id[seq.int(n1)], each = t1),
  t = s + rep.int(seq.int(t1), n1)
)
dt3 <- rbind(df1, df2)
setDT(dt3)

# generate figure 6
setkey(dt1, i)
setkey(dt2, i)
setkey(dt3, i)
dt1[, type := "Pattern 1"]
dt2[, type := "Pattern 2"]
dt3[, type := "Pattern 3"]
dt <- rbind(dt1, dt2, dt3)
ggplot(dt, aes(t, i)) +
  geom_point(shape = 20) +
  facet_grid(~type) +
  theme_classic() +
  xlab("Time Period") +
  ylab("Unit")
ggsave("figures/supp_figure6.pdf", width = 12, height = 5)
