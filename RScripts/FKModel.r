rm(list = ls())
gc()
set.seed(1954)

setwd("~/Desktop/classes/Stan-and-Torsten/Rscripts")
.libPaths("~/Rlib/")
library(rstan)
library(plyr)
library(tidyr)
library(dplyr)
rstan_options(auto_write = TRUE)

data <- read_rdump("data/FKModel.data.r")

# Specify initial conditions: check if needed.
init <- function() {
  list(CL = rlnorm(1, log(10), 0.2),
       Q = rlnorm(1, log(15), 0.2),
       VC = rlnorm(1, log(35), 0.2),
       VP = rlnorm(1, log(105), 0.2),
       ka = rlnorm(1, log(2), 0.2),
       sigma = abs(rcauchy(1, 0, 1)),
       #
       mtt = rlnorm(1, log(125), 0.2),
       circ0 = rlnorm(1, log(5), 0.05),
       gamma = rlnorm(1, log(0.17), 0.1),
       sigmaNeut = abs(rcauchy(1, 0, 1)))
}


# more control over the parameters
fit <- stan(file = "model/FKModel_solution.stan",
            data = data,
            init = init,
            chains = 3, iter = 500,
            cores = 3,
            control = list(adapt_delta = 0.95, stepsize = 0.01))


# Focus on the parameters we care about
pars = c("lp__", "CL", "Q", "VC", "VP", "ka", "sigma",
         "mtt", "circ0", "gamma", "sigmaNeut")
summary(fit, pars)[1]

pdf(file = file.path("deliv", "FKModelPlots%03d.pdf"),
    width = 6, height = 6, onefile = FALSE)

# Let's do some graphical hecks.
traceplot(fit, inc_warmup = TRUE, pars = pars)
traceplot(fit, pars = pars)
stan_dens(fit, separate_chains = TRUE, pars = pars)
pairs(fit, pars = pars)

# Let's do some predictive checks.
data_pred <- data.frame(data$cObs, data$time[data$iObsPK])
names(data_pred) <- c("cObs", "time")

predPK <- as.data.frame(fit, pars = "concentrationObsPred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p1 <- ggplot(predPK, aes(x = time, y = cObs)) +
  geom_point() +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

data_pred2 <- data.frame(data$neutObs, data$time[data$iObsPD])
names(data_pred2) <- c("neutObs", "time")

p1

predPD <- as.data.frame(fit, pars = "neutObsPred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred2)

p2 <- ggplot(predPD, aes(x = time, y = neutObs)) +
  geom_point() +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

p2

dev.off()
