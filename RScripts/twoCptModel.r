rm(list = ls())
gc()
set.seed(1954)

library(rstan)
library(plyr)
library(tidyr)
library(dplyr)
rstan_options(auto_write = TRUE)

data <- read_rdump("data/twoCpt.data.r")

# Specify initial conditions: check if needed.
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(15), 0.2)),
       VC = exp(rnorm(1, log(35), 0.2)),
       VP = exp(rnorm(1, log(105), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       sigma = abs(rcauchy(1, 0, 1)))
}

# Run Stan
fit <- stan(file = "model/solutions/twoCptModel_solution2.stan",
            data = data, 
            init = init,
            iter = 1000, chains = 3, cores = 3)

summary(fit)[1]

# Focus on the parameters we care about
pars = c("lp__", "CL", "Q", "VC", "VP", "ka", "sigma")
summary(fit, pars)[1]

# Let's do some graphical hecks.
traceplot(fit, inc_warmup = TRUE, pars = pars)
traceplot(fit, pars = pars)
stan_dens(fit, separate_chains = TRUE, pars = pars)
pairs(fit, pars = pars)

# This makes us confident we are correctly exploring the
# posterior distribution.
#
# Let's do some predictive checks.

data_pred <- data.frame(data$cObs, data$time[-1])
names(data_pred) <- c("cObs", "time")

pred <- as.data.frame(fit, pars = "concentrationObsPred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p1 <- ggplot(pred, aes(x = time, y = cObs)) +
  geom_point() +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

