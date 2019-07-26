rm(list = ls())
gc()

.libPaths("~/Rlib/")

library(rstan)  # should load version 2.18.2
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)

rstan_options(auto_write = TRUE)
set.seed(1954)

## read data
data <- read_rdump("data/linear.data.r")

## create initial estimates
init <- function() {
  list(sigma = rgamma(1, 1),
       beta = rnorm(1, mean = 1, sd = 1))
}

## run Stan
fit <- stan(file = "model/linear.stan",
            data = data,
            init = init)

# Remark: the result is stored in the fit object.
summary(fit)[1]

# Let's start with some graphical checks to get a 
# sense of what happened.
pars = c("beta", "sigma", "lp__")
summary(fit, pars)[1]
traceplot(fit, inc_warmup = TRUE, pars = pars)
traceplot(fit, pars = pars)

stan_dens(fit, separate_chains = TRUE, pars = pars)
pairs(fit, pars = pars)
# pairs(fit, pars = c("beta", "sigma", "lp__"))

# This makes us confident we are correctly exploring
# the posterior distribution.

# Let us now do some predictive checks.
# (Add generated quantities block)

data_pred <- data.frame(data$x, data$y)
names(data_pred) <- c("x", "y")


pred <- as.data.frame(fit, pars = "y_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p1 <- ggplot(pred, aes(x = x, y = y))
p1 <- p1 + geom_point()
p1 + geom_line(aes(x = x, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

# Clearly, our model is not accurately describing our data.
# How can we improve our model?

fit2 <- stan(file = "model/linear2.stan",
             data = data,
             init = init)

# Remark: the result is stored in the fit object.
summary(fit2)[1]

# Let's start with some graphical checks to get a 
# sense of what happened.
pars = c("beta", "beta_0", "sigma", "lp__")
traceplot(fit2, pars = pars)

stan_dens(fit2, separate_chains = TRUE, pars = pars)

# This makes us confident we are correctly exploring
# the posterior distribution.

# Let us now do some predictive checks.
# (Add generated quantities block)

pred <- as.data.frame(fit2, pars = "y_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p2 <- ggplot(pred, aes(x = x, y = y))
p2 <- p2 + geom_point()
p2 + geom_line(aes(x = x, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
