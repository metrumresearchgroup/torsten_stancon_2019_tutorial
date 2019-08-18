rm(list = ls())
gc()
set.seed(1954)

.libPaths("~/Rlib/")
library(rstan)
library(plyr)
library(tidyr)
library(dplyr)
rstan_options(auto_write = TRUE)

data <- read_rdump("model/chemical_reactions/chem.data.R")

# Read in data and output of model fit
n_chains <- 3
chains <- 1:n_chains
modelDir <- file.path("model", "chemical_reactions")
fit <- read_stan_csv(file.path(modelDir, paste0("chain", chains, ".csv")))


# Parameters of interest
pars = c("lp__", "y0_mu", "y0_1", "sigma")

# Inspect table
summary(fit, pars)[1]

# Let's do some graphical hecks.
traceplot(fit, pars = pars)
stan_dens(fit, separate_chains = TRUE, pars = pars)

#####################################################################
# Posterior predictive plots

ID <- c()
nSubjects <- data$nsub
# ID <- rep(1:nSubjects, each = data$len)

for (i in 1:nSubjects) ID <- c(ID, rep(i, data$len[i]))

data_pred <- data.frame(data$obs,data$ts, ID)
names(data_pred) <- c("obs", "time", "ID")

# prediction for furture observations in the same patient
pred <- as.data.frame(fit, pars = "obs_pred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%  # fix quantiles
  bind_cols(data_pred)

p1 <- ggplot(pred, aes(x = time, y = obs)) +
  geom_point(size = 0.1) +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  facet_wrap(~ ID)

# prediction for future observations in new patients.
pred2 <- as.data.frame(fit, pars = "obs_new_pred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p2 <- ggplot(pred2, aes(x = time, y = obs)) +
  geom_point(size = 0.1) +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
