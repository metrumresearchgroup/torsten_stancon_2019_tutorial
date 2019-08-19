rm(list = ls())
gc()
set.seed(1954)

# Adjust to your setting
modelName <- "twoCptModel"
scriptDir <- getwd()  # use Rscripts as working directory
modelDir <- file.path(scriptDir, "model")
dataDir <- file.path(scriptDir, "data")
stanDir <- "~/Desktop/Code/torsten_stancon_2019_tutorial/Torsten/cmdstan"

library(rstan)
library(plyr)
library(tidyr)
library(dplyr)
library(parallel)
rstan_options(auto_write = TRUE)

source("tools/cmdStanTools.r")

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

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

nChains <- 4
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

tempDir <- file.path(modelDir, modelName, "temp")
dir.create(tempDir)

chains <- 1:nChains
mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init) {
           tempDir <- file.path(tempDir, chain)
           dir.create(tempDir)
           inits <- init()
           with(inits, stan_rdump(ls(inits), 
                                  file = file.path(tempDir, "init.R")))
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = file.path(tempDir, "init.R"),
                    seed = sample(1:99999, 1),
                    chain = chain, refresh = 100,
                    adapt_delta = 0.95, stepsize = 0.01)
         },
         model = file.path(modelDir, modelName),
         data = file.path(dataDir, "twoCpt.data.r"),
         init = init, iter = 500, warmup = 500, thin = 1,
         mc.cores = min(nChains, detectCores()))

fit <- read_stan_csv(file.path(modelDir, modelName, 
                               paste0(modelName, chains, ".csv")))

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
