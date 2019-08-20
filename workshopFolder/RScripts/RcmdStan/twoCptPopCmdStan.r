rm(list = ls())
gc()
set.seed(1954)

# modelName <- "twoCptPopModel"
modelName <- "twoCptPopulationModel2"
scriptDir <- getwd()  # use Rscripts as working directory
modelDir <- file.path(scriptDir, "model", "solutions")  # FIX ME -- remove solutions
dataDir <- file.path(scriptDir, "data")
stanDir <- "~/Desktop/Code/torsten_stancon_2019_tutorial/Torsten/cmdstan"

.libPaths("~/Rlib/")
library(rstan)
library(plyr)
library(tidyr)
library(dplyr)
source("tools/cmdStanTools.r")
rstan_options(auto_write = TRUE)

data <- read_rdump("data/twoCptPop.data.r")

nIIV <- data$nIIV
nSubjects <- data$nSubjects

# Specify initial conditions: check if needed.
init <- function() {
  list(CL_pop = exp(rnorm(1, log(10), 0.2)),
       Q_pop = exp(rnorm(1, log(15), 0.2)),
       VC_pop = exp(rnorm(1, log(35), 0.2)),
       VP_pop = exp(rnorm(1, log(105), 0.2)),
       ka_pop = exp(rnorm(1, log(2), 0.2)),
       sigma = abs(rcauchy(1, 0, 1)),
       omega = abs(runif(nIIV, 0, 1)),
       # for non-centered param
       alpha = matrix(rnorm(nSubjects * nIIV, 0, 1),
                      nrow = nSubjects, ncol = nIIV))
}

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

nChains <- 3
nSamples <- 100
nWarmUp <- 100
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
                    chain = chain, refresh = 10,
                    adapt_delta = 0.8, stepsize = 0.01)
         },
         model = file.path(modelDir, modelName),
         data = file.path(dataDir, "twoCptPop.data.r"),
         init = init, iter = nSamples, warmup = nWarmUp, thin = 1,
         mc.cores = min(nChains, detectCores()))

fit <- read_stan_csv(file.path(modelDir, modelName, 
                               paste0(modelName, chains, ".csv")))

# Run Stan
# fit <- stan(file = "model/solutions/twoCptPopulationModel.stan",
#             data = data, init = init,
#             iter = 100, chains = 3, cores = 3)  #,
#             # control = list(adapt_delta = 0.95))

# Focus on the parameters we care about
pars = c("lp__", "CL_pop", "Q_pop", "VC_pop", "VP_pop", "ka_pop",
         "sigma", "omega")
pars_select = c("omega[1]", "log_theta[1, 1]", "log_theta[2, 1]", "log_theta[3, 1]")
summary(fit, pars)[1]

# Let's do some graphical hecks.
traceplot(fit, pars = pars)
stan_dens(fit, separate_chains = TRUE, pars = pars)

# Look at parameters which may cause an issue.
# pairs(fit, pars = c("omega[1]", "log_theta[1, 1]", "log_theta[2, 1]", "log_theta[3, 1]"))
# pairs(fit, pars = c("omega[2]", "log_theta[1, 2]", "log_theta[2, 2]", "log_theta[3, 2]"))
# pairs(fit, pars = c("omega[3]", "log_theta[1, 3]", "log_theta[2, 3]", "log_theta[3, 3]"))
# pairs(fit, pars = c("omega[4]", "log_theta[1, 4]", "log_theta[2, 4]", "log_theta[3, 4]"))
# pairs(fit, pars = c("omega[5]", "log_theta[1, 5]", "log_theta[2, 5]", "log_theta[3, 5]"))


#####################################################################
# Posterior predictive plots

ID <- c()
for (i in 1:data$nSubjects) {
  # if (i == data$nSubjects) data$start[i + 1] <- length(data$time) + 1
  newID <- rep(i, data$end[i] - data$start[i] + 1)
  ID <- c(ID, newID)
}
data_pred <- data.frame(data$cObs, data$time[data$iObs], ID[data$iObs])
names(data_pred) <- c("cObs", "time", "ID")

# prediction for furture observations in the same patient
pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p1 <- ggplot(pred, aes(x = time, y = cObs)) +
  geom_point(size = 0.1) +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  facet_wrap(~ ID)

# prediction for future observations in new patients.
pred2 <- as.data.frame(fit, pars = "cObsNewPred") %>%
  gather(factor_key =  TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

p2 <- ggplot(pred2, aes(x = time, y = cObs)) +
  geom_point(size = 0.1) +
  geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
