## Relative Precision of Mean Weight
##
## The purpose of this script is to estimate the relative precision of mean weight
## of lake trout in Lake Louise and Susitna Lake, given the proposed sample sizes
## in the operational plan.
##
## Sampling will consist of three events (April, June, and December), and there
## is potential for differing size selectivity between the three events.
##
## The basic outline of this script is to
## (1) load a comparable weight dataset to sample from
## (2) simulate the relative precision by resampling from this sample
##     - investigate the effect of a difference in mean weights between the two samples
##     - investigate the effect of non-equal fishing mortality (winter vs summer)
##
## NOTE: This script is based heavily on a similar method created for Crosswind
## Lake, which consisted of two events, not three.  It is likely that there is
## non-complete correspondence between the two methods.



# load packages
library(tidyverse)
library(dsftools)
library(MCMCpack)




## (1) Look in the Master Laketrout database to see if there are comparable
## samples of weights to use in sample size calculations

## look for Crosswind samples in master LT database
# morphometry <- read_csv("OP_2024/flat_data/all_lakes/lake_morphometry2.csv", skip=2)
# laketrout_all <- read_csv("OP_2024/flat_data/all_lakes/length_weight2.csv", skip=2) %>%
#   left_join(morphometry)

## current lake trout database as of March 2026
morphometry <- read_csv("https://github.com/ADFG-DSF/lake_trout_lengthweight/raw/refs/heads/main/flat_data/lake_morphometry_25_12_23.csv", skip=1)
laketrout_all <- read_csv("https://github.com/ADFG-DSF/lake_trout_lengthweight/raw/refs/heads/main/flat_data/length_weight_25_12_23.csv", skip=1) %>%
  left_join(morphometry) %>%
  filter(!is.na(Weight_g))

table(laketrout_all$LakeName)

justCrosswind <- laketrout_all %>% filter(LakeName=="Crosswind Lake")
nrow(justCrosswind)  # 269


hist(justCrosswind$Weight_g)
# understandably right-skewed

mean(justCrosswind$Weight_g)  # 2673.541
sd(justCrosswind$Weight_g)  # 1731.328
se(justCrosswind$Weight_g)  # 105.561

tapply(justCrosswind$Weight_g,
       justCrosswind$ProjectTitle,
       mean) %>% as.numeric %>% (\(x) x[1]/x[2]) # 1.268316

table(justCrosswind$ProjectTitle)
boxplot(justCrosswind$Weight_g ~ justCrosswind$ProjectTitle)

# renaming the vector of weights to use later
C_wts <- justCrosswind$Weight_g



## (2) Simulate the relative precision of estimating mean weight from stratified sample

# sample sizes
n1 <- 100
n2 <- 100
n3 <- 100

## - relative weights for each stratum
# from email:
# "My guess is that overall, harvest is highest in summer, followed by April
# (mid-late winter, and lowest in December (early winter)."
## note: these are normalized
strat_wts <- c(2, 3, 1)  # 7.4%
# strat_wts <- c(1, 1, 1)  # 7.8%
# strat_wts <- c(1, 0.5, 1.5)  # 8.8%
# strat_wts <- rev(c(1, 0.5, 1.5))  # 8.3%

## - (multiplicative) weight adjustment for each stratum
# including this will hopefully account for the additional uncertainty due to
# differences in mean weight by sampling event (April vs June vs December)
# from email:
# "I think that the length distribution will be largest during the December
# (early winter) sampling event, followed by the April (late winter), and
# smallest in June."
# wt_adj <- c(1, 0.5, 1.5)    # 7.4% with strat_wts = 2, 3, 1
wt_adj <- c(1, 1, 1)    # 7.8% with strat_wts = 2, 3, 1

nsim <- 10000
av_wt1 <- av_wt2 <- av_wt3 <- est_mn <- rep(NA, nsim)  # initializing vectors
for(i_sim in 1:nsim) {
  av_wt1[i_sim] <- mean(sample(C_wts*wt_adj[1], size=n1, replace=TRUE))
  av_wt2[i_sim] <- mean(sample(C_wts*wt_adj[2], size=n2, replace=TRUE))
  av_wt3[i_sim] <- mean(sample(C_wts*wt_adj[3], size=n3, replace=TRUE))
  est_mn[i_sim] <- sum(c(av_wt1[i_sim], av_wt2[i_sim], av_wt3[i_sim]) *
                         strat_wts/sum(strat_wts))
}
hist(est_mn)

mn_true <- sum(strat_wts*wt_adj*mean(C_wts)/sum(strat_wts))
abline(v=mn_true, col=2, lty=2, lwd=2)


# relative accuracy in this case
quantile(abs(est_mn-mn_true)/mn_true, .95)


## let's try a bunch of versions of weight adjustment..
adj_amount <- seq(0, .99, by=.01)  # multiplicative adjustment amounts to consider
quantile_vec <- rep(NA, length(adj_amount))  # corresponding vector of 95th percentile values of relative accuracy
strat_wts_list <- list(c(1, 1, 1),
                       c(1.5, 2, 1),
                       c(2, 3, 1),
                       c(3, 5, 1))
par(mfrow=c(2,2))
for(i_strat in seq_along(strat_wts_list)) {
  strat_wts <- strat_wts_list[[i_strat]]
for(i_amount in seq_along(adj_amount)) {
  # strat_wts <- c(1, 1, 1)  ## relative weights for each stratum
  # strat_wts <- c(2, 3, 1)  ## relative weights for each stratum
  # strat_wts <- c(3, 5, 1)  ## relative weights for each stratum
  wt_adj <- 1 + c(0, -1, 1)*adj_amount[i_amount]  ## (multiplicative) weight adjustment for each stratum?

  nsim <- 10000
  av_wt1 <- av_wt2 <- av_wt3 <- est_mn <- rep(NA, nsim)  # initializing vectors
  for(i_sim in 1:nsim) {
    av_wt1[i_sim] <- mean(sample(C_wts*wt_adj[1], size=n1, replace=TRUE))
    av_wt2[i_sim] <- mean(sample(C_wts*wt_adj[2], size=n2, replace=TRUE))
    av_wt3[i_sim] <- mean(sample(C_wts*wt_adj[3], size=n3, replace=TRUE))
    est_mn[i_sim] <- sum(c(av_wt1[i_sim], av_wt2[i_sim], av_wt3[i_sim]) *
                           strat_wts/sum(strat_wts))
  }
  mn_true <- sum(strat_wts*wt_adj*mean(C_wts)/sum(strat_wts))
  quantile_vec[i_amount] <- quantile(abs(est_mn-mn_true)/mn_true, .95)
}

# slight increase in uncertainty (decrease in precision) when the adjustment
# amount goes up (that is, samples have greater difference in means)
# par(mfrow=c(1,1))
plot(adj_amount, quantile_vec)
# abline(v=0.25)
# abline(h=0.077)
}

## now let's try a bunch of versions of TRUE strat_wts...
# this is done by drawing random samples from a Beta distribution centered at 50%
# Which Beta distribution you ask?  Let's plot a few.  Beta(20,20) looks conservative enough.

### UPDATED to Dirichlet

# curve(dbeta(x,50,50))
# curve(dbeta(x,20,20), add=TRUE, col=2)
# curve(dbeta(x,10,10), add=TRUE, col=3)
# curve(dbeta(x,5,5), add=TRUE, col=4)
# bp <- c(50,20,10,5)
# legend("topleft", col=1:4, lty=1,
#        legend=paste0("beta(", bp, ",", bp, ")"))

par(mfrow=c(2,2))
base_wts <- c(2, 3, 1)/2
dirmult <- c(50, 20, 10, 5)
for(i in seq_along(dirmult)) {
  dirparms <- dirmult[i]*base_wts
  dirsim <- rdirichlet(100000, dirparms)
  plot(density(dirsim[,3]), main=c("Dirichlet", paste(dirparms, collapse = ", ")), xlim=0:1, col=4)
  lines(density(dirsim[,1]), col=2)
  lines(density(dirsim[,2]), col=3)
  legend("topright", legend=c("April","June","December"), col=2:4, lwd=1)
  abline(v=base_wts/sum(base_wts), lty=3)
}

# # further interpretation of a beta(20,20)
# pbeta(.6, 20, 20) - pbeta(.4, 20, 20)  # 0.7958827
# pbeta(.65, 20, 20) - pbeta(.35, 20, 20) # 0.9464266

apply(rdirichlet(100000, 20*base_wts),
      2,
      quantile, c(0.025, 0.1, 0.5, 0.9, 0.975)) %>%
  round(3)
#        [,1]  [,2]  [,3]
# 2.5%  0.220 0.374 0.085
# 10%   0.257 0.417 0.108
# 50%   0.332 0.500 0.163
# 90%   0.413 0.582 0.230
# 97.5% 0.457 0.625 0.270


## ok let's actually do it now
betap <- 5:50 # candidate values for beta parameters (Dirichlet in this case)
base_wts <- c(2, 3, 1)/2
quantile_vec <- rep(NA, length(betap))
for(i_amount in seq_along(betap)) {
  # strat_wts <- c(1, 1, 1)  ## relative weights for each stratum
  strat_wts <- c(2, 3, 1)  ## relative weights for each stratum
  # strat_wts <- c(3, 5, 1)  ## relative weights for each stratum
  wt_adj <- 1 + c(0, -1, 1)*0.5  ## (multiplicative) weight adjustment for each stratum?

  nsim <- 10000
  av_wt1 <- av_wt2 <- av_wt3 <- est_mn <- rep(NA, nsim)  # initializing vectors
  for(i_sim in 1:nsim) {
    av_wt1[i_sim] <- mean(sample(C_wts*wt_adj[1], size=n1, replace=TRUE))
    av_wt2[i_sim] <- mean(sample(C_wts*wt_adj[2], size=n2, replace=TRUE))
    av_wt3[i_sim] <- mean(sample(C_wts*wt_adj[3], size=n3, replace=TRUE))
    est_mn[i_sim] <- sum(c(av_wt1[i_sim], av_wt2[i_sim], av_wt3[i_sim]) *
                           strat_wts/sum(strat_wts))
  }
  # # strat_wts_true <- rep(NA, 2)
  # # strat_wts_true[1] <- rbeta(1, betap[i_amount], betap[i_amount])
  # # strat_wts_true[2] <- 1-strat_wts_true[1]
  # beta1 <- rbeta(nsim, betap[i_amount], betap[i_amount])
  # beta2 <- 1-beta1
  dirsim <- rdirichlet(nsim, betap[i_amount]*base_wts)
  mn_true <- rep(NA, nsim)
  for(i in 1:nsim) {
    strat_wts_true <- dirsim[i,]
    mn_true[i] <- sum(strat_wts_true*wt_adj*mean(C_wts)/sum(strat_wts_true))
  }
  # mn_true <- sum(strat_wts_true*wt_adj*mean(C_wts)/sum(strat_wts_true))
  quantile_vec[i_amount] <- quantile(abs(est_mn-mn_true)/mn_true, .95)
}
plot(betap, quantile_vec)
## just as expected, precision gets better as Dirichlet parameters get bigger


#### FINAL WORST-CASE SCENARIO TO PUT IN THE OP PLAN:
#### - weight adjustment of 67% & 133%  (mean weight of March sample is twice as big as June)
#### - beta parameters of 20 (80% chance that March harvest is between 40% and 60% of total)

strat_wts <- c(2, 3, 1)   ## relative weights for each stratum
wt_adj <- 1 + c(0, -1, 1)*0.5  ## (multiplicative) weight adjustment for each stratum?
betap_fixed <- 20
nsim <- 10000
av_wt1 <- av_wt2 <- av_wt3 <- est_mn <- rep(NA, nsim)  # initializing vectors
for(i_sim in 1:nsim) {
  av_wt1[i_sim] <- mean(sample(C_wts*wt_adj[1], size=n1, replace=TRUE))
  av_wt2[i_sim] <- mean(sample(C_wts*wt_adj[2], size=n2, replace=TRUE))
  av_wt3[i_sim] <- mean(sample(C_wts*wt_adj[3], size=n3, replace=TRUE))
  est_mn[i_sim] <- sum(c(av_wt1[i_sim], av_wt2[i_sim], av_wt3[i_sim]) *
                         strat_wts/sum(strat_wts))
}
# strat_wts_true <- rep(NA, 2)
# strat_wts_true[1] <- rbeta(1, betap[i_amount], betap[i_amount])
# strat_wts_true[2] <- 1-strat_wts_true[1]
dirsim <- rdirichlet(nsim, betap_fixed*base_wts)
mn_true <- rep(NA, nsim)
for(i in 1:nsim) {
  strat_wts_true <- dirsim[i,]
  mn_true[i] <- sum(strat_wts_true*wt_adj*mean(C_wts)/sum(strat_wts_true))
}
# mn_true <- sum(strat_wts_true*wt_adj*mean(C_wts)/sum(strat_wts_true))

quantile(abs(est_mn-mn_true)/mn_true, .95)   # 15%
