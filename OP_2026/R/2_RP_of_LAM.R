## Relative Precision of Lake Area Model

## The purpose of this script is to simulate the estimation (prediction) variance
## associated with applying the Lake Area Model to a new lake.

## First, the residual standard deviation of the LAM reported by Evans, et al. was
## estimated by trial and error (trying values until the model R^2 matched what
## was reported)

# load packages
library(tidyverse)
library(dsftools)
library(MCMCpack)


## (3) approximating the variance within the Lake Area Model
# newlake_la <- 3717  # Crosswind

lakenames <- c("Louise", "Susitna")
lakeareas <- c(5913, 3635)
# ilake <- 1
ilake <- 2

newlake_la <- lakeareas[ilake]

## all lake areas reported in Evans et al (LA model)
la <- c(3582.1,97.8,71.7,633.8,126.3,2014.6,195.1,2226.3,3798.5,39.1,87.4,124,169,365.2,196.2,1174.5,6479.2,172.2,204,101.9,195.1,42.9,67.7,62.4,28.6,23.5,716.4,36.2,143.3,189.3,244.1,291.4,221,278.2,167.9,144.1,1388.5,2519.1,115.7,22.4,167.3,1730.5,137.6,90.3,441.9,483.2,93.8,189.7,26.1,26.4,402.3,1061.5,165.9,80.9,2079.9,975.4,229.5,152.8,1002.7,234.7,945.8,197.1,615.1,65.6,711.5,113.4,336.1,256,951,228.7,194.8,317.6,417.9,87.4,271.7,226.4,74.9,3060.3,93.6,162,995.8,81.6,1012.9,543.2,178.3,289.8,73.7,224.8,1788,46.1,42.1,168.8,158,155.6,113.6,2197.5,817.9,1172,2821,102,1374.6,212.5,36.6,181.7,409.1,1374.6,639.5,178.4,85,2427.6,744.9,52.9,69,21.9,325.6,344.3,229.6,25.3,56.7,95.5,63.3,200.3,31.5,34.5,738.2,489.3,767,16.3,145.6,157.8,5310,324.9,175.6,230.2,1066.3,183.1,35.6,620.3,151,157.8,263.8,507.5,654.5,76.6,289.7,137.6,25,12215.3,383.4,83.1,50.1,115.7,97.2,5154.2,305.6,1077,562.8,754.8,1077,864.8,1197.1,75.8,1320.6,132.8,56.7,136.8,45.2,560.7,210.5,140.3,383.9,287.6,105.2,421.2,6374.4,2071.1,192.2,29.7,171.9,267.5,52.9,252,118.3,246.1,25.1,1256.1,51.7,607.5,96.7,46.9,384.4,118.3,234.7,63.7,48.3,78.1,150,111.7,196,1123,167.9,27.2,157.6,137.2,336.8,1363,44.6,199.1,1955.5,574.8,37.9,363.3,109.4)
hist(la)
abline(v=newlake_la)
hist(log10(la))
abline(v=log10(newlake_la))
length(la)

## these could POSSIBLY have been used in the LA model
whichreg <- c(4,6,8,9,11,10,13,14,17,18,27,28,35,37,38,52,55,56,59,61,63,64,66,69,73,75,81,83,84,94,97,99,101,102,105,108,110,120,123,138,141,143,147,148,152,154,175,182,183,188,198,201,205,210,212)
lareg <- la[whichreg]
hist(log10(lareg))
abline(v=log10(newlake_la))


## This was a trial and error routine to figure out what the error standard
## deviation would have been (arriving at epssd = 0.305)
## The reason lakes are sampled without replacement is because I didn't know
## which 43 of the 50-some lakes were used to inform the LAM
epssd <- 0.305
n <- 10000
b0 <- b1 <- r2 <- predvar <- rep(NA,n)

## would really like to visualize this relationship
all_la <- all_log10_la <- all_log10_yp <- all_yp <- matrix(nrow=n, ncol=43)

for(i in 1:n) {
  x <- sample(lareg, 43)
  # x <- sample(la, 43)
  log10x <- log10(x)
  log10YPpred <- 0.6 + 0.72*log10(x)
  eps <- rnorm(43, 0, epssd)
  log10YPsim <- log10YPpred + eps
  mod <- lm(log10YPsim~log10x)
  b0[i] <- unname(mod$coefficients[1])
  b1[i] <- unname(mod$coefficients[2])
  r2[i] <- summary(mod)$r.squared
  predvar[i] <- epssd^2 + predict(mod,newdata=data.frame(log10x=log10(newlake_la)),se.fit=T)$se.fit^2

  # storing stuff to plot
  all_la[i,] <- x
  all_log10_la[i,] <- log10(x)
  all_log10_yp[i,] <- log10YPsim
  all_yp[i,] <- 10^log10YPsim
}

# plotting regression coefficients and R^2 values, to verify that they match
# the values reported by Evans, et al.
par(mfrow=c(2,2))
hist(b0)
abline(v=.6)
abline(v=median(b0), lty=2, col=2, lwd=2)
hist(b1)
abline(v=.72)
abline(v=median(b1), lty=2, col=2, lwd=2)
hist(r2)
abline(v=.69)
abline(v=median(r2), lty=2, col=2, lwd=2)
hist(predvar)
abline(v=median(predvar), lty=2, col=2, lwd=2)


#### this is the best estimate of the prediction variance for a NEW LAKE
predvarsim <- median(predvar)   # [1] 0.09958763

## vector of simulated YP of new lake, on the log scale
predlog10YP <- rnorm(n, 0.6 + 0.72*log10(newlake_la), sqrt(predvarsim))
hist(predlog10YP)

# quick function to return relative precision from a vector of simulations
rp <- function(vec, trueval, q=0.95) {
  quantile(abs(vec-trueval)/trueval, q)
}

rp(vec=predlog10YP, trueval=median(predlog10YP))  # 20%

## back to natural scale
predYP <- 10^predlog10YP
hist(predYP)

rp(vec=predYP, trueval=median(predYP))  # 240%, yikes

quantile(predYP, c(0.025,0.975))/median(predYP)
#        2.5%     97.5%
#   0.2382935 4.1333496




# visualizing by doing a prediction interval on the log scale and backtransforming

# CI on log scale & natural scale
CIlog <- 0.6 + 0.72*log10(newlake_la) + c(-1,1)*qnorm(0.95)*sqrt(predvarsim)
CInat <- 10^CIlog


par(mfrow=c(1,2))
plot(as.numeric(all_log10_la)[1:200], as.numeric(all_log10_yp)[1:200],
     xlab="log10(lake area)", ylab="log10(YP)",
     main=paste("Pred Intvl for", lakenames[ilake], "(log scale)"))
for(i in 1:100) abline(b0[i], b1[i], col=adjustcolor(4,alpha.f=.2))
abline(.6,.72, lwd=2)
# abline(v=log10(newlake_la))
lines(x=rep(log10(newlake_la),2), y=CIlog, lwd=4, col=4, lend=2)

plot(as.numeric(all_la)[1:200], as.numeric(all_yp)[1:200],
     xlab="lake area", ylab="YP",
     main=paste("Pred Intvl for", lakenames[ilake], "(log scale)"))
for(i in 1:100) curve(10^(b0[i]+b1[i]*log10(x)), add=TRUE, col=adjustcolor(4,alpha.f=.2))
curve(10^(.6+.72*log10(x)), add=TRUE, lwd=2)
# abline(v=newlake_la)
lines(x=rep(newlake_la,2), y=CInat, lwd=4, col=4, lend=2)



## (4) Propegating the LAM variance to RP of YP in number of fish

# # sampling (actually resampling) from weights at Fielding Lake
# laketrout_all <- read_csv("OP_2024/flat_data/all_lakes/length_weight2.csv", skip=2)
#
# justFielding <- laketrout_all %>% filter(LakeName=="Fielding Lake" & !is.na(Weight_g))
# # nrow(justFielding)  # 625
# F_wts <- justFielding$Weight_g

## current lake trout database as of March 2026
morphometry <- read_csv("https://github.com/ADFG-DSF/lake_trout_lengthweight/raw/refs/heads/main/flat_data/lake_morphometry_25_12_23.csv", skip=1)
laketrout_all <- read_csv("https://github.com/ADFG-DSF/lake_trout_lengthweight/raw/refs/heads/main/flat_data/length_weight_25_12_23.csv", skip=1) %>%
  left_join(morphometry) %>%
  filter(!is.na(Weight_g))

justCrosswind <- laketrout_all %>% filter(LakeName=="Crosswind Lake")

C_wts <- justCrosswind$Weight_g


# simulating a vector of estimated mean weights, as well as a vector of true mean weights
# note: these will be applied to the vector of predicted YP simulated earlier

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

hist(predYP/est_mn)
hist(log(predYP/est_mn))

# relative precision incorporating variability from LAM as well as weight sampling
rp(vec = predYP/est_mn,
   trueval = median(predYP)/mn_true)  # 235%!!
quantile((predYP/est_mn)/(median(predYP)/mn_true), c(0.025, 0.975))

# relative precision treating LAM as constant
rp(vec = median(predYP)/est_mn,
   trueval = median(predYP)/mn_true)  # 15% is a little more friendly
#      2.5%     97.5%
# 0.2329519 4.1649215

# treating Lester as constant: MSY should be some fraction of YP but unknown
# will it behave the same way?  Looks like it.

lakenames[ilake]
rp(vec = 0.5*median(predYP)/est_mn,
   trueval = 0.5*median(predYP)/mn_true)  # 15% - yes
