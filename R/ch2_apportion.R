source("FrontMatter.R")

## Chapter 2
plotDIRch2 <- paste(plotDIR, "chapter2", "figures", sep="/")
##########################################################################
############# Apportionment - P(TL|TS) - in Stan for BEESRStan  ###########
##########################################################################
#use this when generating figures for the book
#source("FrontMatter.R")


## load packages
packages(rstan)
##packages(rstantools)
packages(withr)
packages(rv)
##library(rstudioapi)
packages(arm)



### the set up
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))
nchains <- min(c(parallel::detectCores(), 8))
niters <- 5000 #50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## ## data
## ## 1. hydro-acoustic data Data file is too big for GitHub
## data <- read.csv(paste(dataDIR, "HA_targ_18.csv", sep="/"),
##                  header=T, sep=",")
## data <- subset(data, Grid == 4)

## data <- data[order(data$Grid,data$Region_name),]
## head(data)
## dim(data)

## save(data, file="Ch2_HAtarg.RData")
## preparing data

load("Ch2_HAtarg.RData")
y <- log10(data$BS)
par(mfrow=c(1,2))
hist(log10(data$BS))
hist(data$TS_comp)
G <- as.numeric(as.factor(data$uID))
max(G)
##############################################

##############################################
## TS mean model --
##    generating mean target strength from
##    acoustic readings
##############################################
pTL_TS <-"
  data {
    int<lower = 1> N;   // Total number of trials
    int<lower = 1> N_G;
    int<lower=1,upper=N_G> G[N];     // TS groups - fish
    vector[N] y;        // Score in each trial
    }

  parameters {
    vector[N_G] mu;
    real<lower = 0> sigma;
    real<lower = 0> sigma_mu;
    }

  model {
    sigma_mu ~ normal(0,25);
    sigma ~ normal(0, 5);
    mu ~ normal(0, sigma_mu);
    y ~normal(mu[G], sigma);
    }
  "
## The data and initial values
##############################################
pTL_TS_in <- function(y_in, G, n.chains=nchains){
  N <- length(y[!is.na(y)])
  N_G <- max(G[!is.na(y_in)])
  G <- G[!is.na(y_in)]
  y <- y[!is.na(y_in)]

  data <- list(N=N, N_G=N_G, G=G, y=y)
  inits <- list()
  for (i in 1:n.chains)
      inits[[i]] <- list(mu=(rnorm(N_G,-5,2)),
                         sigma=(runif(1,0,7)),
                         sigma_mu=(runif(1,0,2)))
  paras <- c("mu","sigma","sigma_mu")
  return(list(data=data, init=inits,
              para=paras, nchains=n.chains))
}


##############################################
## Compiling, input data, and running the model
##############################################
fit <- stan_model(model_code = pTL_TS)

input.to.stan <- pTL_TS_in(y,G)

keep1 <- sampling(fit, data=input.to.stan$data,
                  init=input.to.stan$init,
                  pars=input.to.stan$para,
                  iter=niters,thin=nthin,
                  chains=input.to.stan$nchains)
save(keep1,file="pTL_TS.RData")
#############################################

##############################################
## Processing Stan output and summaries - for book
##############################################
#temp <- rstan::extract(keep, permuted=T)
# stan_sim <- matrix(c(temp$zbeta0,
#                      temp$zbeta1,
#                      temp$beta0,
#                      temp$beta1,
#                      temp$beta2,
#                      temp$beta3,
#                      temp$phi), ncol=6, byrow=F)

#summary(stan_sim)
##keep1 <- load("pTL_TS.RData")
# use this when you get above code working...
fitcoef1 <- rvsims(as.matrix(as.data.frame(
  extract(keep1, pars="mu",permute=T))))


##options(digits=2)
##summary(fitcoef1)


TSm <- 10*fitcoef1
TLm <- 10^((TSm + 64.23941176 )/19.2)

TSm_sub <- subset(TSm, rvmean(TSm) > -40)
TLm_sub <- 10^((TSm_sub + 64.23941176 )/19.2)
oo <- rev(order(summary(TSm_sub)[,2]))

names(TSm_sub) <- rep("-", length(TSm_sub))
names(TLm_sub) <- rep("-", length(TLm_sub))


## coef plot
tikz(file=paste(plotDIRch2, "pTL_TS_ests.tex", sep="/"),
    width=5,height=3, standAlone=F)

par(mfrow=c(1,2), mar=c(3,1,1,1), oma=c(1,3,1,1), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
mlplot(TSm_sub[oo], xlab="Mean TS (dB)", cex=0.5)
mlplot(TLm_sub[oo], xlab="Mean TL (cm)", cex=0.5)
abline(v=35, lty=3, lwd=2)
mtext("Individual Fish Tracks", 2, 0, outer=T, las=0, cex=1.5)

dev.off()
##############################################

##############################################
## Logistic regression model
##############################################
## data:

dataWL <- read.csv(paste(dataDIR, "GN_length.csv", sep="/"),
                   header=T, sep=",") ## survey data
dataWL <- dataWL[order(dataWL$Grid),]
head(dataWL)

y <- dataWL$walleye
x <- dataWL$LENGTH/10

pW_TL <-"
data{
  int<lower=1> n;
  int<lower=0> y[n];
  real<lower=0> x[n];
}
parameters{
  real alpha;
  real beta;
}
transformed parameters{
  vector[n] theta;
  for(i in 1:n){
    theta[i] = alpha + beta*x[i];
  }
}
model{
  target += bernoulli_logit_lpmf(y | theta);
}
"

##############################################
## The data and initial values
##############################################
pW_TL_in <- function(y_in,x, n.chains=nchains){
  n <- length(y[!is.na(y)])
  y <- y[!is.na(y_in)]
  x <- x[!is.na(y_in)]

  data <- list(n=n, y=y, x=x)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(alpha=(rnorm(1,-30,50)),
                       beta=(rnorm(1,0,30)))
  paras <- c("alpha","beta")
  return(list(data=data, init=inits,
              para=paras, nchains=n.chains))
}

fit <- stan_model(model_code = pW_TL)

input.to.stan <- pW_TL_in(y,x)
keep2 <- sampling(fit, data=input.to.stan$data,
                  init=input.to.stan$init,
                  pars=input.to.stan$para,
                  iter=niters,thin=nthin,
                  chains=input.to.stan$nchains)

save(keep2, file="pW_TLlogit.RData")

fitcoef2 <- (rvsims(as.matrix(as.data.frame(
  extract(keep2, permute=T)))))

options(digits=2)
summary(fitcoef2)


## histogram of paramters
png("pW_TL_params.png",
    width=8,height=4, units = 'in', res = 450, pointsize=8, bg="white")

par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,2,3,1))
rvhist(fitcoef2[1], xlab="Intercept", main="", col="gray")
abline(v=rvmean(fitcoef2[1]), lwd=3)
rvhist(fitcoef2[2], xlab="Slope", main="", col="gray")
abline(v=rvmean(fitcoef2[2]), lwd=3)
mtext("Logistic model parameter estimates", 3, 0, outer=T, cex=1.5)

dev.off()

## plot logistic curve
x <- seq(10,70, by=0.1)

tikz(file=paste(plotDIRch2, "pW_TL_ests.tex", sep="/"),
   width=3.5,height=2.5,standAlone=F)

par(mar=c(3,3,1,1),mgp=c(1.25,0.125,0), las=1,tck=0.01)
plot(dataWL$LENGTH/10,dataWL$walleye, ylab="",col=gray(0,0.2), xlab="",
     col.axis="black",col.lab="black",cex=1, cex.axis=1.5, cex.lab=1.5,
     bty="n", fg="black",xlim=c(10,70), las=1)
curve(invlogit(rvmean(fitcoef2[1]) + rvmean(fitcoef2[2])*x), add=T, lwd=3)

mtext("$\\Pr(W | TL, \\hat{\\beta})$",2,2,cex=1.5, las=0)
mtext("Mean length (cm)",1,2,cex=1.5)
dev.off()

##############################################


############################################################################
## Step 3: Integrate out TS and logistic model parameters - p(W|TS,Betas)
############################################################################
pTL_TS <- TLm
pW_TL <- fitcoef2

pW <- invlogit(pW_TL[1]+pW_TL[2]*pTL_TS)
len <- (logit(rvmean(pW))-rvmean(pW_TL[1]))/rvmean(pW_TL[2]) ## length associate with expected pSC

meanPrW <- rvmean(pW)
likely <- sum(rvmean(pW)>=0.8)
uncertain <- sum(rvmean(pW)<0.8&rvmean(pW)>0.2)
unlikely <- sum(rvmean(pW)<=0.2)

##png("pW_TS.png",
##    width=6,height=4, units = 'in', res = 450, pointsize=8, bg="white")
tikz(file=paste(plotDIRch2, "pW_TL.tex", sep="/"),
     height=2.5, width=3.25, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), las=1,tck=0.01)
x <- seq(0,50,by=0.5)
plot(x,invlogit(rvmean(pW_TL[1])+rvmean(pW_TL[2])*x), typ="l",
     ylab="$\\Pr(\\mathrm{Walleye} | TS,\\beta)$",
     xlab="Total length (cm)")
points(len,meanPrW, pch=16, col=gray(0,0.2), cex=0.5)
abline(h=0.8, lty=3, col=gray(0,0.2))
abline(h=0.2, lty=3, col=gray(0,0.2))
text(10,0.9,paste(likely,"highly probable"), cex=0.75)
text(10,0.5,paste(uncertain,"uncertain"), cex=0.75)
text(10,0.1,paste(unlikely,"unlikely"), cex=0.75)
dev.off()


###########################################################################
## Step 4: Generate uncertainty in Walleye counts
###########################################################################
## summing the expected probabilities produces the
## estimated mean number of Walleye in the hydroacoustic sample
sum(meanPrW) ## out of 15 fish tracks

## But, there is binomial uncertainty associated with each fish track
## For example, a fish with 80% pSC also has 20% !pSC
## Over 1,000 trials that fish would be SC 800 times and !SC 200 times
## We preform 1,000 trials for each fish track and sum within trial
## to generate uncertainty in the expected count

## A function to carry out the Monte Carlo trials
cW <- rvbinom(n=1,size=1,prob=pW) ## Take 1000 random draws (trials) fro

mean(cW) ## should be ~ the same as sum(rvmean(pW))

##png("Walleye_counts.png", width=8,height=4,
##units = 'in', res = 450, pointsize=8, bg="white")
tikz(file=paste(plotDIRch2, "walleyetotalN.tex", sep="/"),
     height=2.25, width=3, standAlone=F)
par(mar=c(3, 2, 1, 1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
rvhist(mean(cW)*length(cW), xlab="Estimated Walleye Count",col="gray",
       ylab="", main="", freq=F)

dev.off()

###########################################################################

