source("FrontMatter.R")

## Chapter 4

plotDIR4 <- paste(plotDIR, "chapter4", "figures", sep="/")

packages(rv)
packages(rstan)
packages(car)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## GOM Hypoxia
## This example uses three model formulations
## The original model was presented in Qian et al (2009)
## The alternative models 1 and 2 are proposed to resolve
##    numerical issues

benthic.data <- read.csv(paste(dataDIR, "BenthicData.csv", sep="/"),
                         header=T)
benthic.data$area2 <-
    ordered(benthic.data$area2,
            levels=levels(ordered(benthic.data$area2))[c(2,1,3)])
benthic.data$Area <-
    ordered(benthic.data$Area,
            levels=levels(ordered(benthic.data$Area))[c(2,1,3)])

head(benthic.data)

## Station: core,
## Area: zone
## The original model
stan_gom <- " /* Stan Code - GOM model*/
data{
  int K;  //total sample size
  int J;  //number of sediment cores
  int I;  //number of zones
  real y[K]; //observed response
  int core[K]; //core index
  int zone[K]; //zone index
  int core_zone[J]; //zone
}
parameters{
  real mu[J];
  real theta[I];
  real mu_hyp;
  real<lower=0> sigma_y[I];
//  real<lower=0> sigma_y;
//  real<lower=0> sigma_i[I];
  real<lower=0> sigma_i;
  real<lower=0> sigma_hyp;
}
model{
  sigma_hyp ~ normal(0,1);
  for (i in 1:I){
    theta[i] ~ normal(mu_hyp, sigma_hyp);
  }
  for (j in 1:J){
    mu[j] ~ normal(theta[core_zone[j]], sigma_i);
  }
  for (k in 1:K){
    y[k] ~ normal(mu[core[k]], sigma_y[zone[k]]);
  }
}
generated quantities{
  real delta1;
  real delta2;
  real delta3;
  delta1 = theta[2]-theta[1];
  delta2 = theta[2]-theta[3];
  delta3 = theta[1]-theta[3];
}
"

stan.in3 <- function(data = benthic.data, y.col=5,
                     chains=nchains){ ## no slope
    n <- dim(data)[1]
    y <- log(data[,y.col])
    core <- as.numeric(ordered(data$Station))
    n.core <- max(core)
    zone <- as.numeric(ordered(data$Area))
    n.zone <- max(zone)
    oo <- order(core)
    ind <- cumsum(table(core[oo]))
    Core.zone <- zone[oo][ind] ## each core belongs to which zone
    stan.dat <- list(K=n, J=n.core, I=n.zone, y=y, core=core,
                     zone=zone, core_zone=Core.zone)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list( mu = rnorm(n.core), theta=rnorm(n.zone),
                           mu_hyp=rnorm(1), sigma_i=runif(1),
                           sigma_y=runif(n.zone), sigma_hyp=runif(1))
    parameters <- c("mu","theta","mu_hyp", "sigma_y", "sigma_i",
                    "sigma_hyp",  "delta1", "delta2","delta3")
    return(list(para=parameters, data=stan.dat, inits=inits,
                n.chains=chains))
}

input.to.stan <- stan.in3()  ## long-runs -- results without sigma_hyp prior saved
stan.fit <- stan_model(model_code=stan_gom)
fit2keep <- sampling(stan.fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99, max_treedepth=20))
print(fit2keep)
rich_stan <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))
save(rich_stan, file="ch4_rich_stanN01.RData")

input.to.stan <- stan.in3(y.col=6)
fit2keep <- sampling(stan.fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99, max_treedepth=20))
print(fit2keep)
abun_stan <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))
save(abun_stan, file="ch4_abun_stanN01.RData")

## estimated core means from the default model
tempA <- abun_stan[1:15]
names(tempA) <- paste("$\\mu_{", 1:15, "}$", sep="")
tempR <- rich_stan[1:15]
names(tempR) <- paste("$\\mu_{", 1:15, "}$", sep="")

delta0R <- rich_stan[25:27]
delta0A <- abun_stan[25:27]
names(delta0R) <- c("H-I","H-O","I-O")
names(delta0A) <- c("H-I","H-O","I-O")

tikz(file=paste(plotDIR4, "ModelCores.tex", sep="/"),
     height=4, width=5.75, standAlone=F)
par(mfrow=c(1,2), mar=c(3, 3, 2, 1), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
mlplot(tempA, xlab="log abundance")
mlplot(tempR, xlab="log richness")
dev.off()

#### Posterior simulations
## replicate data:
## saved posterior rv objects
## "ch4_abun_stan.RData" -- defalt prior for sigma_hyp
## "ch4_abun_stanN01.RData" -- prior for sigma_hyp is N(0,1)
load("ch4_abun_stanN01.RData")
load("ch4_rich_stanN01.RData")
rich_stan01 <- rich_stan
abun_stan01 <- rich_stan
load("ch4_abun_stan.RData")
load("ch4_rich_stan.RData")
tikz(paste(plotDIRch4, "funnel_rich.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mfrow=c(1,2),mar=c(3,3,1,1),
    mgp=c(1.25,0.125,0), las=1,tck=0.01)
plot(rich_stan$mu_hyp, log(rich_stan$sigma_hyp), cex=0.5,
     ylim=c(-3.5,log(150)), xlim=c(-100,100),
     xlab="$\\mu_{hyp}$", ylab="$\\log(\\sigma_{hyp})$")
plot(rich_stan01$mu_hyp, log(rich_stan01$sigma_hyp), cex=0.5,
##     ylim=c(-3.5,log(150)), xlim=c(-100,100),
     xlab="$\\mu_{hyp}$", ylab="$\\log(\\sigma_{hyp})$")
dev.off()

## 1. zone means from hyper-distribution
nms <- names(abun_stan)
## zone means
Abun_sim1 <- rvnorm(3, abun_stan[nms=="mu_hyp"], abun_stan[nms=="sigma_hyp"])
## core means
Abun_sim2 <- c(rvnorm(4, Abun_sim1[1], abun_stan[nms=="sigma_i"]),
               rvnorm(5, Abun_sim1[2], abun_stan[nms=="sigma_i"]),
               rvnorm(6, Abun_sim1[3], abun_stan[nms=="sigma_i"]))
## No need to move on -- the hyper distribution is not defined!

benthic.data$Zone <- with(benthic.data,
                          factor(as.factor(Area):as.factor(Station)))

## reparameterizing hyper-parameters

stan_gom_eta <- " /* Stan Code - GOM model*/
data{
  int K;  //total sample size
  int J;  //number of sediment cores
  int I;  //number of zones
  real y[K]; //observed response
  int core[K]; //core index
  int zone[K]; //zone index
  int core_zone[J]; //zone
}
parameters{
  real mu[J];
  real mu_hyp;
  real eta[I];
  real<lower=0> sigma_y[I];
  real<lower=0> sigma_i;
  real<lower=0> sigma_hyp;
}
transformed parameters{
  real theta[I];
  for (i in 1:I)
    theta[i] = mu_hyp+sigma_hyp*eta[i];
}
model{
  eta ~ normal(0,1);
  for (j in 1:J){
    mu[j] ~ normal(theta[core_zone[j]], sigma_i);
  }
  for (k in 1:K){
    y[k] ~ normal(mu[core[k]], sigma_y[zone[k]]);
  }
}
generated quantities{
  real delta1;
  real delta2;
  real delta3;
  delta1 = theta[2]-theta[1];
  delta2 = theta[2]-theta[3];
  delta3 = theta[1]-theta[3];
}
"
stan.in3_5 <- function(data = benthic.data, y.col=5,
                       chains=nchains){ ## no slope
    n <- dim(data)[1]
    y <- log(data[,y.col])
    core <- as.numeric(ordered(data$Station))
    n.core <- max(core)
    zone <- as.numeric(ordered(data$Area))
    n.zone <- max(zone)
    oo <- order(core)
    ind <- cumsum(table(core[oo]))
    Core.zone <- zone[oo][ind] ## each core belongs to which zone
    stan.dat <- list(K=n, J=n.core, I=n.zone, y=y, core=core,
                     zone=zone, core_zone=Core.zone)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(mu = rnorm(n.core), eta=rnorm(n.zone),
                           mu_hyp=rnorm(1), sigma_i=runif(1),
                           sigma_y=runif(n.zone), sigma_hyp=runif(1))
    parameters <- c("mu","theta","mu_hyp", "sigma_y", "sigma_i",
                    "sigma_hyp",  "delta1", "delta2","delta3")
    return(list(para=parameters, data=stan.dat, inits=inits,
                n.chains=chains))
}

input.to.stan <- stan.in3_5()  ##
stan.fit <- stan_model(model_code=stan_gom_eta)
fit2keep <- sampling(stan.fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99, max_treedepth=20))
print(fit2keep)
rich_stan <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))
save(rich_stan, file="ch4_rich_stanN01.RData")

input.to.stan <- stan.in3(y.col=6)
fit2keep <- sampling(stan.fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99, max_treedepth=20))
print(fit2keep)


###
### First alternative (One-way ANOVA)
###
stan_gom2 <- " /* Stan Code - GOM model*/
data{
  int K;  //total sample size
  int J;  //number of sediment cores
  real y[K]; //observed response
  int core[K]; //core index
}
parameters{
  real mu[J];
  real mu_hyp;
  real<lower=0> sigma_y;
  real<lower=0> sigma_hyp;
}
model{
//  mu_hyp ~ normal(0, 0.5);
  sigma_hyp ~ normal(0,0.5);
  sigma_y ~ normal(0,0.5);
  for (j in 1:J){
    mu[j] ~ normal(mu_hyp, sigma_hyp);
  }
  for (k in 1:K){
    y[k] ~ normal(mu[core[k]], sigma_y);
  }
}
"
stan.in4 <- function(data = benthic.data, y.col=5,
                     chains=nchains){ ## no slope
  n <- dim(data)[1]
  y <- log(data[,y.col])
  y_ave <- mean(y)
  y_sd <- sd(y)

  y <- (y-y_ave)/y_sd
  core <- as.numeric(ordered(data$Station))
  n.core <- max(core)

  stan.dat <- list(K=n, J=n.core, y=y, core=core)
  inits <- list()
  for (i in 1:chains)
    inits[[i]] <- list( mu = rnorm(n.core),
                        mu_hyp=rnorm(1),
                        sigma_y=runif(1), sigma_hyp=runif(1))
  parameters <- c("mu","mu_hyp", "sigma_y", "sigma_hyp")
  return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains, y_cen=y_ave, y_spd=y_sd))
}

input.to.stan <- stan.in4()
stan.fit2 <- stan_model(model_code=stan_gom2)
fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)
rich_stan2 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

## processing output
core <- as.numeric(ordered(benthic.data$Station))
n.core <- max(core)
zone <- as.numeric(ordered(benthic.data$Area))
n.zone <- max(zone)
oo <- order(core)
ind <- cumsum(table(core[oo]))
Core.zone <- zone[oo][ind] ## each core belongs to which zone
input_sd <- input.to.stan$y_spd
input_mu <- input.to.stan$y_cen

core1.musR <- input_mu + input_sd*rich_stan2[1:15]
zone11.Rmu <- mean(core1.musR[Core.zone==1])
zone12.Rmu <- mean(core1.musR[Core.zone==2])
zone13.Rmu <- mean(core1.musR[Core.zone==3])

deltaR11 = zone12.Rmu-zone11.Rmu
deltaR12 = zone12.Rmu-zone13.Rmu
deltaR13 = zone11.Rmu-zone13.Rmu

delta1R <- c(deltaR11, deltaR12, deltaR13)
names(delta1R) <- c("H-I","H-O","I-O")

input.to.stan <- stan.in4(y.col=6)
fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)
abun_stan2 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

input_sd <- input.to.stan$y_spd
input_mu <- input.to.stan$y_cen
core1.musA <- input_mu + input_sd * abun_stan2[1:15]
zone11.Amu <- mean(core1.musA[Core.zone==1])
zone12.Amu <- mean(core1.musA[Core.zone==2])
zone13.Amu <- mean(core1.musA[Core.zone==3])

deltaA11 = zone12.Amu-zone11.Amu
deltaA12 = zone12.Amu-zone13.Amu
deltaA13 = zone11.Amu-zone13.Amu

delta1A <- c(deltaA11, deltaA12, deltaA13)
names(delta1A) <- c("H-I","H-O","I-O")

save(rich_stan2, abun_stan2, "ch4_alt2_stan.RData")

## The second alternative model (nested ANOVA)
## Referencing the lme4 book

stan_gom3 <- " /* Stan Code - GOM model*/
data{
  int K;  //total sample size
  int J;  //number of sediment cores
  int I;  //number of zones
  real y[K]; //observed response
  int core[K]; //core index
  int zone[K]; //zone index
}
parameters{
  real muK[J];
  real muZ[I];
  real mu0;
  real<lower=0> sigmaY;
  real<lower=0> sigmaK;
  real<lower=0> sigmaZ;
}
model{
  sigmaK ~ normal(0,0.5);
  sigmaZ ~ normal(0,0.5);
  sigmaY ~ normal(0,0.5);
  for (i in 1:I){
    muZ[i] ~ normal(0, sigmaZ);
  }
  for (j in 1:J){
    muK[j] ~ normal(0, sigmaK);
  }
  for (k in 1:K){
    y[k] ~ normal(mu0+muK[core[k]]+muZ[zone[k]], sigmaY);
  }
}
generated quantities{
  real delta1;
  real delta2;
  real delta3;
  delta1 = muZ[2]-muZ[1];
  delta2 = muZ[2]-muZ[3];
  delta3 = muZ[1]-muZ[3];
}
"
stan.in5 <- function(data = benthic.data, y.col=5, chains=nchains){
  n <- dim(data)[1]
  y <- log(data[,y.col])
  y_ave <- mean(y)
  y_sd <- sd(y)
  y <- (y-y_ave)/y_sd
  core <- as.numeric(ordered(data$Station))
  n.core <- max(core)
  zone <- as.numeric(ordered(data$Area))
  n.zone <- max(zone)

  stan.dat <- list(K=n, J=n.core, I=n.zone, y=y, core=core, zone=zone)
  inits <- list()
  for (i in 1:chains)
    inits[[i]] <- list( muK = rnorm(n.core),
                        muZ=rnorm(n.zone),
                        mu0=rnorm(1),
                       sigmaY=runif(1), sigmaK=runif(1),
                       sigmaZ=runif(1))
  parameters <- c("mu0","muK", "muZ", "sigmaY", "sigmaK",
                  "sigmaZ", "delta1", "delta2", "delta3")
  return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains, y_cen=y_ave, y_spd=y_sd))
}

input.to.stan <- stan.in5()
stan.fit3 <- stan_model(model_code=stan_gom3)
fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99, max_treedepth=15))
print(fit2keep)
rich_stan3 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

### posterior simulation (Chapter 8) ###
#### replicating input data
input_sd <- input.to.stan$y_spd
input_mu <- input.to.stan$y_cen

nms3 <- names(rich_stan3)
## 1. core means
muK <- rich_stan3[substring(nms3, 1, 3)=="muK"]
## zone 1
tmp <- unique(input.to.stan$data$core[input.to.stan$data$zone==1])
zn1_core <- rich_stan3[nms3=="mu0"]+muK[tmp]+rich_stan3[nms3=="muZ.1"]
zn1_rep <- rvnorm(3, zn1_core, rich_stan3[nms3=="sigmaY"])
zn1_mean <- input_mu+input_sd*mean(zn1_rep)
zn1_range <- input_mu+input_sd*range(zn1_core)

## zone 2
tmp <- unique(input.to.stan$data$core[input.to.stan$data$zone==2])
zn2_core <- rich_stan3[nms3=="mu0"]+muK[tmp]+rich_stan3[nms3=="muZ.2"]
zn2_rep <- rvnorm(3, zn2_core, rich_stan3[nms3=="sigmaY"])
zn2_mean <- input_mu+input_sd*mean(zn2_rep)
zn2_range <- input_mu+input_sd*range(zn2_core)
## zone 3
tmp <- unique(input.to.stan$data$core[input.to.stan$data$zone==3])
zn3_core <- rich_stan3[nms3=="mu0"]+muK[tmp]+rich_stan3[nms3=="muZ.3"]
zn3_rep <- rvnorm(3, zn3_core, rich_stan3[nms3=="sigmaY"])
zn3_mean <- input_mu+input_sd*mean(zn3_rep)
zn3_range <- input_mu+input_sd*range(zn3_core)

## simulated differences
delta1 <- zn2_mean-zn1_mean
delta2 <- zn2_mean-zn3_mean
delta3 <- zn1_mean-zn3_mean

## data means
core.zone <- cbind(input.to.stan$data$zone, input.to.stan$data$core)
core.zone <- core.zone[order(core.zone[,2]),]
ttt <- cumsum(table(core.zone[,2]))
czind <- core.zone[ttt,1]
names(czind) <- 1:length(czind)
coreMean <- input_mu+input_sd*
    tapply(input.to.stan$data$y, input.to.stan$data$core, mean)
zoneMean <- tapply(coreMean, czind, mean)
zoneRange <- tapply(coreMean, czind, range)

tikz(paste(plotDIRch4, "gompostSim1_rich.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01)
plot(c(zn1_mean, zn2_mean, zn3_mean), xlab="Zone", ylab="zone means",
     axes=F)
axis(1, at=1:3)
axis(2)
box()
points(1:3, zoneMean, pch=4)

plot(c(delta1, delta2, delta3), xlab="$\\delta$",
     ylab="zone mean differences", axes=F)
axis(1, at=1:3)
axis(2, at=1:3, labels=c("$\\delta_1$","$\\delta_3$","$\\delta_3$"))
box()
points(1:3, c(zoneMean[2]-zoneMean[1],
              zoneMean[2]-zoneMean[3],
              zoneMean[1]-zoneMean[3]),
       pch=4)
dev.off()

tikz(paste(plotDIRch4, "gompostSim2_rich.tex", sep="/"),
     height=2.5, width=2.5, standAlone=F)
par( mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01)
plot(c(zn1_range[2]-zn1_range[1],
       zn2_range[2]-zn2_range[1],
       zn3_range[2]-zn3_range[1]), xlab="Zone",
     ylab="zone mean range",
     axes=F)
axis(1, at=1:3)
axis(2)
box()
points(1:3, unlist(lapply(zoneRange, diff)), pch=4)
dev.off()

###### end of post sim ######

input_sd <- input.to.stan$y_spd
input_mu <- input.to.stan$y_cen
zone2.musR <- input_sd*rich_stan3[17:19]
core2.musR <- input_mu+input_sd*rich_stan3[2:16]

zone21.Rmu <- mean(core2.musR[Core.zone==1]) + zone2.musR[1]
zone22.Rmu <- mean(core2.musR[Core.zone==2]) + zone2.musR[2]
zone23.Rmu <- mean(core2.musR[Core.zone==3]) + zone2.musR[3]

deltaR21 = zone22.Rmu-zone21.Rmu
deltaR22 = zone22.Rmu-zone23.Rmu
deltaR23 = zone21.Rmu-zone23.Rmu

delta2R <- c(deltaR21, deltaR22, deltaR23)
names(delta2R) <- c("H-I","H-O","I-O")

input.to.stan <- stan.in5(y.col=6)
fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta=0.99,
                                  max_treedepth=15))
print(fit2keep)
abun_stan3 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

input_sd <- input.to.stan$y_spd
input_mu <- input.to.stan$y_cen
zone2.musA <- input_sd*abun_stan3[17:19]
core2.musA <- input_mu+input_sd*abun_stan3[2:16]


### saving all model outputs for later use
save(rich_stan2, abun_stan2, rich_stan3, abun_stan3,
     core1.musR, core1.musA, zone2.musR, zone2.musA,
     core2.musR, core2.musA,
     file="GOMalt23.RData")

zone21.Amu <- mean(core2.musA[Core.zone==1]) + zone2.musA[1]
zone22.Amu <- mean(core2.musA[Core.zone==2]) + zone2.musA[2]
zone23.Amu <- mean(core2.musA[Core.zone==3]) + zone2.musA[3]

deltaA21 = zone22.Amu-zone21.Amu
deltaA22 = zone22.Amu-zone23.Amu
deltaA23 = zone21.Amu-zone23.Amu

delta2A <- c(deltaA21, deltaA22, deltaA23)
names(delta2A) <- c("H-I","H-O","I-O")

## expected
delta2A.Exp <- zone22.Amu-0.5*(zone21.Amu+zone23.Amu)
delta2R.Exp <- zone22.Rmu-0.5*(zone21.Rmu+zone23.Rmu)


## Figures
tikz(file=paste(plotDIR4, "ModelCoresR.tex", sep="/"),
     height=5, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(tempR, xlab="Original model", top.axis=F)
mlplot(core1.musR, xlab="Alternative model 1",
       main="log Richness", axes=F)
axis(1)
mlplot(core2.musR, xlab="Alternative model 2", axes=F)
axis(1)
dev.off()

tikz(file=paste(plotDIR4, "ModelCoresA.tex", sep="/"),
     height=5, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(tempA, xlab="Original model", top.axis=F)
mlplot(core1.musA, xlab="Alternative model 1",
       main="log Abundance", axes=F)
axis(1)
mlplot(core2.musA, xlab="Alternative model 2", axes=F)
axis(1)
dev.off()

zone0R.thetas <- rich_stan[16:18]
names(zone0R.thetas) <- c("Inshore", "Hypoxic","Offshore")
zone0A.thetas <- abun_stan[16:18]
names(zone0A.thetas) <- c("Inshore", "Hypoxic","Offshore")

zone1R.thetas <- c(zone11.Rmu, zone12.Rmu, zone13.Rmu)
names(zone1R.thetas) <- c("Inshore", "Hypoxic","Offshore")

zone1A.thetas <- c(zone11.Amu, zone12.Amu, zone13.Amu)
names(zone1A.thetas) <- c("Inshore", "Hypoxic","Offshore")

zone2R.thetas <- c(zone21.Rmu, zone22.Rmu, zone23.Rmu)
names(zone2R.thetas) <- c("Inshore", "Hypoxic","Offshore")

zone2A.thetas <- c(zone21.Amu, zone22.Amu, zone23.Amu)
names(zone2A.thetas) <- c("Inshore", "Hypoxic","Offshore")

tikz(file=paste(plotDIR4, "ModelZonesR.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(zone0R.thetas, xlab="Original model", top.axis=F, xlim=c(2,3.5))
abline(v=0)
mlplot(zone1R.thetas, xlab="Alternative model 1", main="log Richness",
       axes=F, xlim=c(2,3.5))
axis(1)
abline(v=0)
mlplot(zone1R.thetas, xlab="Alternative model 2", axes=F, xlim=c(2,3.5))
axis(1)
abline(v=0)
dev.off()

tikz(file=paste(plotDIR4, "ModelZonesA.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(zone0A.thetas, xlab="Original model", top.axis=F, xlim=c(3.6,6))
mlplot(zone1A.thetas, xlab="Alternative model 1", main="log Abundance",
       axes=F, xlim=c(3.6,6))
axis(1)
mlplot(zone1A.thetas, xlab="Alternative model 2", axes=F, xlim=c(3.6,6))
axis(1)
dev.off()

tikz(file=paste(plotDIR4, "ModelZonesDelA.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(delta0A, xlab="Original model 1", top.axis=F, xlim=c(-2,1.5))
abline(v=0)
mlplot(delta1A, xlab="Alternative model 1",
       main="log Abundance", axes=F,
       xlim=c(-2,1.5))
axis(1)
abline(v=0)
mlplot(delta2A, xlab="Alternative model 2", axes=F, xlim=c(-2,1.5))
axis(1)
abline(v=0)
dev.off()

tikz(file=paste(plotDIR4, "ModelZonesDelR.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(delta0R, xlab="Original model", top.axis=F, xlim=c(-1.5,0.5))
abline(v=0)
mlplot(delta1R, xlab="Alternative model 1",
       main="log Richness", axes=F, xlim=c(-1.5,0.5))
abline(v=0)
axis(1)
mlplot(delta2R, xlab="Alternative model 2", axes=F, xlim=c(-1.5,0.5))
axis(1)
abline(v=0)
dev.off()


