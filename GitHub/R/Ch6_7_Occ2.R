source("FrontMatter.R")

## Chapter 6
plotDIRch6 <- paste(plotDIR, "chapter6", "figures", sep="/")
packages(tidyverse)
## simulation
packages(rv)
packages(rstan)
packages(readxl)
packages(car)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 20000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Atrizine occuurence in US drinking water
#### Large data -- runs very slow
atrazine <- read.csv(file=paste(dataDIR, "atrazine.csv", sep="/"))
names(atrazine)

## Detection limit units -- missing and mislabeled units
unique(as.character(atrazine$Detection.Limit.Unit))
table(atrazine$Detection.Limit.Value, atrazine$Detection.Limit.Unit)

atrazine$Adj.Detect.Limit <- atrazine$Detection.Limit.Value

temp <- atrazine$Detection.Limit.Value < 0.005 & !is.na(atrazine$Detection.Limit.Value)
atrazine$Adj.Detect.Limit[temp] <- 1000*atrazine$Adj.Detect.Limit[temp]
atrazine$Adj.Detect.Limit[atrazine$Adj.Detect.Limit==0 & !is.na(atrazine$Adj.Detect.Limit)] <- 0.1
## If the detection limit value is >= 0.005 the unit should be ug/L, otherwise mg/L
## If the detection limit value is 0, it is most likely 0.1 ug/L

## separating source type: finished water versus raw water (and unknown)
table(atrazine$Source.Type.Code)
atr.RAW <- atrazine[atrazine$Source.Type.Code=="RW",]
atr.FIN <- atrazine[atrazine$Source.Type.Code=="FN",]
atr.OTH <- atrazine[atrazine$Source.Type.Code!="FN"&atrazine$Source.Type.Code!="RW",]

rnpop <- range(atr.FIN$Adjusted.Total.Population.Served)
pop.served <- cut(atr.FIN$Adjusted.Total.Population.Served,
                  breaks=c(rnpop[1],500, 3300, 10000, 50000,rnpop[2]), labels=1:5)
atr.FIN$group <- paste(substring(atr.FIN$Source.Water.Type, 1,1), pop.served)
atr.FIN$group[is.na(pop.served)] <- NA
atr.FIN <- atr.FIN[!is.na(atr.FIN$group),]

rnpop <- range(atr.RAW$Adjusted.Total.Population.Served)
pop.served <- cut(atr.RAW$Adjusted.Total.Population.Served,
                  breaks=c(rnpop[1],500, 3300, 10000, 50000,rnpop[2]),
                  labels=1:5)
atr.RAW$group <- paste(substring(atr.RAW$Source.Water.Type, 1,1), pop.served)
atr.RAW$group[is.na(pop.served)] <- NA
atr.RAW <- atr.RAW[!is.na(atr.RAW$group),]

rnpop <- range(atr.OTH$Adjusted.Total.Population.Served)
pop.served <- cut(atr.OTH$Adjusted.Total.Population.Served,
                  breaks=c(rnpop[1],500, 3300, 10000, 50000,rnpop[2]), labels=1:5)
atr.OTH$group <- paste(substring(atr.OTH$Source.Water.Type, 1,1), pop.served)
atr.OTH$group[is.na(pop.served)] <- NA
atr.OTH <- atr.OTH[!is.na(atr.OTH$group),]

## all data
rnpop <- range(atrazine$Adjusted.Total.Population.Served)
pop.served <- cut(atrazine$Adjusted.Total.Population.Served,
                  breaks=c(rnpop[1],500, 3300, 10000, 50000,rnpop[2]), labels=1:5)
atrazine$group <- paste(substring(atrazine$Source.Water.Type, 1,1), pop.served)
atrazine$group[is.na(pop.served)] <- NA
atrazinALL <- atrazine[!is.na(atrazine$group),]


## stan model
occ_stan1 <- "
data {
  int<lower=0> N;
  int<lower=0> Ncens;
  int<lower=0> Ngrp;
  int<lower=0> Nsys;
  int<lower=1,upper=Nsys> sid[N];
  int<lower=1,upper=Nsys> sidcens[Ncens];
  int<lower=1,upper=Ngrp> sys_grp[Nsys];
  vector[N] y;
  vector[Ncens] ycens;
}
parameters {
  vector[Nsys] theta;
  vector[Ngrp] gamma;
  real mu;
  real<lower=0,upper=10> sigma_1;
  real<lower=0,upper=10> sigma_2;
  real<lower=0,upper=10> sigma_y;
}
model {
  mu ~ normal(0, 10);
  gamma ~ normal(mu, sigma_2);
  for (j in 1:Nsys)
    theta[j] ~ normal(gamma[sys_grp[j]], sigma_1);
  for (i in 1:N)
    target += normal_lpdf(y[i] | theta[sid[i]], sigma_y);
  for (i in 1:Ncens)
    target += normal_lcdf(ycens[i] | theta[sidcens[i]], sigma_y);

}
"
fit1 <- stan_model(model_code=occ_stan1)
stan_in1 <- function(data=atr.FIN, chains=nchains){
    yAll <- data$Value
    detect <- data$Detect
    y <- yAll[detect==1]
    ycens <- data$Adj.Detect.Limit[detect==0]
    ycens[is.na(ycens)] <- 0.1
    n <- length(y)
    ncens <- length(ycens)
    sysAll <- as.numeric(ordered(data$PWSID))
    sid <- sysAll[detect==1]
    sidcens <- sysAll[detect==0]
    nsys <- max(sysAll)
    oo <- order(sysAll)
    sysAll <- sysAll[oo]
    grpAll <- as.numeric(ordered(data$group))[oo]
    ngrp <- max(grpAll)
    sys_grp <- grpAll[cumsum(table(sysAll))]
    stan_data <- list(N=n, Ncens=ncens, Ngrp=ngrp, Nsys=nsys,
                      y=log(y), ycens=log(ycens),
                      sys_grp=sys_grp, sid=sid, sidcens=sidcens)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(theta=rnorm(nsys,-5),
                                gamma=rnorm(ngrp,-5), mu=rnorm(1,-5),
                                sigma_1=runif(1),sigma_2=runif(1),
                                sigma_y=runif(1))
    stan_pars <- c("gamma", "mu", "sigma_1","sigma_2","sigma_y")
    return(list(data=stan_data, inits=stan_inits, pars=stan_pars,
                n.chains=chains))
}

input.to.stan <- stan_in1()
fit2keepF <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
#                     control=list(max_treedepth=25))

print(fit2keepF)

input.to.stan <- stan_in1(data=atr.RAW)
fit2keepR <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
#                     control=list(max_treedepth=25))
print(fit2keepR)

input.to.stan <- stan_in1(data=atr.OTH)
fit2keepO <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
#                     control=list(max_treedepth=25))

print(fit2keepO)

input.to.stan <- stan_in1(data=atrazinALL)
fit2keepA <- sampling(fit1, data=input.to.stan$data,
#                      init=input.to.stan$inits,
                      pars=input.to.stan$pars,
                      iter=niters,thin=nthin,
                      chains=input.to.stan$n.chains)#,
#                     control=list(max_treedepth=25))

print(fit2keepA)
save(fit2keepF, fit2keepR, fit2keepO, file="Occruns.RData")
load("Occruns.RData")

print(fit2keepR)
print(fit2keepF)
print(fit2keepO)

FINstan <- rstan::extract(fit2keepF)
FINtheta <- rvsims(FINstan$theta)
FINmu <- rvsims(FINstan$mu)
FINsigy <- rvsims(FINstan$sigma_y)
FINsig0 <- rvsims(FINstan$sigma_0)

RAWstan <- rstan::extract(fit2keepR)
RAWtheta <- rvsims(RAWstan$theta)
RAWmu <- rvsims(RAWstan$mu)
RAWsigy <- rvsims(RAWstan$sigma_y)
RAWsig0 <- rvsims(RAWstan$sigma_0)

OTHstan <- rstan::extract(fit2keepO)
OTHtheta <- rvsims(OTHstan$theta)
OTHmu <- rvsims(OTHstan$mu)
OTHsigy <- rvsims(OTHstan$sigma_y)
OTHsig0 <- rvsims(OTHstan$sigma_0)
