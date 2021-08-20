source("FrontMatter.R")

## Chapter 5
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
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## EUSE example -- non-nested
rtol2 <- read.csv(file=paste(dataDIR, "rtolforMS.csv", sep="/"), header=T)

## environmental data
euse.env <- read.csv(paste(dataDIR, "EUSE_NAT_ENV.csv", sep="/"), header=T)
names(euse.env)[13]<-"MAX.ELEV"
names(euse.env)[12]<-"MIN.ELEV"

AvePrec <- tapply(euse.env$AnnMeanP, euse.env$CITY, mean)
AveTemp <- tapply(euse.env$AnnMeanT, euse.env$CITY, mean)
AveElev <- tapply(euse.env$MEANELEV, euse.env$CITY, mean)
AveMaxT <- tapply(euse.env$AnnMaxT, euse.env$CITY, mean)
AveMaxP <- tapply(euse.env$AnnMaxP, euse.env$CITY, mean)
AveMinP <- tapply(euse.env$AnnMinP, euse.env$CITY, mean)
AvePdif <- tapply(euse.env$AnnMaxP-euse.env$AnnMinP, euse.env$CITY, mean)

city_ag<-read.csv(paste(dataDIR, "City_AG_Grassland.csv", sep="/"),
                  header=T, na.strings = ".")
city_ag<-cbind(city_ag[,1:2],city_ag[,3:5]/100)
city_ag[order(city_ag[ ,1]), ]
ag<-city_ag[order(city_ag[ ,1]), ][,5]
ag.cat <- as.numeric(ag>0.5)
site <- as.numeric(ordered(rtol2$city))
ag.full <- as.vector(ag.cat[site])
temp.full <- as.vector(AveTemp[site])
euse.lmer1 <- lmer(richtol ~ nuii + (1+nuii|city), data=rtol2)
summary(euse.lmer1)

euse.lmer <- lmer(richtol ~ nuii+temp.full+nuii:temp.full+
                       ag.full+ag.full:nuii+
                       ag.full:temp.full+
                       ag.full:temp.full:nuii+
                       (1+nuii|city), data=rtol2)
summary(euse.lmer)
## convergence failed

euse.lmer2 <-  lmer(richtol ~ nuii+temp.full+nuii:temp.full+
                      ag.full+ag.full:nuii +
                      ag.full:temp.full:nuii+
                      (1+nuii|site), data=rtol2)
summary(euse.lmer2)

euse_stan1 <- "
data {
  int<lower=0> N;
  int<lower=0> Nreg;
  int<lower=1> K;
  int<lower=1,upper=Nreg> region[N];
  matrix[N,K] x;
  vector[N] y;
}
parameters {
  vector[K] beta[Nreg];
  corr_matrix[K] Omega;
  row_vector[K] mu_b;
  vector<lower=0>[K] tau;
  real<lower=0,upper=10> sigma_y;
}
transformed parameters {
  vector[N] mu_y;
  for (i in 1:N)
    mu_y[i] = x[i] * beta[region[i]];
}
model {
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  beta ~ multi_normal(mu_b, quad_form_diag(Omega, tau));
  y ~ normal(mu_y, sigma_y);
}
"
fit1 <- stan_model(model_code=euse_stan1)
stan_in1 <- function(data=rtol2, chains=nchains){
    y <- data$richtol
    x <- as.matrix(cbind(1, data$nuii))
    n <- dim(data)[1]
    k <- dim(x)[2]
    reg <- as.numeric(ordered(data$city))
    nreg <- max(reg)
    stan_data <- list(N=n, Nreg=nreg, K=k, x=x, y=y, region=reg)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(beta=matrix(rnorm(nreg*k),
                                            nrow=nreg),
                                mu_b=rnorm(k), tau=runif(k),
                                sigma_y=runif(1))
    stan_pars <- c("beta","mu_b","tau","sigma_y")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in1()
fit2keep <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
#                     control=list(adapt_delta = 0.9))
print(fit2keep)
## EUSE no group-level predictor
euse_multilevel1 <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep, pars=c("beta","mu_b",
                   "tau", "sigma_y")))))
## Cholesky factorization
euse_stan15 <- "
data {
  int<lower=0> N;
  int<lower=0> Nreg;
  int<lower=1> K;
  int<lower=1,upper=Nreg> region[N];
  matrix[N,K] x;
  vector[N] y;
}
parameters {
  matrix[K,Nreg] z;
  cholesky_factor_corr[K] L_Omega;
  row_vector[K] mu_b;
  vector<lower=0>[K] tau;
  real<lower=0,upper=10> sigma_y;
}
transformed parameters {
  matrix[Nreg,K] beta;
//  vector<lower=0>[K] tau;
//  for (k in 1:K) tau[k]=2.5*tan(tau_unif[k]);
  beta = rep_matrix(mu_b, Nreg) +
               (diag_pre_multiply(tau,L_Omega)*z)';
}
model {
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  y ~ normal(rows_dot_product(beta[region], x), sigma_y);
}
"
fit15 <- stan_model(model_code=euse_stan15)
stan_in15 <- function(data=rtol2, chains=nchains){
    y <- data$richtol
    x <- as.matrix(cbind(1, data$nuii))
    n <- dim(data)[1]
    k <- dim(x)[2]
    reg <- as.numeric(ordered(data$city))
    nreg <- max(reg)
    stan_data <- list(N=n, Nreg=nreg, K=k, x=x, y=y, region=reg)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(mu_b=rnorm(k), tau=runif(k),
                                sigma_y=runif(1))
    stan_pars <- c("beta","mu_b","tau","sigma_y")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in15()
fit2keep <- sampling(fit15, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
##                     control=list(adapt_delta = 0.9))
print(fit2keep)
euse_multilevel15 <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep,
                   pars=c("beta","mu_b","tau","sigma_y")))))

## With group-level predictors
euse_stan2 <- "
data {
  int<lower=0> N;
  int<lower=0> Nreg;
  int<lower=1> K; // individual predictor
  int<lower=1> J; // group-level predictor
  int<lower=1,upper=Nreg> region[N];
  matrix[N,K] x;
  vector[N] y;
  row_vector[J] gr[Nreg];
}
parameters {
  vector[K] beta[Nreg];
  matrix[J,K] gamma;
  corr_matrix[K] Omega;
  vector<lower=0>[K] tau;
  real<lower=0,upper=10> sigma_y;
}
transformed parameters {
  vector[N] mu_y;
  row_vector[K] u_gamma[Nreg];
  for (j in 1:Nreg)
    u_gamma[j] = gr[j] * gamma;
  for (i in 1:N)
    mu_y[i] = x[i] * beta[region[i]];
}
model {
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  beta ~ multi_normal(u_gamma, quad_form_diag(Omega, tau));
  y ~ normal(mu_y, sigma_y);
}
"

fit2 <- stan_model(model_code=euse_stan2)

stan_in2 <- function(data=rtol2, Ag=ag.cat, Temp=AveTemp,
                     y="richtol", x="nuii", regn="city",
                     Log=F, int=F, chains=nchains){
    if (Log){
        y <- log(data[,y])
        x <- log(data[,x])
    } else {
        y <- data[,y]
        x <- data[,x]
    }
    x <- as.matrix(cbind(1, x))
    n <- dim(data)[1]
    k <- dim(x)[2]
    reg <- as.numeric(ordered(data[,regn]))
    nreg <- max(reg)
    if (int) regPred <- cbind(1, 2*(Ag-0.5), Temp-mean(Temp),
                              2*(Ag-0.5)*Temp-mean(Temp))
    else
        regPred <- cbind(1, 2*(Ag-0.5), Temp-mean(Temp))
    nregpred <- dim(regPred)[2]
    stan_data <- list(N=n, Nreg=nreg, K=k, J=nregpred, x=x, y=y,
                      region=reg, gr = regPred)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(beta=matrix(rnorm(nreg*k),
                                            ncol=k),
                                gamma=matrix(rnorm(nregpred*k),
                                             ncol=k),
                                tau=runif(k), sigma_y=runif(1))
    stan_pars <- c("beta","gamma","tau","sigma_y")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in2()
fit2keep <- sampling(fit2, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta = 0.95,
                                  max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("gamma","tau"))

euse_multilevel2 <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep, pars=c("beta","gamma",
                   "tau", "sigma_y")))))

## Cholesky decomposition
euse_stan3 <- "
data {
  int<lower=0> N;
  int<lower=0> Nreg;
  int<lower=1> K; // individual predictor
  int<lower=1> J; // group-level predictor
  int<lower=1,upper=Nreg> region[N];
  matrix[N,K] x;
  vector[N] y;
  matrix[Nreg,J] gr;
}
parameters {
  matrix[K,Nreg] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0,upper=pi()/2>[K] tau_unif;
  matrix[J,K] gamma;
  real<lower=0,upper=10> sigma_y;
}
transformed parameters {
  matrix[Nreg,K] beta;
  vector<lower=0>[K] tau;

  for (k in 1:K)
    tau[k] = 2.5*tan(tau_unif[k]);
  beta = gr*gamma + (diag_pre_multiply(tau,L_Omega)*z)';
}
model {
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0,5);
  y ~ normal(rows_dot_product(beta[region], x), sigma_y);
}
"
fit3 <- stan_model(model_code=euse_stan3)
stan_in3 <- function(data=rtol2, Ag=ag.cat, Temp=AveTemp,
                     y="richtol", x="nuii", regn="city",
                     Log=F, int=F, chains=nchains){
    if (Log){
        y <- log(data[,y])
        x <- log(data[,x])
    } else {
        y <- data[,y]
        x <- data[,x]
    }
    x <- as.matrix(cbind(1, x))
    n <- dim(data)[1]
    k <- dim(x)[2]
    reg <- as.numeric(ordered(data[,regn]))
    nreg <- max(reg)
    if (int)
        regPred <- cbind(1, 2*(Ag-0.5), Temp-mean(Temp),
                         2*(Ag-0.5)*(Temp-mean(Temp)))
    else regPred <- cbind(1, 2*(Ag-0.5), Temp-mean(Temp))
    nregpred <- dim(regPred)[2]
    stan_data <- list(N=n, Nreg=nreg, K=k, J=nregpred, x=x, y=y,
                      region=reg, gr = regPred)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(z=matrix(rnorm(k*nreg), nrow=k),
                                gamma=matrix(rnorm(nregpred*k),
                                             ncol=k),
                                sigma_y=runif(1))
    stan_pars <- c("beta","gamma","tau","sigma_y")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in3()
fit2keep <- sampling(fit3, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains,
                     control=list(adapt_delta = 0.9))

print(fit2keep)

euse_multilevel3 <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep, pars=c("beta","gamma",
                   "tau", "sigma_y")))))

save(euse_multilevel1, euse_multilevel2, euse_multilevel3,
     file="euse_mult.RData")

## figure
beta_rv <- summary(rvmatrix(euse_multilevel3[1:18], nrow=9))
gamma_rv <- summary(rvmatrix(euse_multilevel3[19:24], nrow=3))
gr <- input.to.stan$data$gr

tikz(file=paste(plotDIRch6, "euseMLM.tex", sep="/"),
     height=3.25, width=6, standAlone=T)
par(mfrow=c(1,2), mar=c(3, 3, 1, 0.25), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(AveTemp,beta_rv[1:9, 1], xlab="average temperature",
     ylab="regression intercept", ylim=c(4,7.25))
a <- sum(gamma_rv[1:2,1]*c(1,-1)) - gamma_rv[3,1]*mean(AveTemp)
abline(a, gamma_rv[3,1])
a <- sum(gamma_rv[1:2,1]*c(1,1)) - gamma_rv[3,1]*mean(AveTemp)
abline(a, gamma_rv[3,1])
segments(y0=beta_rv[1:9, 4], y1=beta_rv[1:9, 8], x0=AveTemp, x1=AveTemp)

plot(AveTemp,beta_rv[10:18, 1], xlab="average temperature",
     ylab="regression slope", ylim=c(-0.01, 0.05))
a <- sum(gamma_rv[4:5,1]*c(1,-1)) - gamma_rv[6,1]*mean(AveTemp)
abline(a, gamma_rv[6,1])
a <- sum(gamma_rv[4:5,1]*c(1,1)) - gamma_rv[6,1]*mean(AveTemp)
abline(a, gamma_rv[6,1])
segments(y0=beta_rv[10:18, 4], y1=beta_rv[10:18, 8],
         x0=AveTemp, x1=AveTemp)
dev.off()

## group level coef
gamma_rv[,1]
gamma_rv[,2]

for (i in 1:dim(gr)[1]){

}


gr%*%gamma_rv
## processing output in `RV`
## without group-level predictor -- comparing to `lmer` results

## fixed effects:
summary(euse_multilevel1[19:20])
fixef(euse.lmer1)
se.fixef(euse.lmer1)

## random effects:
summary(euse_multilevel1[1:9]-summary(euse_multilevel1[19])$mean)
summary(euse_multilevel1[10:18]-summary(euse_multilevel1[20])$mean)
ranef(euse.lmer1)
se.ranef(euse.lmer1)

## make a table

euse.lmer2 <- lmer(richtol ~ nuii+temp.full+nuii:temp.full+
                      ag.full + ag.full:nuii +
                      (1+nuii|site), data=rtol2)
summary(euse.lmer)

## With group-level predictor ag and temp

## Lakes in China (Tang et al, 2019) -- nested example
ChinaLakes <- read.csv(paste(dataDIR, "Tang_etal_2019.csv",
                             sep="/"))

grps <- as.numeric(ordered(substring(ChinaLakes$group, 2,2)))
lake_res <- as.numeric(ordered(ChinaLakes$type))-1
## nested
lake_res2 <- paste(grps, ChinaLakes$type)
oo <- order(grps)
##grps <- grps[oo]
##lake_res <- lake_res[oo]
Iron <- ChinaLakes$iron[oo]
Iron <- Iron[cumsum(table(grps[oo]))]

iron_mu <- mean(ChinaLakes$iron)
iron.c <- ChinaLakes$iron-iron_mu
lr.c <- lake_res-mean(lake_res)

chlk_lmer1 <- lmer(log(TP)~log(prec)+lr.c*iron.c+lr.c:log(prec)+
                      iron.c:log(prec)+log(prec):iron.c:lr.c+
                      (1+log(prec)|grps), data=ChinaLakes)
summary(chlk_lmer1)

chlk_lmer2 <- lmer(log(Chla)~log(TP)+lr.c*iron.c+lr.c:log(TP)+
                      iron.c:log(TP)+log(TP):iron.c:lr.c+
                      (1+log(TP)|grps), data=ChinaLakes)
summary(chlk_lmer2)

## nested region:lr
chlk_lmer3 <- lmer(log(Chla)~log(TP)+(1+log(TP)|grps)+(1+log(TP)|lake_res2),
                   data=ChinaLakes)
summary(chlk_lmer3)
ranef(chlk_lmer3)

nested_stan2 <- "
data {
  int<lower=0> N;
  int<lower=0> Nreg;
  int<lower=0> region[N];
  vector[N] x;
  vector[N] y;
  vector[N] typ;
  vector[Nreg] Temp;
}
parameters {
  matrix[2,Nreg] z;
  cholesky_factor_corr[2] L_Omega;
  vector<lower=0,upper=pi()/2>[2] tau_unif;
  real a[4];
  real b[4];
  real<lower=0,upper=10> sigma_y;
}
transformed parameters {
  vector[N] mu_y;
  vector[Nreg] be0;
  vector[Nreg] be0Ag;
  vector[Nreg] be1;
  vector[Nreg] be1Ag;
  matrix[Nreg,2] beta;
  vector<lower=0>[2] tau;

  for (k in 1:2)
    tau[k] = 2.5*tan(tau_unif[k]);
  beta = (diag_pre_multiply(tau,L_Omega)*z)';

  for (j in 1:Nreg){
    be0[j] = a[1] + a[2]*Temp[j];
    be0Ag[j] = a[3] + a[4]*Temp[j];
    be1[j] = b[1] + b[2]*Temp[j];
    be1Ag[j] = b[3] + b[4]*Temp[j];
  }
  for (i in 1:N)
    mu_y[i] = (beta[region[i],1]+be0[region[i]]+
               typ[i]*be0Ag[region[i]]) +
              (beta[region[i],2]+be1[region[i]]+
               typ[i]*be1Ag[region[i]])*x[i];
}
model {
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  y ~ normal(mu_y, sigma_y);
}
"

fit2n <- stan_model(model_code=nested_stan2)
stan_in2n <- function(data=ChinaLakes, LR="type", temp=Iron,
                      y="TP", x="prec", regn=grps, Log=T,
                      chains=nchains){
    if (Log){
        y <- log(data[,y])
        x <- log(data[,x])
    } else {
        y <- data[,y]
        x <- data[,x]
    }
    n <- dim(data)[1]
    typ <- as.numeric(ordered(data[,LR]))
    region <- regn
    nreg <- max(region)
    stan_data <- list(N=n, Nreg=nreg, x=x-mean(x),
                      y=y, region=region,
                      typ=2*(typ-1.5), Temp=temp-mean(temp))
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(beta=matrix(rnorm(nreg*2),
                                            nrow=nreg),
                                a=rnorm(4), b=rnorm(4),
                                sigma_y=runif(1))
    stan_pars <- c("a","b", "beta", "sigma_y")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in2n()
fit2keepM1 <- sampling(fit2n, data=input.to.stan$data,
                       init=input.to.stan$inits,
                       pars=input.to.stan$pars,
                       iter=niters,thin=nthin,
                       chains=input.to.stan$n.chains,
                       control=list(adapt_delta = 0.99,
                                    max_treedepth=15))

print(fit2keepM1)
input.to.stan <- stan_in2n(y="Chla", x="TP", chains=nchains)
fit2keepM2 <- sampling(fit2n, data=input.to.stan$data,
                       init=input.to.stan$inits,
                       pars=input.to.stan$pars,
                       iter=niters,thin=nthin,
                       chains=input.to.stan$n.chains,
                       control=list(adapt_delta = 0.99,
                                    max_treedepth=15))
print(fit2keepM2)

save(fit2keepM1, fit2keepM2, file="ChinaLkesStan.RData")

## model 1 TP-precipitation
CLR_m1 <- rvsims(as.matrix(
    as.data.frame(extract(fit2keepM1, pars=c("a","b", "beta", "sigma_y")))))
## figure
beta1_rv <- summary(rvmatrix(CLR_m1[9:22], nrow=7))
a1_rv <- summary(CLR_m1[1:4])
b1_rv <- summary(CLR_m1[5:8])

## coefficients
## retransforming the centered group-level predictor
be0L1 = (a1_rv[1,2]-a1_rv[3,2])+(a1_rv[2,2]-a1_rv[4,2])*mean(Iron) -
    (a1_rv[2,2]-a1_rv[4,2])*Iron
be0R1 = (a1_rv[1,2]-a1_rv[3,2])-(a1_rv[2,2]-a1_rv[4,2])*mean(Iron) +
    (a1_rv[2,2]-a1_rv[4,2])*Iron
be1L1 = (b1_rv[1,2]-b1_rv[3,2])+(b1_rv[2,2]-b1_rv[4,2])*mean(Iron) -
    (b1_rv[2,2]-b1_rv[4,2])*Iron
be1R1 = (b1_rv[1,2]-b1_rv[3,2])-(b1_rv[2,2]-b1_rv[4,2])*mean(Iron) +
    (b1_rv[2,2]-b1_rv[4,2])*Iron

## intercepts
plot(Iron, be0L1, ylim=c(3,6))
lines(Iron, be0R1)
## slopes
plot(Iron, be1L1, ylim=c(-0.5,0.5))
lines(Iron, be1R1)

## model 2 chla -- tp

CLR_m2 <- rvsims(as.matrix(
    as.data.frame(extract(fit2keepM2, pars=c("a","b", "beta", "sigma_y")))))
beta2_rv <- rvmatrix(CLR_m2[9:22], nrow=7)
a2_rv <- summary(CLR_m2[1:4])
b2_rv <- summary(CLR_m2[5:8])

be0L2 = (a2_rv[1,2]-a2_rv[3,2])+(a2_rv[2,2]-a2_rv[4,2])*mean(Iron) -
    (a2_rv[2,2]-a2_rv[4,2])*Iron
be0R2 = (a2_rv[1,2]+a2_rv[3,2])-(a2_rv[2,2]-a2_rv[4,2])*mean(Iron) +
    (a2_rv[2,2]-a2_rv[4,2])*Iron
be1L2 = (b2_rv[1,2]-b2_rv[3,2])+(b2_rv[2,2]-b2_rv[4,2])*mean(Iron) -
    (b2_rv[2,2]-b2_rv[4,2])*Iron
be1R2 = (b2_rv[1,2]-b2_rv[3,2])-(b2_rv[2,2]-b2_rv[4,2])*mean(Iron) +
    (b2_rv[2,2]-b2_rv[4,2])*Iron


## intercepts

tikz(file=paste(plotDIRch6, "TangetalMLM.tex", sep="/"),
     height=3.25, width=6, standAlone=T)
par(mfrow=c(1,2), mar=c(3, 3, 1, 0.25), mgp=c(1.5,0.125,0), las=1, tck=0.01)

be0Lerr <- summary(be0L2 + beta2_rv[1:7])
be0Rerr <- summary(be0R2 + beta2_rv[1:7])

plot(Iron, be0L2, ylim=c(1.5,3.5), type="l",lty=2,
     xlab="mean soil iron (\\%)", ylab="regression intercept")
lines(Iron, be0R2)
segments(x0=Iron+0.01, x1=Iron+0.01, y0=be0Lerr[,4], y1=be0Lerr[,8], lty=2)
segments(x0=Iron-0.01, x1=Iron-0.01, y0=be0Rerr[,4], y1=be0Rerr[,8])
points(x=Iron, y=be0Lerr[,1])
points(x=Iron, y=be0Rerr[,1], pch=16)


## slopes
be1Lerr <- summary(be1L2 + beta2_rv[8:14])
be1Rerr <- summary(be1R2 + beta2_rv[8:14])
plot(Iron, be1L2, ylim=c(-1,2), type="l", lty=2,
          xlab="mean soil iron (\\%)", ylab="regression slope")
lines(Iron, be1R2)
segments(x0=Iron+0.01, x1=Iron+0.01, y0=be1Lerr[,4], y1=be1Lerr[,8], lty=2)
segments(x0=Iron-0.01, x1=Iron-0.01, y0=be1Rerr[,4], y1=be1Rerr[,8])
points(x=Iron, y=be1Lerr[,1])
points(x=Iron, y=be1Rerr[,1], pch=16)
dev.off()
