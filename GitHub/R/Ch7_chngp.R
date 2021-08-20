source("FrontMatter.R")

## Chapter 7
plotDIRch7 <- paste(plotDIR, "chapter7", "figures", sep="/")
packages(tidyverse)
## simulation
packages(rv)
packages(rstan)
packages(readxl)
packages(car)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Nuese River TP time series -- change point example
### Reading and processing data

upper.neuse <- read.table(paste(dataDIR,
                                "Data_upp_20080509_220930_RegResults.txt",
                                sep="/"),
                          sep="~", header=T)

J417 <- upper.neuse[upper.neuse$Station.ID=="J4170000       " , ]

J417.TP <- J417[J417$Character=="Phosphorus as P",]
J417.TP$Date <- as.Date(J417.TP$Activity.Start,
                        format="%Y-%m-%d %H:%M:%S")
names(J417.TP)[24] <- "Value"

plot(Value ~ Date, data=J417.TP)

## calculating monthly means
TP<- rep(NA, 12*(2001-1970))
k <- 0
for (i in 1971:2001){ ## year
  for (j in 1:12){ ## month
     k <- k+1
     temp <- as.numeric(format(J417.TP$Date, "%m"))==j &
         as.numeric(format(J417.TP$Date, "%Y"))==i
     if (sum(temp)>0) TP[k] <- mean(J417.TP$Value[temp])
  }
}

## Time series plot
TP.ts <- ts(TP, start=c(1971,1), end=c(2001,12), freq=12)
tikz(file=paste(plotDIRch7, "neuseTPchngp.tex", sep="/"),
     height=3, width=4.75, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01,las=1)
plot(TP.ts, type="p", ylab="TP Concentration (mg/L)")
abline(v=1988)
dev.off()

### US Lilac first bloom dates
## reading and processing data
USLilac <- read.csv(paste(dataDIR, "NAmlilac.csv", sep="/"))

USLilac$type <- "Syringa chinensis clone"
USLilac$type[USLilac$Ptype==2] <-"Syringa vulgaris"

USLilac$FirstLeaf[USLilac$FirstLeaf==999] <- NA
USLilac$FirstBloom[USLilac$FirstBloom==999] <- NA

## keep stations with at least 30 years of data

keep <- !is.na(match(USLilac$STID,
                     as.numeric(names(table(USLilac$STID))
                                [table(USLilac$STID)>=30])))
uslilacs <- USLilac[keep,]

temp <- USLilac[USLilac$STID==354147,]
lilactemp <- loess.smooth(y=temp$FirstBloom, x=temp$Year,
                          degree=1, span=0.75)

tikz(file=paste(plotDIRch7, "lilacFB.tex", sep="/"),
     height=2.5, width=4.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.75,0.125,0), tck=0.01,las=1)
plot(FirstBloom~Year, data=temp)
lines(lilactemp$y~lilactemp$x)
dev.off()


### Everglade dosing study -- Utricularia (bladderwort) stem density
utric <- read.csv(paste(dataDIR, "UtriGM.csv", sep="/"))

tikz(file=paste(plotDIRch7, "UtricChngp.tex", sep="/"),
     height=2.5, width=4.5, standAlone=T)
par(mar=c(3,3,1,1), mgp=c(1.75,0.125,0), tck=0.01,las=1)
plot(Utot~GM.UTP, data=utric, subset=DATE=="25-Aug-95",
     xlab="TP Concentration ($\\mu$g/L)", ylab="Stems per m$^2$")
dev.off()

#### Change point models
### Neuse River TP concentration -- Normal response
TPdata <- data.frame(TP=TP.ts, Month=1:length(TP.ts))
TPdata$mn <- as.numeric(time(TP.ts))
TPdata <- TPdata[!is.na(TPdata[,1]),]

neuseTP_step <- "
data{
  int N;
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real beta0;
  real delta;
  real<lower=minX, upper=maxX> phi;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0 + delta*step(x[i]-phi);
}
model{
  beta0 ~ normal(mu0, 10);
  delta ~ normal(0, 5);
  sigma ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
"
fit1 <- stan_model(model_code=neuseTP_step)

neuseStep_input <- function(data=TPdata, y=1, x=2, Log=T,
                            n.chains=nchains){
    n <- dim(data)[1]
    if (Log)
        y <- log(data[,y])
    else
        y <- data[,y]
    x <- data[,x]
    standata <- list(N=n, y=y, x=x, minX=min(x), maxX=max(x),
                     mu0=mean(y[x<1988]))
    staninits <- list()
    for (i in 1:n.chains)
        staninits[[i]] <- list(beta0=rnorm(1, 1, 1),
                               delta=rnorm(1),
                               phi=min(x)+runif(1)*(max(x)-min(x)),
                               sigma=runif(1))
    stanpara=c("beta0","delta","sigma","phi")
    return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains))
}

input.to.stan <- neuseStep_input(x="mn", Log=T)

fit2keep <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
                     ##control=list(max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("beta0","delta","phi", "sigma"))

neuseStan <- rstan::extract(fit2keep,
                            pars=c("beta0","delta","phi","sigma"),
                            permute=T)
save(neuseStan, file="neuseStep.RData")
hist(neuseStan$phi)

## using four-parameter logit function
neuseTP_fpl <- "
data{
  int N;
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real A;
  real<upper=0> delta;
  real<lower=minX, upper=maxX> phi;
  real scal;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = A + delta/(1+exp((phi-x[i])/scal));
}
model{
  A ~ normal(mu0, 10);
  delta ~ normal(0, 5);
  sigma ~ normal(0, 5);
  scal ~ normal(0,2.5);
  y ~ normal(mu, sigma);
}
"
fitfpl <- stan_model(model_code=neuseTP_fpl)

neuseFPL_input <- function(data=TPdata, y=1, x=2, Log=T,
                           n.chains=nchains){
    n <- dim(data)[1]
    if (Log)
        y <- log(data[,y])
    else
        y <- data[,y]
    x <- data[,x]
    standata <- list(N=n, y=y, x=x, minX=min(x), maxX=max(x),
                     mu0=mean(y[x<1988]))
    staninits <- list()
    for (i in 1:n.chains)
        staninits[[i]] <- list(A=rnorm(1, 1, 1), delta=-runif(1),
                               phi=min(x)+runif(1)*(max(x)-min(x)),
                               scal=runif(1), sigma=runif(1))
    stanpara=c("A","delta","phi","scal","sigma")
    return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains))
}

input.to.stan <- neuseFPL_input(x="mn", Log=T)

fit2keep <- sampling(fitfpl, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("A", "delta", "phi", "scal", "sigma"))

neuseStanFPL <- rstan::extract(fit2keep,
                               pars=c("A","delta","phi","scal","sigma"),
                               permute=T)
save(neuseStanFPL, file="neuseFPL.RData")
hist(neuseStanFPL$phi)


## plotting the fitted model
fplsumm_rv <- rvsims(as.matrix(as.data.frame(neuseStanFPL)))
xx <- seq(1972,2002,,100)
fplpred <- fplsumm_rv[1] +
    fplsumm_rv[2]/(1+exp((fplsumm_rv[3]-xx)/fplsumm_rv[4]))
fplsumm <- summary(fplpred)

tikz(file=paste(plotDIRch7, "NeuseChngpFPL.tex", sep="/"),
     height=2.5, width=4, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(log(TP.ts), type="n", xlab="Year",
     ylab="TP Concentration (mg/L)", axes=F)
axis(1)
axis(2, at=log(c(0.05,0.1,0.25,0.5,1,2)),
     label=c(0.05,0.1,0.25,0.5,1,2))
box()
polygon(x=c(xx, rev(xx)), y=c(fplsumm[,4], rev(fplsumm[,8])),
        border=grey(0.5), col=grey(0.5))
polygon(x=c(xx, rev(xx)), y=c(fplsumm[,5], rev(fplsumm[,7])),
        border=grey(0.2), col=grey(0.2))
lines(xx, fplsumm$mean, col="white")
points(log(TP.ts), cex=0.5)
dev.off()

#### Utricularia stem count -- Poisson response
utric_fpl <- "
data{
  int N;
  int y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real beta0;
  real<upper=0> delta;
  real<lower=minX, upper=maxX> phi;
  real scal;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0 + delta/(1+exp((phi-x[i])/scal));
}
model{
  beta0 ~ normal(mu0, 10);
  delta ~ normal(0, 5);
  scal ~ normal(0, 2.5);
  y ~ poisson(exp(mu));
}
"

fit2 <- stan_model(model_code=utric_fpl)

utricFPL_input <- function(data=TPdata, y=1, x=2, Log=T,
                           n.chains=nchains){
    n <- dim(data)[1]
    if (Log)
        x <- log(data[,x])
    else
        x <- data[,x]
    y <- data[,y]
    standata <- list(N=n, y=y, x=x, minX=min(x), maxX=max(x),
                     mu0=log(mean(y[x<20], na.rm=T)))
    staninits <- list()
    for (i in 1:n.chains)
        staninits[[i]] <- list(beta0=rnorm(1), delta=-runif(1),
                               phi=min(x)+runif(1)*(max(x)-min(x)),
                               scal=runif(1))
    stanpara=c("beta0","delta","phi", "scal")
    return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains))
}

input.to.stan <- utricFPL_input(data=utric[utric$DATE=="25-Aug-95",],
                                 y="Utot", x="GM.UTP", Log=F)

fit2keep <- sampling(fit2, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("beta0","delta","phi", "scal"))
plot(input.to.stan$data$x, input.to.stan$data$y)
utricStan <- rstan::extract(fit2keep,
                            pars=c("beta0","delta","phi", "scal"),
                            permute=T)

save(utricStan, file="utricStan.RData")

## plotting the fitted model

plotdata <- utric[utric$DATE=="25-Aug-95",]
fplutric_rv <- rvsims(as.matrix(as.data.frame(utricStan)))
xx <- seq(5,65,,100)
fplpred <- fplutric_rv[1] +
    fplutric_rv[2]/(1+exp((fplutric_rv[3]-xx)/fplutric_rv[4]))
fplsumm <- summary(exp(fplpred))

tikz(file=paste(plotDIRch7, "UtricChngpFPL.tex", sep="/"),
     height=2.5, width=4, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(Utot~GM.UTP, data=plotdata, type="n", xlab="TP ($\\mu$g/L)",
     ylab="Total Stem Count", axes=F)
axis(1)
axis(2)
box()
polygon(x=c(xx, rev(xx)), y=c(fplsumm[,4], rev(fplsumm[,8])),
        border=grey(0.5), col=grey(0.5))
polygon(x=c(xx, rev(xx)), y=c(fplsumm[,5], rev(fplsumm[,7])),
        border=grey(0.2), col=grey(0.2))
lines(xx, fplsumm$mean, col="white")
points(Utot~GM.UTP, data=plotdata, cex=0.5)
dev.off()

### Plot in log-scale:
fplsumm_log <- summary(fplpred)

tikz(file=paste(plotDIRch7, "UtricChngpFPL_log.tex", sep="/"),
     height=2.5, width=4, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(log(Utot+0.01)~GM.UTP, data=plotdata, type="n",
     xlab="TP ($\\mu$g/L)",
     ylab="Log Total Stem Count", ylim=c(-6,3),axes=F)
axis(1)
axis(2)
box()
##polygon(x=c(xx, rev(xx)), y=c(fplsumm_log[,4], rev(fplsumm_log[,8])),
##        border=grey(0.5), col=grey(0.5))
polygon(x=c(xx, rev(xx)), y=c(fplsumm_log[,5], rev(fplsumm_log[,7])),
        border=grey(0.7), col=grey(0.7))
lines(xx, fplsumm_log$mean, col="white")
points(log(Utot+0.1)~GM.UTP, data=plotdata, cex=0.5)
abline(v=summary(fplutric_rv)[3,2])
dev.off()


### Utricularia hierarchical model
utric_fplHier <- "
data{
  int N;
  int<lower=2> Ngrp;
  int<lower=1,upper=Ngrp> grps[N];
  int y[N];
  real x[N];
  real minX;
  real maxX;
  real mub0;
  real phiInit;
}
parameters{
  matrix[4,Ngrp] z;
  cholesky_factor_corr[4] L_Omega;
  vector<lower=0,upper=pi()/2>[4] tau_unif;
  row_vector[4] mu_b;
}
transformed parameters{
  matrix[Ngrp,4] beta;
  real mu[N];
  vector<lower=0>[4] tau;
  for (i in 1:4)
    tau[i] = 2.5*tan(tau_unif[i]);
  beta = rep_matrix(mu_b, Ngrp)+
          (diag_pre_multiply(tau,L_Omega)*z)';
  for (i in 1:N)
    mu[i] = beta[grps[i],1] +
            beta[grps[i],2]/(1+exp((beta[grps[i],3]-x[i])/beta[grps[i],4]));
}
model{
  mu_b[1] ~ normal(mub0, 1.5);
  mu_b[2] ~ normal(0, 2.5);
  mu_b[3] ~ normal(phiInit, 2.5);
  mu_b[4] ~ normal(0, 1);
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  y ~ poisson(exp(mu));
}
"
fitHier_utric <- stan_model(model_code=utric_fplHier)

## without multivariate normal hyper-distribution
utric_fplHier2 <- "
data{
  int N;
  int<lower=2> Ngrp;
  int<lower=1,upper=Ngrp> grps[N];
  int y[N];
  real x[N];
  real minX;
  real maxX;
  real mub0;
  real phiInit;
}
parameters{
  real beta0[Ngrp];
  real<lower=minX,upper=maxX> phi[Ngrp];
  real<upper=0> delta[Ngrp];
  real<lower=0> scal[Ngrp];
  real mu_b0;
  real mu_phi;
  real<upper=0> mu_delta;
  real<lower=0> mu_scal;
  real<lower=0> sigma_b0;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_scal;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0[grps[i]] +
            delta[grps[i]]/(1+exp((phi[grps[i]]-x[i])/scal[grps[i]]));
}
model{
  mu_b0 ~ normal(mub0, 5);
  mu_delta ~ normal(0, 5);
  mu_phi ~ normal(phiInit, 5);
  mu_scal ~ normal(0, 2.5);
  beta0 ~ normal(mu_b0, sigma_b0);
  delta ~ normal(mu_delta, sigma_delta);
  phi ~ normal(mu_phi, sigma_phi);
  scal ~ normal(mu_scal, sigma_scal);
  y ~ poisson(exp(mu));
}
"

fitHier_utric2 <- stan_model(model_code=utric_fplHier2)

temp <- utric$DATE=="25-Aug-95" |
    utric$DATE=="27-Mar-96" |
    utric$DATE=="16-Apr-98" |
    utric$DATE=="4-Aug-98"
utric_hier <- function(data=utric[temp,], y="Utot",
                        x="GM.UTP", n.chains=nchains){
    y <- data[,y]
    nas <- is.na(y)
    data <- data[!nas,]
    n <- dim(data)[1]
    x <- data[,x]
    y <- y[!nas]
    grps <- as.numeric(ordered(data$DATE))
    ngrps <- max(grps)
    standata <- list(N=n, Ngrp=ngrps, y=y, x=x,
                     minX=min(x), maxX=max(x),
                     grps=grps, mub0=log(20),
                     phiInit=20)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(mu_b=c(rnorm(1,log(20)),
                                    -runif(1),
                                    runif(1, min(x),max(x)),
                                    runif(1)),
                             tau=runif(4))
    }
    stanpara=c("beta","mu_b","tau")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

input.to.stan <- utric_hier()

fit2keep <- sampling(fitHier_utric, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=20,
                                  adapt_delta=0.95))

print(fit2keep)
pairs(fit2keep, pars=c("mu_b"))
utricStanHier <- rstan::extract(fit2keep, permute=T,
                                  pars=c("beta","mu_b","tau"))
save(utricStanHier, file="utricStanHier.RData")
HierUtric_rv <- rvsims(as.matrix(as.data.frame(utricStanHier)))

utric_hier2 <- function(data=utric[temp,], y="Utot",
                        x="GM.UTP", n.chains=nchains){
    y <- data[,y]
    nas <- is.na(y)
    data <- data[!nas,]
    n <- dim(data)[1]
    x <- data[,x]
    y <- y[!nas]
    grps <- as.numeric(ordered(data$DATE))
    ngrps <- max(grps)
    standata <- list(N=n, Ngrp=ngrps, y=y, x=x,
                     minX=min(x), maxX=max(x),
                     grps=grps, mub0=log(20),
                     phiInit=20)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(beta0=rnorm(ngrps, log(20)),
                             phi=runif(ngrps, min(x),max(x)),
                             delta= -runif(ngrps), scal=runif(ngrps),
                             sigma_b0=runif(1), sigma_phi=runif(1),
                             sigma_delta=runif(1), sigma_scal=runif(1),
                             mu_b0=rnorm(1, log(20)),
                             mu_phi=runif(1,min(x),max(x)),
                             mu_delta= -runif(1), mu_scal=runif(1))
    }
    stanpara=c("beta0","phi","delta","scal",
               "mu_b0","mu_delta","mu_phi","mu_scal",
               "sigma_b0","sigma_delta","sigma_phi","sigma_scal")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

input.to.stan <- utric_hier2()

fit2keep <- sampling(fitHier_utric2,
                     data=input.to.stan$data,
                     ##init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters*10,thin=nthin*10,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=25,
                                  adapt_delta=0.95))

print(fit2keep)
pairs(fit2keep, pars=c("mu_b0","mu_delta","mu_phi","mu_scal"))

### Lilac first bloom date -- piecewise linear model
lilac_hockey <- "
data{
  int N;
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real beta0;
  real beta1;
  real<upper=0> delta;
  real<lower=minX, upper=maxX> phi;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0 + (beta1+delta*int_step(x[i]-phi))*(x[i]-phi);
}
model{
  beta0 ~ normal(mu0, 5);
  beta1 ~ normal(0, 0.125);
  delta ~ normal(0, 2);
  sigma ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
"
fit3 <- stan_model(model_code=lilac_hockey)

## beta1=0
lilac_hockey2 <- "
data{
  int N;
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real beta0;
  real<upper=0> delta;
  real<lower=minX, upper=maxX> phi;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0 + delta*int_step(x[i]-phi)*(x[i]-phi);
}
model{
  beta0 ~ normal(mu0, 5);
  delta ~ normal(0, 2);
  sigma ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
"
fit4 <- stan_model(model_code=lilac_hockey2)
lilac_input <- function(data=temp, y="FirstBloom",
                        x="Year", n.chains=nchains,
                        beta1=T){
    n <- dim(data)[1]
    x <- data[,x]
    y <- data[,y]
    standata <- list(N=n, y=y, x=x, minX=min(x),
                     maxX=max(x), mu0=115)
    staninits <- list()
    for (i in 1:n.chains){
        if (beta1)
            staninits[[i]] <- list(beta0=rnorm(1, 115),
                                   beta1=rnorm(1),
                                   delta=-abs(rnorm(1)),
                                   sigma=runif(1),
                                   phi=min(x)+runif(1)*(max(x)-min(x)))
        else
            staninits[[i]] <- list(beta0=rnorm(1, 115),
                                   delta=-abs(rnorm(1)),
                                   sigma=runif(1),
                                   phi=min(x)+runif(1)*(max(x)-min(x)))
    }
    if (beta1)
        stanpara=c("beta0","beta1","delta", "phi","sigma")
    else
        stanpara=c("beta0","delta", "phi","sigma")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

hockey_smooth <- "
data{
  int N;
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mu0;
}
parameters{
  real beta0;
  real beta1;
  real<upper=0> delta;
  real<lower=0> scal;
  real<lower=minX, upper=maxX> xmid;
  real<lower=0> sigma;
}
transformed parameters{
  real mu[N];
  for (i in 1:N)
    mu[i] = beta0 + beta1 *(x[i]-xmid) +
            delta*scal*log1p_exp((x[i]-xmid)/scal);
}
model{
  beta0 ~ normal(mu0, 5);
  beta1 ~ normal(0, 0.5);
  scal ~ normal(0, 1);
  delta ~ normal(0, 2);
  sigma ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
"
fit5 <- stan_model(model_code=hockey_smooth)
smooth_input <- function(data=temp, y="FirstBloom", x="Year", n.chains=nchains,
                         beta1=T){
    n <- dim(data)[1]
    x <- data[,x]
    y <- data[,y]
    standata <- list(N=n, y=y, x=x, minX=min(x),
                     maxX=max(x), mu0=115)
    staninits <- list()
    for (i in 1:n.chains)
        staninits[[i]] <- list(beta0=rnorm(1, 115),
                               beta1=rnorm(1),
                               delta=-abs(rnorm(1)),
                               scal=runif(1),
                               sigma=runif(1),
                               phi=min(x)+runif(1)*(max(x)-min(x)))
    stanpara=c("beta0","beta1","delta", "xmid",
               "scal","sigma")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

input.to.stan <- lilac_input()

fit2keep <- sampling(fit3, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
##                     control=list(max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("beta0","beta1","delta","phi", "sigma"))

lilacStan<-rstan::extract(fit2keep,
                          pars=c("beta0","beta1",
                                 "delta","phi","sigma"),
                          permute=T)

input.to.stan <- lilac_input(beta1=F)

fit2keep <- sampling(fit4, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
##                     control=list(max_treedepth=15))

lilacStan2 <-rstan::extract(fit2keep,
                          pars=c("beta0","delta","phi","sigma"),
                          permute=T)


save(lilacStan, lilacStan2, file="lilacStan.RData")
load("lilacStan.RData")
input.to.stan <- smooth_input()

fit2keep <- sampling(fit5, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
##                     control=list(max_treedepth=15))

print(fit2keep)
pairs(fit2keep, pars=c("beta0","beta1","delta",
                       "xmid","scal","sigma"))
lilacStanSM<-rstan::extract(fit2keep,
                           pars=c("beta0","beta1",
                                  "delta","xmid","scal",
                                  "sigma"),
                           permute=T)

save(lilacStanSM, file="lilacStanSmooth.RData")
write.table(rbind(
    summary(rvsims(as.matrix(as.data.frame(lilacStan))))[,c(1:2,5:9)],
    summary(rvsims(as.matrix(as.data.frame(lilacStanSM))))[,c(1:2,5:9)]),
    file="lilacHKvSM.txt", sep="&")

### Hierarchical Hockey Stick model
#### using hockey2 -- beta1 = 0

lilac_hockeyHier <- "
data{
  int N;
  int<lower=2> Ngrp;
  int<lower=1,upper=Ngrp> grps[N];
  real y[N];
  real x[N];
  real minX;
  real maxX;
  real mub0;
  real phiInit;
}
parameters{
  matrix[3,Ngrp] z;
  cholesky_factor_corr[3] L_Omega;
  vector<lower=0,upper=pi()/2>[3] tau_unif;
  row_vector[3] mu_b;
  real<lower=0> sigma;
}
transformed parameters{
  matrix[Ngrp,3] beta;
  real mu[N];
  vector<lower=0>[3] tau;
  for (i in 1:3) tau[i] = 2.5*tan(tau_unif[i]);
  beta = rep_matrix(mu_b, Ngrp)+
          (diag_pre_multiply(tau,L_Omega)*z)';
  for (i in 1:N)
    mu[i] = beta[grps[i],1] +
            beta[grps[i],2]*int_step(x[i]-beta[grps[i],3])*
            (x[i]-beta[grps[i],3]);
}
model{
  mu_b[1] ~ normal(mub0,20);
  mu_b[2] ~ normal(0,5);
  mu_b[3] ~ normal(phiInit, 20);
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  sigma ~ normal(0, 5);
  y ~ normal(mu, sigma);
}
"
fitHier <- stan_model(model_code=lilac_hockeyHier)

lilac_hier <- function(data=uslilacs, y="FirstBloom",
                        x="Year", n.chains=nchains){
    y <- data[,y]
    nas <- is.na(y)
    data <- data[!nas,]
    n <- dim(data)[1]
    x <- data[,x]
    y <- y[!nas]
    grps <- as.numeric(ordered(data$STID))
    ngrps <- max(grps)
    standata <- list(N=n, Ngrp=ngrps, y=y, x=x,
                     minX=min(x), maxX=max(x),
                     grps=grps, mub0=115, phiInit=1980)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(mu_b=c(rnorm(1,115),
                                    -runif(1),
                                    runif(1, min(x),max(x))),
                             tau=runif(3),
                             sigma=runif(1))
    }
    stanpara=c("beta","mu_b","tau","sigma")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

input.to.stan <- lilac_hier()

fit2keep <- sampling(fitHier, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=15,
                                  adapt_delta=0.9))

print(fit2keep)
pairs(fit2keep, pars=c("mu_b", "tau"))
lilacHier<-rstan::extract(fit2keep,
                           pars=c("beta","mu_b","tau","sigma"),
                           permute=T)

save(lilacHier, file="lilacStanHier.RData")
load("lilacStanHier.RData")
lilacHier_rv <- rvsims(as.matrix(as.data.frame(lilacHier)))

nsites <- max(as.numeric(ordered(uslilacs$STID)))
stIDcomp <- (1:53)[levels(ordered(uslilacs$STID))=="354147"]

lilacStan2_rv <- rvsims(as.matrix(as.data.frame(lilacStan2)))
lilacHier_rv[c(stIDcomp, stIDcomp+nsites, stIDcomp+nsites*2, nsites*3+(1:3))]

tikz(file=paste(plotDIRch7,"lilacHierComp.tex", sep="/"),
     standAlone=F, height=2.75, width=5.5)
par(mfrow=c(1,3), oma=c(1, 1, 1, 1), mar=c(1,1,1,1), mgp=c(1.75,0.125,0),
    las=1, tck= 0.01)
temp <- c(lilacStan2_rv[1],lilacHier_rv[c(stIDcomp, nsites*3+1)])
names(temp) <- c("site-only","site-hier","hyper")
mlplot(temp, main="$\\beta_0$", top.axis=F)
temp <- c(lilacStan2_rv[2],lilacHier_rv[c(stIDcomp+nsites, nsites*3+2)])
names(temp) <- c("","","")
mlplot(temp, main="$\\delta$", top.axis=F, xlim=c(-3,0))
temp <- c(lilacStan2_rv[3],lilacHier_rv[c(stIDcomp+nsites*2, nsites*3+3)])
names(temp) <- c("","","")
mlplot(temp, main="$\\phi$", top.axis=F, xlim=c(1965,1985))
dev.off()

## Threshold-changepoint model (Alameddine et al 2011)

## reading data
neuse_flowN <- read.csv(paste(dataDIR, "neuse_flow_data.csv", sep="/"))
head(neuse_flowN)
neuse_flowN$yr <- ordered(format(as.Date(neuse_flowN[,"datetime"],
                                         format="%m/%d/%Y"),
                                 format="%Y"))
neuse_flowN$Rdate <- as.Date(neuse_flowN$datetime, format="%m/%d/%Y")
neuse_flowN$year <- as.numeric(format(neuse_flowN$Rdate, format="%Y"))
## data plot
tikz(file=paste(plotDIRch7, "neuseChThdata.tex", sep="/"),
     height=3, width=3.5, standAlone=T)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(conc~I(flow/86400), log="xy",
     type="n", data=neuse_flowN,
     xlim=range(neuse_flowN$flow/86400),
     ylim=range(neuse_flowN$conc),
     xlab="Flow (m$^3$/s)", ylab="Concentration (mg/L)")
points(conc~I(flow/86400), data=neuse_flowN,
       subset=year>=2000, pch=16, cex=0.75)
points(conc~I(flow/86400), data=neuse_flowN,
       subset=year<2000,col=grey(0.5))
legend(x=750, y=2, pch=c(1, 16), #cex=c(0.75,1),
       legend=c("before","after"), bty="n")
dev.off()

neuse_FLWN <- "
data{
   int N;
   int Nyr;
   real y[N];
   real x[N];
   real wtrTemp[N];
   int yr[N];
   real phi_low;
   real phi_up;
   real max_yr;
}
parameters{
   real beta_temp;
   real<lower=phi_low,upper=phi_up> phi;
   real<lower=1,upper=max_yr> chngp;
   real<lower=0> sigma;
   real<lower=0> scal;
   real a[3];
   real b[3];
   real<lower=0> d[3];
}
transformed parameters{
   real mu_hat[N];
   real beta0[Nyr];
   real beta1[Nyr];
   real delta[Nyr];
   for (j in 1:Nyr){
     beta0[j] = a[1] + (b[1]-a[1])/(1+exp((chngp-j)/d[1]));
     beta1[j] = a[2] + (b[2]-a[2])/(1+exp((chngp-j)/d[2]));
     delta[j] = a[3] + (b[3]-a[3])/(1+exp((chngp-j)/d[3]));
   }
   for (i in 1:N){
     mu_hat[i] = beta0[yr[i]] +
                 beta1[yr[i]]*(x[i]-phi)+
                 delta[yr[i]]*scal*log1p_exp((x[i]-phi)/scal)+
                 beta_temp*wtrTemp[i];
   }
}
model{
   d ~ normal(0,1);
   sigma ~ cauchy(0, 2.5);
   y ~ normal (mu_hat, sigma);
}
"
fitPhi_Chngp <- stan_model(model_code=neuse_FLWN)

neuse_chngp <- function(data=neuse_flowN,  y="conc",
                        x="flow", n.chains=nchains){
    y <- log(data[,y])
    x <- log(data[,x])
    yr <- as.numeric(data$yr)
    n <- length(y)
    maxYr <- max(yr)
    lowF <- min(x)
    uppF <- max(x)

    standata <- list(N=n, y=y, x=x, yr=yr, Nyr=max(yr), wtrTemp=data$waterT,
                     phi_low=lowF, phi_up=uppF, max_yr=maxYr)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(a=runif(3),b=runif(3),d=runif(3),
                             sigma=runif(1),phi=runif(1,lowF,uppF),
                             chngp=runif(1,1,maxYr))
    }
    stanpara=c("beta0","beta1","delta","phi","chngp","sigma", "a","b","d")
    return(list(data=standata, inits=staninits,
                pars=stanpara, nchains=n.chains))
}

input.to.stan <- neuse_chngp ()

fit_neuseThCh <- sampling(fitPhi_Chngp, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=15,
                                  adapt_delta=0.95))

print(fit_neuseThCh)
pairs(fit_neuseThCh, pars=c("a","b","d"))

## Using Alameddine et al 2011 formulation:
##\begin{array}{rcl}
##y_i &= &(\beta_0+\delta_0*I(yr_i-\phi_{yr})) + \beta_TTemp +\\
##    &&  (\beta_1+\delta_1*I(yr_i-\phi_{yr}) +
##         \delta_q*I(\log(Q_i)-\phi_q)*I(yr_i-\phi_{yr}))*\\
##    &&   (\log(Q_i)-\phi_q*I(yr_i-\phi_{yr}))
##    \end{array}

neuse_step <- "
data{
   int N;
   int Nyr;
   real y[N];
   real x[N];
   real wtrTemp[N];
   real yr[N];
   real phi_low;
   real phi_up;
   real max_yr;
}
parameters{
   real<lower=phi_low,upper=phi_up> phi_q;
   real<lower=1,upper=max_yr> phi_yr;
   real<lower=0> sigma;
   real beta0;
   real<upper=0> beta1;
   real betaT;
   real<lower=0> delta0;
   real delta1;
   real<upper=0> deltaq;
}
transformed parameters{
   real mu_hat[N];
   real phiYR[N];
   real phiQ[N];
   for (i in 1:N){
     phiYR[i] = step(yr[i]-phi_yr);
     phiQ[i] = step(x[i]-phi_q);
     mu_hat[i] = (beta0+delta0*phiYR[i])+betaT*wtrTemp[i] +
                 (beta1+delta1*phiYR[i]+
                  deltaq*phiQ[i]*phiYR[i])*(x[i]-phi_q*phiYR[i]);
   }
}
model{
   beta0 ~ normal(2,1);
   beta1 ~ normal(-1,1);
   delta1 ~ normal(2,1);
   deltaq ~ normal(-2,1);
   sigma ~ cauchy(0, 2.5);
   y ~ normal (mu_hat, sigma);
}
"
fitphiphi_step <- stan_model(model_code=neuse_step)


neuse_step <- function(data=neuse_flowN,  y="conc",
                        x="flow", n.chains=nchains, centerX=F){
    y <- log(data[,y])
    x <- log(data[,x])
    wT <- data$waterT
    if (centerX){
        avgX <- mean(x)
        x <- x-avgX
        avgT <- mean(wT)
        wT <- data$waterT-avgT
    }
    yr <- as.numeric(data$yr)
    n <- length(y)
    maxYr <- max(yr)
    lowF <- min(x)
    uppF <- max(x)

    standata <- list(N=n, y=y, x=x, yr=yr, Nyr=max(yr), wtrTemp=wT,
                     phi_low=lowF, phi_up=uppF, max_yr=maxYr)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(beta0=rnorm(1), beta1=-runif(1), betaT=rnorm(1),
                             delta0=runif(1), delta1=rnorm(1), deltaq=-runif(1),
                             sigma=runif(1),phi_q=runif(1,lowF,uppF),
                             phi_yr=runif(1,1,maxYr))
    }
stanpara=c("beta0","beta1","betaT","delta0","delta1","deltaq",
           "phi_q","phi_yr","sigma")
    return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains, meanX=ifelse(centerX, c(avgX, avgT), NA)))
}

input.to.stan <- neuse_step ()
fit_neuseStep <- sampling(fitphiphi_step, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control=list(max_treedepth=15,
                                  adapt_delta=0.95))

## centering did not help
print(fit_neuseStep)
pairs(fit_neuseStep, pars=c("beta0","beta1","betaT","phi_yr","phi_q",
                            "delta0","delta1","deltaq"))

## highly correlated beta0 and beta1, delta1 and deltaq
## revising the model
## 1. centering the predictor (log(Q))
input.to.stan <- neuse_step (centerX=T)
fit_neuseStep_cX <- sampling(fitphiphi_step, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)#,
#                     control=list(max_treedepth=15,
#                                  adapt_delta=0.95))
print(fit_neuseStep_cX)
pairs(fit_neuseStep_cX, pars=c("beta0","beta1","betaT","phi_yr","phi_q",
                            "delta0","delta1","deltaq"))
save(fit_neuseStep_cX, file="neuseThCh_cX.RData")

## Using Alameddine et al 2011 formulation, alternative:
## \begin{array}{rcl}
##     y_i &= &(\alpha_0+\beta_TTemp +\alpha_1\log(Q_i))*I(yr_i-\phi_{yr})\\

##         &&   \left [ (\beta_0 +\beta_TTemp+\beta_1*(\delta+I(\log(Q_i)-\phi_q))*(\log(Q_i)-\phi_q) \right ]*(1-I(yr_i-\phi_{yr}))
##         \end{array}

## beta0 + beta1*(x[i]-phi)+
##         delta*scal*log1p_exp((x[i]-phi)/scal)+
##     beta_temp*wtrTemp[i];

Neuse_step2 <- "
data{
   int N;
   int Nyr;
   real y[N];
   real x[N];
   real wtrTemp[N];
   real yr[N];
   real phi_low;
   real phi_up;
   real max_yr;
   real<lower=0> scal;
}
parameters{
   real<lower=phi_low,upper=phi_up> phi_q;
   real<lower=1,upper=max_yr> phi_yr;
   real<lower=0> sigma;
   real alpha0;
   real beta0;
   real alpha1;
   real beta1;
   real betaT;
   real<upper=0> delta;
//   real delta1;
//   real<upper=0> deltaq;
}
transformed parameters{
   real mu_hat[N];
   real phiYR[N];
   real phiQ[N];
   real mu1[N];
   real mu2[N];
   for (i in 1:N){
     phiYR[i] = 1-step(yr[i]-phi_yr);
     mu1[i] = alpha0+betaT*wtrTemp[i]+alpha1*x[i];
     mu2[i] = beta0 + beta1*(x[i]-phi_q)+
         delta*scal*log1p_exp((x[i]-phi_q)/scal)+
         betaT*wtrTemp[i];
     mu_hat[i] = phiYR[i]*mu1[i]+(1-phiYR[i])*mu2[i];
   }
}
model{
   phi_q~normal(0,1);
   sigma ~ cauchy(0, 2.5);
   y ~ normal (mu_hat, sigma);
}
"
fitphiphi_step2 <- stan_model(model_code=Neuse_step2)


neuse_step2 <- function(data=neuse_flowN,  y="conc",
                        x="flow", n.chains=nchains, centerX=F){
    y <- log(data[,y])
    x <- log(data[,x])
    wT <- data$waterT
    if (centerX){
        avgX <- mean(x)
        x <- x-avgX
        avgT <- mean(wT)
        wT <- data$waterT-avgT
    }
    print(ifelse(centerX, c(avgX, avgT), NA))
    yr <- as.numeric(data$yr)
    n <- length(y)
    maxYr <- max(yr)
    lowF <- min(x)
    uppF <- max(x)

    standata <- list(N=n, y=y, x=x, yr=yr, Nyr=max(yr), wtrTemp=wT,
                     phi_low=lowF, phi_up=uppF, max_yr=maxYr,
                     scal=0.1*(uppF-lowF))
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(alpha0=rnorm(1), beta0=rnorm(1),
                             alpha1=rnorm(1), beta1=rnorm(1), betaT=rnorm(1),
                             delta=-runif(1), sigma=runif(1),
                             phi_q=runif(1,lowF,uppF),
                             phi_yr=runif(1,1,maxYr))
    }
    stanpara=c("alpha0", "alpha1",
               "beta0","beta1","betaT","delta",
               "phi_q","phi_yr","sigma")
    if (centerX)
        return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains, meanX=c(avgX, avgT)))
    else
        return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains))
}

## 1. centering the predictor (log(Q))
input.to.stan <- neuse_step2 (centerX=T)
fit_neuseStep_cX2 <- sampling(fitphiphi_step2,
                             data=input.to.stan$data,
                             init=input.to.stan$inits,
                             pars=input.to.stan$pars,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains)#,
#                     control=list(max_treedepth=15,
#                                  adapt_delta=0.95))
print(fit_neuseStep_cX2)
pairs(fit_neuseStep_cX2, pars=c("alpha0","beta0","alpha1","beta1","betaT",
                               "phi_yr","phi_q","delta"))
save(fit_neuseStep_cX2, file="neuseThCh_cX2.RData")

Neuse_step3 <- "
data{
   int N;
   int Nyr;
   real y[N];
   real x[N];
   real wtrTemp[N];
   real yr[N];
   real phi_low;
   real phi_up;
   real max_yr;
}
parameters{
   real<lower=phi_low,upper=phi_up> phi_q;
   real<lower=1,upper=max_yr> phi_yr;
   real<lower=0> sigma;
   real alpha0;
   real beta0;
   real alpha1;
   real beta1;
   real betaT;
   real<upper=0> delta;
   real<lower=0> scal;
}
transformed parameters{
   real mu_hat[N];
   real phiYR[N];
   real phiQ[N];
   real mu1[N];
   real mu2[N];
   for (i in 1:N){
     phiYR[i] = 1-step(yr[i]-phi_yr);
     mu1[i] = alpha0+betaT*wtrTemp[i]+alpha1*x[i];
     mu2[i] = beta0 + beta1*(x[i]-phi_q)+
         delta*scal*log1p_exp((x[i]-phi_q)/scal)+
         betaT*wtrTemp[i];
     mu_hat[i] = phiYR[i]*mu1[i]+(1-phiYR[i])*mu2[i];
   }
}
model{
   phi_q~normal(0,1);
   scal ~ normal(0,1);
   sigma ~ cauchy(0, 2.5);
   y ~ normal (mu_hat, sigma);
}
"
fitphiphi_step3 <- stan_model(model_code=Neuse_step3)


neuse_step3 <- function(data=neuse_flowN,  y="conc",
                        x="flow", n.chains=nchains, centerX=F){
    y <- log(data[,y])
    x <- log(data[,x])
    wT <- data$waterT
    if (centerX){
        avgX <- mean(x)
        x <- x-avgX
        avgT <- mean(wT)
        wT <- data$waterT-avgT
    }
    print(ifelse(centerX, c(avgX, avgT), NA))
    yr <- as.numeric(data$yr)
    n <- length(y)
    maxYr <- max(yr)
    lowF <- min(x)
    uppF <- max(x)

    standata <- list(N=n, y=y, x=x, yr=yr, Nyr=max(yr), wtrTemp=wT,
                     phi_low=lowF, phi_up=uppF, max_yr=maxYr)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(alpha0=rnorm(1), beta0=rnorm(1),
                             alpha1=rnorm(1), beta1=rnorm(1), betaT=rnorm(1),
                             delta=-runif(1), sigma=runif(1),
                             phi_q=runif(1,lowF,uppF), scal=runif(1),
                             phi_yr=runif(1,1,maxYr))
    }
    stanpara=c("alpha0", "alpha1",
               "beta0","beta1","betaT","delta",
               "phi_q","phi_yr","sigma", "scal")
    if (centerX)
        return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains, meanX=c(avgX, avgT)))
    else
        return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains))
}

## 1. centering the predictor (log(Q))
input.to.stan <- neuse_step3 (centerX=T)
fit_neuseStep_cX3 <- sampling(fitphiphi_step3,
                             data=input.to.stan$data,
                             init=input.to.stan$inits,
                             pars=input.to.stan$pars,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains)#,
#                     control=list(max_treedepth=15,
#                                  adapt_delta=0.95))
print(fit_neuseStep_cX3)
pairs(fit_neuseStep_cX3, pars=c("alpha0","beta0","alpha1","beta1","betaT",
                               "phi_yr","phi_q","delta","scal"))
save(fit_neuseStep_cX3, file="neuseThCh_cX3.RData")

## posterior plots
neuseThCh2 <- rvsims(as.matrix(as.data.frame(extract(fit_neuseStep_cX2))))
(neuseThCh2)

tmp <- sample(1:length(neuseThCh2$beta1), size=500)
tikz(paste(plotDIRch7,"neusePostPairs.tex", sep="/"),
     height=2, width=5., standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(neuseThCh2$beta1[tmp], neuseThCh2$delta[tmp],
     xlab="$\\beta_1$", ylab="$\\delta$")
plot(neuseThCh2$beta1[tmp], neuseThCh2$phi_q[tmp],
     xlab="$\\beta_1$", ylab="$\\phi_q$")
plot(neuseThCh2$phi_q[tmp], neuseThCh2$delta[tmp],
     xlab="$\\phi_q$",  ylab="$\\delta$")
dev.off()

stanin2 <- neuse_step2 (centerX=T)
stanin3 <- neuse_step3 (centerX=T)

### packages(Rmpfr) ## for function log1pexp
tmp <- range(stanin2$data$x)
X <- seq(tmp[1], tmp[2],,50)
plotX <- X+stanin2$meanX[1]-log(86400)

YR <- as.numeric(neuse_flowN$yr)
phiYR = 1-as.numeric(YR > neuseThCh2[8])
mu1 <- neuseThCh2[1]+##neuseThCh2[5]*stanin2$data$wtrTemp+
    neuseThCh2[2]*X
mu2 <- neuseThCh2[3]+neuseThCh2[4]*(X-neuseThCh2[7])+
    neuseThCh2[6]*stanin2$data$scal*
    log1p(exp((X-neuseThCh2[7])/stanin2$data$scal))#+
    ##neuseThCh2[5]*stanin2$data$wtrTemp
mu1_summ <- summary(mu1)
mu2_summ <- summary(mu2)

th <- summary(neuseThCh2[7]+stanin2$meanX[1]-log(86400))

phi_yr <- summary(phiYR)$mean==0
tikz(file=paste(plotDIRch7, "neuseChThmodel.tex", sep="/"),
     height=3, width=3.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(log(conc)~log(flow/86400), data=neuse_flowN, col="gray", type="n",
     xlab="log flow (m$^3$/s)", ylab="log concentration (mg/L)")
points(log(conc)~log(flow/86400), data=neuse_flowN, subset=!phi_yr, pch=16, col=gray(0.7), cex=0.5)
points(log(conc)~log(flow/86400), data=neuse_flowN, subset=phi_yr, pch=16, col=gray(0.3), cex=0.5)
polygon(x=c(plotX,rev(plotX)),
        y=c(mu1_summ$"2.5%", rev(mu1_summ$"97.5%")), col=gray(0.7), border=F)
polygon(x=c(plotX,rev(plotX)),
        y=c(mu2_summ$"2.5%", rev(mu2_summ$"97.5%")), col=gray(0.3), border=F)
lines(plotX,mu1_summ$mean)
lines(plotX,mu2_summ$mean)
segments(x0=th$"2.5%", x1=th$"97.5", y0=-4.5, y1=-4.5)
segments(x0=th$"25%", x1=th$"75", y0=-4.5, y1=-4.5, lwd=3)
points(x=th$mean, y=-4.5)
dev.off()

### the gamma model ###
Neuse_step4 <- "
data{
   int N;
   int Nyr;
   real y[N];
   real x[N];
   real wtrTemp[N];
   real yr[N];
   real max_yr;
}
parameters{
   real<lower=1,upper=max_yr> phi_yr;
   real<lower=0> sigma;
   real alpha0;
   real beta0;
   real alpha1;
   real<upper=0> beta1;
   real betaT;
   real<upper=-1> beta2;
}
transformed parameters{
   real mu_hat[N];
   real phiYR[N];
   real phiQ[N];
   real mu1[N];
   real mu2[N];
   for (i in 1:N){
     phiYR[i] = 1-step(yr[i]-phi_yr);
     mu1[i] = alpha0+betaT*wtrTemp[i]+alpha1*x[i];
     mu2[i] = beta0 + beta1*x[i]+ beta2*log(x[i])+betaT*wtrTemp[i];
     mu_hat[i] = phiYR[i]*mu1[i]+(1-phiYR[i])*mu2[i];
   }
}
model{
   sigma ~ cauchy(0, 2.5);
   y ~ normal (mu_hat, sigma);
}
"
fitphiphi_step4 <- stan_model(model_code=Neuse_step4)


neuse_step4 <- function(data=neuse_flowN,  y="conc",
                        x="flow", n.chains=nchains){
    y <- log(data[,y])
    x <- log(data[,x])
    wT <- data$waterT
    avgT <- mean(wT)
    wT <- data$waterT-avgT
    yr <- as.numeric(data$yr)
    n <- length(y)
    maxYr <- max(yr)
    lowF <- min(x)
    uppF <- max(x)

    standata <- list(N=n, y=y, x=x, yr=yr, Nyr=max(yr), wtrTemp=wT,
                     max_yr=maxYr)
    staninits <- list()
    for (i in 1:n.chains){
        staninits[[i]]<-list(alpha0=rnorm(1), beta0=rnorm(1),
                             alpha1=rnorm(1), beta1=-runif(1), betaT=rnorm(1),
                             beta2=-1-runif(1), sigma=runif(1),
                             phi_yr=runif(1,1,maxYr))
    }
    stanpara=c("alpha0", "alpha1","beta0","beta1","beta2","betaT",
               "phi_yr","sigma")
        return(list(data=standata, inits=staninits, pars=stanpara,
                nchains=n.chains, meanX=avgT))
}

## 1. centering the predictor (log(Q))
input.to.stan <- neuse_step4 ()
fit_neuseStep_cX4 <- sampling(fitphiphi_step4,
                             data=input.to.stan$data,
                             init=input.to.stan$inits,
                             pars=input.to.stan$pars,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                     control=list(max_treedepth=15))#,
#                                  adapt_delta=0.95))
print(fit_neuseStep_cX4)
pairs(fit_neuseStep_cX4, pars=c("alpha0","beta0","alpha1","beta1","beta1",
                                "betaT","phi_yr"))
save(fit_neuseStep_cX4, file="neuseThCh_cX4.RData")

neuseThCh4 <- rvsims(as.matrix(as.data.frame(extract(fit_neuseStep_cX4))))


stanin4 <- neuse_step4 ()

### packages(Rmpfr) ## for function log1pexp
tmp <- range(stanin4$data$x)
X <- seq(tmp[1], tmp[2],,50)
plotX <- X-log(86400)

YR <- as.numeric(neuse_flowN$yr)
phiYR = 1-as.numeric(YR > neuseThCh4[7])
mu1 <- neuseThCh4[1]+##neuseThCh2[5]*stanin2$data$wtrTemp+
    neuseThCh4[2]*X
mu2 <- neuseThCh4[3]+neuseThCh4[4]*(X)+
    neuseThCh4[5]*log(X)
mu1_summ <- summary(mu1)
mu2_summ <- summary(mu2)

phi_yr <- summary(phiYR)$mean==0
tikz(file=paste(plotDIRch7, "neuseChThmodel.tex", sep="/"),
     height=3, width=3.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(log(conc)~log(flow/86400), data=neuse_flowN, col="gray", type="n",
     xlab="log flow (m$^3$/s)", ylab="log concentration (mg/L)")
points(log(conc)~log(flow/86400), data=neuse_flowN, subset=!phi_yr, pch=16, col=gray(0.7), cex=0.5)
points(log(conc)~log(flow/86400), data=neuse_flowN, subset=phi_yr, pch=16, col=gray(0.3), cex=0.5)
polygon(x=c(plotX,rev(plotX)),
        y=c(mu1_summ$"2.5%", rev(mu1_summ$"97.5%")), col=gray(0.7), border=F)
polygon(x=c(plotX,rev(plotX)),
        y=c(mu2_summ$"2.5%", rev(mu2_summ$"97.5%")), col=gray(0.3), border=F)
lines(plotX,mu1_summ$mean)
lines(plotX,mu2_summ$mean)
dev.off()
