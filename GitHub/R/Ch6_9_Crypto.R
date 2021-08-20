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
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Cryptosporidium survey data
icr <- read.table(paste(dataDIR,"CryptoData2.txt",sep="/"), header=T,
        na.string = ".", sep = "\t")
## example of Poisson multilevel model

levels(ordered(icr$MSrcCat)) -> Lmsrc
levels(ordered(icr$M.WTP.Type))->Lwtptype

dcts.data <- icr[icr$MSrcCat==Lmsrc[3] | icr$MSrcCat==Lmsrc[4],]
dcts.data <- dcts.data[dcts.data$M.WTP.Type=="Y"|dcts.data$M.WTP.Type=="N",]
dcts.data$M.WTP.Type <- ordered(as.vector( dcts.data$M.WTP.Type ))
dcts.data$MSrcCat <- ordered(as.vector( dcts.data$MSrcCat ))
dcts.data$ICR.PWSID <- ordered(dcts.data$ICR.PWSID)

## method of moments for alpha and beta
r_bar <- 0.44
sr2 <- (0.33/4)^2
alpha=(r_bar^2*(1-r_bar)-r_bar*(sr2))/(sr2)
beta=alpha*(1-r_bar)/r_bar

icr.lmer1 <- glmer(n.cT ~ 1+(1|ICR.PWSID),
                   data=dcts.data, family="poisson",
                   offset=log(volume*0.44))

icr.size <- as.vector(table(dcts.data$ICR.PWSID))
size <- icr.size + runif(length(icr.size), -0.1,0.1)

## Figure 10.29
tikz(file=paste(plotDIRch10, "cryptolmer1.tex",
                sep="/"),
     width=3.5, height=2.5, standAlone=F)
lmer.mean <- fixef(icr.lmer1)[1] +
    ranef(icr.lmer1)[[1]][,1]
lmer.se <- sqrt(se.fixef(icr.lmer1)[1]^2 +
                se.ranef(icr.lmer1)[[1]][,1]^2)
lower <- lmer.mean-lmer.se
upper <- lmer.mean+lmer.se

par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0),
    tck=0.01, las=1)
plot(size, rnorm(length(size)),type="n",
     xlab="Sample size",
     ylab="Log Mean", log="x", ylim=range(lower,upper))
abline(h=fixef(icr.lmer1)[[1]][1])
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, lmer.mean, pch=16,cex=0.5)
dev.off()

cs <- coef(icr.lmer1)[[1]][,1]
n.cs <- length(cs)

## Figure 10.30
tikz(file=paste(plotDIRch10, "cryptocdf.tex", sep="/"),
     width=3.5, height=2.5, standAlone=F)
par(mar=c(3,3,0.5,0.5), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
plot(sort(cs), ((1:n.cs)-0.5)/n.cs, axes=F,
     xlab="Crypto Concentration (oocyst/L)", ylab="CDF",
     type="n")
points(sort(cs), ((1:n.cs)-0.5)/n.cs, pch=1,
       cex=0.5, col="gray")
curve(pnorm(x, -5.384, 2.08), add=T)
axis(1, at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
     labels=c("0.0001","0.001","0.01","0.1","1","10","100"))
axis(2)
box()
dev.off()

##
## stan model
Crypto_stan1 <- "
data {
  int<lower=0> Nobs;
  int<lower=0> Nmax;
  int<lower=0> Npwsid;
  int<lower=0,upper=Npwsid> pwsid[Nobs];
  int y[Nobs];
  vector[Nobs] vol;
  real alpha;
  real beta;
}
parameters {
  vector[Npwsid] lconc;
  real mu;
  real<lower=0,upper=10> sigma;
  real<lower=0,upper=1> r;
}
transformed parameters {
  real lambda[Nobs];
  for (i in 1:Nobs)
    lambda[i] = exp(lconc[pwsid[i]]) * vol[i];
}
model{
  int k;
  r ~ beta(alpha, beta);
  lconc ~ normal(mu, sigma);
  for (i in 1:Nobs){
    vector[Nmax-y[i]+1] temp;
    for (j in y[i]:Nmax){
      k = j-y[i]+1;
      temp[k] = binomial_lpmf(y[i] | j, r) +
                poisson_lpmf(j | lambda[i]);
    }
    target += log_sum_exp(temp);
  }
}
"
fit1 <- stan_model(model_code=Crypto_stan1)
stan_in1 <- function(data=dcts.data, a=alpha, b=beta, chains=nchains,
                     Nmax=100){
    data <- data[data$n.cT<100,]
    y <- data$n.cT
    n <- length(y)
    pwsid <- as.numeric(ordered(data$ICR.PWSID))
    npwsid <- max(pwsid)
    stan_data <- list(Nobs=n, Nmax=Nmax, Npwsid=npwsid,
                      pwsid=pwsid, vol=data$volume,
                      y=y, alpha=a, beta=b)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(lconc=rnorm(npwsid, -2),
                                mu=rnorm(1, -2),
                                sigma=runif(1), r=runif(1))
    stan_pars <- c("lconc", "mu", "sigma", "r")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

## full data
input.to.stan <- stan_in1(Nmax=200)
fit2keepFull <- sampling(fit1, data=input.to.stan$data,
                         init=input.to.stan$inits,
                         pars=input.to.stan$pars,
                         iter=niters,thin=nthin,
                         chains=input.to.stan$n.chains)
##                     control=list(max_treedepth=25))

pwsid <- table(dcts.data$ICR.PWSID)
subdata <- dcts.data$ICR.PWSID %in% names(pwsid)[pwsid>20 & pwsid<25]
tem.data <- dcts.data[subdata, ]
input.to.stan <- stan_in1(data=tem.data, Nmax=75)
fit2keep75 <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)
##                     control=list(max_treedepth=25))
input.to.stan <- stan_in1(data=tem.data, Nmax=150)
fit2keep150 <- sampling(fit1, data=input.to.stan$data,
                        init=input.to.stan$inits,
                        pars=input.to.stan$pars,
                        iter=niters,thin=nthin,
                        chains=input.to.stan$n.chains)
##                     control=list(max_treedepth=25))

save(fit2keep75, fit2keep150, file="Crypto75v150r.RData")
save(fit2keep150, file="Crypto150r.RData")
load("Crypto75v100r.RData")
print(fit2keep150)
pairs(fit2keep150, pars=c("mu","sigma", "r"))

## 3600 v 4600 seconds
stan75 <- rstan::extract(fit2keep75)
lconc75 <- rvsims(stan75$lconc)
sigma75 <- rvsims(stan75$sigma)
mu75 <- rvsims(stan75$mu)
r75 <- rvsims(stan75$r)

stan150 <- rstan::extract(fit2keep150)
lconc150 <- rvsims(stan150$lconc)
sigma150 <- rvsims(stan150$sigma)
mu150 <- rvsims(stan150$mu)
r150 <- rvsims(stan150$r)

write.table((round(rbind(summary(sigma75), summary(sigma150),
                         summary(mu75), summary(mu150),
                         summary(r75), summary(r150))[,c(1,2,4:8)], 3)),
      file=paste(plotDIRch6, "CryptoNmaxTabr.tex", sep="/"), sep="&")

save(stan75, stan150, file="CryptoPois.RData")

tikz(paste(plotDIRch6, "CryptoNmaxComp.tex", sep="/"),
     height=4, width=4, standAlone=F)
par(mfrow=c(2,1), mar=c(0,1,0,1), oma=c(3,3,3,3),
    mgp=c(1.25,0.125,0), las=0, tck=0.01)
plot(sort(lconc75), ylim=c(-20,0), axes=F)
axis(2, las=1, cex.axis=0.75)
axis(3, cex.axis=0.75)
box()
plot(sort(lconc150), ylim=c(-20,0), axes=F)
axis(1, cex.axis=0.75)
axis(4, las=1, cex.axis=0.75)
box()
mtext("lab number", side=1, outer=T, line=1.)
mtext("log concentration", side=2, outer=T, line=1.)
dev.off()

### Binary formulation
Crypto_stan2 <- "
data {
  int<lower=0> Nobs;
  int<lower=0> Npwsid;
  int<lower=0,upper=Npwsid> pwsid[Nobs];
  int y[Nobs];
  vector[Nobs] vol;
  real a;
  real b;
}
parameters {
  vector[Npwsid] lconc;
  real<lower=0,upper=1> r;
  real mu;
  real<lower=0,upper=10> sigma;
}
transformed parameters{
  vector[Npwsid] conc;
  real lambda[Nobs];
  for (i in 1:Nobs)
    lambda[i] = exp(lconc[pwsid[i]]) * vol[i];
  conc = exp(lconc);
}
model {
  real temp2[2];
  r ~ beta(a, b);
  lconc ~ normal(mu, sigma);
  for (i in 1:Nobs){
    temp2[1] = -conc[pwsid[i]]*vol[i];
    temp2[2] = log1m_exp(-conc[pwsid[i]]*vol[i])+log(1-r);
    target += y[i]*(log1m_exp(-conc[pwsid[i]]*vol[i]) + log(r)) +
              (1-y[i])*log_sum_exp(temp2);
  }
}
"
fit2<- stan_model(model_code=Crypto_stan2)


stan_in2 <- function(data=dcts.data, a=alpha, b=beta,
                     chains=nchains){
    data <- data[data$n.cT<100,]
    y <- as.numeric(data$n.cT>0)
    n <- length(y)
    pwsid <- as.numeric(ordered(data$ICR.PWSID))
    npwsid <- max(pwsid)
    stan_data <- list(Nobs=n, Npwsid=npwsid,
                      pwsid=pwsid, vol=data$volume,
                      y=y, a=a, b=b)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(lconc=rnorm(npwsid, -2),
                                mu=rnorm(1, -2),
                                sigma=runif(1))
    stan_pars <- c("lconc", "mu", "sigma", "r")
    return(list(data=stan_data, inits=stan_inits,
                pars=stan_pars, n.chains=chains))
}

input.to.stan <- stan_in2()
fit2keepBin <- sampling(fit2, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)
##                     control=list(max_treedepth=25))
print(fit2keepBin)
save(fit2keepBin,file= "CryptoBin.RData")
load("CryptoBin.RData")
pairs(fit2keepBin, pars=c("mu","sigma","r"))
## 5050 seconds
stanBin <- rstan::extract(fit2keepBin)
lconcB <- rvsims(stanBin$lconc)
sigmaB <- rvsims(stanBin$sigma)
muB <- rvsims(stanBin$mu)
rB <- rvsims(stanBin$r)

plot(sort(lconcB))
## Running the selected data

input.to.stan <- stan_in2(data=tem.data)
fit2keepBinSMP <- sampling(fit2, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)#,
##                     control=list(adapt_delta=0.9,max_treedepth=25))
print(fit2keepBinSMP)
pairs(fit2keepBinSMP, pars=c("mu","sigma","r"))
## 50 seconds

stanBinSMP <- rstan::extract(fit2keepBinSMP)
lconcBSMP <- rvsims(stanBinSMP$lconc)
sigmaBSMP <- rvsims(stanBinSMP$sigma)
muBSMP <- rvsims(stanBinSMP$mu)
rBSMP <- rvsims(stanBinSMP$r)

write.table(round(rbind(summary(sigma75),summary(sigma150),
                         summary(sigmaBSMP),summary(sigmaB),
                      summary(mu75),summary(mu150),summary(muBSMP),summary(muB),
                         summary(r75),summary(r150),summary(rBSMP),summary(rB))[,c(1,2,4:8)], 3),
      file=paste(plotDIRch6, "CryptoNmaxTab2.tex", sep="/"), sep="&")

tikz(paste(plotDIRch6, "CryptoNmaxComp2.tex", sep="/"),
     height=6, width=4, standAlone=F)
par(mfrow=c(4,1), mar=c(0,1,0,1), oma=c(3,3,3,3),
    mgp=c(1.25,0.125,0), las=0, tck=0.01)
plot(sort(lconc75), (1:length(lconc75))/length(lconc75),
     xlim=c(-25,5), axes=F)
axis(2, las=1, cex.axis=0.75)
axis(3, cex.axis=0.75)
box()
text(x=-20, y=0.8, "(a) Poisson mixture, subset, \\texttt{Nmax=75}")
plot(sort(lconc150), (1:length(lconc150))/length(lconc150),
     xlim=c(-25,5), axes=F)
axis(4, las=1, cex.axis=0.75)
box()
text(x=-20, y=0.8, "(b) Poisson mixture, subset, \\texttt{Nmax=150}")
plot(sort(lconcBSMP), (1:length(lconcBSMP))/length(lconcBSMP),
     xlim=c(-25,5), axes=F)
axis(2, las=1, cex.axis=0.75)
box()
text(x=-20, y=0.8, "(c) Binary mixture, sub-set")
plot(sort(lconcB), (1:length(lconcB))/length(lconcB),
     xlim=c(-25,5), axes=F)
axis(1, cex.axis=0.75)
axis(4, las=1, cex.axis=0.75)
box()
text(x=-20, y=0.8, "(d) Binary mixture, full data")
mtext("Cumulative Mean Distribution", side=2, outer=T, line=1.)
mtext("log concentration", side=1, outer=T, line=1.)
dev.off()

