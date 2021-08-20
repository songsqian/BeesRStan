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

## Floride occuurence in Chinese drinking water source waters
#### Large data -- runs very slow
chinaDW <- read.csv(file=paste(dataDIR, "Source_Data.csv", sep="/"))
names(chinaDW)

##
summary(chinaDW$F)
table(chinaDW$Type)
length(table(chinaDW$Province))

##
## stan model

CNdw_stan1 <- "
data {
  int<lower=0> N;
  int<lower=0> Ncens;
  int<lower=0> Nprov;
  int<lower=0> Ncity;
  int<lower=1,upper=Ncity> cid[N];
  int<lower=1,upper=Ncity> cidcens[Ncens];
  int<lower=1,upper=Nprov> city_prov[Ncity];
  vector[N] y;
  vector[Ncens] ycens;
}
parameters {
  vector[Ncity] theta;
  vector[Nprov] gamma;
  real mu;
  real<lower=0,upper=10> sigma_1;
  real<lower=0,upper=10> sigma_2;
  real<lower=0,upper=10> sigma_y;
}
model {
  mu ~ normal(0, 10);
  gamma ~ normal(mu, sigma_2);
  for (j in 1:Ncity)
    theta[j] ~ normal(gamma[city_prov[j]], sigma_1);
  for (i in 1:N)
    target += normal_lpdf(y[i] | theta[cid[i]], sigma_y);
  for (i in 1:Ncens)
    target += normal_lcdf(ycens[i] | theta[cidcens[i]], sigma_y);
}
"
fit1 <- stan_model(model_code=CNdw_stan1)
stan_in1 <- function(data=chinaDW, chains=nchains){
    data <- data[!is.na(data$F),]
    yAll <- data$F
    detect <- data$Note4
    y <- yAll[detect!="<"]
    ycens <- data$F[detect=="<"]
    n <- length(y)
    ncens <- length(ycens)
    cityAll <- as.numeric(ordered(data$City.num))
    cid <- cityAll[detect!="<"]
    cidcens <- cityAll[detect=="<"]
    ncity <- max(cityAll)
    oo <- order(cityAll)
    cityAll <- cityAll[oo]
    provAll <- as.numeric(ordered(data$Province))[oo]
    nprov <- max(provAll)
    city_prov <- provAll[cumsum(table(cityAll))]
    stan_data <- list(N=n, Ncens=ncens, Nprov=nprov, Ncity=ncity,
                      y=log(y), ycens=log(ycens),
                      city_prov=city_prov, cid=cid, cidcens=cidcens)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(theta=rnorm(ncity, -2),
                                gamma=rnorm(nprov, -2), mu=rnorm(1, -2),
                                sigma_1=runif(1),
                                sigma_2=runif(1), sigma_y=runif(1))
    stan_pars <- c("theta", "gamma", "mu", "sigma_1","sigma_2","sigma_y")
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

Fstan <- rstan::extract(fit2keepF)
Ftheta <- rvsims(Fstan$theta)
Fgamma <- rvsims(Fstan$gamma)
Fmu <- rvsims(Fstan$mu)
Fsigy <- rvsims(Fstan$sigma_y)
Fsig1 <- rvsims(Fstan$sigma_1)
Fsig2 <- rvsims(Fstan$sigma_2)

Ftheta[input.to.stan$data$city_prov==19]
Fgamma[19]

dens_plot <- function(theta = Ftheta,
                      input = input.to.stan$data$city_prov,
                      gamma=Fgamma, prov=26,
                      Xlab="Fluoride Concentration (mg/L)",
                      Ylab="", mn=NULL, Axis1=T, Axis2=T,
                      Prvs=levels(chinaDW$Province)){
    theta <- theta[input==prov]
    gamma <- gamma[prov]
    summ <- summary(theta)
    upper <- max(c(round(max(summ$"99%")), log(4)))
    lower <- min(c(round(min(summ$"1%")), log(0.7)))
    if(is.null(mn)) mn <- Prvs[prov]
    plot(density(sims(gamma)[,1], from=lower, to=upper),
         xlab=Xlab, ylab=Ylab, main=mn, type="n",axes=F)
    for (i in 1:length(theta))
        lines(density(sims(theta[i])[,1]), col="grey")
    abline(v=log(c(4,2,1.5,1,0.7)), lty=2)
    lines(density(sims(gamma)[,1], from=lower, to=upper),
          lwd=3)
    box()
    if (Axis2)
        axis(2, outer=T)
    else
        axis(4, outer=T)
    if (Axis1)
        axis(1, at=log(c(0.01,0.05,0.1,0.5,1,1.5,2,5)),
             label=c(0.01,0.05,0.1,0.5,1,1.5,2,5), outer=T)
    invisible()
}

cdf_plot <- function(theta = Ftheta, sig=Fsig1,
                     input = input.to.stan$data$city_prov,
                     gamma=Fgamma, prov=26,
                     Xlab="Fluoride Concentration (mg/L)",
                     Ylab="cdf", mn=NULL, Axis1=T, Axis2=T,
                     Prvs=levels(chinaDW$Province)){
    theta <- theta[input==prov]
    gamma <- gamma[prov]
    summ <- summary(theta)
    upper <- max(c(round(max(summ$"99%")), log(4)))
    lower <- min(c(round(min(summ$"1%")), log(0.7)))
    cut_x <- seq(lower,upper,length=100)
    gm <- rvnorm(1, gamma, sig)
    if(is.null(mn)) mn <- Prvs[prov]
    plot(cut_x, 1-Pr(gm>cut_x), type="n", lwd=3,
         xlab=Xlab, ylab=Ylab, main=mn, axes=F)
    for (i in 1:length(theta))
        lines(cut_x, 1-Pr(theta[i]>cut_x), col="grey")
    abline(v=log(c(4,2,1.5,1,0.7)), lty=2)
    lines(cut_x, 1-Pr(gm>cut_x), lwd=3)
    box()
    if (Axis2)
        axis(2, outer=T)
    else
        axis(4, outer=T)
    if (Axis1)
        axis(1, at=log(c(0.01,0.05,0.1,0.5,1,1.5,2,5)),
             label=c(0.01,0.05,0.1,0.5,1,1.5,2,5), outer=T)
    invisible()
}

## Sichuan
cdf_plot(prov=26)
dens_plot()

## Anhui
dens_plot(prov=1)
cdf_plot(prov=1)

## Inner Mongolia
dens_plot(prov=19, mn="Inner Mongolia")
cdf_plot(prov=19, mn="Inner Mongolia")

dens_plot(prov=10)
cdf_plot(prov=10)

tikz(file=paste(plotDIRch6, "ChinaDW.tex", sep="/"),
     width=5, height=5, standAlone=F)
par(mfrow=c(2,2), mar=c(rep(0,4)), mgp=c(1.25,0.125,0),
    oma=c(rep(3,4)), las=0, tck=0.01)
dens_plot(prov=26, Axis1=F, mn="")
mtext("pdf", side=2, line=1.5)
mtext(levels(chinaDW$Province)[26],
      side = 3, line = 1, outer = F)
dens_plot(prov=19, Axis1=F, Axis2=F, mn="")
mtext("pdf", side=4, line=1.5)
mtext("Inner Mongolia",
      side = 3, line = 1, outer = F)

cdf_plot(prov=26, mn="")
mtext("cdf", side=2, line=1.5)
cdf_plot(prov=19, Axis2=F, mn="")
mtext("cdf", side=4, line=1.5)
mtext("Fluoride Concentration (mg/L)", side=1, line=1.5, outer=T)
dev.off()

