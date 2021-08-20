source("FrontMatter.R")

## Chapter 4

plotDIRch4 <- paste(plotDIR, "chapter4", "figures", sep="/")

packages(rv)
packages(rstan)
packages(car)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Snake fungal disease
## Initial model
sfd_stan <- "
  data{
  int n;
  int x;
  real alpha;
  real beta;
  real ap;
  real an;
  real bp;
  real bn;
}
parameters{
  real theta;
  real fp;
  real fn;
}
transformed parameters{
  real ppos;
  ppos = theta*(1-fn)+(1-theta)*fp;
}
model{
  theta ~ beta(alpha, beta);
  fp ~ beta(ap, bp);
  fn ~ beta(an, bn);
  target += x*log(ppos) + (n-x)*log1p(-ppos);
}
"
fit <- stan_model(model_code = sfd_stan)
sfd_input <- function(alpha = 1, beta = 1, ap = 5, bp = 20,
                      an = 2, bn = 25, n=20, x=5, n.chains=nchains){
    data  <- list(alpha=alpha,beta=beta,ap=ap,bp=bp,an=an,bn=bn,
                  n = n, x = x)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(fp=rbeta(1, 5,20), fn=rbeta(1,2,25),
                           theta=runif(1))
    pars <- c("theta", "fp", "fn")
    return (list(data=data, inits=inits, pars=pars))
}

input.to.stan <- sfd_input()
fit2keep <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=nchains)
print(fit2keep)
pairs(fit2keep, pars=input.to.stan$pars)
firstrun <- extract(fit2keep)
div <- get_divergent_iterations(fit2keep)

tikz(file=paste(plotDIRch4,"sfd_pairs1.tex", sep="/"),
     height=2, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(firstrun$theta, firstrun$fp, xlab="$\\theta$", ylab="$f_p$")
points(firstrun$theta[div], firstrun$fp[div], pch=16, cex=0.75, col="gray")

plot(firstrun$theta, firstrun$fn, xlab="$\\theta$", ylab="$f_n$")
points(firstrun$theta[div], firstrun$fn[div], pch=16, cex=0.75, col="gray")

plot(firstrun$fn, firstrun$fp, xlab="$f_n$", ylab="$f_p$")
points(firstrun$fn[div], firstrun$fp[div], pch=16, cex=0.75, col="gray")
dev.off()

## posterior simulation (see Chapter 8)
p_posit <- firstrun$theta*(1-firstrun$fn)+
           (1-firstrun$theta)*firstrun$fp
post_sim <- rbinom(length(p_posit), size=20, p_posit)
postsims <- table(post_sim)

tikz(file=paste(plotDIRch4,"sfd_postsim1.tex", sep="/"),
     height=2.5, width=2.5, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.5,0.125,0), las=1, tck=0.01)
barplot(postsims, xlab="$x$")
lines(1:length(postsims),
      max(postsims)*cumsum(table(post_sim))/length(p_posit))
dev.off()

summary(post_sim)

## increasing sample size
input.to.stan <- sfd_input(n=2000, x=500)
fit2keep <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=nchains)
print(fit2keep)
pairs(fit2keep, pars=input.to.stan$pars)
firstrun <- extract(fit2keep)
div <- get_divergent_iterations(fit2keep)

tikz(file=paste(plotDIRch4,"sfd_pairs2.tex", sep="/"),
     height=2, width=5.5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(firstrun$theta, firstrun$fp, xlab="$\\theta$", ylab="$f_p$")
points(firstrun$theta[div], firstrun$fp[div], pch=16, cex=0.75, col="gray")

plot(firstrun$theta, firstrun$fn, xlab="$\\theta$", ylab="$f_n$")
points(firstrun$theta[div], firstrun$fn[div], pch=16, cex=0.75, col="gray")

plot(firstrun$fn, firstrun$fp, xlab="$f_n$", ylab="$f_p$")
points(firstrun$fn[div], firstrun$fp[div], pch=16, cex=0.75, col="gray")
dev.off()

## setting f_p=0
sfd_stan0 <- "
  data{
  int n;
  int x;
  real alpha;
  real beta;
  real an;
  real bn;
}
parameters{
  real theta;
  real fn;
}
transformed parameters{
  real ppos;
  ppos = theta*(1-fn);
}
model{
  theta ~ beta(alpha, beta);
  fn ~ beta(an, bn);
  target += x*log(ppos) + (n-x)*log1p(-ppos);
}
"
fit <- stan_model(model_code = sfd_stan0)

sfd_input0 <- function(alpha = 1, beta = 1, an = 2, bn = 25,
                       n=20, x=5, n.chains=nchains){
    data  <- list(alpha=alpha,beta=beta,an=an,bn=bn,
                  n = n, x = x)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(fn=rbeta(1,2,25),
                           theta=runif(1))
    pars <- c("theta", "fn")
    return (list(data=data, inits=inits, pars=pars))
}

##increase sample size

input.to.stan <- sfd_input0(an=2,bn=25,n=2000, x=500)
fit2keep1 <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=nchains)
input.to.stan <- sfd_input0(an=20,bn=250,n=2000, x=500)
fit2keep2 <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=nchains)

firstrun1 <- extract(fit2keep1)
div1 <- get_divergent_iterations(fit2keep1)
firstrun2 <- extract(fit2keep2)
div2 <- get_divergent_iterations(fit2keep2)

tikz(file=paste(plotDIRch4,"sfd_pairs3.tex", sep="/"),
     height=2, width=3.75, standAlone=F)
par(mfrow=c(1,2), mar=c(3,3,1,0.5), mgp=c(1.5,0.125,0), las=1, tck=0.01)
plot(firstrun1$fn, firstrun1$theta, ylab="$\\theta$", xlab="$f_n$",
     xlim=c(0, 0.35), ylim=c(0.2,0.4))
points(firstrun1$fn[div], firstrun1$theta[div], pch=16, cex=0.75, col="gray")

plot(firstrun2$fn, firstrun2$theta, ylab="$\\theta$", xlab="$f_n$",
     xlim=c(0, 0.35), ylim=c(0.2,0.4))
points(firstrun2$fn[div], firstrun2$theta[div], pch=16, cex=0.75, col="gray")
dev.off()

## Seaweed grazer
seaweed <- read.csv(paste(dataDIR, "seaweed.csv", sep="/"))
seaweed$TREAT <- ordered(seaweed$TREAT,
                         levels=c("CONTROL","f","fF","L","Lf","LfF"))
seaweed$y <- logit(seaweed$COVER/100)
names(seaweed)[2:3] <- c("Block","Treatment")

## Stan code (ANOVA-1)
stan_aov1 <- "
data{
  int N;
  int J;
  real y[N];
  int treatment[N];
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> n0;
  real<lower=0> mu0;
}
parameters{
  real<lower=0> sigma1sq;
  real<lower=0> sigma2sq;
  real theta;
  real mu_j[J];
}
model{
  sigma1sq ~ inv_gamma(alpha, beta);
  sigma2sq ~ inv_gamma(a, b);
  theta ~ normal(mu0, sqrt(sigma2sq/n0));
  mu_j ~ normal(theta, sqrt(sigma2sq));
  for (i in 1:N){
    y[i] ~ normal(mu_j[treatment[i]], sqrt(sigma1sq));
  }
}
generated quantities{
  real grand_mu;
  real trt[J];
  grand_mu = mean(mu_j[]);
  for (j in 1:J){
    trt[j] = mu_j[j];
  }
}
"
fit <- stan_model(model_code = stan_aov1)

## input data
input_Stan <- function(data=seaweed, chains=nchains,
                       a, b, alpha, beta, mu0, n0){
  y <- car::logit(data$COVER)
  trt <- as.numeric(ordered(data$Treatment))
  J <- max(trt)
  N <- length(y)
  stan_data <- list(N=N, J=J, y=y, treatment=trt,
                    a=a, b=b, alpha=alpha,
                    beta=beta, mu0=mu0, n0=n0)
  stan_inits <- list()
  for (i in 1:chains)
    stan_inits[[i]] <- list(sigma1sq=runif(1), sigma2sq=runif(1),
                            mu_j=rnorm(J), theta=rnorm(1))
  para=c("mu_j","sigma1sq","sigma2sq", "theta", "grand_mu", "trt")
  return(list(data=stan_data, inits=stan_inits, paras=para,
              n.chains=chains))
}

input.to.stan <- input_Stan(a=0.1, b=0.1,
                            alpha=0.1, beta=0.1,
                            mu0=0, n0=1)

## running stan (model output saved)
fit2keep <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)

## processing output

mcmcdraws <- as.data.frame(extract(fit2keep))
save(mcmcdraws, file="seaweed1w.RData")

## Comparisons
## 1. estimated means
swd.aov <- aov(car::logit(COVER) ~ Treatment, data=seaweed)
model.tables(swd.aov, "means", se=T)  -> aovtables
summary(swd.aov)

mcmcdraws_summ <- summary(rvsims(as.matrix(mcmcdraws)))
tikz(file=paste(plotDIRch4, "swdAOVcomp.tex", sep="/"), width=3.5,
     height=3, standAlone=F)
par(mar=c(3, 6, 0.25,0.25), mgp=c(1.25, 0.125, 0), tck=0.01)
line.plots.compare(est1=unlist(aovtables$tables$Treatment),
                   se1 =rep(unlist(aovtables$se), 6),
                   est2=mcmcdraws_summ$mean[1:6],
                   se2 =mcmcdraws_summ$sd[1:6],
                   Xlab="Means",
                   Ylabel=levels(seaweed$Treatment),
                   V=unlist(aovtables$tables[[1]]))
dev.off()

## 2. variance components

outputnames <- names(mcmcdraws)
mujs <- mcmcdraws[substring(outputnames, 1, 2)=="mu"]
muj_var <- rvsims(apply(mujs, 1, var))
sigma1var <- rvsims(unlist(mcmcdraws[outputnames=="sigma1sq"]))
var_comp <- sqrt(c(muj_var, sigma1var))
names(var_comp) <- c("Treatment", "Residuals")

tikz(file=paste(plotDIRch4, "swdVARcomp.tex", sep="/"),
     height=3, width=3, standAlone=T)
mlplot(var_comp, xlab="standard deviation")
dev.off()

## Two-way ANOVA
## Stan code (ANOVA-2)
stan_aov2 <- "
data{
  int N;
  int J;
  int K;
  real y[N];
  int treatment[N];
  int blck[N];
}
parameters{
  real<lower=0> sigma1sq;
  real<lower=0> sigma2sq;
  real<lower=0> sigma3sq;
  real theta;
  real delta_j[J];
  real delta_k[K];
}
model{
  delta_j ~ normal(0, sqrt(sigma2sq));
  delta_k ~ normal(0, sqrt(sigma3sq));
  for (i in 1:N){
    y[i] ~ normal(theta+delta_j[treatment[i]]+delta_k[blck[i]],
                  sqrt(sigma1sq));
  }
}
generated quantities{
  real mu;
  real trt[J];
  real blk[K];
  mu = theta+mean(delta_j[])+mean(delta_k[]);
  for (j in 1:J){
    trt[j] = delta_j[j]-mean(delta_j[]);
  }
  for (k in 1:K){
    blk[k] = delta_k[k]-mean(delta_k[]);
  }
}
"
fit <- stan_model(model_code = stan_aov2)

## input data
input_Stan <- function(data=seaweed, chains=nchains){
  y <- car::logit(data$COVER)
  trt <- as.numeric(ordered(data$Treatment))
  blk <- as.numeric(ordered(data$Block))
  J <- max(trt)
  K <- max(blk)
  N <- length(y)
  stan_data <- list(N=N, J=J, K=K, y=y, treatment=trt, blck=blk)
  stan_inits <- list()
  for (i in 1:chains)
    stan_inits[[i]] <- list(sigma1sq=runif(1), sigma2sq=runif(1), sigma3sq=runif(1),
                            delta_j=rnorm(J), delta_k=rnorm(K), theta=rnorm(1))
  para=c("mu", "trt", "blk", "sigma1sq","sigma2sq", "sigma3sq")
  return(list(data=stan_data, inits=stan_inits, paras=para,
              n.chains=chains))
}

input.to.stan <- input_Stan()

## running stan (output saved, no need to rerun)
fit2keep <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
	             iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)

## processing output
mcmcdraws2 <- as.data.frame(extract(fit2keep))
save(mcmcdraws2, file="seaweed2w.RData")
mcmcdraws2_summ <- summary(rvsims(as.matrix(mcmcdraws2)))
## Comparisons
## 1. estimated means
swd.aov2 <- aov(car::logit(COVER) ~ Treatment+Block, data=seaweed)
model.tables(swd.aov2, se=T)  -> aovtables2
summary(swd.aov2)

tikz(file=paste(plotDIRch4, "swdAOVcomp2w.tex", sep="/"), width=7,
     height=3, standAlone=F)
par(mfrow=c(1,2), mar=c(3, 6, 0.25,0.25), mgp=c(1.25, 0.125, 0), tck=0.01)
line.plots.compare(est1=unlist(aovtables2$tables$Treatment),
                   se1 =rep(unlist(aovtables2$se[1]), 6),
                   est2=mcmcdraws2_summ$mean[2:7],
                   se2 =mcmcdraws2_summ$sd[2:7],
                   Xlab="Means",
                   Ylabel=levels(seaweed$Treatment))

line.plots.compare(est1=unlist(aovtables2$tables$Block),
                   se1 =rep(unlist(aovtables2$se[2]), 8),
                   est2=mcmcdraws2_summ$mean[8:15],
                   se2 =mcmcdraws2_summ$sd[8:15],
                   Xlab="Means",
                   Ylabel=levels(seaweed$Block))
dev.off()
## 2. variance components

drawnames <- names(mcmcdraws2)
deltajs <- mcmcdraws2[substring(drawnames, 1, 3)=="trt"]
deltaks <- mcmcdraws2[substring(drawnames, 1, 3)=="blk"]

deltaj_var <- rvsims(apply(deltajs, 1, var))
deltak_var <- rvsims(apply(deltaks, 1, var))

var_comp <- sqrt(c(deltaj_var, deltak_var, mcmcdraws2[16]))
names(var_comp) <- c("Treatment", "Block", "Residuals")

tikz(file=paste(plotDIRch4, "swdVARcomp2.tex", sep="/"),
     height=3, width=3, standAlone=F)
mlplot(var_comp, xlab="standard deviation")
dev.off()

## GOM Hypoxia
## This example uses three model formulations
## The original model was presented in Qian et al (2009)
## The alternative models 1 and 2 are proposed to resolve
##    numerical issues

########### See ch4_GOM.R
## benthic.data <- read.csv(paste(dataDIR, "BenthicData.csv", sep="/"),
##                          header=T)
## benthic.data$area2 <- ordered(benthic.data$area2,
##                               levels=unique(benthic.data$area2)[c(2,1,3)])
## benthic.data$Area <- ordered(benthic.data$Area,
##                              levels=unique(benthic.data$Area)[c(2,1,3)])

## head(benthic.data)

## ## Station: core,
## ## Area: zone
## ## The original model
## stan_gom <- " /* Stan Code - GOM model*/
## data{
##   int K;  //total sample size
##   int J;  //number of sediment cores
##   int I;  //number of zones
##   real y[K]; //observed response
##   int core[K]; //core index
##   int zone[K]; //zone index
##   int core_zone[J]; //zone
## }
## parameters{
##   real mu[J];
##   real theta[I];
##   real mu_hyp;
##   real<lower=0> sigma_y[I];
## //  real<lower=0> sigma_y;
## //  real<lower=0> sigma_i[I];
##   real<lower=0> sigma_i;
##   real<lower=0> sigma_hyp;
## }
## model{
##   for (i in 1:I){
##     theta[i] ~ normal(mu_hyp, sigma_hyp);
##   }
##   for (j in 1:J){
##     mu[j] ~ normal(theta[core_zone[j]], sigma_i);
##   }
##   for (k in 1:K){
##     y[k] ~ normal(mu[core[k]], sigma_y[zone[k]]);
##   }
## }
## generated quantities{
##   real delta1;
##   real delta2;
##   real delta3;
##   delta1 = theta[2]-theta[1];
##   delta2 = theta[2]-theta[3];
##   delta3 = theta[1]-theta[3];
## }
## "

## stan.in3 <- function(data = benthic.data, y.col=5,
##                      chains=nchains){ ## no slope
##     n <- dim(data)[1]
##     y <- log(data[,y.col])
##     core <- as.numeric(ordered(data$Station))
##     n.core <- max(core)
##     zone <- as.numeric(ordered(data$Area))
##     n.zone <- max(zone)
##     oo <- order(core)
##     ind <- cumsum(table(core[oo]))
##     Core.zone <- zone[oo][ind] ## each core belongs to which zone
##     stan.dat <- list(K=n, J=n.core, I=n.zone, y=y, core=core,
##                      zone=zone, core_zone=Core.zone)
##     inits <- list()
##     for (i in 1:chains)
##         inits[[i]] <- list( mu = rnorm(n.core), theta=rnorm(n.zone),
##                            mu_hyp=rnorm(1), sigma_i=runif(1),
##                            sigma_y=runif(n.zone), sigma_hyp=runif(1))
##     parameters <- c("mu","theta","mu_hyp", "sigma_y", "sigma_i",
##                     "sigma_hyp",  "delta1", "delta2","delta3")
##     return(list(para=parameters, data=stan.dat, inits=inits,
##                 n.chains=chains))
## }

## ##niters <- 1500000
## ##nthin <- ceiling((niters/2)*nchains/nkeep)

## input.to.stan <- stan.in3()  ## long-runs -- results saved
## stan.fit <- stan_model(model_code=stan_gom)
## fit2keep <- sampling(stan.fit, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains,
##                      control=list(adapt_delta=0.99, max_treedepth=20))
## print(fit2keep)
## rich_stan <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))
## save(rich_stan, file="ch4_rich_stan.RData")

## input.to.stan <- stan.in3(y.col=6)
## fit2keep <- sampling(stan.fit, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains,
##                      control=list(adapt_delta=0.99, max_treedepth=20))
## print(fit2keep)
## abun_stan <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))
## save(abun_stan, file="ch4_abun_stan.RData")

## ## estimated core means from the default model
## tempA <- abun_stan[1:15]
## names(tempA) <- paste("$\\mu_{", 1:15, "}$", sep="")
## tempR <- rich_stan[1:15]
## names(tempR) <- paste("$\\mu_{", 1:15, "}$", sep="")

## delta0R <- rich_stan[25:27]
## delta0A <- abun_stan[25:27]
## names(delta0R) <- c("H-I","H-O","I-O")
## names(delta0A) <- c("H-I","H-O","I-O")

## tikz(file=paste(plotDIRch4, "ModelCores.tex", sep="/"),
##      height=4, width=5.75, standAlone=F)
## par(mfrow=c(1,2), mar=c(3, 3, 2, 1), mgp=c(1.25,0.125,0),
##     las=1, tck=0.01)
## mlplot(tempA, xlab="log abundance")
## mlplot(tempR, xlab="log richness")
## dev.off()

## ## Alternative model 1
## benthic.data$Zone <- with(benthic.data,
##                           factor(as.factor(Area):as.factor(Station)))
## ###
## ### First alternative (One-way ANOVA)
## ###
## stan_gom2 <- " /* Stan Code - GOM model*/
## data{
##   int K;  //total sample size
##   int J;  //number of sediment cores
##   real y[K]; //observed response
##   int core[K]; //core index
## }
## parameters{
##   real mu[J];
##   real mu_hyp;
##   real<lower=0> sigma_y;
##   real<lower=0> sigma_hyp;
## }
## model{
## //  mu_hyp ~ normal(0, 0.5);
##   sigma_hyp ~ normal(0,0.5);
##   sigma_y ~ normal(0,0.5);
##   for (j in 1:J){
##     mu[j] ~ normal(mu_hyp, sigma_hyp);
##   }
##   for (k in 1:K){
##     y[k] ~ normal(mu[core[k]], sigma_y);
##   }
## }
## "
## stan.in4 <- function(data = benthic.data, y.col=5,
##                      chains=nchains){ ## no slope
##   n <- dim(data)[1]
##   y <- log(data[,y.col])
##   y_ave <- mean(y)
##   y_sd <- sd(y)

##   y <- (y-y_ave)/y_sd
##   core <- as.numeric(ordered(data$Station))
##   n.core <- max(core)

##   stan.dat <- list(K=n, J=n.core, y=y, core=core)
##   inits <- list()
##   for (i in 1:chains)
##     inits[[i]] <- list( mu = rnorm(n.core),
##                         mu_hyp=rnorm(1),
##                         sigma_y=runif(1), sigma_hyp=runif(1))
##   parameters <- c("mu","mu_hyp", "sigma_y", "sigma_hyp")
##   return(list(para=parameters, data=stan.dat, inits=inits,
##               n.chains=chains, y_cen=y_ave, y_spd=y_sd))
## }

## input.to.stan <- stan.in4()
## stan.fit2 <- stan_model(model_code=stan_gom2)
## fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains)
## print(fit2keep)
## rich_stan2 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

## ## processing output
## core <- as.numeric(ordered(benthic.data$Station))
## n.core <- max(core)
## zone <- as.numeric(ordered(benthic.data$Area))
## n.zone <- max(zone)
## oo <- order(core)
## ind <- cumsum(table(core[oo]))
## Core.zone <- zone[oo][ind] ## each core belongs to which zone
## input_sd <- input.to.stan$y_spd
## input_mu <- input.to.stan$y_cen

## core1.musR <- input_mu + input_sd*rich_stan2[1:15]
## zone11.Rmu <- mean(core1.musR[Core.zone==1])
## zone12.Rmu <- mean(core1.musR[Core.zone==2])
## zone13.Rmu <- mean(core1.musR[Core.zone==3])

## deltaR11 = zone12.Rmu-zone11.Rmu
## deltaR12 = zone12.Rmu-zone13.Rmu
## deltaR13 = zone11.Rmu-zone13.Rmu

## delta1R <- c(deltaR11, deltaR12, deltaR13)
## names(delta1R) <- c("H-I","H-O","I-O")

## input.to.stan <- stan.in4(y.col=6)
## fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains)
## print(fit2keep)
## abun_stan2 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

## input_sd <- input.to.stan$y_spd
## input_mu <- input.to.stan$y_cen
## core1.musA <- input_mu + input_sd * abun_stan2[1:15]
## zone11.Amu <- mean(core1.musA[Core.zone==1])
## zone12.Amu <- mean(core1.musA[Core.zone==2])
## zone13.Amu <- mean(core1.musA[Core.zone==3])

## deltaA11 = zone12.Amu-zone11.Amu
## deltaA12 = zone12.Amu-zone13.Amu
## deltaA13 = zone11.Amu-zone13.Amu

## delta1A <- c(deltaA11, deltaA12, deltaA13)
## names(delta1A) <- c("H-I","H-O","I-O")

## ## The second alternative model (nested ANOVA)
## ## Referencing the lme4 book

## stan_gom3 <- " /* Stan Code - GOM model*/
## data{
##   int K;  //total sample size
##   int J;  //number of sediment cores
##   int I;  //number of zones
##   real y[K]; //observed response
##   int core[K]; //core index
##   int zone[K]; //zone index
## }
## parameters{
##   real muK[J];
##   real muZ[I];
##   real mu0;
##   real<lower=0> sigmaY;
##   real<lower=0> sigmaK;
##   real<lower=0> sigmaZ;
## }
## model{
##   sigmaK ~ normal(0,0.5);
##   sigmaZ ~ normal(0,0.5);
##   sigmaY ~ normal(0,0.5);
##   for (i in 1:I){
##     muZ[i] ~ normal(0, sigmaZ);
##   }
##   for (j in 1:J){
##     muK[j] ~ normal(0, sigmaK);
##   }
##   for (k in 1:K){
##     y[k] ~ normal(mu0+muK[core[k]]+muZ[zone[k]], sigmaY);
##   }
## }
## generated quantities{
##   real delta1;
##   real delta2;
##   real delta3;
##   delta1 = muZ[2]-muZ[1];
##   delta2 = muZ[2]-muZ[3];
##   delta3 = muZ[1]-muZ[3];
## }
## "
## stan.in5 <- function(data = benthic.data, y.col=5, chains=nchains){
##   n <- dim(data)[1]
##   y <- log(data[,y.col])
##   y_ave <- mean(y)
##   y_sd <- sd(y)

##   y <- (y-y_ave)/y_sd
##   core <- as.numeric(ordered(data$Station))
##   n.core <- max(core)
##   zone <- as.numeric(ordered(data$Area))
##   n.zone <- max(zone)

##   stan.dat <- list(K=n, J=n.core, I=n.zone, y=y, core=core, zone=zone)
##   inits <- list()
##   for (i in 1:chains)
##     inits[[i]] <- list( muK = rnorm(n.core),
##                         muZ=rnorm(n.zone),
##                         mu0=rnorm(1),
##                        sigmaY=runif(1), sigmaK=runif(1),
##                        sigmaZ=runif(1))
##   parameters <- c("mu0","muK", "muZ", "sigmaY", "sigmaK",
##                   "sigmaZ", "delta1", "delta2", "delta3")
##   return(list(para=parameters, data=stan.dat, inits=inits,
##               n.chains=chains, y_cen=y_ave, y_spd=y_sd))
## }

## input.to.stan <- stan.in5()
## stan.fit3 <- stan_model(model_code=stan_gom3)
## fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains,
##                      control=list(adapt_delta=0.99, max_treedepth=15))
## print(fit2keep)
## rich_stan3 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

## input_sd <- input.to.stan$y_spd
## input_mu <- input.to.stan$y_cen
## zone2.musR <- input_sd*rich_stan3[17:19]
## core2.musR <- input_mu+input_sd*rich_stan3[2:16]

## zone21.Rmu <- mean(core2.musR[Core.zone==1]) + zone2.musR[1]
## zone22.Rmu <- mean(core2.musR[Core.zone==2]) + zone2.musR[2]
## zone23.Rmu <- mean(core2.musR[Core.zone==3]) + zone2.musR[3]

## deltaR21 = zone22.Rmu-zone21.Rmu
## deltaR22 = zone22.Rmu-zone23.Rmu
## deltaR23 = zone21.Rmu-zone23.Rmu

## delta2R <- c(deltaR21, deltaR22, deltaR23)
## names(delta2R) <- c("H-I","H-O","I-O")

## input.to.stan <- stan.in5(y.col=6)
## fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
##                      init=input.to.stan$inits,
##                      pars = input.to.stan$para,
##                      iter=niters, thin=nthin,
##                      chains=input.to.stan$n.chains,
##                      control=list(adapt_delta=0.99, max_treedepth=15))
## print(fit2keep)
## abun_stan3 <- rvsims(as.matrix(as.data.frame(extract(fit2keep))))

## input_sd <- input.to.stan$y_spd
## input_mu <- input.to.stan$y_cen
## zone2.musA <- input_sd*abun_stan3[17:19]
## core2.musA <- input_mu+input_sd*abun_stan3[2:16]


## ### saving all model outputs for later use
## save(rich_stan2, abun_stan2, rich_stan3, abun_stan3,
##      core1.musR, core1.musA, zone2.musR, zone2.musA,
##      core2.musR, core2.musA,
##      file="GOMalt23.RData")

## zone21.Amu <- mean(core2.musA[Core.zone==1]) + zone2.musA[1]
## zone22.Amu <- mean(core2.musA[Core.zone==2]) + zone2.musA[2]
## zone23.Amu <- mean(core2.musA[Core.zone==3]) + zone2.musA[3]

## deltaA21 = zone22.Amu-zone21.Amu
## deltaA22 = zone22.Amu-zone23.Amu
## deltaA23 = zone21.Amu-zone23.Amu

## delta2A <- c(deltaA21, deltaA22, deltaA23)
## names(delta2A) <- c("H-I","H-O","I-O")

## ## Figures
## tikz(file=paste(plotDIRch4, "ModelCoresR.tex", sep="/"),
##      height=5, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(tempR, xlab="Original model", top.axis=F)
## mlplot(core1.musR, xlab="Alternative model 1",
##        main="log Richness", axes=F)
## axis(1)
## mlplot(core2.musR, xlab="Alternative model 2", axes=F)
## axis(1)
## dev.off()

## tikz(file=paste(plotDIRch4, "ModelCoresA.tex", sep="/"),
##      height=5, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(tempA, xlab="Original model", top.axis=F)
## mlplot(core1.musA, xlab="Alternative model 1",
##        main="log Abundance", axes=F)
## axis(1)
## mlplot(core2.musA, xlab="Alternative model 2", axes=F)
## axis(1)
## dev.off()

## zone0R.thetas <- rich_stan[16:18]
## names(zone0R.thetas) <- c("Inshore", "Hypoxic","Offshore")
## zone0A.thetas <- abun_stan[16:18]
## names(zone0A.thetas) <- c("Inshore", "Hypoxic","Offshore")

## zone1R.thetas <- c(zone11.Rmu, zone12.Rmu, zone13.Rmu)
## names(zone1R.thetas) <- c("Inshore", "Hypoxic","Offshore")

## zone1A.thetas <- c(zone11.Amu, zone12.Amu, zone13.Amu)
## names(zone1A.thetas) <- c("Inshore", "Hypoxic","Offshore")

## zone2R.thetas <- c(zone21.Rmu, zone22.Rmu, zone23.Rmu)
## names(zone2R.thetas) <- c("Inshore", "Hypoxic","Offshore")

## zone2A.thetas <- c(zone21.Amu, zone22.Amu, zone23.Amu)
## names(zone2A.thetas) <- c("Inshore", "Hypoxic","Offshore")

## tikz(file=paste(plotDIRch4, "ModelZonesR.tex", sep="/"),
##      height=3, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(zone0R.thetas, xlab="Original model", top.axis=F, xlim=c(2,3.5))
## abline(v=0)
## mlplot(zone1R.thetas, xlab="Alternative model 1", main="log Richness",
##        axes=F, xlim=c(2,3.5))
## axis(1)
## abline(v=0)
## mlplot(zone1R.thetas, xlab="Alternative model 2", axes=F, xlim=c(2,3.5))
## axis(1)
## abline(v=0)
## dev.off()

## tikz(file=paste(plotDIRch4, "ModelZonesA.tex", sep="/"),
##      height=3, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(zone0A.thetas, xlab="Original model", top.axis=F, xlim=c(3.6,6))
## mlplot(zone1A.thetas, xlab="Alternative model 1", main="log Abundance",
##        axes=F, xlim=c(3.6,6))
## axis(1)
## mlplot(zone1A.thetas, xlab="Alternative model 2", axes=F, xlim=c(3.6,6))
## axis(1)
## dev.off()

## tikz(file=paste(plotDIRch4, "ModelZonesDelA.tex", sep="/"),
##      height=3, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(delta0A, xlab="Original model 1", top.axis=F, xlim=c(-2,1.5))
## abline(v=0)
## mlplot(delta1A, xlab="Alternative model 1",
##        main="log Abundance", axes=F,
##        xlim=c(-2,1.5))
## axis(1)
## abline(v=0)
## mlplot(delta2A, xlab="Alternative model 2", axes=F, xlim=c(-2,1.5))
## axis(1)
## abline(v=0)
## dev.off()

## tikz(file=paste(plotDIRch4, "ModelZonesDelR.tex", sep="/"),
##      height=3, width=5.5, standAlone=F)
## par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
## mlplot(delta0R, xlab="Original model", top.axis=F, xlim=c(-1.5,0.5))
## abline(v=0)
## mlplot(delta1R, xlab="Alternative model 1",
##        main="log Richness", axes=F, xlim=c(-1.5,0.5))
## abline(v=0)
## axis(1)
## mlplot(delta2R, xlab="Alternative model 2", axes=F, xlim=c(-1.5,0.5))
## axis(1)
## abline(v=0)
## dev.off()

## Neuse River Estuary
## Neuse River Estuary Example (long-term mean distribution) ##
## prior parameters
## Model 1 TN (from chapter 3)
pst <- list(alpha =32,
            beta = 28,
            nn = 45,
            mu = 7.6)

NGpost <- function(x, alpha, beta, n0, mu0){
    x_bar <- mean(x)
    n <- length(x)
    s2 <- sd(x)^2
    return(list(nn=n+n0, mu=(n*x_bar+n0*mu0)/(n+n0),
           alpha = alpha+n/2,
           beta = beta+0.5*(n*s2 + (n*n0)*(x_bar-mu0)^2/(n+n0))))
}

tmp <- neuse$Year>=1992 & neuse$Year <= 2000 & ##neuse$SECTION=="UPPER")|
        neuse$SECTION=="MIDDLE"
neuse2 <- neuse[tmp,]

tikz(file=paste(plotDIRch3, "neuseTN.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$TOTN), prob=T, nclass=20,
     axes=F, xlab="Total Nitrogen ($\\mu$g/L)",
     ylim=c(0,1.5), main="", ylab="")
axis(1, at = log(c(10, 100, 250, 500, 1000, 2500, 5000)),
                 labels=c(10, 100,250, 500, 1000, 2500, 5000))
pst <- list(alpha =32,
            beta = 28,
            nn = 45,
            mu = 7.6)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2
curve(dt((x-pr_mut)/pr_sigmat, df=pr_df)/pr_sigmat, add=T, lwd=3)

for (i in 1992:2000){
    tmp <- neuse2$Year==i & !is.na(neuse2$TOTN)
    pst <- NGpost(log(neuse2$TOTN[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
    curve(dt((x-mu_t)/sigma_t, df=df_t)/sigma_t, lty=i-1990, add=T)
}
dev.off()


tikz(file=paste(plotDIRch3, "neuseCHLA.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$SURFCHLA), prob=T, nclass=20,
     axes=F, xlab="Chlorophyll a ($\\mu$g/L)",
     ylim=c(0,1), main="", ylab="")
axis(1, at = log(c(0.1,1, 5, 10, 100, 250, 500)),
                 labels=c(0.1,1, 5, 10, 100,250, 500))
pst <- list(alpha =100,
            beta = 40,
            nn = 40,
            mu = 3.4)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2
curve(dt((x-pr_mut)/pr_sigmat, df=pr_df)/pr_sigmat, add=T, lwd=3)

for (i in 1992:2000){
    tmp <- neuse2$Year==i & !is.na(neuse2$SURFCHLA)
    pst <- NGpost(log(neuse2$SURFCHLA[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
    curve(dt((x-mu_t)/sigma_t, df=df_t)/sigma_t, lty=i-1990, add=T)
}
dev.off()


neuse <- read.csv(paste(dataDIR, "NeuseEstuary.csv", sep="/"))
neuse[neuse<0] <- NA

neuse$RDate <- as.Date(as.character(neuse$DATE), format="%m/%d/%Y")
neuse$Year <- as.numeric(format(neuse$RDate, "%Y"))

## stan model (normal with conjugate prior for hyper-parameters)
stan.mod <- "
  data{
    int N; // observation count
    real y[N];
    real mu0; // prior parameters
    real n0;
    real alpha;
    real beta;
  }
  parameters{
    real<lower=0> sigma_sq;
    real mu;
  }
  model{
    sigma_sq ~ inv_gamma(alpha, beta);
    mu ~ normal(mu0, sqrt(sigma_sq/n0));
    y ~ normal(mu, sqrt(sigma_sq));
  }
"

stan_in <- function(data=neuse, yVAR="SURFCHLA", Sub=NULL,
                    chains=nchains,
                    mu0=0,n0=1, alpha=0.1,beta=0.1){
    if (is.null(Sub))
        y <- log(data[,yVAR])
    else y <- log(data[Sub, yVAR])
    y <- y[!is.na(y)]
    stan_data <- list(N=length(y), y=y, mu0=mu0,n0=n0,
                      alpha=alpha, beta=beta)
    stan_inits <- list()
    for (i in 1:chains)
        stan_inits[[i]] <- list(sigma_sq=runif(1), mu=rnorm(1))
    return (list(data=stan_data, inits=stan_inits,
                 pars=c("mu", "sigma_sq"), n.chains=chains))
}

uniqueYr <- as.numeric(levels(ordered(neuse$Year)))
input.to.stan <- stan_in(Sub=neuse$Year==uniqueYr[1])
fit <- stan_model(model_code = stan.mod)
fit2keep <- sampling(fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$pars,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$n.chains)

fitcoef79 <- rvsims(as.matrix(as.data.frame(extract(fit2keep,
                                                    permute=T))))

priors_normal <- function(fit_coef=summary(fitcoef79)){
    mu <- fit_coef$name=="mu"
    sig2 <- fit_coef$name=="sigma_sq"
    ex <- fit_coef$mean[mu]
    vx <- fit_coef$sd[mu]^2
    esig2 <- fit_coef$mean[sig2]
    vsig2 <- fit_coef$sd[sig2]^2
    return(list(mu0 = ex, n0 = ex/vx, alpha=(alpha <- 2+esig2^2/vsig2),
                beta=(alpha-1)*esig2))
}


## now automate the updating process
StanFit <- fitcoef79
StanSave <- list()
StanSave[[1]] <- StanFit
j <- 1
for (i in uniqueYr[-1]){
    j <- j+1
    pr <- priors_normal(summary(StanFit))
    input.to.stan <- stan_in(Sub=neuse$Year==i, mu0=pr$mu0,
                             n0=pr$n0, alpha=pr$alpha,
                             beta=pr$beta)
    fit2keep <- sampling(fit, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$pars,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    StanFit <- rvsims(as.matrix(as.data.frame(extract(fit2keep,
                                                      permute=T))))
    StanSave[[j]] <- StanFit
}

nn <- length(StanSave)
Stan_mu <- StanSave[[1]][1]
Stan_sig <- StanSave[[1]][2]
for (i in 2:nn){
    Stan_mu <- c(Stan_mu, StanSave[[i]][1])
    Stan_sig <- c(Stan_sig, sqrt(StanSave[[i]][2]))
}

tikz(file=paste(plotDIRch4, "neuse_old.tex", sep="/"),
     height=4, width=5, standAlone=F)
par(mfrow=c(2,1), oma=c(3,3,3,3), mar=c(0,3,0,0),
    mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(Stan_mu, axes=F, ylab="$\\mu$")
axis(2)
axis(3, labels=NA)
axis(4, labels=NA)
box()
plot(Stan_sig, axes=F, ylab="$\\sigma$")
axis(1, at=1:22, labels=uniqueYr)
axis(2, labels=NA)
axis(4)
box()
dev.off()

### Linear Regression
## Lake Trout (linear model)
laketrout <- read.csv(paste(dataDIR, "laketrout2.csv", sep="/"))
laketrout$size_cm <- round(laketrout$length * 2.54, 1)

nsims <- 5000
lake.lm1 <- lm(log(pcb)~size_cm, data=laketrout)
summ <- summary(lake.lm1)
beta_hat <- coef(lake.lm1)
V_beta <- summ$cov.unscaled
n <- sum(summ$df[1:2])
k <- summ$df[1]
sigma <- summ$sigma

## random numbers of sigma
rsigma <- sigma * sqrt((n - k)/rchisq(nsims, n-k))
hist(rsigma)
## random numbers for beta
beta <- array(NA, c(nsims, k))
dimnames(beta) <- list(NULL, rownames(beta_hat))
for (i in 1:nsims){
    beta[i, ] <- MASS::mvrnorm(1, beta_hat, V_beta *
                                 rsigma[i]^2)
}

## predictive distribution
setnsims(5000)
tildeX <- c(50, 60, 70, 80, 90)
lakeM1_rv <- posterior(lake.lm1)
lakeM1_pred <- rvnorm(1, lakeM1_rv$beta[1]+lakeM1_rv$beta[2]*tildeX, lakeM1_rv$sigma)
Pr(exp(lakeM1_pred) > 1)

lm1_pred <- predict(lake.lm1, newdata=data.frame(size_cm=tildeX),
                    se.fit=T, interval="prediction")

