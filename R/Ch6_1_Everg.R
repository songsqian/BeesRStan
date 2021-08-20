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
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

everg_stan <- "
data {
  int<lower=0> J; // number of schools
  real y[J]; // estimated treatment effects
  real<lower=0> sigma[J]; // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
"

y.hat <- c(19.2,  13,  13, 12.4, 8.2, 19.9, 18.3,
           15.6, 14.8, 14.5, 23.5, 15.3)
sigma.hat <- c( 1.6, 6.4, 9.6, 16.2, 0.8,  4.2,  1.2,
                0.9,  2.1, 10.2, 26.3, 10.3)/4

metrics <- c("Mat Cover","BCD","% Tolerant Species","% Sensitive Species",
             "% Predators","% Crustacean","% Oligchaeta","Total Utricularia",
             "Utricularia purpurea","% diatom (stem)","% diatom (plexi)",
             "% diatom (mat)")
metricsTEX <- c("Mat Cover","BCD","\\% Tol Sp",
                "\\% Sen Sp", "\\% Pred","\\% Crust",
                "\\% Oligchaeta","Tot Utr","Utr P.",
                "\\% diatom (stem)","\\% diatom (plexi)",
                "\\% diatom (mat)")


everg_in <- function(y=y.hat, sig=sigma.hat, n.chains=nchains){
  J <- length(y)
  data <- list(y=y, sigma=sig, J=J)
  inits<-list()
  for (i in 1:n.chains)
    inits[[i]] <- list(eta=rnorm(J), mu=rnorm(1), tau=runif(1))
  pars <- c("theta", "mu", "eta", "tau")
  return(list(data=data, inits=inits, pars=pars, chains=n.chains))
}

input.to.stan <- everg_in()
fit1 <- stan_model(model_code = everg_stan)
fit2keep <- sampling(fit1, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(max_treedepth=25))

print(fit2keep)

everg_fit1 <- rvsims(as.matrix(
    as.data.frame(rstan::extract(fit2keep, permuted=T))))

## raw data plot
tikz(file=paste(plotDIRch6, "chngpRaw.tex", sep="/"),
     width=3.5, height=3,
     standAlone=F)
par(mar=c(3, 7, 1, 0.5), mgp=c(1.25,0.25,0),tck=0.01)
plot(range(y.hat-2*sigma.hat, y.hat+2*sigma.hat),
     c(1,length(y.hat)), type="n",
     xlab="TP Threshold", ylab=" ", axes=F)
axis(1, cex.axis=0.75)
axis(2, at=1:length(y.hat), labels=metricsTEX, las=1, cex.axis=0.75)
segments(x0=y.hat+sigma.hat, x1=y.hat-sigma.hat,
         y0=1:length(y.hat), y1=1:length(y.hat), lwd=3)#,
#         col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4))
segments(x0=y.hat+2*sigma.hat, x1=y.hat-2*sigma.hat,
         y0=1:length(y.hat), y1=1:length(y.hat), lwd=1)#,
#         col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4))
points(x=y.hat, y=1:length(y.hat))
dev.off()

## shrinkage effect
everg_theta <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep,
                                                             permuted=T,
                                                             pars="theta"))))
everg_mu <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep,
                                                             permuted=T,
                                                             pars="mu"))))
everg_tau <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep,
                                                             permuted=T,
                                                             pars="tau"))))
theta <- summary(everg_theta)
mu <- summary(everg_mu)
tau <- summary(everg_tau)

tikz(file=paste(plotDIRch6, "chngpShrink.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3, 7, 1, 0.5), mgp=c(1.25,0.125,0), tck=0.01)
plot(range(y.hat-1*sigma.hat, y.hat+1*sigma.hat),
     c(1,length(y.hat)), type="n",
     xlab="TP Threshold", ylab=" ", axes=F)
axis(1)
axis(2, at=seq(1,length(y.hat)), labels=metricsTEX, las=1)
segments(x0=y.hat+sigma.hat, x1=y.hat-sigma.hat,
         y0=seq(1,length(y.hat))-0.125,
         y1=seq(1,length(y.hat))-0.125,
         lwd=1, lty=2)
## col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4),
segments(x0=theta$"25%", x1=theta$"75%",
         y0=(seq(1,length(y.hat)))+0.125,
         y1=(seq(1,length(y.hat)))+0.125)
##         col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4))
points(x=y.hat, y=seq(1,length(y.hat))-0.125, cex=0.5)
##       col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4),
points(x=theta$mean, y=0.125+(seq(1,length(y.hat))), cex=0.5)
##       pch=16,col=c(1,1, 2,2,2,2,2, 3,3, 4, 4,4))
abline(v=mu$mean)
dev.off()

## mu versus N(mu, tau)
mu_tau <- rvnorm(1, everg_mu, everg_tau)
p1 <- hist(sims(everg_mu)[,1], freq=F)
p2 <- hist(sims(mu_tau)[,1], nclass=35)

tikz(file=paste(plotDIRch6, "evergHist.tex", sep="/"),
     width=3, height=3, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0), tck=0.01)
plot(p1, col=rgb(0.1,0.1,.1,1/4),
     xlim=c(0,30), ylim=c(0,0.35), freq=F,
     xlab="TP Threshold ($\\mu$g/L)", main="")  # first histogram
plot(p2, col=rgb(.7,.7,.7,1/4),
     xlim=c(0,30), ylim=c(0,0.35), freq=F, add=T)  # second
box()
dev.off()
