source("FrontMatter.R")

## Chapter 4

plotDIRch4 <- paste(plotDIR, "chapter4", "figures", sep="/")

packages(rv)
packages(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

LCAAPsoil.data <- read.csv(file=paste(dataDIR, "LCAAPsoil.csv", sep="/"))

lead.data <- LCAAPsoil.data[LCAAPsoil.data$chemical.name=="Lead",]
lead.data$final.detection.limit
table(lead.data$bkg.char.group)
oo <- order(lead.data$bkg.char.group)
cbind(lead.data$bkg.char.group, lead.data$soil.type)[oo,]

mercury.data <- LCAAPsoil.data[LCAAPsoil.data$chemical.name=="Mercury",]
sqrt(var(log(mercury.data$final.result.value), na.rm=T))
table(mercury.data$bkg.char.group)
oo <- order(mercury.data$bkg.char.group)
cbind(mercury.data$bkg.char.group, mercury.data$soil.type)[oo,]

######################
### The Stan Model ###
######################

mixture1 <- "
data{
    int N;
    vector[N] x;
}
parameters{
    real mu1;
    real mu2;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
    real<lower=0,upper=1> p;
}
model {
    sigma1 ~ normal(0,5);
    sigma2 ~ normal(0,5);
    mu1 ~ normal(0,5);
    mu2 ~ normal(0,5);
    target += log_mix(p, normal_lpdf(x|mu1, sigma1),
                      normal_lpdf(x|mu2, sigma2));
}
"

fit <- stan_model(model_code=mixture1)

mixture_input <- function(data=lead.data, n.chains=nchains)
{
    y <- data$final.result.value
    tmp <- is.na(y)
    y <- y[!tmp]
    n <- length(y)
    data_in <- list(N=n,x=log(y))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(p=runif(1,0.5,1),mu1=rnorm(1),mu2=rnorm(1),
                           sigma1=runif(1), sigma2=runif(1))
    pars <- c("mu1","mu2","sigma1","sigma2", "p")
    return(list(data=data_in, inits=inits, pars=pars,nchains=n.chains))
}

set.seed(201)
input.to.stan <- mixture_input()
fit2keep1 <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     chains=input.to.stan$nchains,
                     iter=niters,
                     thin=nthin)
print(fit2keep1)
pairs(fit2keep1, pars=c("mu1","mu2", "p"))
run1_p <- as.data.frame(extract(fit2keep1, pars="p", permuted=F))
run1_mu1 <- as.data.frame(extract(fit2keep1, pars="mu1", permuted=F))
run1_mu2 <- as.data.frame(extract(fit2keep1, pars="mu2", permuted=F))

tikz(file=paste(plotDIRch4, "Lead_mixture1.tex", sep="/"),
     height=5, width=4, standAlone=F)
par(mfrow=c(3,2),mar=c(3,2,1.5,1),mgp=c(1.25,0.125,0),las=1,tck=0.01)
hist(run1_p[,1], main="chain 1", xlab="$p$", ylab="")
hist(run1_p[,2], main="chain 2", xlab="$p$", ylab="")
hist(run1_mu1[,1], xlim=c(-15,15), main="", xlab="$\\mu_1$", ylab="")
hist(run1_mu1[,2], xlim=c(-15,15), main="", xlab="$\\mu_2$", ylab="")
hist(run1_mu2[,1], xlim=c(-15,15), main="", xlab="$\\mu_1$", ylab="")
hist(run1_mu2[,2], xlim=c(-15,15), main="", xlab="$\\mu_2$", ylab="")
dev.off()


mixture2 <- "
data{
    int N;
    vector[N] x;
}
parameters{
    real mu1;
    real<lower=0> delta;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
    real<lower=0,upper=1> p;
}
model {
    sigma1 ~ normal(0,5);
    sigma2 ~ normal(0,5);
    mu1 ~ normal(0,5);
    delta ~ normal(0,2);
    target += log_mix(p, normal_lpdf(x|mu1, sigma1),
                      normal_lpdf(x|mu1+delta, sigma2));
}
"

fit2 <- stan_model(model_code=mixture2)

mixture_input2 <- function(data=lead.data, n.chains=nchains)
{
    y <- data$final.result.value
    tmp <- is.na(y)
    y <- y[!tmp]
    n <- length(y)
    data_in <- list(N=n,x=log(y))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(p=runif(1,0.5,1),mu1=rnorm(1),delta=abs(rnorm(1)),
                           sigma1=runif(1), sigma2=runif(1))
    pars <- c("mu1","delta","sigma1","sigma2", "p")
    return(list(data=data_in, inits=inits, pars=pars,nchains=n.chains))
}

input.to.stan <- mixture_input2()
fit2keep2 <- sampling(fit2, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     chains=input.to.stan$nchains,
                     iter=niters,
                     thin=nthin)
print(fit2keep2)
pairs(fit2keep2, pars=c("mu1","delta", "p"))

run2_p <- as.data.frame(extract(fit2keep2, pars="p", permuted=F))
run2_mu1 <- as.data.frame(extract(fit2keep2, pars="mu1", permuted=F))
run2_delta <- as.data.frame(extract(fit2keep2, pars="delta", permuted=F))

tikz(file=paste(plotDIRch4, "Lead_mixture1.tex", sep="/"),
     height=5, width=4, standAlone=F)
par(mfrow=c(3,2),mar=c(3,2,1.5,1),mgp=c(1.25,0.125,0),las=1,tck=0.01)
hist(run1_p[,1], main="chain 1", xlab="$p$", ylab="")
hist(run1_p[,2], main="chain 2", xlab="$p$", ylab="")
hist(run1_mu1[,1], xlim=c(-15,15), main="", xlab="$\\mu_1$", ylab="")
hist(run1_mu1[,2], xlim=c(-15,15), main="", xlab="$\\mu_2$", ylab="")
hist(run1_mu2[,1], xlim=c(-15,15), main="", xlab="$\\mu_1$", ylab="")
hist(run1_mu2[,2], xlim=c(-15,15), main="", xlab="$\\mu_2$", ylab="")
dev.off()


tikz(file=paste(plotDIRch4, "Lead_mixture2.tex", sep="/"),
     height=5, width=4, standAlone=F)
par(mfrow=c(3,2),mar=c(3,2,1.5,1),mgp=c(1.25,0.125,0),las=1,tck=0.01)
hist(run2_p[,1], main="chain 1", xlab="$p$", ylab="")
hist(run2_p[,2], main="chain 2", xlab="$p$", ylab="")
hist(run2_mu1[,1], xlim=c(1,4), main="", xlab="$\\mu_1$", ylab="")
hist(run2_mu1[,2], xlim=c(1,4), main="", xlab="$\\mu_1$", ylab="")
hist(run2_delta[,1], main="", xlab="$\\delta$", ylab="")
hist(run2_delta[,2], main="", xlab="$\\delta$", ylab="")
dev.off()

mixture2_5 <- "
data{
    int N;
    vector[N] x;
}
parameters{
    real mu1;
    real<lower=0> delta;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
    real<lower=0.5,upper=1> p;
}
model {
    sigma1 ~ normal(0,5);
    sigma2 ~ normal(0,5);
    mu1 ~ normal(0,5);
    delta ~ normal(0,2);
    target += log_mix(p, normal_lpdf(x|mu1, sigma1),
                      normal_lpdf(x|mu1+delta, sigma2));
}
"

fit2_5 <- stan_model(model_code=mixture2_5)

input.to.stan <- mixture_input2()
fit2keep2_5 <- sampling(fit2_5, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     chains=input.to.stan$nchains,
                     iter=niters,
                     thin=nthin)
print(fit2keep2_5)
pairs(fit2keep2_5, pars=c("mu1","delta", "p"))

run25_p <- as.data.frame(extract(fit2keep2_5, pars="p", permuted=F))
run25_mu1 <- as.data.frame(extract(fit2keep2_5, pars="mu1", permuted=F))
run25_delta <- as.data.frame(extract(fit2keep2_5, pars="delta", permuted=F))

tikz(file=paste(plotDIRch4, "Lead_mixture25.tex", sep="/"),
     height=5, width=4, standAlone=F)
par(mfrow=c(3,2),mar=c(3,2,1.5,1),mgp=c(1.25,0.125,0),las=1,tck=0.01)
hist(run25_p[,1], main="chain 1", xlab="$p$", ylab="")
hist(run25_p[,2], main="chain 2", xlab="$p$", ylab="")
hist(run25_mu1[,1], xlim=c(1,4), main="", xlab="$\\mu_1$", ylab="")
hist(run25_mu1[,2], xlim=c(1,4), main="", xlab="$\\mu_1$", ylab="")
hist(run25_delta[,1], main="", xlab="$\\delta$", ylab="")
hist(run25_delta[,2], main="", xlab="$\\delta$", ylab="")
dev.off()

run2_mu1RV <- rvsims(as.matrix(as.data.frame(extract(fit2keep2_5,
                                                     pars="mu1",
                                                     permuted=T))))
run2_deltaRV <- rvsims(as.matrix(as.data.frame(extract(fit2keep2_5,
                                                       pars="delta",
                                                       permuted=T))))
mixture3 <- "
data{
    int N;
    vector[N] x;
    int sites[N];
    int Ns;
}
parameters{
    real mu1;
    real<lower=0> delta[Ns];
    real<lower=0> sigma1;
    real<lower=0> sigma2;
    real<lower=0> mu_delta;
    real<lower=0> sigma_delta;
    real<lower=0.5,upper=1> p;
}
model {
    sigma1 ~ normal(0,5);
    sigma2 ~ normal(0,5);
    sigma_delta ~ normal(0,5);
    mu1 ~ normal(0,5);
    mu_delta ~ normal(0,2);
    delta ~ normal(mu_delta,sigma_delta);
    for (i in 1:N)
      target += log_mix(p, normal_lpdf(x[i]|mu1, sigma1),
                        normal_lpdf(x[i]|mu1+delta[sites[i]], sigma2));
}
"

fit3 <- stan_model(model_code=mixture3)

mixture_input3 <- function(data=lead.data, n.chains=nchains)
{
    y <- data$final.result.value
    tmp <- is.na(y)
    y <- y[!tmp]
    sites <- as.numeric(ordered(data$bkg.char.group[!tmp]))
    n <- length(y)
    ns <- max(sites)
    data_in <- list(N=n,Ns=ns, x=log(y),sites=sites)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(p=runif(1,0.5,1),mu1=rnorm(1),delta=abs(rnorm(ns)),
                           sigma_delta=runif(1), mu_delta=abs(rnorm(1)),
                           sigma1=runif(1), sigma2=runif(1))
    pars <- c("mu1","delta","sigma1","sigma2","p","mu_delta","sigma_delta")
    return(list(data=data_in, inits=inits, pars=pars,nchains=n.chains))
}

input.to.stan <- mixture_input3()
fit2keep3 <- sampling(fit3, data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     chains=input.to.stan$nchains,
                     iter=niters*10,
                     thin=nthin*10)
print(fit2keep3)
pairs(fit2keep3, pars=c("mu1","mu_delta", "sigma_delta","p"))

run3_pars <-rvsims(as.matrix(as.data.frame(extract(fit2keep3,permuted=T))))
run3_mu1RV <- run3_pars[1]

tmp=c("$\\mu_{b,m_2}$", "$\\mu_{c,m_2}$","$\\mu_{b,m_3}$",
      paste("$\\mu_{c,", 1:19, "}$", sep=""))
tikz(file=paste(plotDIRch4, "lead_mixture4.tex", sep="/"),
     height=4.5, width=3, standAlone=F)
par(mar=c(3, 5, 1, 1), mgp=c(1.25, 0.125,0), las=1, tck=0.01)
mlplot(c(run2_mu1RV, run2_mu1RV+run2_deltaRV, run3_pars[1], run3_pars[1]+run3_pars[2:20]), xlim=c(2,10),
       xlab="log Pb concentrations", axes=F, cex=0.75)
axis(1)
axis(2, at=1:22, labels=tmp)
box()
dev.off()

