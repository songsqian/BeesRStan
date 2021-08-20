source("FrontMatter.R")

## Chapter 5
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")
## simulation
packages(rv)
packages(rstan)
packages(MASS)
packages(pscl)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

####################
## Zero Inflation ##
####################

## simulated data
theta <- 0.25
lambda <- 2.5
n <- 100

y <- rpois(n, lambda)
zr <- rbinom(n, 1, theta)
## the bare basic ZIP model

y <- y*zr


zip1 <- "
data{
  int<lower=1> n0;
  int<lower=1> np;
  int<lower=1> yp[np];
}
parameters{
  real<lower=0,upper=1> theta;
  real<lower=0> lambda;
}
model{
  theta ~ beta(1,1);
  lambda ~ normal(0,2);
  target += n0*log_sum_exp(log(theta), log1m(theta)-lambda);
  target += np*log1m(theta)+poisson_lpmf(yp|lambda);
}
"
ZIP_in <- function(y, n.chains=nchains){
    n0 <- sum(y==0)
    np <- sum(y!=0)
    yp <- y[y!=0]
    data <- list(n0=n0, np=np, yp=yp)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(theta=runif(1), lambda=abs(rnorm(1)))
    paras <- c("theta", "lambda")
    return(list(data=data, init=inits, para=paras, nchains=n.chains))
}

fit <- stan_model(model_code=zip1)
input.to.stan <- ZIP_in(y)
fit2keep <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)

pairs(fit2keep, pars=c("theta","lambda"))
print(fit2keep)

##theta <- rstan::extract(fit2keep, pars="theta", permuted=T)
##lambda <- rstan::extract(fit2keep, pars="lambda", permuted=T)
##cor(lambda[[1]], theta[[1]])

nsims <- 1000
simulated <- matrix(rpois(n*nsims, lambda)*rbinom(n*nsims, 1, theta),
                    ncol=n, nrow=nsims)

mle_sim <- t(apply(simulated, 1, FUN=function(x){
    temp <- zeroinfl(x~1|1)
    return(c(invlogit(temp$coef[[2]]), exp(temp$coef[[1]])))
}))

stan_sim <- matrix(0, nrow=nsims, ncol=2)
for (i in 1:nsims){
    print(paste(i, "of ", nsims))
    x <- rpois(n, lambda) * rbinom(n, 1, theta)
    input.to.stan <- ZIP_in(x)
    keep <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)
    temp <- rstan::extract(keep, pars=c("theta", "lambda"), permuted=T)
    stan_sim[i,] <- c(mean(temp$theta), mean(temp$lambda))
}

sim_data <- as.data.frame(rbind(stan_sim, mle_sim))
names(sim_data) <- c("$\\theta$", "$\\lambda$")
save(sim_data, file="ZIP_sim.RData")

sim_data$method <- c(rep("Bayesian", 1000), rep("MLE", 1000))

tikz(file=paste(plotDIRch5, "simzip.tex", sep="/"),
     height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(3, 3, 2, 2),
    mgp=c(1.25,0.125,0), las=1, tck=0.01)
hist(stan_sim[,1], col="gray", border="white", main="", xlab="",
     ylab="Bayesian",axes=F)
box()
axis(3)
abline(v=0.75, lwd=2)

hist(stan_sim[,2], col="gray", border="white", main="", xlab="", ylab="",
     axes=F, xlim=c(1,3.5))
box()
axis(3)
abline(v=2.5, lwd=2)

hist(mle_sim[,1], col="gray", border="white", main="", xlab="$\\theta$",
     ylab="MLE", axes=F)
abline(v=0.75, lwd=2)
axis(1)
box()

hist(mle_sim[,2], col="gray", border="white", main="", xlab="$\\lambda$",
     ylab="", axes=F, xlim=c(1,3.5))
abline(v=2.5, lwd=2)
axis(1)
box()
dev.off()

### EUSE example ###
##  reading data ##

rthE <- read.csv(paste(dataDIR, "rthE.csv",sep="/"))
rthP <- read.csv(paste(dataDIR, "rthP.csv",sep="/"))
rthT <- read.csv(paste(dataDIR, "rthT.csv",sep="/"))
## the following two should be both 0
sum(rthE$SampleID-rthT$SampleID)
sum(rthE$SampleID-rthP$SampleID)

tikz(file=paste(plotDIRch5, "euseNLCD2.tex", sep="/"),
     height=2.5, width=4.25, standAlone=F)
bwplot(NLCD2~City, data=rthE, ylab="\\% developed land")
dev.off()

NNe <- (1:dim(rthE)[2])[names(rthE)=="MANUII"] -1
NNp <- (1:dim(rthP)[2])[names(rthP)=="MANUII"] -1
NNt <- (1:dim(rthT)[2])[names(rthT)=="MANUII"] -1
## "MANUII" is the first non-taxon name column
## number of taxa

city.data <- function(City="BOS"){
    e.data <- rthE[rthE$City==City,]
    p.data <- rthP[rthP$City==City,]
    t.data <- rthT[rthT$City==City,]

    names(e.data)[1:NNe] <- paste("e", 1:NNe, sep="")
    names(p.data)[1:NNp] <- paste("p", 1:NNp, sep="")
    names(t.data)[1:NNt] <- paste("t", 1:NNt, sep="")

## keeping taxa with at least 2 occurrences (to avoid complete separation?)
    non0colE <- apply(e.data[,1:NNe], 2,
                      FUN=function(x) return(sum(x>0)>1))
    e.data <- cbind(e.data[,non0colE],
                    e.data[,(NNe+1):dim(e.data)[2]])

    non0colP <- apply(p.data[,1:NNp], 2,
                      FUN=function(x) return(sum(x>0)>1))
    p.data <- cbind(p.data[,non0colP],
                    p.data[,(NNp+1):dim(p.data)[2]])

    non0colT <- apply(t.data[,1:NNt], 2,
                      FUN=function(x) return(sum(x>0)>1))
    t.data <- cbind(t.data[,non0colT],
                    t.data[,(NNt+1):dim(t.data)[2]])
    Ne <- sum(substring(names(e.data),1,1)=="e")
    Np <- sum(substring(names(p.data),1,1)=="p")
    Nt <- sum(substring(names(t.data),1,1)=="t")

    return(list(e.data, p.data, t.data, c(Ne,Np,Nt)))
}

## three taxa example in Qian and Cuffney
temp <- city.data(City="BOS")
e.data <- temp[[1]]
p.data <- temp[[2]]
t.data <- temp[[3]]

plot(e15~NLCD2, data=e.data)
points(e.data$NLCD2, e.data$e25, col="red")
points(e.data$NLCD2, e.data$e29, col="blue")

## individual taxon

zip2 <- "
data{
  int<lower=1> n0;
  int<lower=1> np;
  int<lower=0> yp[np];
  vector[np] xp;
  vector[n0] x0;
}
parameters{
  real alpha0;
  real<lower=0> alpha1;
  real<upper=10> beta0;
  real<lower=0,upper=5> beta1;
  real<lower=-1,upper=5> beta2;
}
transformed parameters{
  vector[n0] theta0;
  vector[n0] lambda0;
  vector[np] theta1;
  vector[np] lambda1;
  theta0 = inv_logit(alpha0+alpha1*logit(x0/100));
  lambda0 = exp(beta0-beta1*x0 + beta2*log(x0));
  theta1 = inv_logit(alpha0+alpha1*logit(xp/100));
  lambda1 = exp(beta0-beta1*xp + beta2*log(xp));
}
model{
  for (i in 1:n0){
    target += log_sum_exp(log(theta0[i]),
               log1m(theta0[i])-lambda0[i]);
  }
  target += log1m(theta1);
  target += poisson_lpmf(yp|lambda1);
}
"
fit <- stan_model(model_code=zip2)

ZIP2_in <- function(y, x, n.chains=nchains){
     n0 <- sum(y==0)
     np <- sum(y!=0)
     yp <- y[y!=0]
     x0 <- x[y==0]
     xp <- x[y!=0]
   data <- list(n0=n0, np=np, yp=yp, x0=x0, xp=xp)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(alpha0=runif(1),alpha1=runif(1),
                       beta0=runif(1),beta1=runif(1),
                       beta2= 1+runif(1))
  paras <- c("alpha0","alpha1","beta0","beta1","beta2")
  return(list(data=data, init=inits,
              para=paras, nchains=n.chains))
}

## e15:
input.to.stan  <- ZIP2_in(y=e.data$e15, x=e.data$NLCD2)
fit2keep_e15 <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)
##pairs(fit2keep_e15)
## beta0, beta1, beta2 are highly correlated
## identification problem due to low number of non-zero observations

print(fit2keep_e15)

## e25:
input.to.stan  <- ZIP2_in(y=e.data$e25, x=e.data$NLCD2)
##fit <- stan_model(model_code=zip2)
fit2keep_e25 <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control = list(adapt_delta = 0.999))
##pairs(fit2keep_e25)
print(fit2keep_e25)
tikz(file=paste(plotDIRch5, "zipBOS_pairs.tex", sep="/"),
     height=3.5,width=3.5, standAlone=F)
pairs(fit2keep_e25, pars=c("beta0","beta1","beta2","alpha0","alpha1"))
dev.off()
## e29:
input.to.stan  <- ZIP2_in(y=e.data$e29, x=e.data$NLCD2)
##fit <- stan_model(model_code=zip2)
fit2keep_e29 <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains,
                     control = list(adapt_delta = 0.99))
##pairs(fit2keep_e29)
print(fit2keep_e29)
save(fit2keep_e15, fit2keep_e25, fit2keep_e29, file="euseINDe3sp.RData")

## extracting coefficients
coef_e15 <- summary(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e15,
                                      permuted=T)))))
coef_e25 <- summary(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e25,
                                                          permuted=T)))))
coef_e29 <- summary(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e29,
                                                          permuted=T)))))
save(coef_e15, coef_e25, coef_e29, file="EPT_BOS_single.RData")

load("EPT_BOS_single.RData")

tikz(file=paste(plotDIRch5, "etaxa3.tex", sep="/"), height=5.25, width=4,
     standAlone=F)
par(mfrow=c(3,2), mar=c(0,0,0,0), oma=c(3.5,3.5,2,3.5),
    mgp=c(1.25,0.125,0), tck=0.01)
plot(e.data$NLCD2, e.data$e15, xlab="", ylab="", axes=F)
curve(exp(coef_e15$mean[3]-coef_e15$mean[4]*x+coef_e15$mean[5]*log(x))*
    (1-invlogit(coef_e15$mean[1]+coef_e15$mean[2]*logit(x))), add=T)
curve(exp(coef_e15$mean[3]-coef_e15$mean[4]*x+coef_e15$mean[5]*log(x)),
      add=T, lty=2)
axis(2, outer=T)
axis(3, outer=T)
box()
legend("topright", "(a)", bty="n")

curve(invlogit(coef_e15$mean[1]+coef_e15$mean[2]*logit(x)),
      xlim=range(e.data$NLCD2), ylim=c(0,1),
      ylab="", xlab="", axes=F)
axis(3, outer=T)
axis(4, outer=T)
box()
points(e.data$NLCD2,jitter.binary(e.data$e15==0))

plot(e.data$NLCD2, e.data$e25, xlab="", ylab="", axes=F)
curve(exp(coef_e25$mean[3]-coef_e25$mean[4]*x+coef_e25$mean[5]*log(x))*
    (1-invlogit(coef_e25$mean[1]+coef_e25$mean[2]*logit(x))), add=T)
curve(exp(coef_e25$mean[3]-coef_e25$mean[4]*x+coef_e25$mean[5]*log(x)),
      lty=2, add=T)
axis(2, outer=T)
box()
legend("topright", "(b)", bty="n")

curve(invlogit(coef_e25$mean[1]+coef_e25$mean[2]*logit(x)),
      xlim=range(e.data$NLCD2), ylim=c(0,1),
      ylab="", xlab="", axes=F)
axis(4, outer=T)
box()
points(e.data$NLCD2,jitter.binary(e.data$e25==0))

plot(e.data$NLCD2, e.data$e29, xlab="", ylab="", axes=F)
curve(exp(coef_e29$mean[3]-coef_e29$mean[4]*x+coef_e29$mean[5]*log(x))*
      (1-invlogit(coef_e29$mean[1]+coef_e29$mean[2]*logit(x))), add=T)
curve(exp(coef_e29$mean[3]-coef_e29$mean[4]*x+coef_e29$mean[5]*log(x)),
      lty=2, add=T)
axis(1, outer=T)
axis(2, outer=T)
box()
legend("topright", "(c)", bty="n")

curve(invlogit(coef_e29$mean[1]+coef_e29$mean[2]*logit(x)),
      xlim=range(e.data$NLCD2), ylim=c(0,1),
      ylab="", xlab="", axes=F)
axis(1)
axis(4, outer=T)
box()
points(e.data$NLCD2,jitter.binary(e.data$e29==0))
mtext("\\% Developed Land", side=1, outer=T, line=1.75)
mtext("Observed Abundance", side=2, outer=T, line=1.75)
mtext("Prob of Extinction", side=4, outer=T, line=1.75)
dev.off()

## Hierarchical ZIP
## combining EPT data:

## E.P.T.
## number of taxa
## multilevel model -- E.P.T.
mult_data <- function(temp=city.data(City="DEN")){
    e.data <- temp[[1]]
    p.data <- temp[[2]]
    t.data <- temp[[3]]
    n.ept <- temp[[4]]

    if (n.ept[1]==0) {
        Edata.mult <- NULL
        print("No e species")
    } else
        Edata.mult <-
            data.frame(y=unlist(e.data[,1:n.ept[1]]),
                       NUII = rep(e.data$NUII, n.ept[1]),
                       taxa = rep(names(e.data)[1:n.ept[1]],
                                  each=dim(e.data)[1]),
                       NLCD2= rep(e.data$NLCD2, n.ept[1]))
    if (n.ept[2]==0){
        Pdata.mult <- NULL
        print("No p species")
    } else
        Pdata.mult <- data.frame(y=unlist(p.data[,1:n.ept[2]]),
                                 NUII = rep(p.data$NUII, n.ept[2]),
                                 taxa = rep(names(p.data)[1:n.ept[2]],
                                            each=dim(p.data)[1]),
                                 NLCD2= rep(p.data$NLCD2, n.ept[2]))
    if (n.ept[3]==0){
        Tdata.mult <- NULL
        print("No t species")
    } else
        Tdata.mult <- data.frame(y=unlist(t.data[,1:n.ept[3]]),
                                 NUII = rep(t.data$NUII, sum(n.ept[3])),
                                 taxa = rep(names(t.data)[1:n.ept[3]],
                                            each=dim(e.data)[1]),
                                 NLCD2= rep(t.data$NLCD2, sum(n.ept[3])))

    return(rbind(Edata.mult,Pdata.mult,Tdata.mult))
}

zip3 <- "
data{
  int<lower=1> nsp;
  int<lower=1> n0;
  int<lower=1> np;
  int<lower=0> yp[np];
  vector[np] xp;
  vector[n0] x0;
  int<lower=0> TX0[n0];
  int<lower=0> TXp[np];
}
parameters{
  vector[nsp] alpha0;
  vector<lower=0>[nsp] alpha1;
  real a0;
  real a1;
  vector[nsp] beta0;
  real b0;
  vector<lower=0>[nsp] beta1;
  real<lower=0> b1;
  vector<lower=-1>[nsp] beta2;
  real<lower=-1> b2;
  real<lower=0> sigma[5];
}
transformed parameters{
  vector[n0] theta0;
  vector[n0] lambda0;
  vector[np] theta1;
  vector[np] lambda1;
  for (i in 1:n0){
    theta0[i] = inv_logit(alpha0[TX0[i]]+alpha1[TX0[i]]*logit(x0[i]/100));
    lambda0[i] = exp(beta0[TX0[i]]-beta1[TX0[i]]*x0[i] +
                 beta2[TX0[i]]*log(x0[i]));
  }
  for (i in 1:np){
    theta1[i] = inv_logit(alpha0[TXp[i]]+alpha1[TXp[i]]*logit(xp[i]/100));
    lambda1[i] = exp(beta0[TXp[i]]-beta1[TXp[i]]*xp[i] +
                 beta2[TXp[i]]*log(xp[i]));
  }
}
model{
  beta0 ~ normal(b0, sigma[1]);
  beta1 ~ normal(b1, sigma[2]);
  beta2 ~ normal(b2, sigma[3]);
  alpha0 ~ normal(a0, sigma[4]);
  alpha1 ~ normal(a1, sigma[5]);
  for (i in 1:n0){
    target += log_sum_exp(log(theta0[i]),
               log1m(theta0[i])-lambda0[i]);
  }
  target += log1m(theta1);
  target += poisson_lpmf(yp|lambda1);
}
generated quantities{
  vector[nsp] logitEC50;
  for (i in 1:nsp)
    logitEC50[i] = -alpha0[i]/alpha1[i];
}
"
fit3 <- stan_model(model_code=zip3)
Stan.inMult <- function(y=EPT.data$y, x=EPT.data$NUII,
                        sp=EPT.data$taxa, n.chains=nchains,
                        adjust=T, a=0.025){
  n <- length(y)
  iszero <- y==0
  n0 <- sum(iszero)
  np <- n-n0
  if (adjust & (range(x)[1]==0 | range(x)[2]==100))
      x <- 50 + (1-2*a) * (x - 50)
  ## avoiding logit/log (0 or 100)
  y0 <- y[iszero]
  x0 <- x[iszero]
  yp <- y[!iszero]
  xp <- x[!iszero]
  sp_names <- levels(ordered(sp))
  sp <- as.numeric(ordered(sp))
  nsp <- max(sp)
  taxa0 <- sp[iszero]
  taxap <- sp[!iszero]
  bugs.ini <- list()
  bugs.data <- list(n0=n0, np=np, nsp=nsp, yp=yp, x0=x0, xp=xp, TXp=taxap,
                    TX0=taxa0)
  for (i in 1:n.chains)
    bugs.ini [[i]] <- list(alpha0=runif(nsp,0,0.1),alpha1=runif(nsp,0,0.1),
                           a0=runif(1,0,0.1),a1=runif(1,0,0.1),
                           beta0=runif(nsp,0,0.1),beta1=runif(nsp,0,0.1),
                           beta2=runif(nsp,0,0.1),b0=runif(1,0,0.1),
                           b1=runif(1,0,0.1),b2=runif(1,0,0.1),
                           sigma=runif(5, 1,2))
  paras <- c("alpha0","alpha1","beta0","beta1","beta2",
             "a0","a1","b0","b1","b2", "sigma", "logitEC50")
  return(list(data=bugs.data, init=bugs.ini, para=paras,
              nchains=n.chains, sp=sp_names))
}

## BOS

BOS_mult <- mult_data(temp=city.data(City="BOS"))
input.to.stan  <- Stan.inMult(y=BOS_mult$y, x=BOS_mult$NUII,
                     sp=BOS_mult$taxa, n.chains=nchains)
fit2keep_multBOS <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_mult)
save(fit2keep_multBOS, file="stanout_multBOS_zip.RData")
BOS_taxa <- input.to.stan$sp

load("stanout_multBOS_zip.RData")
pairs(fit2keep_multBOS, pars=c("a0","a1","b0","b1","b2"))
pairs(fit2keep_multBOS, pars=c("alpha0[2]","alpha1[2]",
                            "beta0[2]","beta1[2]","beta2[2]"))

## ATL
ATL_mult <- mult_data(temp=city.data(City="ATL"))
input.to.stan  <- Stan.inMult(y=ATL_mult$y, x=ATL_mult$NUII,
                     sp=ATL_mult$taxa, n.chains=nchains)
fit2keep_multATL <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
##pairs(fit2keep_e25)
print(fit2keep_multATL)
save(fit2keep_multATL, file="stanout_multATL_zip.RData")
ATL_taxa <- input.to.stan$sp
load("stanout_multATL_zip.RData")
## BIR
BIR_mult <- mult_data(temp=city.data(City="BIR"))
input.to.stan  <- Stan.inMult(y=BIR_mult$y, x=BIR_mult$NUII,
                     sp=BIR_mult$taxa, n.chains=nchains)
fit2keep_multBIR <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
##pairs(fit2keep_e25)
print(fit2keep_multBIR)
save(fit2keep_multBIR, file="stanout_multBIR_zip.RData")
BIR_taxa <- input.to.stan$sp
load("stanout_multBIR_zip.RData")

## DEN
DEN_mult <- mult_data(temp=city.data(City="DEN"))
input.to.stan  <- Stan.inMult(y=DEN_mult$y, x=DEN_mult$NUII,
                     sp=DEN_mult$taxa, n.chains=nchains)
fit2keep_multDEN <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_multDEN)
save(fit2keep_multDEN, file="stanout_multDEN_zip.RData")
DEN_taxa <- input.to.stan$sp
load("stanout_multDEN_zip.RData")

## DFW
DFW_mult <- mult_data(temp=city.data(City="DFW"))
input.to.stan  <- Stan.inMult(y=DFW_mult$y, x=DFW_mult$NUII,
                     sp=DFW_mult$taxa, n.chains=nchains)
fit2keep_multDFW <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.99,
                                            max_treedepth = 20))
print(fit2keep_multDFW)
save(fit2keep_multDFW, file="stanout_multDFW_zip.RData")
DFW_taxa <- input.to.stan$sp
load("stanout_multDFW_zip.RData")

## MGB
MGB_mult <- mult_data(temp=city.data(City="MGB"))
input.to.stan  <- Stan.inMult(y=MGB_mult$y, x=MGB_mult$NUII,
                     sp=MGB_mult$taxa, n.chains=nchains)
fit2keep_multMGB <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_multMGB)
save(fit2keep_multMGB, file="stanout_multMGB_zip.RData")
MGB_taxa <- input.to.stan$sp
load("stanout_multMGB_zip.RData")

## POR
POR_mult <- mult_data(temp=city.data(City="POR"))
input.to.stan  <- Stan.inMult(y=POR_mult$y, x=POR_mult$NUII,
                     sp=POR_mult$taxa, n.chains=nchains)
fit2keep_multPOR <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_multPOR)
save(fit2keep_multPOR, file="stanout_multPOR_zip.RData")
POR_taxa <- input.to.stan$sp
load("stanout_multPOR_zip.RData")

## RAL
RAL_mult <- mult_data(temp=city.data(City="RAL"))
input.to.stan  <- Stan.inMult(y=RAL_mult$y, x=RAL_mult$NUII,
                     sp=RAL_mult$taxa, n.chains=nchains)
fit2keep_multRAL <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_multRAL)
save(fit2keep_multRAL, file="stanout_multRAL_zip.RData")
RAL_taxa <- input.to.stan$sp
load("stanout_multRAL_zip.RData")

## SLC
SLC_mult <- mult_data(temp=city.data(City="SLC"))
input.to.stan  <- Stan.inMult(y=SLC_mult$y, x=SLC_mult$NUII,
                     sp=SLC_mult$taxa, n.chains=nchains)
fit2keep_multSLC <- sampling(fit3, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains,
                             control = list(adapt_delta = 0.95))
print(fit2keep_multSLC)
save(fit2keep_multSLC, file="stanout_multSLC_zip.RData")
SLC_taxa <- input.to.stan$sp
load("stanout_multSLC_zip.RData")

pairs(fit2keep_multSLC, pars=c("a0","a1","b0","b1","b2"))

### Processing output
##### 1. Figures of ST50 comparing individual species model to multilevel model
##### 2. Hyper-parameter ST50 for 9 regions
###

### figure: comparing model coefficients of e15, e25, and e29
elist <- BOS_taxa %in% c("e15","e25","e29")
enumber <- (1:length(elist))[elist]
BOS_alpha0 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars=paste("alpha0[", 1:length(elist), "]", sep="")))))
BOS_alpha1 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars=paste("alpha1[", 1:length(elist), "]", sep="")))))
BOS_beta0 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars=paste("beta0[", 1:length(elist), "]", sep="")))))
BOS_beta1 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars=paste("beta1[", 1:length(elist), "]", sep="")))))
BOS_beta2 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars=paste("beta2[", 1:length(elist), "]", sep="")))))

st50_e15M <- summary(invlogit(-BOS_alpha0[enumber[1]]/BOS_alpha1[enumber[1]]))
st50_e25M <- summary(invlogit(-BOS_alpha0[enumber[2]]/BOS_alpha1[enumber[2]]))
st50_e29M <- summary(invlogit(-BOS_alpha0[enumber[3]]/BOS_alpha1[enumber[3]]))

coef_e15rv <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e15,
                                      permuted=T))))
coef_e25rv <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e25,
                                                          permuted=T))))
coef_e29rv <- rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep_e29,
                                                          permuted=T))))
st50_e15 <- summary(invlogit(-coef_e15rv[1]/coef_e15rv[2]))
st50_e25 <- summary(invlogit(-coef_e25rv[1]/coef_e25rv[2]))
st50_e29 <- summary(invlogit(-coef_e29rv[1]/coef_e29rv[2]))

BOS_a0 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars="a0"))))
BOS_a1 <- rvsims(as.matrix(as.data.frame(extract(fit2keep_multBOS, pars="a1"))))

st50_BOShyp <- summary(invlogit(-BOS_a0/BOS_a1))

tikz(file=paste(plotDIRch5, "st50shrink.tex", sep="/"), height=3, width=4,
     standAlone=F)  ## zipShrink
par(mar=c(3,7,1,0.5), mgp=c(1.25,0.125,0),las=1,tck=0.01)
plot(c(0,1), c(1,4), type="n", xlab="ST50", ylab="", axes=F)
points(x=c(st50_BOShyp$mean, st50_e15M$mean,  st50_e25M$mean,  st50_e29M$mean),y=c(1, (2:4)-0.25/4))
segments(x0=c(st50_BOShyp$"2.5%", st50_e15M$"2.5%",
              st50_e25M$"2.5%",  st50_e29M$"2.5%"),
         x1=c(st50_BOShyp$"97.5%", st50_e15M$"97.5%",
              st50_e25M$"97.5%",  st50_e29M$"97.5%"),
         y0=c(1, (2:4)-0.125/2),y1=c(1, (2:4)-0.25/4))
segments(x0=c(st50_BOShyp$"25%", st50_e15M$"25%",
              st50_e25M$"25%",  st50_e29M$"25%"),
         x1=c(st50_BOShyp$"75%", st50_e15M$"75%",
              st50_e25M$"75%",  st50_e29M$"75%"),
         y0=c(1, (2:4)-0.125/2),y1=c(1, (2:4)-0.25/4), lwd=3)
points(x=c(st50_e15$mean, st50_e25$mean, st50_e29$mean),y=(2:4)+0.125/2)
segments(x0=c(st50_e15$"2.5%",
              st50_e25$"2.5%",  st50_e29$"2.5%"),
         x1=c(st50_e15$"97.5%",
              st50_e25$"97.5%",  st50_e29$"97.5%"),
         y0=(2:4)+0.125/2,y1=(2:4)+0.125/2, lty=2)
segments(x0=c(st50_e15$"25%",
              st50_e25$"25%",  st50_e29$"25%"),
         x1=c(st50_e15$"75%",
              st50_e25$"75%",  st50_e29$"75%"),
         y0=(2:4)+0.125/2,y1=(2:4)+0.125/2, lwd=3, lty=2)
abline(v=st50_BOShyp$mean)
axis(1)
axis(2, at=1:4, labels=c("Hyper", "\\emph{S. serrata}","\\emph{A. turbida}",
                         "\\emph{B. flavistriga}"))
dev.off()

## figure 4 of Qian and Cuffney 2014, only the right column
plotFUN <- function(fit=fit2keep_multATL, mult_data=ATL_mult){
    input.to.stan  <- Stan.inMult(y=mult_data$y, x=mult_data$NUII,
                                  sp=mult_data$taxa, n.chains=nchains)
    FITout <- extract(fit, permuted=T)
    FITst50_hyp <- invlogit(-rvsims(FITout$a0)/rvsims(FITout$a1))
    FITout_st50 <- invlogit(rvsims(FITout$logitEC50))
    oo <- order(summary(FITout_st50)$mean)
    return(list(st50=FITout_st50, st50_hyp=FITst50_hyp, oo=oo,
                taxa=input.to.stan$sp))
}

tmp <- plotFUN(fit=fit2keep_multATL, mult_data=ATL_mult)
ATLst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multBIR, mult_data=BIR_mult)
BIRst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multBOS, mult_data=BOS_mult)
BOSst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multPOR, mult_data=POR_mult)
PORst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multRAL, mult_data=RAL_mult)
RALst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multSLC, mult_data=SLC_mult)
SLCst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multDEN, mult_data=DEN_mult)
DENst50 <- tmp$st50_hyp

tmp <- plotFUN(fit=fit2keep_multDFW, mult_data=DFW_mult)
DFWst50 <- tmp$st50_hyp
tmp <- plotFUN(fit=fit2keep_multMGB, mult_data=MGB_mult)
MGBst50 <- tmp$st50_hyp

ttt <- summary(c(ATLst50, BIRst50, BOSst50, PORst50,
         RALst50, SLCst50, DENst50, DFWst50, MGBst50))

tikz(file=paste(plotDIRch5, "eusest50s.tex", sep="/"), height=3, width=4,
     standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(c(0,0.65), c(1,9), type="n", axes=F,
     xlab="Regional ST50", ylab="")
points(ttt$"50%", 1:9)
segments(x0=ttt$"2.5%", x1=ttt$"97.5%", y0=1:9, y1=1:9)
segments(x0=ttt$"25%", x1=ttt$"75%", y0=1:9, y1=1:9, lwd=3)
axis(2, at=1:9, labels=c("ATL","BIR","BOS","POR","RAL","SLC",
                         "DEN","DFW","MGB"), las=1)
axis(1)
dev.off()



### Stop here.
## using multivariate normal prior
### 1.
#zip3mvn <- "
#data{
#  int<lower=1> nsp;
#  int<lower=1> n0;
#  int<lower=1> np;
#  int<lower=0> yp[np];
#  vector[np] xp;
#  vector[n0] x0;
#  int<lower=0> TX0[n0];
#  int<lower=0> TXp[np];
#}
#parameters{
#  vector[nsp] alpha0;
#  vector<lower=0>[nsp] alpha1;
#  real a0;
#  real a1;
#  vector[nsp] beta0;
#  real b0;
#  vector<lower=0>[nsp] beta1;
#  real<lower=0> b1;
#  vector<lower=-1>[nsp] beta2;
#  real<lower=-1> b2;
#  vector<lower=0>[3] tau_b;
#  vector<lower=0>[2] tau_a;
#  corr_matrix[3] Omega_b;
#  corr_matrix[2] Omega_a;
#}
#transformed parameters{
#  vector[2] alpha[nsp];
#  vector[3] beta[nsp];
#  row_vector[3] mu_b;
#  row_vector[2] mu_a;
#  vector[n0] theta0;
#  vector[n0] lambda0;
#  vector[np] theta1;
#  vector[np] lambda1;
#  for (i in 1:nsp){
#    alpha[i][1]=alpha0[i];
#    alpha[i][2]=alpha1[i];
#    beta[i][1]=beta0[i];
#    beta[i][2]=beta1[i];
#    beta[i][3]=beta2[i];
#  }
#  mu_a[1]=a0;
#  mu_a[2]=a1;
#  mu_b[1]=b0;
#  mu_b[2]=b1;
#  mu_b[3]=b2;
#  for (i in 1:n0){
#    theta0[i] = inv_logit(alpha[TX0[i]][1]+
#                          alpha[TX0[i]][2]*logit(x0[i]));
#    lambda0[i] = exp(beta[TX0[i]][1]-beta[TX0[i]][2]*x0[i] +
#                 beta[TX0[i]][3]*log(x0[i]));
#  }
#  for (i in 1:np){
#    theta1[i] = inv_logit(alpha[TXp[i]][1]+
#                          alpha[TXp[i]][2]*logit(xp[i]/100));
#    lambda1[i] = exp(beta[TXp[i]][1]-beta[TXp[i]][2]*xp[i] +
#                     beta[TXp[i]][3]*log(xp[i]));
#  }
#}
#model{
#  tau_a ~ normal(0,2.5);
#  tau_b ~ normal(0,2.5);
#  Omega_a ~ lkj_corr(2);
#  Omega_b ~ lkj_corr(2);
#  alpha ~ multi_normal(mu_a, quad_form_diag(Omega_a, tau_a));
#  beta ~ multi_normal(mu_b, quad_form_diag(Omega_b, tau_b));
#  for (i in 1:n0){
#    target += log_sum_exp(log(theta0[i]),
#               log1m(theta0[i])-lambda0[i]);
#  }
#  target += log1m(theta1);
#  target += poisson_lpmf(yp|lambda1);
#}
#"
#fit3mvn <- stan_model(model_code=zip3mvn)


