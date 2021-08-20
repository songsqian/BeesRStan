source("FrontMatter.R")

## Chapter 5 -- multinomial regression

plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")

packages(rv)
packages(rstan)
packages(MASS)
packages(nnet)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

##data
euseSPdata <- read.csv(paste(dataDIR, "usgsmultinomial.csv", sep="/"))
envpara <- read.csv(paste(dataDIR, "envpara.csv", sep="/"), header=T)
names(euseSPdata)[30:33] <- c("Sensitive","Mod Tol","Tol", "Unknown")

# fit a multinomial logit model where PID is the response variable,
# and the predictor is nuii:

multinom.BOS <- multinom(as.matrix(euseSPdata[,30:33])~NLCD2,data=euseSPdata,
                         subset=City=="BOS")
summary(multinom.BOS, corr=FALSE)

beta <- coef(multinom.BOS)
X <- cbind(1,0:100)

n.taxa <- 4
Xb <- matrix(0, nrow=dim(X)[1], ncol=n.taxa-1)
for (i in 1:(n.taxa-1)) Xb[,i] <- X %*% beta[i,]

denomsum <- apply(exp(Xb), 1, sum)

pp <-matrix(0, nrow=dim(X)[1], ncol=n.taxa)
pp[,1] <- 1/(1+denomsum)
for (i in 2:n.taxa) pp[,i] <- exp(Xb[,i-1])/(1+denomsum)

## alternative:
Xb1 <- X %*% beta[1,]
Xb2 <- X %*% beta[2,]
Xb3 <- X %*% beta[3,]

denomsum <- exp(Xb1) + exp(Xb2) + exp(Xb3)

p0 <- 1/(1+denomsum)
p1 <- exp(Xb1)/(1+denomsum)
p2 <- exp(Xb2)/(1+denomsum)
p3 <- exp(Xb3)/(1+denomsum)

dataP <- t(apply(euseSPdata[,30:33], 1, function(x) return(x/sum(x))))
## figure 1
rows <- euseSPdata$City=="BOS"

tikz(file=paste(plotDIR5, "multN_BOS1.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
plot(0:100, p0, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
lines(0:100, p0)
text(50, 0.9, "Sensitive", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(0:100, p1, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
lines(0:100, p1)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)
plot(0:100, p2, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
lines(0:100, p2)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)
plot(0:100, p3, xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
lines(0:100, p3)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "relative abundance", outer=T, line=1.75, las=0)
dev.off()

## Using Stan for uncertainty analysis

multN <- "
data {
  int<lower = 2> K;
  int<lower = 0> N;
  int<lower = 1> D;
  int y[N, K];
  matrix[N, D] x;
}

transformed data {
  row_vector[D] zeros = rep_row_vector(0, D);
}

parameters {
  matrix[K - 1, D] beta_raw;
}

transformed parameters {
  matrix[K,D] beta;
  beta = append_row(beta_raw, zeros);
}
model {
  matrix[N, K] x_beta = x * beta';
  to_vector(beta) ~ normal(0, 5);
  for (n in 1:N)
    y[n, ] ~ multinomial(softmax(to_vector(x_beta[n,])));
}
"

fit <- stan_model(model_code=multN)

multN_in <- function(data=euseSPdata, subset="BOS",
                     ycol=30:33, xcol="NLCD2", chains=nchains)
{
    data <- data[data$City==subset, ]
    y <- as.matrix(data[,ycol])
    x <- cbind(1, data[,xcol])
    N <- dim(data)[1]
    K <- dim(y)[2]
    D <- dim(x)[2]
    data_in <- list(N=N,K=K,D=D,y=y,x=x)
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(beta_raw=matrix(rnorm((K-1)*D), nrow=K-1,ncol=D))
    paras <- c("beta_raw")
    return(list(data=data_in, init=inits, para=paras, nchains=chains))
}

input.to.stan <- multN_in()

fit2keep <- sampling(fit, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)
pairs(fit2keep, pars=c("beta_raw"))
print(fit2keep)
traceplot(fit2keep)

mulN_stan <-
    rvmatrix(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep,
                                                           pars="beta_raw",
                                                           permuted=T)))),
             nrow=3)

X <- cbind(1,1:100)
n.taxa <- 4
Xb_rv <- '%*%.rv' (X , t(mulN_stan))

temp <- 1+exp(Xb_rv)[,1]+exp(Xb_rv)[,2]+exp(Xb_rv)[,3]

pp <-rvmatrix(0, nrow=dim(X)[1], ncol=n.taxa)
pp[,4] <- 1/(temp)
for (i in 1:(n.taxa-1)) pp[,i] <- exp(Xb_rv[,i])/(temp)

tikz(file=paste(plotDIR5, "multN_BOS2.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
plot(pp[,1], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
text(50, 0.9, "Sensitive", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(pp[,2], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)

plot(pp[,3], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)

plot(pp[,4], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "relative abundance", outer=T, line=1.75, las=0)
dev.off()

## mlplot(pp)

## Poisson-Multinomial transformation
PoimultN <-"
data {
  int<lower = 2> K;
  int<lower = 0> N;
  int<lower = 1> D;
  int<lower=1> Np;
  int y[N,K];
  matrix[N,D] x;
  matrix[Np,D] x_pred;
}
parameters {
  matrix[K,D] beta;
//  vector[N] lambda;
}
transformed parameters{
  matrix[N,K] x_beta;
  matrix[Np,K] mu_pred;
  x_beta = x*beta';
//  for (i in 1:N)
//    mu[i,] = lambda[i]+x_beta[i,];
  mu_pred = x_pred*beta';
}
model {
  for (k in 1:K){
    y[,k] ~ poisson(exp(x_beta[,k]));
  }
}
generated quantities{
  vector[K] p[Np];
  for (n in 1:Np){
    p[n] = softmax(to_vector(mu_pred[n,]));
  }
}
"

fitPoiMn <- stan_model(model_code=PoimultN)

PoimultN_in <- function(data=euseSPdata, subset="BOS",
                        ycol=30:33, xcol="NLCD2",
                        x.pred=1:100, chains=nchains)
{
    data <- data[data$City==subset, ]
    y <- as.matrix(data[,ycol])
    x <- cbind(1, data[,xcol], log(data[,xcol]))
    N <- dim(data)[1]
    K <- dim(y)[2]
    D <- dim(x)[2]
    data_in <- list(N=N,K=K,D=D,Np=length(x.pred),y=y,x=x,
                    x_pred=cbind(1,x.pred,log(x.pred)))
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(beta=matrix(rnorm(K*D), nrow=K,ncol=D))
    paras <- c("beta","p")
    return(list(data=data_in, init=inits, para=paras, nchains=chains))
}

PoimultN_in2 <- function(data=euseSPdata, subset="BOS",
                        ycol=30:33, xcol="NLCD2",
                        x.pred=1:100, chains=nchains)
{
    data <- data[data$City==subset, ]
    y <- as.matrix(data[,ycol])
    x <- cbind(1, data[,xcol])
    N <- dim(data)[1]
    K <- dim(y)[2]
    D <- dim(x)[2]
    data_in <- list(N=N,K=K,D=D,Np=length(x.pred),y=y,x=x,
                    x_pred=cbind(1,x.pred))
    inits <- list()
    for (i in 1:chains)
        inits[[i]] <- list(beta=matrix(rnorm(K*D), nrow=K,ncol=D))
    paras <- c("beta","p")
    return(list(data=data_in, init=inits, para=paras, nchains=chains))
}

input.to.stan <- PoimultN_in()

fit2keep <- sampling(fitPoiMn, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)
pairs(fit2keep, pars="beta")
print(fit2keep)
traceplot(fit2keep, pars="beta")

input.to.stan <- PoimultN_in2()

fit2keep2 <- sampling(fitPoiMn, data=input.to.stan$data,
                     init=input.to.stan$init,
                     pars=input.to.stan$para,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$nchains)

print(fit2keep2)

PoimulN_beta <-
    rvmatrix(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep2,
                                                           pars="beta",
                                                           permuted=T)))),
             nrow=4)

PoimulN_prob <-
    rvmatrix(rvsims(as.matrix(as.data.frame(rstan::extract(fit2keep2,
                                                           pars="p",
                                                           permuted=T)))),
             nrow=100)

tikz(file=paste(plotDIR5, "multN_BOS3.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3), mgp=c(1.25,0.25,0),
    tck=-0.015, las=1)
plot(PoimulN_prob[,1], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,1], col="gray", cex=0.5)
text(50, 0.9, "Sensitive", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(PoimulN_prob[,2], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,2], col="gray", cex=0.5)
text(50, 0.9, "Moderate Tol", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)

plot(PoimulN_prob[,3], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,3], col="gray", cex=0.5)
text(50, 0.9, "Tolerant", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)

plot(PoimulN_prob[,4], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(euseSPdata$NLCD2[rows], dataP[rows,4], col="gray", cex=0.5)
text(50, 0.9, "Unknown", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "\\% developed land", outer=T, line=1.5)
mtext(side=2, "relative abundance", outer=T, line=1.75, las=0)
dev.off()

