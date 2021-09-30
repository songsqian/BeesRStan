source("FrontMatter.R")

## Chapter 5
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")
## simulation
packages(rv)
packages(rstan)
packages(MASS)
packages(nnet)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

plotDIRch5 <- paste(plotDIR, "chapter5","figures", sep="/")
## Data -- the insect oviposition behavior (Ecological Detective)
oviposition <- data.frame(complement=4:23,
                          one=c(0,0,1,1,0,0,0,0,0,1,rep(0,10)),
                          two=c(2,5,11,5,2,1,4,3,1,2,0,3,2,2,0,0,0,1,0,0),
                          three=c(1,1,3,1,1,0,3,4,6,4,3,4,6,4,2,6,2,1,0,1),
                          four=c(rep(0,5),1,rep(0,4),1,rep(0,9))
                          )
oviposition <- oviposition[-19,]

## multinomial model (glm)
ovi_M1 <- multinom(as.matrix(oviposition[,-1]) ~ complement,
                      data=oviposition)
summary(ovi_M1, corr=F)
beta <- coef(ovi_M1)
X  <- cbind(1, 4:23)
n.grps <- 4
Xb <- matrix(0, nrow=dim(X)[1], ncol=n.grps-1)
for (i in 1:(n.grps-1)) Xb[,i] <- X%*% beta[i,]

denomsum <- apply(exp(Xb), 1, sum)
PP <- matrix(0, nrow=dim(X)[1], ncol=n.grps)
PP[,1] <- 1/(1+denomsum)
for (i in 2:n.grps) PP[,i] <- exp(Xb[,i-1])/(1+denomsum)


dataP <- t(apply(oviposition[,-1], 1,
                 function(x) return(x/sum(x))))

tikz(file=paste(plotDIRch5, "ovi_glm.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3),
    mgp=c(1.25,0.25,0), tck=-0.015, las=1)
plot(4:23, PP[,1], xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(oviposition$complement, dataP[,1], cex=0.5)
lines(4:23, PP[,1])
text(10, 0.9, "One", cex=0.75)
axis(3, outer=T)
axis(2,labels=FALSE)
box()
plot(4:23, PP[,2], xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(oviposition$complement, dataP[,2], cex=0.5)
lines(4:23,PP[,2])
text(10, 0.9, "Two", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)
plot(4:23, PP[,3], xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(oviposition$complement, dataP[,3], cex=0.5)
lines(4:23, PP[,3])
text(10, 0.9, "Three", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)
plot(4:23, PP[,4], xlab="", ylab="", ylim=c(0,1),
     type="n", col=1, pch=16, axes=F)
points(oviposition$complement, dataP[,4], cex=0.5)
lines(4:23, PP[,4])
text(10, 0.9, "Four", cex=0.75)
box()
axis(1, outer=T)
axis(4, labels=FALSE)
mtext(side=1, "Complement", outer=T, line=1.5)
mtext(side=2, "Pr(clutch size)", outer=T, line=1.75, las=0)
dev.off()

## using Stan

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

multN_in <- function(data=oviposition,ycol=2:5, xcol="complement",
                     chains=nchains)
{
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
    rvmatrix(rvsims(as.matrix(as.data.frame(
        rstan::extract(fit2keep, pars="beta_raw",
                       permuted=T)))), nrow=3)

X <- cbind(1,4:23)
n.grps <- 4
Xb_rv <- '%*%.rv' (X , t(mulN_stan))

temp <- 1+exp(Xb_rv)[,1]+exp(Xb_rv)[,2]+exp(Xb_rv)[,3]

pp <-rvmatrix(0, nrow=dim(X)[1], ncol=n.grps)
pp[,4] <- 1/(temp)
for (i in 1:(n.grps-1)) pp[,i] <- exp(Xb_rv[,i])/(temp)

tikz(file=paste(plotDIRch5, "ovi_stan.tex", sep="/"),
           height=4, width=4, standAlone=F)
par(mfrow=c(2,2), mar=c(0.,0.,0.,0.),oma=c(3,3,3,3),
    mgp=c(1.25,0.25,0), tck=-0.015, las=1)
plot(pp[,1], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(oviposition$complement-3, dataP[,1], pch=16, cex=0.5)
text(12.5-3, 0.9, "One", cex=0.75)
axis(3, at=seq(5,20,5)-3,labels=seq(5,20,5),outer=T)
axis(2,labels=FALSE)
box()
plot(pp[,2], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(oviposition$complement-3, dataP[,2], pch=16, cex=0.5)
text(12.5-3, 0.9, "Two", cex=0.75)
box()
axis(4, outer=T)
axis(3, labels=FALSE)

plot(pp[,3], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(oviposition$complement-3, dataP[,3], pch=16, cex=0.5)
text(12.5-3, 0.9, "Three", cex=0.75)
box()
axis(2, outer=T)
axis(1, labels=FALSE)

plot(pp[,4], xlab="", ylab="", ylim=c(0,1),
     axes=F)
points(oviposition$complement-3, dataP[,4], pch=16, cex=0.5)
text(12.5-3, 0.9, "Four", cex=0.75)
box()
axis(1, at=seq(5,20,5)-3,labels=seq(5,20,5), outer=T)
axis(4, labels=FALSE)
mtext(side=1, "Complement", outer=T, line=1.5)
mtext(side=2, "Pr(clutch size)", outer=T, line=1.75, las=0)
dev.off()

