source("FrontMatter.R")

## Chapter 5
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")
## simulation
packages(rv)
packages(rstan)
packages(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

###################
## Mixture model ##
###################

## Beta method of moment
beta_mm <- function(x){
    mux <- mean(x)
    sigmax <- sd(x)
    alpha <- (mux*(mux*(1-mux)-sigmax^2))/sigmax^2
    beta <- alpha*(1-mux)/mux
    return(c(alpha, beta))
}

## Binomial -- imperfect tests
## Gibbs sampler
latent.gibbs <- function(Y1, Y2, alpha=2, beta=2, a1=2, b1=2,
                         a2=2, b2=2, nsims=50000, warmup=T, thin=10){
    ## Data: Y1 (positive), Y2 (negative), N=Y1+Y2
    ## Latent variable: Z1 (present & positive), Z2 (present & negative)
    N <- Y1+Y2
    Z1 <- ceiling(runif(1)*Y1)
    Z2 <- ceiling(runif(1)*Y2)
    nwarmup <- round(nsims/2)
    if (warmup){
        for (i in 1:nwarmup){
          theta <- rbeta(1, alpha+Z1+Z2, beta+N-Z1-Z2)
          fp <- rbeta(1, a1+Y1-Z1, b1+Y2-Z2)
          fn <- rbeta(1, a2+Z2, b2+Z1)
          phi1 <- theta*(1-fn)/(theta*(1-fn)+(1-theta)*fp)
          Z1 <- rbinom(1, Y1, phi1)
          phi2 <- theta*fn/(theta*fn+(1-theta)*(1-fp))
          Z2 <- rbinom(1, Y2, phi2)
        }
        nn <- nwarmup+1
      } else {
        theta <- rbeta(1, alpha+Z1+Z2, beta+N-Z1-Z2)
        fp <- rbeta(1, a1+Y1-Z1, b1+Y2-Z2)
        fn <- rbeta(1, a2+Z2, b2+Z1)
        phi1 <- theta*(1-fn)/(theta*(1-fn)+(1-theta)*fp)
        Z1 <- rbinom(1, Y1, phi1)
        phi2 <- theta*fn/(theta*fn+(1-theta)*(1-fp))
        Z2 <- rbinom(1, Y2, phi2)
        nn <- 2
      }
    th <- pfp <- pfn <- pZ1 <- pZ2 <- numeric()
    i <- 0
    for (j in nn:nsims){
        theta <- rbeta(1, alpha+Z1+Z2, beta+N-Z1-Z2)
        fp <- rbeta(1, a1+Y1-Z1, b1+Y2-Z2)
        fn <- rbeta(1, a2+Z2, b2+Z1)
        phi1 <- theta*(1-fn)/(theta*(1-fn)+(1-theta)*fp)
        Z1 <- rbinom(1, Y1, phi1)
        phi2 <- theta*fn/(theta*fn+(1-theta)*(1-fp))
        Z2 <- rbinom(1, Y2, phi2)
        if (j%%thin==0){
            i <- i+1
            th[i] <- theta; pfp[i] <- fp; pfn[i] <- fn
            pZ1[i] <- Z1; pZ2[i] <- Z2
        }
    }
    return(data.frame(theta=th, fp=pfp, fn=pfn, Z1=pZ1, Z2=pZ2))
}

try1 <- latent.gibbs(Y1=6, Y2=14, alpha=4, beta=10,
                     a1=1, b1=10, a2=2, b2=10,
                     nsims=500000, warmup=T, thin=100)

try2 <- latent.gibbs(Y1=6, Y2=14, alpha=5, beta=10,
                     a1=1, b1=10, a2=2, b2=10,
                     nsims=1000, warmup=F, thin=1)
par(mfrow=c(1,3))
plot(try1$theta, type="l")
plot(try1$fp, type="l")
plot(try1$fn, type="l")

hist(try1$theta)
hist(try1$fp)
hist(try1$fn)


## random number generator
rbin_imperfect <- function(n=10, N=100, theta=0.2, fp=0.15, fn=0.4){
## n: number of random samples
## N: total number tested
    InF <- rbinom(n=n, size=N, prob=theta)
    InP <- rbinom(n, size=InF, prob=1-fn)
    NinF <- N-InF
    NiP <- rbinom(n, size=NinF, fp)
    y1 <- InP+NiP
    return(data.frame(Y1=y1, Y2=N-y1))
}

## simulation
sim_data50 <- rbin_imperfect(n=50, theta=0.2, fp=0.1, fn=0.3,N=100)
sim_sum <- NULL
for (i in 1:dim(sim_data50)[1]){
    if (i%%25==0) print(paste(i, "of ", dim(sim_data50)[1]))
    sim_sum <- rbind(sim_sum,
                     apply(latent.gibbs(Y1=sim_data50[i,1],
                                        Y2=sim_data50[i,2],
                                        alpha=4, beta=10, a1=6, b1=20,
                                        a2=1, b2=9, nsims=100000,
                                        warmup=T, thin=50),
                           2, mean))
}

## simulation 2: sequential updating
sim_sum2 <- NULL
al=2; be=2; a1=3; b1=23; a2=2; b2=22;
nn <- dim(sim_data50)[1]
for (i in 1:nn){
    if (i%%round(nn/10)==0) print(paste(i, "of ", dim(sim_data50)[1]))
    y1 <- sim_data50$Y1[i]
    y2 <- sim_data50$Y2[i]
    temp <- latent.gibbs(Y1=y1, Y2=y2, alpha=al, beta=be,
                         a1=a1, b1=b1, a2=a2, b2=b2,
                         nsims=100000, warmup=T, thin=500)
   ## theta <- try(fitdistr(temp$theta, densfun="beta",
   ##                       start=list(shape1=al, shape2=be))[[1]],
   ##              silent=T)
   ## if (!is.null(attributes(theta)$class))
        theta <- beta_mm(temp$theta)
   ## fp <- try(fitdistr(temp$fp, densfun="beta",
   ##                   start=list(shape1=a1, shape2=b1))[[1]], silent=T)
   ## if (!is.null(attributes(fp)$class))
        fp <- beta_mm(temp$fp)
   ## fn <- try(fitdistr(temp$fn, densfun="beta",
   ##                   start=list(shape1=a2, shape2=b2))[[1]], silent=T)
   ## if (!is.null(attributes(fn)$class))
        fn <- beta_mm(temp$fn)

    sim_sum2 <- rbind(sim_sum2, rvsims(as.matrix(temp)))
    al <- theta[1]
    be <- theta[2]
    ##if (al+be>100){
    ##    al <- al/10
    ##    be <- be/10
    ##}
    a1 <- fp[1]
    b1 <- fp[2]
    ##if (a1+b1>100){
    ##    a1 <- a1/10
    ##    b1 <- b1/10
    ##}
    a2 <- fn[1]
    b2 <- fn[2]
    ##if (a2+b2>100){
    ##    a2 <- a2/10
    ##    b2 <- b2/10
   ## }
}

## prior
al=2; be=2; a1=3; b1=23; a2=2; b2=22;

tikz(file=paste(plotDIRch5, "simbinlatent.tex", sep="/"),
     width=4.75, height=2., standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
hist(sim_sum[,1], prob=T, xlab="$\\theta$", main="",
     xlim=range(c(sim_sum[,1], 0.1)))
abline(v=0.2, col="red", lwd=2) ## true rate
curve(dbeta(x,al,be), add=T, col=grey(0.4)) ## prior density
hist(sim_sum[,2], prob=T, xlab="$f_p$", main="",
     xlim=range(c(sim_sum[,2], 0.2)))
abline(v=0.1,col="red", lwd=2)
curve(dbeta(x, a1, b1), add=T, col=grey(0.4))
hist(sim_sum[,3], prob=T, xlab="$f_n$", main="",
     xlim=range(c(sim_sum[,3], 0.3)))
abline(v=0.3,col="red", lwd=2)
curve(dbeta(x, a2,b2), add=T, col=grey(0.4))
dev.off()

tikz(file=paste(plotDIRch5, "simbinupdate.tex", sep="/"),
     width=4.75, height=2., standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
mlplot(sim_sum2[,1], xlab="$\\theta$", main="")
abline(v=0.1, col="red", lwd=2) ## true rate
mlplot(sim_sum2[,2], xlab="$f_p$", main="")
abline(v=0.2,col="red", lwd=2)
mlplot(sim_sum2[,3], xlab="$f_n$", main="")
abline(v=0.3,col="red", lwd=2)
dev.off()

## Ohio COVID-19 test data
OhioCovid <- read.csv(paste(dataDIR, "Ohio_COVID_daily.csv", sep="/"))
OhioCovid$RDate <- as.Date(as.character(OhioCovid$date), format="%Y%m%d")

OhioCovid <- OhioCovid[, c("date", "RDate", "positive","totalTestResults")]
OhioCovid$positiveIncrease <-  c(-diff(OhioCovid$positive, 1), NA)
OhioCovid$negativeIncrease <- c(-diff(OhioCovid$totalTestResults-OhioCovid$positive, 1),NA)

early <- OhioCovid$date>=20200309 & OhioCovid$date <= 20200315
OhioInit <- data.frame(Y1=OhioCovid[early, "positiveIncrease"],
                       Y2=OhioCovid[early, "negativeIncrease"],
                       Date=OhioCovid[early, "RDate"])
OhioInit <- OhioInit[order(OhioInit$Date), ]
OhioInit[1,2] <- 11
OhioInit <- OhioInit[!(OhioInit$Y1==0 & OhioInit$Y2==0),]

shelter <- OhioCovid$date>=20200326 & OhioCovid$date <= 20200530
OhioShelter <- data.frame(Y1=OhioCovid[shelter, "positiveIncrease"],
                          Y2=OhioCovid[shelter, "negativeIncrease"],
                          Date=OhioCovid[shelter, "RDate"])
OhioShelter <- OhioShelter[order(OhioShelter$Date), ]
OhioShelter <- OhioShelter[OhioShelter$Y2>0&OhioShelter$Y1>0,]

reopen <- OhioCovid$date>=20200531 & OhioCovid$date <= 20200630
OhioReopen <- data.frame(Y1=OhioCovid[reopen, "positiveIncrease"],
                          Y2=OhioCovid[reopen, "negativeIncrease"],
                          Date=OhioCovid[reopen, "RDate"])
OhioReopen <- OhioReopen[order(OhioReopen$Date), ]
OhioReopen <- OhioReopen[OhioReopen$Y2>0&OhioReopen$Y1>0,]

OhioSim <- NULL
for (i in 1:dim(OhioInit)[1]){
    y1 <- OhioInit$Y1[i]
    y2 <- OhioInit$Y2[i]
    temp <- latent.gibbs(Y1=y1, Y2=y2, alpha=al, beta=be,
                         a1=a1, b1=b1, a2=a2, b2=b2,
                         nsims=10000, warmup=T, thin=5)
    theta <- fitdistr(temp$theta, densfun="beta",
                      start=list(shape1=2, shape2=2))
    fp <- fitdistr(temp$fp, densfun="beta",
                      start=list(shape1=2, shape2=2))
    fn <- fitdistr(temp$fn, densfun="beta",
                      start=list(shape1=2, shape2=2))

    OhioSim <- rbind(OhioSim, rvsims(as.matrix(temp)))
    al <- theta[[1]][1]
    be <- theta[[1]][2]
    a1 <- fp[[1]][1]
    b1 <- fp[[1]][2]
    a2 <- fn[[1]][1]
    b2 <- fn[[1]][2]
}

al=2; be=2; a1=3; b1=23; a2=2; b2=22;
temp <- OhioCovid[OhioCovid$date>=20200309 & OhioCovid$date <= 20200625,
                  c("RDate",  "positiveIncrease",
                    "negativeIncrease")]

OhioUp <- data.frame(Y1=temp$positiveIncrease, Y2=temp$negativeIncrease,
                     Date=temp$RDate)
OhioUp <- OhioUp[order(OhioUp$Date), ]
OhioUp <- OhioUp[OhioUp$Y1>0 & OhioUp$Y2>0, ]

OhioSim2 <- NULL
for (i in 1:dim(OhioUp)[1]){
    y1 <- OhioUp$Y1[i]
    y2 <- OhioUp$Y2[i]
    temp <- latent.gibbs(Y1=y1, Y2=y2, alpha=al, beta=be,
                         a1=a1, b1=b1, a2=a2, b2=b2,
                         nsims=10000, warmup=T, thin=5)
    theta <- try(fitdistr(temp$theta, densfun="beta",
                          start=list(shape1=al, shape2=be))[[1]], silent=T)
    if (!is.null(attributes(theta)$class))
        theta <- beta_mm(temp$theta)
##    fp <- try(fitdistr(temp$fp, densfun="beta",
  ##                    start=list(shape1=a1, shape2=b1))[[1]], silent=T)
    ##if (!is.null(attributes(fp)$class))
        fp <- beta_mm(temp$fp)
    ##fn <- try(fitdistr(temp$fn, densfun="beta",
      ##                start=list(shape1=a2, shape2=b2))[[1]], silent=T)
##    if (!is.null(attributes(fn)$class))
        fn <- beta_mm(temp$fn)

    OhioSim2 <- rbind(OhioSim2, rvsims(as.matrix(temp)))
    al <- theta[1]
    be <- theta[2]
  ##  if (al+be>1000){
  ##      al <- al/10
  ##      be <- be/10
  ##  }
    a1 <- fp[1]
    b1 <- fp[2]
  ##  if (a1+b1>1000){
  ##      a1 <- a1/10
  ##      b1 <- b1/10
  ##  }
    a2 <- fn[1]
    b2 <- fn[2]
  ##  if (a2+b2>1000){
  ##      a2 <- a2/10
  ##      b2 <- b2/10
  ##  }
}

OutDate <- format(OhioUp$Date, format="%b-%d")
plotTheta <- OhioSim2[,1]
names(plotTheta) <- OutDate
plotFp <- OhioSim2[,2]
names(plotFp) <- OutDate
plotFn <- OhioSim2[,3]
names(plotFn) <- OutDate
tikz(file=paste(plotDIRch5, "OhioCovidGibbs.tex", sep="/"),
     height=4, width=6, standAlone=T)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(plotTheta, xlab="$\\theta$", cex=0.75)
mlplot(plotFp, xlab="$f_p$", cex=0.75)
mlplot(plotFn, xlab="$f_n$", cex=0.75)
dev.off()

## Using Stan ##
## simulation 2: sequential updating
sim_sum2 <- NULL
al=2; be=2; a1=3; b1=23; a2=2; b2=22;
nn <- dim(sim_data50)[1]

## Stan model
latent_Bin <- "
data{
  int<lower=1> n; // sample size
  int<lower=0> y1[n]; // number of positives
  int<lower=1> Nt[n]; // number of trials
  real alpha;
  real beta;
  real a1;
  real b1;
  real a2;
  real b2;
}
parameters{
  real<lower=0,upper=0.75> theta; //prevalence
  real<lower=0,upper=0.75> fp; //false positive probability
  real<lower=0,upper=0.75> fn; //false negative probability
}
model{
  theta ~ beta(alpha, beta);
  fp ~ beta(a1, b1);
  fn ~ beta(a2, b2);
  for (i in 1:n){
    real temp[Nt[i]+1];
    real PP;
    for (j in 1:(Nt[i]+1)){
      PP = ((1-fn)*(j-1)+fp*(Nt[i]-(j-1)))/Nt[i];
      temp[j] = (j-1)*log(theta)+
                (Nt[i]-(j-1))*log(1-theta)+
                y1[i]*log(PP)+
                (Nt[i]-y1[i])*log(1-PP);
    }
    target+=log_sum_exp(temp);
  }
}
"
## Stan model2 ordered fn and fp
latent_Bin2 <- "
data{
  int<lower=1> n; // sample size
  int<lower=0> y1[n]; // number of positives
  int<lower=1> Nt[n]; // number of trials
  real alpha;
  real beta;
  real a1;
  real b1;
  real a2;
  real b2;
}
parameters{
  real<lower=0,upper=0.5> theta; //prevalence
  ordered[2] ff; //false positive probability
}
model{
  theta ~ beta(alpha, beta);
  ff[1] ~ beta(a1, b1);
  ff[2] ~ beta(a2, b2);
  for (i in 1:n){
    real temp[Nt[i]+1];
    real PP;
    for (j in 1:(Nt[i]+1)){
      PP = ((1-ff[2])*(j-1)+ff[1]*(Nt[i]-(j-1)))/Nt[i];
      temp[j] = (j-1)*log(theta)+
                (Nt[i]-(j-1))*log(1-theta)+
                y1[i]*log(PP)+
                (Nt[i]-y1[i])*log(1-PP);
    }
    target+=log_sum_exp(temp);
  }
}
"

## Stan model 3 through Bernoulli process
latent_Bern <- "
data{
  int<lower=1> n; // sample size
  int<lower=0> y1[n]; // number of positives
  int<lower=1> y2[n]; // number of negatives
  real alpha;
  real beta;
  real a1;
  real b1;
  real a2;
  real b2;
}
parameters{
  real theta; //prevalence
  real fp; //false positive probability
  real fn; //false negative probability
}
model{
  real temp1[2];
  real temp2[2];
  theta ~ beta(alpha, beta);
  fp ~ beta(a1, b1);
  fn ~ beta(a2, b2);
  temp1[1] = log(1-theta)+log(fp);
  temp1[2] = log(theta)+log(1-fn);
  temp2[1] = log(1-theta)+log(1-fp);
  temp2[2] = log(theta)+log(fn);
  for (i in 1:n){
    target+=y1[i]*log_sum_exp(temp1)+y2[i]*log_sum_exp(temp2);
  }
}
"
Bin_input <- function(Data=sim_data50, al=3, be=27,
                      a1=4, b1=16, a2=3, b2=7, n.chains=nchains){
    n <- dim(Data)[1]
    y1 <- Data$Y1
    Nt <- Data$Y2+y1
    data_in <- list(n=n, y1=y1, Nt=Nt, alpha=al, beta=be, a1=a1, b1=b1,
                    a2=a2, b2=b2)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(theta=runif(1,0,0.5), fp=runif(1,0,0.5),
                           fn=runif(1,0,0.5))
    paras <- c("theta","fp","fn")
    return(list(data=data_in, inits=inits, paras=paras,
                chains=n.chains))
}

Bin_input2 <- function(Data=sim_data50, al=3, be=27,
                      a1=4, b1=16, a2=3, b2=7, n.chains=nchains){
    n <- dim(Data)[1]
    y1 <- Data$Y1
    Nt <- Data$Y2+y1
    data_in <- list(n=n, y1=y1, Nt=Nt, alpha=al, beta=be, a1=a1, b1=b1,
                    a2=a2, b2=b2)
    inits <- list()
    for (i in 1:n.chains){
        fi <- runif(1,0,0.5)
        inits[[i]] <- list(theta=runif(1,0.05,0.5), ff=c(fi,fi+0.1))
    }
    paras <- c("theta","ff")
    return(list(data=data_in, inits=inits, paras=paras,
                chains=n.chains))
}
Bern_input <- function(Data=sim_data50, al=3, be=27,
                      a1=4, b1=16, a2=3, b2=7, n.chains=nchains){
    n <- dim(Data)[1]
    y1 <- Data$Y1
    y2 <- Data$Y2
    data_in <- list(n=n, y1=y1, y2=y2, alpha=al, beta=be, a1=a1, b1=b1,
                    a2=a2, b2=b2)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(theta=runif(1,0,0.5), fp=runif(1,0.05,0.5),
                           fn=runif(1,0,0.5))
    paras <- c("theta","fp","fn")
    return(list(data=data_in, inits=inits, paras=paras,
                chains=n.chains))
}

## simulation -- works well when priors for fp and fn are apropriate
stan_fit <- stan_model(model_code=latent_Bern)
sim_data50 <- rbin_imperfect(n=50, theta=0.3, fp=0.05, fn=0.1,N=100)
input.to.stan <- Bern_input(Data=sim_data50, al=1, be=1,
                            a1=1, b1=19, a2=1, b2=9)
fit2keep_sim1 <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains)#,
##                     control=list(adapt_delta=0.99))
print(fit2keep_sim1)
pairs(fit2keep_sim1, pars=c("theta","fp","fn"))
traceplot(fit2keep_sim1, pars=c("theta","fp","fn"))

input.to.stan <- Bern_input(Data=sim_data50, al=30, be=70,
                            a1=1, b1=19, a2=1, b2=1)
fit2keep_sim2 <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains)#,
#                     control=list(adapt_delta=0.99))
print(fit2keep_sim2)
pairs(fit2keep_sim2, pars=c("theta","fp","fn"))
traceplot(fit2keep_sim2, pars=c("theta","fp","fn"))

input.to.stan <- Bin_input(Data=OhioInit, al=1, be=1,
                           a1=3, b1=23, a2=2, b2=22)
stan_fit <- stan_model(model_code=latent_Bin)
fit2keep <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains)#,
##                     control=list(adapt_delta=0.9))
print(fit2keep)
pairs(fit2keep, pars=c("theta","fp","fn"))
traceplot(fit2keep, pars=c("theta","fp","fn"))

theta_prior <- beta_mm(rstan::extract(fit2keep, pars="theta",
                                permuted=T)[[1]])
fp_prior <- beta_mm(rstan::extract(fit2keep, pars="fp",
                                permuted=T)[[1]])
fn_prior <- beta_mm(rstan::extract(fit2keep, pars="fn",
                                permuted=T)[[1]])

input.to.stan <- Bin_input(Data=OhioShelter[3:4,], al=theta_prior[1],
                           be=theta_prior[2], a1=fp_prior[1],
                           b1=fp_prior[2], a2=fn_prior[1],
                           b2=fn_prior[2])
stan_fit <- stan_model(model_code=latent_Bin)
fit2keep <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains)
print(fit2keep)
pairs(fit2keep, pars=c("theta","fp","fn"))
traceplot(fit2keep, pars=c("theta","fp","fn"))

input.to.stan <- Bin_input2(Data=sim_data50[1:5,], al=2, be=2,
                           a1=3, b1=23, a2=2, b2=22)
stan_fit2 <- stan_model(model_code=latent_Bin2)
fit2keep2 <- sampling(stan_fit2, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.9))
print(fit2keep2)
pairs(fit2keep2, pars=c("theta","ff"))
traceplot(fit2keep2, pars=c("theta","ff"))

## Bernoulli formulation
stan_fit <- stan_model(model_code=latent_Bern)

## Initial period
input.to.stan <- Bern_input(Data=OhioInit, al=1, be=1,
                           a1=3, b1=23, a2=2, b2=22)
fit2keepI <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepI)
pairs(fit2keepI, pars=c("theta","fp","fn"))
traceplot(fit2keepI, pars=c("theta","fp","fn"))

## Shelter-in-place
theta_prior <- beta_mm(rstan::extract(fit2keepI, pars="theta",
                                permuted=T)[[1]])
fp_prior <- beta_mm(rstan::extract(fit2keepI, pars="fp",
                                permuted=T)[[1]])
fn_prior <- beta_mm(rstan::extract(fit2keepI, pars="fn",
                                permuted=T)[[1]])

input.to.stan <- Bern_input(Data=OhioShelter, al=1, be=1,
                           a1=3, b1=23, a2=2, b2=22)
input.to.stan <- Bern_input(Data=OhioShelter, al=theta_prior[1],
                           be=theta_prior[2], a1=fp_prior[1],
                           b1=fp_prior[2], a2=fn_prior[1],
                           b2=fn_prior[2])
fit2keepS <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepS)
pairs(fit2keepS, pars=c("theta","fp","fn"))
traceplot(fit2keepS, pars=c("theta","fp","fn"))

## Reopening
theta_prior <- beta_mm(rstan::extract(fit2keepS, pars="theta",
                                permuted=T)[[1]])
fp_prior <- beta_mm(rstan::extract(fit2keepS, pars="fp",
                                permuted=T)[[1]])
fn_prior <- beta_mm(rstan::extract(fit2keepS, pars="fn",
                                permuted=T)[[1]])

input.to.stan <- Bern_input(Data=OhioReopen, al=1, be=1,
                           a1=3, b1=23, a2=2, b2=22)
input.to.stan <- Bern_input(Data=OhioReopen, al=theta_prior[1],
                           be=theta_prior[2], a1=fp_prior[1],
                           b1=fp_prior[2], a2=fn_prior[1],
                           b2=fn_prior[2])
fit2keepR <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepR)
pairs(fit2keepR, pars=c("theta","fp","fn"))
traceplot(fit2keepR, pars=c("theta","fp","fn"))

## redo with same difused priors
## Initial period
input.to.stan <- Bern_input(Data=OhioInit, al=1, be=1,
                           a1=3, b1=23, a2=2, b2=22)
fit2keepI <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepI)
pairs(fit2keepI, pars=c("theta","fp","fn"))
traceplot(fit2keepI, pars=c("theta","fp","fn"))
thetaI <- rstan::extract(fit2keepI, pars="theta",permuted=T)[[1]]
fpI <- rstan::extract(fit2keepI, pars="fp",permuted=T)[[1]]
fnI <- rstan::extract(fit2keepI, pars="fn",permuted=T)[[1]]
subb <- sample(1:2504, size=100)
tikz(file=paste(plotDIRch5, "OHcovidIpairs.tex", sep="/"),
     height=1.75, width=5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01)
plot(thetaI[subb], fpI[subb], xlab="$\\theta$", ylab="$f_p$")
plot(thetaI[subb], fnI[subb], xlab="$\\theta$", ylab="$f_n$")
plot(fnI[subb], fpI[subb], xlab="$f_n$", ylab="$f_p$")
dev.off()

## Shelter-in-place
theta_prior <- beta_mm(rstan::extract(fit2keepI, pars="theta",
                                permuted=T)[[1]])
fp_prior <- beta_mm(rstan::extract(fit2keepI, pars="fp",
                                permuted=T)[[1]])
fn_prior <- beta_mm(rstan::extract(fit2keepI, pars="fn",
                                permuted=T)[[1]])

input.to.stan <- Bern_input(Data=OhioShelter, al=theta_prior[1],
                            be=,theta_prior[2],
                            a1=fp_prior[1], b1=fp_prior[2],
                            a2=fn_prior[1], b2=fn_prior[2])
fit2keepS2 <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepS2)
pairs(fit2keepS2, pars=c("theta","fp","fn"))
traceplot(fit2keepS2, pars=c("theta","fp","fn"))

thetaS <- rstan::extract(fit2keepS2, pars="theta",permuted=T)[[1]]
fpS <- rstan::extract(fit2keepS2, pars="fp",permuted=T)[[1]]
fnS <- rstan::extract(fit2keepS2, pars="fn",permuted=T)[[1]]

tikz(file=paste(plotDIRch5, "OHcovidSpairs.tex", sep="/"),
     height=1.75, width=5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01)
plot(thetaS[subb], fpS[subb], xlab="$\\theta$", ylab="$f_p$")
plot(thetaS[subb], fnS[subb], xlab="$\\theta$", ylab="$f_n$")
plot(fnS[subb], fpS[subb], xlab="$f_n$", ylab="$f_p$")
dev.off()

## Reopening
theta_prior <- beta_mm(rstan::extract(fit2keepS2, pars="theta",
                                permuted=T)[[1]])
fp_prior <- beta_mm(rstan::extract(fit2keepS2, pars="fp",
                                permuted=T)[[1]])
fn_prior <- beta_mm(rstan::extract(fit2keepS2, pars="fn",
                                permuted=T)[[1]])

input.to.stan <- Bern_input(Data=OhioReopen, al=theta_prior[1],
                            be=,theta_prior[2],
                            a1=fp_prior[1], b1=fp_prior[2],
                            a2=fn_prior[1], b2=fn_prior[2])
fit2keepR2 <- sampling(stan_fit, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$paras,
                     iter=niters, thin=nthin,
                     chains=input.to.stan$chains,
                     control=list(adapt_delta=0.99))
print(fit2keepR2)
pairs(fit2keepR2, pars=c("theta","fp","fn"))
traceplot(fit2keepR2, pars=c("theta","fp","fn"))


thetaR <- rstan::extract(fit2keepR2, pars="theta",permuted=T)[[1]]
fpR <- rstan::extract(fit2keepR2, pars="fp",permuted=T)[[1]]
fnR <- rstan::extract(fit2keepR2, pars="fn",permuted=T)[[1]]

tikz(file=paste(plotDIRch5, "OHcovidRpairs.tex", sep="/"),
     height=1.75, width=5, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,1,1), mgp=c(1.25,0.125,0), tck=0.01)
plot(thetaR[subb], fpR[subb], xlab="$\\theta$", ylab="$f_p$")
plot(thetaR[subb], fnR[subb], xlab="$\\theta$", ylab="$f_n$")
plot(fnR[subb], fpR[subb], xlab="$f_n$", ylab="$f_p$")
dev.off()
