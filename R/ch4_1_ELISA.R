source("FrontMatter.R")

## Chapter 4

plotDIRch4 <- paste(plotDIR, "chapter4", "figures", sep="/")

packages(rv)
packages(rstan)
packages(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 10000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Nonlinear regression and missing data problems

## Self-starter function of a 4-parameter logistic function:
## y = al4+(al1-al4)/(1+(x/al3)^al2)
## the function SSfpl is for a different form of fpl:
## A+(B-A)/(1+exp((xmid-input)/scal))
## these two functions are the same (input = log(x))

## the mean function
fplModel <- function(input, al1, al2, al3, al4){
    .x <- input+0.0001
    .expr1 <- (.x/al3)^al2
    .expr2 <- al1-al4
    .expr3 <- 1 + .expr1
    .expr4 <- .x/al3
    .value <- al4 + .expr2/.expr3
    .grad <- array(0, c(length(.value), 4L),
                   list(NULL, c("al1","al2","al3","al4")))
    .grad[,"al1"] <- 1/.expr3
    .grad[,"al2"] <- -.expr2*.expr1*log(.expr4)/.expr3^2
    .grad[,"al3"] <- .expr1*.expr2*(al2/al3)/.expr3^2
    .grad[,"al4"] <- .expr1/(1+.expr1)
    attr(.value, "gradient") <- .grad
    .value
}

## initial values
fplModelInit <- function(mCall, LHS, data){
    xy <- sortedXyData(mCall[["input"]], LHS, data)
   if (nrow(xy) < 5) {
        stop("too few distinct input values to fit a four-parameter logistic")
    }
    rng <- range(xy$y)
    drng <- diff(rng)
    xy$prop <- (xy$y-rng[1]+0.05*drng)/(1.1*drng)
    xy$logx <- log(xy$x+0.0001)
    ir <- as.vector(coef(lm(I(log(prop/(1-prop))) ~ logx, data=xy)))
    pars <- as.vector(coef(nls(y~cbind(1, 1/(1+(x/exp(lal3))^al2)),
                               data=xy,
                               start=list(al2=-ir[2], lal3=-ir[1]/ir[2]),
                               algorithm="plinear")))
    value <- c(pars[4]+pars[3], pars[1], exp(pars[2]), pars[3])
    names(value) <- mCall[c("al1","al2","al3","al4")]
    value
}

SSfpl2 <- selfStart(fplModel, fplModelInit, c("al1","al2","al3","al4"))


## Toledo Water Crisis
Toledo <- read.csv(paste(dataDIR, "ToledoCrisis.csv", sep="/"), header=T)
names(Toledo)[1] <- "SampleID"

head(Toledo)

ggplot(data=Toledo, aes(y=Absorbance, x=jitter(Concentration), color=Test))+
    geom_point()

### Default ELISA model
### Stan Model (single test kit)
ELISA1 <- "
data {
int<lower=0> N; // sample size (training data, non-zero)
int<lower=0> M; // observed ODs (for estimating concentration)
real y[N]; // response data
real x[N]; // predictor data
real<lower=1> dlf[M]; // dilution factor of testing data
real y0[M]; // test data response
int MM; //unique samples
int ws[M]; //unique water samples
}
parameters {
real<lower=0> delta;
real<lower=0> th2;
real<lower=0> th3;
real<lower=0> th4;
real<lower=0> sigma;
real<lower=0> x0[MM];
}
transformed parameters{
  real mu[N];
  real mu0[M];
  for (i in 1:N){
    mu[i] = th4 - delta/(1+(x[i]/th3)^(-th2));
  }
  for (i in 1:M){
    mu0[i] = th4 - delta/(1+((x0[ws[i]]/dlf[i])/th3)^(-th2));
  }
}
model {
  x0 ~ normal(0,2.5);
  delta~normal(1,1);
  th2~normal(0,2);
  th3~normal(0,2);
  th4~normal(1,1);
  target += normal_lpdf(y| mu, sigma);
  target += normal_lpdf(y0| mu0, sigma);
}
"
stan.fit2 <- stan_model(model_code=ELISA1)

stan.in1 <- function(data = Toledo, chains=nchains){
    n <- dim(data)[1]
    temp  <- is.na(data$Concentration)
    y <- data$Absorbance[!temp]
    x <- data$Concentration[!temp]
    N <- n-sum(temp)
    M <- sum(temp)
    y0 <- data$Absorbance[temp]
    dlf <- data$Dilution[temp]
    ws <- as.numeric(ordered(data$SampleID[temp]))
    stan.dat <- list(N=N, M=M, MM=max(ws), y=y, x=x, minY=min(c(y,y0)),
                     y0=y0, ws=ws, dlf=dlf)
  inits <- list()
  for (i in 1:chains)
      inits[[i]] <- list(delta = abs(rnorm(1,0,0.25)), th2=abs(rnorm(1)),
                         th3 = abs(rnorm(1)), th4=abs(rnorm(1,1,0.25)),
                         x0 = abs(rnorm(max(ws))),
                         sigma = runif(1))
  parameters <- c("x0","delta", "th2","th3","th4","sigma")
  return(list(para=parameters, data=stan.dat, inits=inits,
              n.chains=chains))
}

test_stan <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stan.in1(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter=niters, thin=nthin,
                         chains=input.to.stan$n.chains)
    print(fit2keep)
    test_stan[[i]] <- rstan::extract(fit2keep)
}
save(test_stan, file="ELISA_stan1.RData")

## compare to nls:

fpl <- function(x, th1, th2, th3, th4){
    return(th4+(th1-th4)/(1+(x/th3)^th2))
}

inv.fpl <- function(y, th1, th2, th3, th4){
    return(th3*((th1-y)/(y-th4))^(1/th2))
}

tm <- NULL
tm_pred <- NULL
for (i in 1:6){
    tm <- rbind(tm, coef(nls(Absorbance ~ (al1-al4)/(1+(Concentration/al3)^al2)+al4,
                start=list(al1=1.2, al2=0.64, al3=0.16, al4=.15),
                data=Toledo[Toledo$Test==i & !is.na(Toledo$Concentration),])))
    ##tm <- rbind(tm, coef(nls(Absorbance ~ SSfpl2(Concentration,
    ##                                             al1, al2, al3, al4),
    ##                         data=Toledo[Toledo$Test==i & !is.na(Toledo$Concentration),])))
    temp <- Toledo[Toledo$Test==i & is.na(Toledo$Concentration),]
    temp$Concentration <- inv.fpl(temp$Absorbance, tm[i,1], tm[i,2], tm[i,3], tm[i,4])
    tm_pred <- rbind(tm_pred, temp)
}

## Compare
## Test 1

test_summ1 <- function(test=1){
    tm_predTest <- tm_pred%>%filter(Test==test)%>%arrange(SampleID)
    testTStan <-
        summary(rvsims(as.matrix(as.data.frame(test_stan[[test]]$x0))))
    tm_predTest$StanPred <-
        testTStan$mean[as.numeric(ordered(tm_predTest$SampleID))]
    tm_predTsmM <- tm_predTest%>%group_by(SampleID)%>%
        summarise(ELISA=mean(Concentration*Dilution, na.rm=T),
                  Stan=mean(StanPred, na.rm=T))
    plot(Concentration*Dilution ~ StanPred, data=tm_predTest, log="xy")
    abline(0,1)
    points(ELISA ~ Stan, data=tm_predTsmM, col="grey", pch=16)
    invisible(list(tm_predTest, tm_predTsmM))
}

test1summary <- test_summ1(1)
test2summary <- test_summ1(2)
test3summary <- test_summ1(3)
test4summary <- test_summ1(4)
test5summary <- test_summ1(5)
test6summary <- test_summ1(6)

## plotting all estimated concentrations
tm_predAll <- rbind(test1summary[[1]], test2summary[[1]], test3summary[[1]],
                    test4summary[[1]], test5summary[[1]], test6summary[[1]])
tm_predAllm <- rbind(test1summary[[2]], test2summary[[2]], test3summary[[2]],
                     test4summary[[2]], test5summary[[2]], test6summary[[2]])

tikz(file=paste(plotDIRch4, "ELISA_comp1.tex", sep="/"),
     height=4, width=4, standAlone=F)
par(mar=c(3,4,1,0.5), mgp=c(1.75,0.125,0), tck=0.01, las=1)
plot(Concentration*Dilution ~ StanPred, data=tm_predAll, log="xy", axes=F,
     xlab="Bayesian estimates", ylab="Least squares estimates")
axis(1)
axis(2, at=c(0.01,0.05,0.1,0.5,1,5,10),
     label=c(0.01,0.05,0.1,0.5,1,5,10))
box()
points(ELISA~Stan, data=tm_predAllm, pch=16, col="grey")
abline(0,1, col="grey")
dev.off()

## FPL in log-concentration scale (limit version)
ELISA2 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; // observed ODs (for estimating concentration)
  int<lower=0> n0; // number of 0 standard solutions
  real zeroY[n0]; // observed 0 conc OD
  real y[N]; // response data
  real z[N]; // predictor data
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
}
parameters {
  real<lower=-5,upper=5> az0[MM];
  real<lower=0> xA;
  real<lower=0> xB;
  real<lower=-5,upper=1.5> xmid;
  real<lower=0> scal;
  real<lower=0> sigma;
}
model {
  real mu[N];
  real mu0[M];
  az0 ~ normal(0,2);
  xA ~ normal(2,1);
  xB ~ normal(0,1);
  xmid ~ normal(-1, 1);
  scal ~ normal(1.3, 1);

  zeroY ~ normal(xA, sigma);
  for (i in 1:N){
    mu[i] = xA + (xB-xA)/(1+exp((xmid-z[i])/scal));
  }
  y ~ normal(mu, sigma);
  for (i in 1:M){
    mu0[i] = xA + (xB-xA)/(1+exp((xmid-(az0[ws[i]]-log(dlf[i])))/scal));
  }
  y0 ~ normal(mu0, sigma);
}
"

## compile Stan input file (limit version)
stanInput <- function(data=Toledo, n.chains=nchains){
    ##infile <- infile[infile$conc!=0,]
    ##testset <- infile$training==6
    testset <- is.na(data$Concentration)
    train <- (!testset) & data$Concentration != 0
    dlf <- data$Dilution[testset]
    zerodata <- (!is.na(data$Concentration)) & data$Concentration==0
    M <- sum(testset) #3x2 additional samples for verification
                                        #if (is.null(trainset)){
    N <- sum(train)  #6x2 original standard solutions
    ws <- as.numeric(ordered(data$SampleID[testset]))
    my.data <- list(N = N,M=M, z=log(data$Concentration[train]),
                    y = data$Absorbance[train],
                    y0 = data$Absorbance[testset],
                    ws=ws, MM=max(ws),
                    dlf=dlf, zeroY = data$Absorbance[zerodata],
                    n0 = sum(zerodata))
    my.paras <- c("az0", "xA", "xB", "xmid", "scal", "sigma")
    my.init <- list()
    for (i in 1:n.chains)
        my.init[[i]] <- list(xA=runif(1), xB=runif(1),
                             xmid=-runif(1),
                             scal=runif(1),
                             az0=-runif(max(ws)),
                             sigma=runif(1))
    return(list(data=my.data, inits=my.init, pars=my.paras, n.chains=n.chains))
}

stan.fit2 <- stan_model(model_code=ELISA2)
test_stan2 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stanInput(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter = niters, thin = nthin,
                         chains=input.to.stan$n.chains)#,
##                         control=list(adapt_delta=0.95, max_treedepth=20))
    print(fit2keep)
    test_stan2[[i]] <- rstan::extract(fit2keep)
}
save(test_stan2, file="ELISA_stan2_limit.RData")

test_summ <- function(stanout=test_stan2, test=1){
    tm_predTest <- tm_pred%>%filter(Test==test)%>%arrange(SampleID)
    testTStan <-
        summary(rvsims(as.matrix(as.data.frame(stanout[[test]]$az0))))
    tm_predTest$StanPred <-
        testTStan$mean[as.numeric(ordered(tm_predTest$SampleID))]
    tm_predTsmM <- tm_predTest%>%group_by(SampleID)%>%
        summarise(ELISA=mean(Concentration*Dilution, na.rm=T),
                  Stan=mean(StanPred, na.rm=T))
    plot(Concentration*Dilution ~ exp(StanPred), data=tm_predTest, log="xy")
    abline(0,1)
    points(ELISA ~ exp(Stan), data=tm_predTsmM, col="grey", pch=16)
    invisible(list(tm_predTest, tm_predTsmM))
}

test1summary <- test_summ(stanout=test_stan2, test=1)
test2summary <- test_summ(test=2)
test3summary <- test_summ(test=3)
test4summary <- test_summ(test=4)
test5summary <- test_summ(test=5)
test6summary <- test_summ(test=6)



## FPL in log-concentration scale (no 0 solutions)
ELISA3 <- "
data {
  int<lower=0> N; // sample size (training data, non-zero)
  int<lower=0> M; // observed ODs (for estimating concentration)
  real y[N]; // response data
  real z[N]; // predictor data
  real<lower=1> dlf[M]; // dilution factor of testing data
  real y0[M]; // test data response
  int MM; //unique samples
  int ws[M]; //unique water samples
}
parameters {
  real<lower=-5,upper=5> az0[MM];
  real<lower=0> xA;
  real<lower=0> xB;
  real<lower=-5,upper=1> xmid;
  real<lower=0> scal;
  real<lower=0> sigma;
}
model {
  real mu[N];
  real mu0[M];
  az0 ~ normal(0,2);
  xA ~ normal(2,1);
  xB ~ normal(0,1);
  xmid ~ normal(-1, 1);
  scal ~ normal(1, 1);

  for (i in 1:N){
    mu[i] = xA + (xB-xA)/(1+exp((xmid-z[i])/scal));
  }
  y ~ normal(mu, sigma);
  for (i in 1:M){
    mu0[i] = xA + (xB-xA)/(1+exp((xmid-(az0[ws[i]]-log(dlf[i])))/scal));
    }
  y0 ~ normal(mu0, sigma);
}
"

## input (no 0 standards)
stanInput3 <- function(data=Toledo, n.chains=nchains, zero=0.0001){
    ##infile <- infile[infile$conc!=0,]
    ##testset <- infile$training==6
    testset <- is.na(data$Concentration)
    if (is.na(zero)){
        train <- (!testset) & data$Concentration != 0
        z  <- data$Concentration[train]
    } else {
        train <- (!testset)
        z  <- data$Concentration[train]
        z[z==0] <- zero
    }
    dlf <- data$Dilution[testset]
    M <- sum(testset) #3x2 additional samples for verification
                                        #if (is.null(trainset)){
    N <- sum(train)  #6x2 original standard solutions
    ws <- as.numeric(ordered(data$SampleID[testset]))
    my.data <- list(N = N,M=M, z=log(z),
                    y = data$Absorbance[train],
                    y0 = data$Absorbance[testset],
                    ws=ws, MM=max(ws), dlf=dlf)
    my.paras <- c("az0", "xA", "xB", "xmid", "scal", "sigma")
    my.init <- list()
    for (i in 1:n.chains)
        my.init[[i]] <- list(xA=runif(1), xB=runif(1),
                             xmid=-runif(1),
                             scal=runif(1),
                             az0=-runif(max(ws)),
                             sigma=runif(1))
    return(list(data=my.data, inits=my.init, pars=my.paras,
                n.chains=n.chains))
}

stan.fit3 <- stan_model(model_code=ELISA3)
test_stan3 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stanInput3(data=Toledo[Toledo$Test==i,])
    fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter = niters, thin = nthin,
                         chains=input.to.stan$n.chains,
                         control=list(adapt_delta=0.95, max_treedepth=20))
    print(fit2keep)
    test_stan3[[i]] <- rstan::extract(fit2keep)
}
save(test_stan3, file="ELISA_stan3_0001.RData")

test1summary <- test_summ(stanout=test_stan3, test=1)
test2summary <- test_summ(stanout=test_stan3, test=2)
test3summary <- test_summ(stanout=test_stan3, test=3)
test4summary <- test_summ(stanout=test_stan3, test=4)
test5summary <- test_summ(stanout=test_stan3, test=5)
test6summary <- test_summ(stanout=test_stan3, test=6)

test_stan4 <- list()
for (i in 1:6){
    print(paste("Test:", i))
    input.to.stan <- stanInput3(data=Toledo[Toledo$Test==i,], zero=NA)
    fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
                         init=input.to.stan$inits,
                         pars = input.to.stan$para,
                         iter = niters, thin = nthin,
                         chains=input.to.stan$n.chains)#,
#                         control=list(adapt_delta=0.95, max_treedepth=20))
    print(fit2keep)
    test_stan4[[i]] <- rstan::extract(fit2keep)
}
save(test_stan4, file="ELISA_stan4_nas.RData")

test1summary <- test_summ(stanout=test_stan4, test=1)
test2summary <- test_summ(stanout=test_stan4, test=2)
test3summary <- test_summ(stanout=test_stan4, test=3)
test4summary <- test_summ(stanout=test_stan4, test=4)
test5summary <- test_summ(stanout=test_stan4, test=5)
test6summary <- test_summ(stanout=test_stan4, test=6)

####################################
## comparing alternatives -- test 2#
####################################
input.to.stan <- stanInput(data=Toledo[Toledo$Test==2,])
stan.fit2 <- stan_model(model_code=ELISA2)
fit2keep <- sampling(stan.fit2, data = input.to.stan$data,
                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter = niters, thin = nthin,
                     chains=input.to.stan$n.chains)##,
##                     control=list(adapt_delta=0.9, max_treedepth=20))
print(fit2keep)
test2_limit <- rstan::extract(fit2keep)

input.to.stan <- stanInput3(data=Toledo[Toledo$Test==2,], zero=NA)
stan.fit3 <- stan_model(model_code=ELISA3)
fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
#                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter = niters, thin = nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)
test2_nas <- rstan::extract(fit2keep)

input.to.stan <- stanInput3(data=Toledo[Toledo$Test==2,])
stan.fit3 <- stan_model(model_code=ELISA3)
fit2keep <- sampling(stan.fit3, data = input.to.stan$data,
#                     init=input.to.stan$inits,
                     pars = input.to.stan$para,
                     iter = niters, thin = nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)
test2_0001 <- rstan::extract(fit2keep)
rv0001sum  <- summary(rvsims(as.matrix(as.data.frame(test2_0001))))
rvnassum  <- summary(rvsims(as.matrix(as.data.frame(test2_nas))))
rvlimitsum  <- summary(rvsims(as.matrix(as.data.frame(test2_limit))))
cbind(rv0001sum$mean, rvnassum$mean, rvlimitsum$mean)

oo <- order(rv0001sum$mean[1:11])
rv0001  <- rvsims(as.matrix(as.data.frame(test2_0001)))
rvnas  <- rvsims(as.matrix(as.data.frame(test2_nas)))
rvlimit  <- rvsims(as.matrix(as.data.frame(test2_limit)))

rvest <- c(rv0001[1:11][oo],rvnas[1:11][oo],rvlimit[1:11][oo])
rvcoef <- c(rv0001[12:16],rvnas[12:16],rvlimit[12:16])

tikz(file=paste(plotDIRch4, "ELISA_comp2.tex", sep="/"),
     height=5, width=5.5, standAlone=F)
par(mfrow=c(2,1), mar=c(3,3,3,0.5), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(rvest, axes=F, xlab="Samples", ylab="log MC concentration")
axis(1,at=1:33, labels=rep(1:11, 3))
     abline(v=c(11.5,22.5))
axis(2)
box()
text(x=c(6,18,29), y=rep(-4.75,3), labels=c("0.0001", "NAs","limit"))

plot(rvcoef, axes=F, xlab="Coefficients", ylab="")
abline(v=c(5.5,10.5))
axis(1, at=c(1,3,5,6,8,10,11,13,15),
     labels=rep(c("$A$",  "$x_{mid}$", "$\\sigma$"), 3))
axis(3, at=c(2,4,7,9,12,14),
     labels=rep(c( "$B$",  "$s_{cal}$"), 3))
axis(2)
text(x=c(3,8,13), y=rep(-4,3), labels=c("0.0001", "NAs","limit"))
box()
dev.off()

save(test2_limit, test2_nas, test2_0001, file="test2compare.RData")

######################################################
## comparing alternatives -- numerical control (0.75)#
######################################################
## Selecting NConl
NCol <- numeric()
for (i in 1:6){
    temp  <- ordered(Toledo$SampleID[(Toledo$Test==i)&(is.na(Toledo$Concentration))])
    n <- 1:length(levels(temp))
    NCol[i] <- n[levels(temp)=="NConl"]
}

all_flp1 <- all_limit <- all_nas <- all_0001 <- rv(length=0)
for (i in 1:6){
    all_flp1 <- c(all_flp1, rvsims(as.matrix(as.data.frame(test_stan[[i]]$x0)))[NCol[i]])
    all_limit <- c(all_limit, rvsims(as.matrix(as.data.frame(test_stan2[[i]]$az0)))[NCol[i]])
    all_nas <- c(all_nas, rvsims(as.matrix(as.data.frame(test_stan4[[i]]$az0)))[NCol[i]])
    all_0001 <- c(all_0001, rvsims(as.matrix(as.data.frame(test_stan3[[i]]$az0)))[NCol[i]])
}

all_elisa <- c(log(all_flp1),all_0001,all_nas,all_limit)
tikz(file=paste(plotDIRch4, "ELISA_comp3.tex", sep="/"),
     height=3, width=5.5, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25,0.125,0), las=1,tck=0.01)
plot(all_elisa, xlab="tests", ylab="log MC concentration",
     ylim=c(-3,1),axes=F)
axis(1, at=1:24,labels=rep(1:6, 4))
axis(2)
box()
abline(v=c(6.5,12.5, 18.5))
abline(h=log(0.75))
text(x=c(3.25, 9.75, 15.76, 21.75), y=-2.75,
     labels=c("FLP1", "0.0001", "NAs", "limit"))
dev.off()
