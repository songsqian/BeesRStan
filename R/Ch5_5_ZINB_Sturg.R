############   JASAMP zero-inflated negative binomial model   ###########
##########################################################################
#use this when generating figures for the book
                                        #
source("FrontMatter.R")
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")

## load packages
packages(rstan)
packages(rv)
packages(arm)
packages(tikzDevice)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))
nchains <- min(c(parallel::detectCores(), 8))
niters <- 5000 #50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

zipdata <- read.csv(paste(dataDIR, "sdata_15_sub.csv", sep="/"),
                 header=T, sep=",") ## survey data
zipdata <- zipdata[order(zipdata$YEAR),]
head(zipdata)
##############################################

##############################################
## The model
##############################################
zinb1 <-"
data {
  int<lower=1> Np;     // total number of observations with positive count
  int<lower=1> N0;
  int<lower=1> Nyr;   // number of years (regions)
  int<lower=0> Yp[Np];  // response variable (positive counts only)
  int<lower=1> Ka;     // number of abundance model coefficients
  int<lower=1,upper=Nyr> yearP[Np];
  int<lower=1,upper=Nyr> year0[N0];
  matrix[Np, Ka] Xp;     // abundance model design matrix
  matrix[N0, Ka] X0;     // abundance model design matrix
  vector[Np] offsetP;  // sampling effort (net-soaking time)
  vector[N0] offset0;
  real Zp[Np];  // zero model design matrix
  real Z0[N0];  // zero model design matrix
}
parameters {
  real alpha0;
  vector[Nyr] b0;
  vector[Ka] beta;       // abundance model coefficients
  real<lower=0> phi;         // NB shape parameter
  real alpha[2];   // zero model coefficient
  real mu_b; // hyper mean of b0
  real<lower=0> tau_b;
}
transformed parameters{
  vector[Np] eta_p;// initialize linear predictor term
  vector[N0] eta_0;// initialize linear predictor term
  vector[Np] lambda_p;
  vector[N0] lambda_0;
  vector[Np] zi_p;          // initialize linear predictor term
  vector[N0] zi_0;          // initialize linear predictor term
  vector[Np] theta_p;
  vector[N0] theta_0;
  for (i in 1:Np){
    eta_p[i] = b0[yearP[i]]+ Xp[i] * beta + offsetP[i];
    zi_p[i]  = alpha[1] + Zp[i] * alpha[2];
  }
  for (i in 1:N0){
    eta_0[i] = b0[year0[i]]+ X0[i] * beta + offset0[i];
    zi_0[i]  = alpha[1] + Z0[i] * alpha[2];
  }
  lambda_p = exp(eta_p);
  theta_p = inv_logit(zi_p);
  lambda_0 = exp(eta_0);
  theta_0 = inv_logit(zi_0);
}
model {
  // priors
  tau_b ~ normal(0,5);
  mu_b ~ normal(0,5);
  b0 ~ normal(mu_b, tau_b);
  beta ~ normal(0,5);
  alpha ~ normal(0,5);
  target += gamma_lpdf(phi | 0.01, 0.01);

  // likelihood
  for(i in 1:N0){
    target += log_sum_exp(log(theta_0[i]),
                       log1m(theta_0[i]) +
                       neg_binomial_2_log_lpmf(0 | eta_0[i], phi));
    }
  target += log1m(theta_p) +
            neg_binomial_2_log_lpmf(Yp | eta_p, phi);
}
"

##############################################
fit <- stan_model(model_code = zinb1)


##############################################
## The data and initial values
##############################################
ZINB_in <- function(data=zipdata, n.chains=nchains){
    y <- data$CATCH
    tmp <- !is.na(y)
    data <- data[tmp,]
    N <- dim(data)[1]
    Y <- y[tmp]
    X <- cbind(data$TEMP-mean(data$TEMP),
               data$DTSF2-mean(data$DTSF2),
               log(data$DTSF2)-mean(log(data$DTSF2)))
    Kbeta <- dim(X)[2]
    Z <- data$TEMP-mean(data$TEMP)
    offsets <- log(data$Effort)
    year <- as.numeric(ordered(data$YEAR))
    nyr <- max(year)
    tmp <- y==0
    np <- sum(!tmp)
    n0 <- sum(tmp)
    data <- list(Np=np, N0=n0, Nyr=nyr, Ka=Kbeta,
                 Xp=X[!tmp,], X0=X[tmp,], Yp=y[!tmp],
                 Zp=Z[!tmp], Z0=Z[tmp],
                 offset0=offsets[tmp], offsetP=offsets[!tmp],
                 yearP=year[!tmp], year0=year[tmp])
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(b0=rnorm(nyr), alpha=rnorm(2), phi=(runif(1)),
                       mu_b=rnorm(1), tau_b=runif(1), beta=rnorm(Kbeta))
  paras <- c("b0","beta","alpha","mu_b","phi", "tau_b")
  return(list(data=data, init=inits, nchains=n.chains, para=paras ))
}

input.to.stan <- ZINB_in()
keep <- sampling(fit, data=input.to.stan$data,
                 init=input.to.stan$init,
                 pars=input.to.stan$para,
                 iter=niters,thin=nthin,
                 chains=input.to.stan$nchains)##,
##                 control=list(max_treedepth=20))

save(keep,file="zinb_sturg2.RData")


##############################################
## Processing Stan output and summaries - for book
##############################################
##keep <- load("zinb_sturg.RData")

fitcoef <- rvsims(as.matrix(as.data.frame(
extract(keep, permute=T))))
summary(fitcoef)


## Use this when generating figures for the book...
tikz(file=paste(plotDIRch5, "zinb_sturg_coefs.tex", sep="/"),
      height=2, width=4.75, standAlone=F)
par(mfrow=c(1,3), mar=c(3,3,3,1), mgp=c(1.5, 0.125,0), las=1, tck=0.01)
## Temp effect on zi
tmp <- substring(names(fitcoef), 1, 1)=="a"
zero_coef <- fitcoef[tmp]
x.tmp <- seq(0,15,by=1)
tmp.zi <- zero_coef[1] + zero_coef[2]*x.tmp
plot((invlogit(tmp.zi)),  pch=1, xaxt="n", cex=0.5,
     xlab="Temperature ($^{\\circ}$C)", ylab="Probability of false zero")
axis(1,at=c(1:16),labels=x.tmp)

## TEMP effect on counts
## b0
int <- substring(names(fitcoef), 1, 2)=="b0"
bets <- fitcoef[substring(names(fitcoef), 1,4)=="beta"]
dtsf <- zipdata$DTSF2
tmp.c <- mean(fitcoef[int]) + bets[1]*(x.tmp-mean(zipdata$TEMP)) + mean(log(zipdata$Effort))
plot(exp(tmp.c),  pch=1, xaxt="n", cex=0.5,
     xlab="Temperature ($^{\\circ}$C)", ylab="CPUE")
axis(1,at=c(1:16),labels=x.tmp)

## DTSF effect on counts
x.dtsf <- seq(0.05,1.05,by=0.05)
dtsf.c <- mean(fitcoef[int]) +  bets[2]*x.dtsf + bets[3]*log(x.dtsf) +  mean(log(zipdata$Effort))
plot(exp(dtsf.c), xaxt="n",  pch=1, cex=0.5,
     xlab="Distance to Salt Front (km)", ylab="CPUE")
axis(1,at=c(1:21),labels=round(seq(min(zipdata$DTSF),max(zipdata$DTSF),length.out=21),0))
abline(v=9.5, lty=3, col=gray(0,0.5))
text(9.5,1.5,"Salt-front",pos=2,srt=90)

dev.off()



## Adjusted index
matrix(c(seq(0,1,length.out=117),
  seq(1,117,by=1)-50),ncol=2,byrow=F)

tikz(file=paste(plotDIRch5, "zinbPred.tex", sep="/"),
     height=2.75, width=3.25, standAlone=F)
par(mar=c(3,3,3,1), mgp=c(1.5, 0.125,0), las=1, tck=0.01)
x.dtsf <- seq(0,1,by=0.01)
std_index <- exp((fitcoef[1:12]) )
plot(jitter(as.numeric(as.factor(zipdata$YEAR))),
     zipdata$CATCH/zipdata$Effort,pch=1,col="gray", cex=0.5, xaxt="n",
     xlab="Year", ylab="CPUE")

#abline(h=10,lty=3, col=gray(0,0.2))
#text(4,10,"Hypothetical \nManagement Target", pos=3,cex=0.75)

nom_index <- aggregate(CATCH/Effort~YEAR, data=zipdata,FUN="mean")
points(nom_index$'CATCH/Effort', col="black", typ="b",pch=4)

points((std_index), cex=0.5)
axis(1,at=c(1:12),labels=c(2004:2015))
dev.off()


## 10% percent change per year
# (nom_eval <- ifelse((nom$'y/offst'>=10)=="FALSE",0,1))
# (std_eval <- rvmean(index>=10))

nom_index
p_chng <- NULL
for(i in 2:12) p_chng[i-1] <- (std_index[i]-std_index[i-1])/std_index[i-1]


matrix(p_chng, ncol=4000,byrow=T)

rvsims(p_chng[[11]])
  p_chng <- c((std_index[2]-std_index[1])/std_index[1],
            (std_index[3]-std_index[2])/std_index[2],
            (std_index[4]-std_index[3])/std_index[3],
            (std_index[5]-std_index[4])/std_index[4],
            (std_index[6]-std_index[5])/std_index[5],
            (std_index[7]-std_index[6])/std_index[6],
            (std_index[8]-std_index[7])/std_index[7],
            (std_index[9]-std_index[8])/std_index[8],
            (std_index[10]-std_index[9])/std_index[9],
            (std_index[11]-std_index[10])/std_index[10],
            (std_index[12]-std_index[11])/std_index[11])

mp_chng <- (100*mean(p_chng))


tikz(file=paste(plotDIRch5, "zinb_percentChng.tex", sep="/"),
     height=2, width=3.75, standAlone=F)

par(mfrow=c(1,2), mar=c(3,3,3,1), mgp=c(1.5, 0.125,0),
    las=1, tck=0.01)
plot(100*p_chng, rvcol="lightgray", ylim=c(-100,500), pch=16,
     ylab="\\% change from previous year", xlab="Year",
     cex.lab=0.8,xaxt="n", cex=0.5)
axis(1,at=c(1:11),labels=c(2005:2015))
abline(h=10,col=gray(0,0.2),lty=3,lwd=3)
hist(unlist(mp_chng), xlim=c(0,40),
     xlab="Mean change 2005-2015 (\\%)", main="",col="gray", cex.lab=0.8)
abline(v=10,col=gray(0,0.2),lty=3, lwd=2)

dev.off()

##############################################


