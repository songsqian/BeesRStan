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

GLdata <- read.csv(paste(dataDIR, "rawdata.csv", sep="/"),header=T)
GLdata <- GLdata[!is.na(GLdata$STAID),]

sites <- read.csv(paste(dataDIR, "sites.csv", sep="/"), header=T)
sites <- sites[!is.na(sites$STAID),]
names(GLdata)
names(sites) <- c("n","STAID","SNAME","Year","DRNsqkm","basinclass",
                  "Ecoregion", "diatonEcoGroup","Algal","agN","agP",
                  "MRB3TNclass")
sites <- unique(sites[,c("STAID","Ecoregion")])

## temp <- merge(GLdata, sites, by="STAID", all=T)
## merge won't work
GLdata$Ecoregion <- NA
temp <- unique(GLdata$STAID)
for (i in 1:length(temp)){
  tt <- GLdata$STAID==temp[i]
  ss <- sites$STAID==temp[i]
  if (!is.null(tt) & sum(ss)!=0) {
    if (sum(ss)>1) stop("multiple")
    GLdata$Ecoregion[tt] <- sites$Ecoregion[ss]
  }
}
### STAID 5082625 not listed in file "sites" --
### TURTLE R AT TURTLE R STATE PARK NR ARVILLA, ND -- ecoregion == 48
GLdata$Ecoregion[is.na(GLdata$Ecoregion)] <- 48

## exploratory plots
GLdata$STAID <- as.numeric(ordered(GLdata$STAID))
GLdata$Date <- as.Date(as.character(GLdata$Date),"%m/%d/%Y")
GLdata$Year <- as.numeric(substring(as.character(GLdata$Date), 1,4))
GLdata <- GLdata[!is.na(GLdata$Year),]
GLdata$State <- ordered(as.vector(GLdata$State))
st <- state.abb[c(13, 14, 15, 22, 23, 32, 34, 35, 49)]
GLdata$state <- st[as.numeric(GLdata$State)]
GLdata$STAIDst <- paste(GLdata$state, GLdata$STAID, sep="")
GLdata$STAIDer <- paste(GLdata$Ecoregion, GLdata$STAID, sep="-")
bwplot(log(Tpwu)~state|factor(Year), data=GLdata)
bwplot(log(Total.Nitrogen.calc.)~state|factor(Year), data=GLdata)
bwplot(log(NH3PorgNwu)~state|factor(Year), data=GLdata)
xyplot(log(Tpwu)~jitter(Year)|state,cex=0.5, layout=c(3,3), data=GLdata,
       xlab="Year", ylab="log TP")
xyplot(log(Total.Nitrogen.calc.)~jitter(Year)|state, cex=0.5,layout=c(3,3),
       data=GLdata, xlab="Year", ylab="log TN")

histogram(~(Tpwu)|state, data=GLdata)
##The only site in the state of New York is set aside.
### fit model without NY and use Bayesian updating
## multilevel model for TP

TP.lmer <- lmer(log(Tpwu)~1+(1|Season)+(1|Ecoregion)+(1|Year)+(1|STAIDer),
                data=GLdata, subset=State!="NEW YORK")
TP.lm <- lm(log(Tpwu)~factor(Season)+factor(Year)+factor(STAIDer),
            data=GLdata, subset=State!="NEW YORK")
## effect of Ecoregion cannot be separated from the same of site using lm.

#TP.aov <- aov(log(Tpwu)~factor(Season)+factor(Year)+factor(STAIDer),
#              data=GLdata, subset=State!="NEW YORK")

## sum of squares for each site

## site + ecoregion effects
### in theory, it should be beta_2^2 + beta_3^2:
lm.coef <- coef(TP.lm)
lm.season <- c(0,lm.coef[2:4])
lm.year <- c(0, lm.coef[5:18])
lm.sn <- lm.season-mean(lm.season)
lm.yr <- lm.year-mean(lm.year)

lmer.yr <- ranef(TP.lmer)$Year[,1]
lmer.sn <- ranef(TP.lmer)$Season[,1]


ss.lmer <- ss.lm <- numeric()
k <- 0
for (i in 1:4)
    for (j in 1:15){
        k <- k+1
        ss.lmer[k] <- lmer.yr[j]^2+lmer.sn[i]^2
        ss.lm[k] <- lm.yr[j]^2+lm.sn[i]^2
    }


stdid <-  row.names(ranef(TP.lmer)$STAID)
stdid.frame<-data.frame(STAIDer=stdid, Ecoregion=substring(stdid, 1, 2))
year.season <- data.frame(Year=rep(1991:2002, each=4),
                          Season=rep(c("Fall","Spring","Summer","Winter"), 12))

pred.data <- data.frame(STAIDer=rep(stdid.frame$STAIDer, each=48),
                        Ecoregion=rep(stdid.frame$Ecoregion, each=48),
                        Year=rep(year.season$Year, 59),
                        Season=rep(year.season$Season, 59))
pred.data$Pred.lm <- predict(TP.lm, new=pred.data)
pred.data$Pred.lmer <- fixef(TP.lmer) +
  ranef(TP.lmer)$Season[as.numeric(ordered(pred.data$Season)),1] +
  ranef(TP.lmer)$Year[as.numeric(ordered(pred.data$Year)),1] +
  ranef(TP.lmer)$Ecoregion[as.numeric(ordered(pred.data$Ecoregion)),1] +
  ranef(TP.lmer)$STAIDer[as.numeric(ordered(pred.data$STAIDer)),1]

## sum of squares -- model estimated seasonal means - site overall mean
GL.clean <- GLdata[(GLdata$State!="NEW YORK") & (!is.na(GLdata$Tpwu)),]
data.mean <- unlist(by(log(GL.clean$Tpwu),GL.clean$STAIDer, mean, na.rm=T))
diff.lmer <- pred.data$Pred.lmer -
  as.numeric(data.mean[as.numeric(ordered(pred.data$STAIDer))])
diff.lm <- pred.data$Pred.lm -
  as.numeric(data.mean[as.numeric(ordered(pred.data$STAIDer))])
SS.lmer <- by(diff.lmer^2, pred.data$STAIDer, sum)
SS.lm <- by(diff.lm^2, pred.data$STAIDer, sum)
boxplot(SS.lmer, SS.lm, names=c("Lmer", "Lm"))
mean(SS.lm > SS.lmer)

## The following lines generate the sample size figure.
## sample size plot

packages(reshape2)
packages(plyr)
gldata.m <- melt(GL.clean[,c("STAIDer","Ecoregion","Season","Year","Tpwu")],1:4,5)
smplsz <- dcast(gldata.m, STAIDer+Year+Season ~.)

smsz2 <- dcast(gldata.m, STAIDer+Year+Season ~., subset = .(STAIDer == "55-60"))

smplsz$Ecoregion <- substring(smplsz$STAIDer, 1, 2)
names(smplsz)[4] <- "sz"

smsz2$Ecoregion <- substring(smsz2$STAIDer, 1, 2)
names(smsz2)[4] <- "sz"

pred.lm <- predict(TP.lm, new=smplsz, se.fit=T, interval="prediction")
pred.lm <- data.frame(fit=pred.lm$fit, se=pred.lm$se.fit)

pred.lmer <- data.frame(fit=fixef(TP.lmer) +
        ranef(TP.lmer)$Season[as.numeric(ordered(smplsz$Season)),1] +
        ranef(TP.lmer)$Year[as.numeric(ordered(smplsz$Year)),1] +
        ranef(TP.lmer)$Ecoregion[as.numeric(ordered(smplsz$Ecoregion)),1] +
        ranef(TP.lmer)$STAIDer[as.numeric(ordered(smplsz$STAIDer)),1],
      se = sqrt((se.fixef(TP.lmer))^2 +
        (se.ranef(TP.lmer)$Season[as.numeric(ordered(smplsz$Season)),1])^2 +
        (se.ranef(TP.lmer)$Year[as.numeric(ordered(smplsz$Year)),1])^2 +
        (se.ranef(TP.lmer)$Ecoregion[as.numeric(ordered(smplsz$Ecoregion)),1])^2 +
        (se.ranef(TP.lmer)$STAIDer[as.numeric(ordered(smplsz$STAIDer)),1])^2))


par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.25,0.25,0), las=1, tck=0.01)
plot(xx <- jitter(smplsz$sz), pred.lm$fit.fit,
     ylim=range(pred.lm[,c(2,3)]), pch=16, cex=0.75,
     xlab="sample size", ylab="predicted means")
segments(x0=xx, x1=xx, y0=pred.lm$fit.lwr, y1=pred.lm$fit.upr)
abline(h=mean(log(GL.clean$Tpwu)))
plot(xx, pred.lmer$fit, ylim=range(pred.lm[,c(2,3)]), pch=16, cex=0.75,
     xlab="sample size", ylab="predicted means")
segments(x0=xx, x1=xx, y0=pred.lmer$fit-1.96*pred.lmer$se,
         y1=pred.lmer$fit+1.96*pred.lmer$se)
abline(h=mean(log(GL.clean$Tpwu)))
## site means:
GL.clean$SY <- paste(GL.clean$Year, GL.clean$Season)
TP.lmer2 <- lmer(log(Tpwu)~1+(1|SY),
                data=GL.clean, subset=STAIDer=="55-60")
TP.lm2 <- lm(log(Tpwu)~factor(SY)-1,
            data=GL.clean, subset=STAIDer=="55-60")
summary(TP.lmer2)
summary(TP.lm2)
clean.m <-  melt(GL.clean[,c("STAIDer","Ecoregion","SY","Tpwu")],1:3,4)
smsz <- dcast(clean.m, SY ~., subset=.(STAIDer=="55-60"))

tikz(file=paste(plotDIR, "smsize.tex", sep="/"), height=3, width=6.75, standAlone=F)
#pdf(file=paste(plotDIR, "smsize.pdf", sep="/"), height=3, width=6.75)
par(mfrow=c(1,2), mgp=c(1.75,0.125,0), mar=c(3,2,0.5,0.05), oma=c(1,2,1,2), tck=0.02)
lmCoef <- summary(TP.lm2)$coef
lower <- lmCoef[,1] - lmCoef[,2]
upper <- lmCoef[,1] + lmCoef[,2]
Ylim <- range(lower, upper)
size <- jitter(smsz[,2])
plot(size, rnorm(length(size)), ylim=Ylim, type="n", xlab="Sample size",
  ylab=" ", log="x")
abline(h=mean(lmCoef[,1]))
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, lmCoef[,1], pch=16, cex=0.5)
#title(main="No Pooling", cex=0.75)
mtext("Log Mean (MLE)", 2, line=2)

par(mgp=c(1.75,0.125,0), mar=c(3,0.05,0.5,2), tck=0.02)
lmer.mean <- fixef(TP.lmer2)[1] + ranef(TP.lmer2)[[1]][,1]
lmer.se <- sqrt(se.fixef(TP.lmer2)[1]^2 + se.ranef(TP.lmer2)[[1]][,1]^2)
lower <- lmer.mean-lmer.se
upper <- lmer.mean+lmer.se

plot(size, rnorm(length(size)), ylim=Ylim, type="n", xlab="Sample size",
  ylab=" ", log="x", axes=F, las=1)
abline(h=fixef(TP.lmer2)[[1]][1])
segments(x0=size, x1=size, y0=lower, y1=upper)
points(size, lmer.mean, pch=16,cex=0.5)
box()
axis(1, las=1)
axis(4)
mtext("Log Mean (Shrinkage)", side=4, line=2)
dev.off()

## Bayesian Estimator for single site
## 1. priors
##    beta0 -- fixef(TP.lmer2)[1], se.fixef(TP.lmer2)[1]
## sigma0 -- summary(TP.lmer2)@sigma
## Random effects:
## Groups    Name        Variance  Std.Dev.
## STAIDer   (Intercept) 0.4389319 0.66252  attr(VarCorr(TP.lmer2)$STAIDer, "stddev")
## Year      (Intercept) 0.0035391 0.05949  attr(VarCorr(TP.lmer2)$Year, "stddev")
## Ecoregion (Intercept) 0.4930635 0.70218  attr(VarCorr(TP.lmer2)$Ecoregion, "stddev")
## Season    (Intercept) 0.0598359 0.24461  attr(VarCorr(TP.lmer2)$Season, "stddev")
## Residual              0.6055830 0.77819  attr(VarCorr(TP.lmer2), "sc")

sequp_ch6 <- "
              data{
                int<lower=0> NOBS;
                int<lower=0> ncens;
                real y[NOBS];
                real ycens;
                real alpha1;
                real beta1;
                real alpha2;
                real beta2;
                real mub0;
                real<lower=0> sigmab0;
              }
              parameters{
                real mu;
                real mu0;
                real<lower=0> sigma1sq;
                real<lower=0> sigma2sq;
              }
              transformed parameters{
                real<lower=0> sigma1;
                real<lower=0> sigma2;
                sigma1 = sqrt(sigma1sq);
                sigma2 = sqrt(sigma2sq);
              }
              model{
                sigma1sq ~ inv_gamma(alpha1,beta1);
                sigma2sq ~ inv_gamma(alpha2,beta2);
                mu0 ~ normal(mub0, sigmab0);
                mu ~ normal(mu0, sigma2);
                target += normal_lpdf(y | mu, sigma1);
                target += ncens*normal_lcdf(ycens | mu, sigma1);
              }
              "
stan.in1 <- function(m=TP.lmer, y=GLdata$Tpwu[GLdata$State=="NEW YORK"],
                    y.censor=GLdata$TpwuCensored[GLdata$State=="NEW YORK"],
                    n.chains=nchains){
  n <- length(y)
  Ecoregion <- 11 ## eco83 ## 33 OH data points in the same region
  b0 <- fixef(m)[1]+ranef(m)$Ecoregion[11,1]
  seb0 <- sqrt(se.fixef(m)[1]^2+se.ranef(m)$Ecoregion[11,1]^2)
  sigma2.sq <- VarCorr(m)$STAIDer[1,1]+VarCorr(m)$Year[1,1] +
    VarCorr(m)$Season[1,1]
  n0 <- 30
  alpha2 <- n0/2
  beta2 <- (alpha2-1)*sigma2.sq
  sigma1.sq <- attr(VarCorr(m), "sc")^2
  alpha1 <- n0/2
  beta1 <- (alpha1-1)*sigma1.sq
  inits <- list()
  n.cens <- sum(y.censor=="<")
  data <- list(NOBS=n-n.cens, ncens=n.cens,
                    y=log(y)[y.censor!="<"],
                    mub0=b0, sigmab0=seb0, ycens=log(0.01),
                    alpha1=alpha1, beta1=beta1, alpha2=alpha2, beta2=beta2)
  for (i in 1:n.chains)
    inits [[i]] <- list(sigma1sq=runif(1), sigma2sq=runif(1), mu0=rnorm(1),
                           mu=rnorm(1))
  para <- c("mu","mu0","sigma1","sigma2")
  return(list(data=data, inits=inits, pars=para,
              n.chains=n.chains))
}

input.to.stan <- stan.in1()
tmp.m <- stan_model(model_code=sequp_ch6)
## initial iteration of 1 for model compilation
fit2keep <- sampling(tmp.m,
                     data=input.to.stan$data,
                     init=input.to.stan$inits,
                     pars=input.to.stan$pars,
                     iter=niters,thin=nthin,
                     chains=input.to.stan$n.chains)
print(fit2keep)

## trace plots
traceplot(fit2keep)

## extract output and save into rv object
simCoda <- rvsims(as.matrix(
    as.data.frame(rstan::extract(fit2keep, permuted=T))))
nms <- names(simCoda)

## predictive distribution (TP)
ny.pred <- rvnorm(1, mean=simCoda[1], sd=simCoda[3])

mub0 <- input.to.stan$data$mub0
seb0 <- input.to.stan$data$sigmab0
##postscript(file=paste(plotDIR, "NYpred.eps", sep="/"), height=2.5, width=3.25, horizontal=F)
## col=rgb() does not support postscript
pdf(file=paste(plotDIR, "NYpred.pdf", sep="/"), height=2.5, width=3.25)
par(mar=c(2.5,2.5,0.125,0.125), mgp=c(1.25,0.125,0), las=1, tck=0.01)
rvhist(exp(ny.pred), main="", xlab="TP (mg/L)", ylab="", xlim=c(0,0.22),ylim=c(0,75),
       col=grey(0.5), border=grey(0.85), density=-1)
rvhist(exp(simCoda[1]), add=T,
       col=rgb(0.75,0.75,0.75,0.5), border=grey(0.7), breaks=10)
abline(v=rvquantile(exp(simCoda[1]), prob=0.9))
##box()
segments(x0=exp(mub0-seb0), x1=exp(mub0+seb0), y0=70,y1=70)

points(x=exp(mub0), y=70)
points(x=exp(input.to.stan$data$y), y=rep(0, input.to.stan$data$NOBS), pch=16)
text(x=exp(input.to.stan$data$ycens), y=-1.5,  "$\\leftarrow$", pos=2, offset=0)
points(x=jitter(rep(exp(input.to.stan$data$ycens), input.to.stan$data$ncens)),
       y=seq(from=0,to=15,,input.to.stan$data$ncens), pch=1, cex=0.5)

dev.off()

cr <- 0.02413
Pr(simCoda[1]> log(cr))
Pr(ny.pred>log(cr))

## testing by year
tmp <- GLdata[GLdata$State=="NEW YORK",]
cr <- 0.02413
## 1. binom.test:
tapply(tmp$Tpwu, tmp$Year, FUN=function(x)binom.test(sum(x>cr),
                               length(x), 0.1, alternative="greater")$p.value)
tapply(tmp$Tpwu, tmp$Year, FUN=function(x)c(sum(x>cr), length(x)))

## counting censored data
tapply(tmp$TpwuCensored, tmp$Year, FUN=function(x)c(sum(x=="<"), length(x)))

## The Bayesian estimation process was programmed into two functions:
tmp.m <- stan_model(model_code=sequp_ch6)
fit.stan1 <- function(y=y.tmp, yC=y.tmpC, priorM=TP.lmer, n.chains=5,
                      n.iters=15000, n.keep=5000){
    input.to.stan <- stan.in1 (m=priorM, y=y, y.censor=yC, n.chains=1)
    fit2keep <- sampling(tmp.m,
                         data=input.to.stan$data,
                         init=input.to.stan$inits,
                         pars=input.to.stan$pars,
                         iter=niters,thin=nthin,
                         chains=input.to.stan$n.chains)
    return(fit2keep)
}

pred.dist <- function(stan.fit){
    simCoda <- rvsims(as.matrix(
        as.data.frame(rstan::extract(stan.fit, permuted=T))))
    ny.pred <- rvnorm(1, mean=simCoda[1], sd=simCoda[3])
    prob1 <- Pr(ny.pred > log(cr))
    prob2 <- Pr(simCoda[1] > log(cr))
    return(list(prob1, prob2, simCoda, ny.pred))
}
## We now run the model one year at a time:
## 1996
y.tmp <- GLdata$Tpwu[GLdata$Year==1996 & GLdata$State=="NEW YORK"]
y.tmpC <- GLdata$TpwuCensored[GLdata$Year==1996 & GLdata$State=="NEW YORK"]
n1996 <- paste(1996, "(", paste(length(y.tmp), sum(y.tmpC=="<"), sep="/"), ")", sep="")
y1996 <- y.tmp
y1996C <- sum(y.tmpC=="<")
fit1996 <- fit.stan1()
pred1996 <- pred.dist(fit1996)

## 1997
y.tmp <- GLdata$Tpwu[GLdata$Year==1997 & GLdata$State=="NEW YORK"]
y.tmpC <- GLdata$TpwuCensored[GLdata$Year==1997 & GLdata$State=="NEW YORK"]
n1997 <-  paste(1997, "(", paste(length(y.tmp), sum(y.tmpC=="<"), sep="/"), ")", sep="")
y1997 <- y.tmp
y1997C <- sum(y.tmpC=="<")
fit1997 <- fit.stan1()
pred1997 <- pred.dist(fit1997)

## 1998
y.tmp <- GLdata$Tpwu[GLdata$Year==1998 & GLdata$State=="NEW YORK"]
y.tmpC <- GLdata$TpwuCensored[GLdata$Year==1998 & GLdata$State=="NEW YORK"]
n1998 <-  paste(1998, "(", paste(length(y.tmp), sum(y.tmpC=="<"), sep="/"), ")", sep="")
y1998 <- y.tmp
y1998C <- sum(y.tmpC=="<")
fit1998 <- fit.stan1()
pred1998 <- pred.dist(fit1998)


tikz(file=paste(plotDIRch6, "NYpred2.tex", sep="/"), height=3.25, width=4,
     standAlone=F)
par(mar=c(2.5,2.5,0.125,0.125), mgp=c(1.25,0.125,0), las=1, tck=0.01)
rvhist(exp(ny.pred), main="", xlab="TP (mg/L)", ylab="", xlim=c(0,0.22),ylim=c(0,75),
       col=grey(0.5), border=grey(0.85), density=-1)
rvhist(exp(simCoda[1]), add=T,
       col=rgb(0.75,0.75,0.75,0.5), border=grey(0.7), breaks=10)
abline(v=rvquantile(exp(simCoda[1]), prob=0.9))
##box()
segments(x0=exp(mub0-seb0), x1=exp(mub0+seb0), y0=70,y1=70, lty=2)

points(x=exp(mub0), y=70, pch=2)
points(x=exp(input.to.stan$data$y), y=rep(0, input.to.stan$data$NOBS), pch=16)
text(x=exp(input.to.stan$data$ycens), y=-1.5,  "$\\leftarrow$", pos=2, offset=0)
points(x=jitter(rep(exp(input.to.stan$data$ycens), input.to.stan$data$ncens)),
       y=seq(from=0,to=15,,input.to.stan$data$ncens), pch=1, cex=0.5)

segments(x0=summary(exp(pred1996[[3]]))[1,5],
         x1=summary(exp(pred1996[[3]]))[1,9],
         y0=35, y1=35, col=grey(0.6))
segments(x0=summary(exp(pred1997[[3]]))[1,5],
         x1=summary(exp(pred1997[[3]]))[1,9],
         y0=40, y1=40, col=grey(0.6))
segments(x0=summary(exp(pred1998[[3]]))[1,5],
         x1=summary(exp(pred1998[[3]]))[1,9],
         y0=45, y1=45, col=grey(0.6))
points(x=c(summary(exp(pred1996[[3]]))[1,7],
           summary(exp(pred1997[[3]]))[1,7],
           summary(exp(pred1998[[3]]))[1,7]),
       y=c(35, 40, 45), col=grey(0.6))
points(x=c(y1996, y1997, y1998),
       y=c(rep(35-0.5,length(y1996)), rep(40-0.5,length(y1997)),
           rep(45-0.5,length(y1998))), col=grey(0.5), pch=15, cex=0.5)
points(x=jitter(rep(exp(input.to.stan$data$ycens), input.to.stan$data$ncens)),
       y=seq(from=0,to=20,,input.to.stan$data$ncens), pch=1, cex=0.5)

text(x=c(summary(exp(pred1996[[3]]))[1,9],
         summary(exp(pred1997[[3]]))[1,9],
         summary(exp(pred1998[[3]]))[1,9]),
     y=c(35,40,45), c(n1996, n1997,n1998), adj=c(0,0), cex=0.75)
dev.off()

