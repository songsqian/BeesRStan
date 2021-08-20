source("FrontMatter.R")

## Chapter 5
plotDIRch5 <- paste(plotDIR, "chapter5", "figures", sep="/")
## simulation
packages(rv)
packages(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

## Liverpool moth
moth <- read.table(paste(dataDIR, "moths.txt", sep="/"), header=T)

## The Stan model
moth_var <- "
data{
  int n; //sample size
  int y[n]; //# of removed
  int NN[n]; //# of placed
  int nM;
  real D[n]; //distance
  int<lower=1> morph[n]; //1=dark, 2=light
}
parameters{
  real beta0;
  real beta1[nM];
  real<lower=0, upper=2> sigma1;
  real beta2;
  real beta3[nM];
  real<lower=0, upper=2> sigma3;
}
model{
  real xb[n];
  beta1 ~ normal(0, sigma1);
  beta3 ~ normal(0, sigma3);
  for (i in 1:n){
    xb[i] = beta0+beta1[morph[i]]+beta2*D[i]+beta3[morph[i]]*D[i];
  }
  y ~ binomial_logit(NN, xb);
}
generated quantities{
  real alpha0;
  real alpha1[nM];
  real alpha2;
  real alpha3[nM];
  real Int[n];
  real sM;
  real sD;
  real sInt;

  alpha0 = beta0 + mean(beta1[]);
  alpha2 = beta2 + mean(beta3[]);
  for (i in 1:nM){
    alpha1[i] = beta1[i]-mean(beta1[]);
    alpha3[i] = beta3[i]-mean(beta3[]);
  }
  for (i in 1:n)
    Int[i] = alpha3[morph[i]]*D[i];
  sM = sd(alpha1[]);
  sD = fabs(alpha2)*sd(D[]);
  sInt = sd(Int[]);
}
"
moth_var2 <- "
data{
  int n; //sample size
  int y[n]; //# of removed
  int NN[n]; //# of placed
  int nM;
  real D[n]; //distance
  int<lower=1> morph[n]; //1=dark, 2=light
}
parameters{
  real beta0;
  real beta1[nM];
  real<lower=0, upper=2> sigma1;
  real beta2;
  real beta3[nM];
  real<lower=0, upper=2> sigma3;
}
model{
  real xb[n];
  beta1 ~ normal(0, sigma1);
  beta3 ~ normal(0, sigma3);
  for (i in 1:n){
    xb[i] = beta0+beta1[morph[i]]+beta2*D[i]+beta3[morph[i]]*D[i];
  }
  y ~ binomial_logit(NN, xb);
}
"

Input <- function(infile=moth, n.chains=nchains, vc=T){
    n <- dim(infile)[1]
    y <- infile$removed
    N <- infile$placed
    morph <- as.numeric(ordered(infile$morph))
    nM <- max(morph)
    Dist <- infile$distance
    data <- list(y=y, n=n, NN=N, nM=nM, morph=morph, D=Dist)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(beta0=rnorm(1), beta1=rnorm(nM),
                           sigma1=runif(1), beta2=rnorm(1),
                           beta3=rnorm(nM), sigma3=runif(1))
    if (vc)
        paras <- c("alpha0", "alpha1", "alpha2", "alpha3", "sM", "sD", "sInt",
               "beta0", "beta1", "sigma1", "beta2", "beta3", "sigma3")
    else paras <- c("beta0", "beta1", "sigma1", "beta2", "beta3", "sigma3")
    return(list(data=data, inits=inits, paras=paras, nchains=n.chains))
}

input.to.stan <- Input()
fit <- stan_model(model_code=moth_var)
fit2keep <- sampling(fit, data = input.to.stan$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan$inits,
                     pars=input.to.stan$paras, thin=nthin)
print(fit2keep)
save(fit2keep, file="mothStan1.RData")
mothStanrun1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep, permuted=T))))
mothStanbeta0 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep, par="beta0", permuted=T))))
mothStanbeta1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep, par="beta1", permuted=T))))
mothStanbeta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep, par="beta2", permuted=T))))
mothStanbeta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep, par="beta3", permuted=T))))

input.to.stan <- Input(vc=F)
fit <- stan_model(model_code=moth_var2)
fit2keep2 <- sampling(fit, data = input.to.stan$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan$inits,
                     pars=input.to.stan$paras, thin=nthin)
print(fit2keep2)

## using rv to calculate variance components
alpha0 <- mothStanbeta0 + mean(mothStanbeta1)
alpha2 <- mothStanbeta2 + mean(mothStanbeta3)
alpha1 <- mothStanbeta1[1]-mean(mothStanbeta1)
alpha1 <- c(alpha1, mothStanbeta1[2]-mean(mothStanbeta1))
alpha3 <- mothStanbeta3[1]-mean(mothStanbeta3)
alpha3 <- c(alpha3, mothStanbeta3[2]-mean(mothStanbeta3))

input.to.stan$data$morph
input.to.stan$data$D
Int <- alpha3[input.to.stan$data$morph]*input.to.stan$data$D;
sM <- sd.rv(alpha1);
sD <- abs(alpha2)*sd(input.to.stan$data$D);
sInt <- sd.rv(Int);
vcomp <- c(sM, sD, sInt)
names(vcomp) <- c("Morph", "Distance", "Interaction")

tikz(file=paste(plotDIRch5, "mothvcomp.tex", sep="/"),
     height=3, width=3.5, standAlone=F)
par(mar=c(3, 5, 3, 1), las=1, mgp=c(1.25, 0.125,0), tck=0.01)
mlplot(vcomp, xlab="standard deviation", ylab="")
dev.off()

## slopes
slopes <- alpha2+alpha3
names(slopes) <- c("dark","light")
tikz(file=paste(plotDIRch5, "mothslopes.tex", sep="/"),
     height=3, width=3.5, standAlone=F)
par(mar=c(3, 5, 3, 1), las=1, mgp=c(1.25, 0.125,0), tck=0.01)
mlplot(slopes, xlab="slope", ylab="")
abline(v=0, col="grey")
dev.off()

Int <- alpha0+alpha1
 names(Int) <- c("dark","light")
 Int

## Using GLM
mothGLM <- glm(cbind(removed, placed-removed)~factor(morph)*distance,
               data=moth, family=binomial)
display(mothGLM, 4)
anova(mothGLM)

#############################################
#### seedling recruitment data from Shen ####
#############################################

transect <- read.csv(paste(dataDIR, "transect2.csv", sep="/"), header=T)
transect.desc <- read.csv(paste(dataDIR, "transectX.csv", sep="/"), header=T)
species.code3 <- read.csv(paste(dataDIR, "species3.csv", sep="/"), header=T)

seedling <- transect[order(transect$Plot),]
temp <- transect.desc[unique(seedling$Plot),]
seedling$Position <- rep(temp$Position, table(seedling$Plot))
seedling$Gap <- rep(temp$Gap, table(seedling$Plot))
seedling$Total.C <- rep(temp$Total.C, table(seedling$Plot))
seedling$ABH <- rep(temp$ABH, table(seedling$Plot))
seedling$Gap <- ifelse(seedling$Gap>1, 1, seedling$Gap)
seedling$gap.center <- seedling$Gap-0.5

seedlingGLM <- glm(Numb ~ factor(type)*logit(Gap)+
                       factor(type)*Total.C+
                       factor(Position)*logit(Gap)+
                       factor(Position)*Total.C+
                       factor(Position)*factor(type),
                   family=poisson, data=seedling)
seedlingGLM <- glm(Numb ~ factor(type)*logit(Gap)+
                        factor(type)*Total.C+
                        factor(Position),
                   family=poisson, data=seedling)
seedlingGLM2 <- glm(Numb ~ factor(type)*logit(Gap)+
                        factor(type)*Total.C,
                   family=poisson, data=seedling)
anova(seedlingGLM)
anova(seedlingGLM2, seedlingGLM)
summary(seedlingGLM)
summary(seedlingGLM2)

matrix(table(paste(seedling$type, seedling$Position, sep="_")), 5,5)
matrix(names(table(paste(seedling$type, seedling$Position, sep="_"))), 5,5)
## Stan models
seedStan <- "
data{
  int n; //sample size
  int nt; //number of tree types
  int np; //number of positions
  int type[n];
  int position[n];
  int y[n];
  real soilC[n];
  real gap[n]; //logit of gap
}
parameters{
  real beta0;
  real beta1[nt];
  real beta2[np];
  real beta3;
  real delta3[nt];
  real beta4;
  real delta4[np];
  real<lower=0> sigma[4];
}
model{
  real log_lambda[n];
  sigma ~ normal(0, 2);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, sigma[1]);
  beta2 ~ normal(0, sigma[2]);
  beta3 ~ normal(0, 1);
  beta4 ~ normal(0, 1);
  delta3 ~ normal(0, sigma[3]);
  delta4 ~ normal(0, sigma[4]);
  for (i in 1:n){
    log_lambda[i] = beta0 + beta1[type[i]] + beta2[position[i]] +
      (beta3+delta3[type[i]])*soilC[i] + (beta4+delta4[type[i]])*gap[i];
  }
  y ~ poisson_log(log_lambda);
}
"

## Stan models with type:position interaction
seedStan2 <- "
data{
  int n; //sample size
  int nt; //number of tree types
  int ntp; //number of type-position combinations
  int type[n];
  int typ_pos[n];
  int y[n];
  real soilC[n];
  real gap[n]; //logit of gap
}
parameters{
  real beta0;
  real beta1[ntp];
  real beta2;
  real delta2[nt];
  real beta3;
  real delta3[nt];
  real<lower=0> sigma[3];
}
model{
  real log_lambda[n];
  sigma ~ normal(0, 2);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, sigma[1]);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  delta2 ~ normal(0, sigma[2]);
  delta3 ~ normal(0, sigma[3]);
  for (i in 1:n){
    log_lambda[i] = beta0 + beta1[typ_pos[i]] +
      (beta2+delta2[type[i]])*soilC[i] + (beta3+delta3[type[i]])*gap[i];
  }
  y ~ poisson_log(log_lambda);
}
"
inputSeed <- function(indata=seedling, n.chains=nchains){
    temp <- !is.na(indata$Numb)
    indata <- indata[temp,]
    n <- dim(indata)[1]
    y <- indata$Numb
    type <- as.numeric(ordered(indata$type))
    position <- as.numeric(ordered(indata$Position))
    np <- max(position)
    nt <- max(type)
    data_in <- list(n=n, nt=nt, np=np, type=type, position=position,
                    y=y, soilC=indata$Total.C, gap=logit(indata$Gap))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(beta0=rnorm(1), beta1=rnorm(nt),
                           beta2=rnorm(np), beta3=rnorm(1),
                           delta3=rnorm(nt), beta4=rnorm(1),
                           delta4=rnorm(nt), sigma=runif(4))
    para=c("beta0","beta1","beta2","beta3","delta3","beta4",
           "delta4","sigma")
    return(list(data=data_in, inits=inits, paras=para, nchains=n.chains))
}

input.to.stan <- inputSeed()
fit1 <- stan_model(model_code=seedStan)
fit2keep1 <- sampling(fit1, data = input.to.stan$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan$inits,
                     pars=input.to.stan$paras, thin=nthin,
                     control=list(adapt_delta=0.95))
print(fit2keep1)

inputSeed2 <- function(indata=seedling, n.chains=nchains){
    temp <- !is.na(indata$Numb)
    indata <- indata[temp,]
    n <- dim(indata)[1]
    y <- indata$Numb
    type <- as.numeric(ordered(indata$type))
    type_position <- as.numeric(ordered(paste(indata$type, indata$Position)))
    ntp <- max(type_position)
    nt <- max(type)
    data_in <- list(n=n, nt=nt, ntp=ntp, type=type, typ_pos=type_position,
                    y=y, soilC=indata$Total.C, gap=logit(indata$Gap))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(beta0=rnorm(1), beta1=rnorm(ntp),
                           beta2=rnorm(1),
                           delta2=rnorm(nt), beta3=rnorm(1),
                           delta3=rnorm(nt), sigma=runif(3))
    para=c("beta0","beta1","beta2","delta2","beta3",
           "delta3","sigma")
    return(list(data=data_in, inits=inits, paras=para, nchains=n.chains))
}

input.to.stan2 <- inputSeed2()
fit2 <- stan_model(model_code=seedStan2)
fit2keep2 <- sampling(fit2, data = input.to.stan2$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan2$inits,
                     pars=input.to.stan2$paras, thin=nthin,
                     control=list(adapt_delta=0.95))
print(fit2keep2)

save(fit2keep1, fit2keep2, file="seedfits.RData")
## using rv to calculate variance components (without type-position interaction)
seedStanbeta0 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="beta0", permuted=T))))
seedStanbeta1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="beta1", permuted=T))))
seedStanbeta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="beta2", permuted=T))))
seedStanbeta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="beta3", permuted=T))))
seedStandelta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="delta3", permuted=T))))
seedStanbeta4 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="beta4", permuted=T))))
seedStandelta4 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="delta4", permuted=T))))
seedStansigma <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep1, par="sigma", permuted=T))))

alpha0 <- seedStanbeta0 + mean(seedStanbeta1) + mean(seedStanbeta2)
alpha1 <- seedStanbeta1 - mean(seedStanbeta1)
alpha2 <- seedStanbeta2 - mean(seedStanbeta2)
alpha3 <- seedStanbeta3 + mean(seedStandelta3)
del3 <- seedStandelta3 - mean(seedStandelta3)
alpha4 <- seedStanbeta4 + mean(seedStandelta4)
del4 <- seedStandelta4 - mean(seedStandelta4)

input.to.stan$data$gap
input.to.stan$data$soilC
Int3 <- del3[input.to.stan$data$type]*input.to.stan$data$soilC;
Int4 <- del4[input.to.stan$data$type]*input.to.stan$data$gap;
sType <- sd.rv(alpha1)
sPos <- sd.rv(alpha2)
sSoilC <- abs(alpha3)*sd(input.to.stan$data$soilC)
sGap <- abs(alpha4)*sd(input.to.stan$data$gap)
sInt3 <- sd.rv(Int3)
sInt4 <- sd.rv(Int4)

vcomp1 <- c(sType, sPos, sSoilC, sGap, sInt3, sInt4)
names(vcomp1) <- c("Type", "Position", "Soil C",
                   "Gap", "Type:Soil C", "Type:Gap")
tikz(file=paste(plotDIRch5, "seedVC1.tex", sep="/"),
     height=2.75, width=3.75, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(vcomp1, xlab="standard deviation", xlim=c(0,0.8), cex=0.75)
dev.off()
## using rv to calculate variance components (with type-position interaction)
seedStan2beta0 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="beta0", permuted=T))))
seedStan2beta1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="beta1", permuted=T))))
seedStan2beta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="beta2", permuted=T))))
seedStan2beta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="beta3", permuted=T))))
seedStan2delta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="delta2", permuted=T))))
seedStan2delta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="delta3", permuted=T))))
seedStan2sigma <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep2, par="sigma", permuted=T))))

alpha0 <- seedStan2beta0 + mean(seedStan2beta1)
alpha1 <- seedStan2beta1 - mean(seedStan2beta1)
alpha2 <- seedStan2beta2 + mean(seedStan2delta2)
alpha3 <- seedStan2beta3 + mean(seedStan2delta3)
del2 <- seedStan2delta2 - mean(seedStan2delta2)
del3 <- seedStan2delta3 - mean(seedStan2delta3)

input.to.stan2$data$gap
input.to.stan2$data$soilC
Int2 <- del2[input.to.stan2$data$type]*input.to.stan2$data$soilC;
Int3 <- del3[input.to.stan2$data$type]*input.to.stan2$data$gap;

sType_Pos <- sd.rv(alpha1)
sSoilC <- abs(alpha2)*sd(input.to.stan2$data$soilC)
sGap <- abs(alpha3)*sd(input.to.stan2$data$gap)
sInt2 <- sd.rv(Int2)
sInt3 <- sd.rv(Int3)

vcomp2 <- c(sType_Pos, sSoilC, sGap, sInt2, sInt3)
names(vcomp2) <- c("Type:Position", "Soil C",
                   "Gap", "Type:Soil C", "Type:Gap")
tikz(file=paste(plotDIRch5, "seedVC2.tex", sep="/"),
     height=2.75, width=3.75, standAlone=F)
par(mar=c(3,3,1,0.25), mgp=c(1.25,0.125,0), las=1,tck=0.01)
mlplot(vcomp2, xlab="standard deviation", xlim=c(0,0.5), cex=0.75)
dev.off()

alp1 <- rvmatrix(alpha1, 5)
row.names(alp1) <- c("Pioneer","Dominant","Companion","Evergreen","Tolerant")
colnames(alp1) <- paste("Pos", 1:5)

par(mfrow=c(1,5), oma=c(3, 3, 1, 1), mar=c(3, 0, 0, 0), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1[1,], xlab=row.names(alp1)[1], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
axis(2, at=1:5, labels=colnames(alp1))
mlplot(alp1[2,], xlab=row.names(alp1)[2], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[3,], xlab=row.names(alp1)[3], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[4,], xlab=row.names(alp1)[4], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[5,], xlab=row.names(alp1)[5], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)

tikz(file=paste(plotDIRch5, "seedInt1.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(t(alp1), xlab="", cex=0.5)
dev.off()
tikz(file=paste(plotDIRch5, "seedInt2.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1, xlab="", cex=0.5)
dev.off()

##### CAR prior
## Adjacency matrix W: w[i,i]=0, w[i,j]=1 if i is adjacent to j,
##                               w[i,j]=0 otherwise
#####

### Stan model (using multivariate normal CAR prior)
seedStanCar1 <- "
data{
  int n; //sample size
  int nt; //number of tree types
  int ntp; //number of type-position combinations
  int N; //number of plots
  int type[n];
  int typ_pos[n];
  int plot[n];
  int y[n];
  real soilC[n];
  real gap[n]; //logit of gap
  matrix<lower = 0, upper = 1>[N, N] W;
}
transformed data{
  vector[N] zeros;
  matrix<lower = 0>[N, N] D;
  {
    vector[N] W_rowsums;
    for (i in 1:N) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, N);
}
parameters{
  real beta0;
  real beta1[ntp];
  real beta2;
  real delta2[nt];
  real beta3;
  real delta3[nt];
  real<lower=0> sigma[3];
  vector[N] phi;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
}
model{
  real log_lambda[n];
  sigma ~ normal(0, 2);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, sigma[1]);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  delta2 ~ normal(0, sigma[2]);
  delta3 ~ normal(0, sigma[3]);
  phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  tau ~ gamma(2, 2);
  for (i in 1:n){
    log_lambda[i] = beta0 + beta1[typ_pos[i]] +
      (beta2+delta2[type[i]])*soilC[i] + (beta3+delta3[type[i]])*gap[i] +
       phi[plot[i]];
  }
  y ~ poisson_log(log_lambda);
}
"

inputSeedCar1 <- function(indata=seedling, n.chains=nchains){
    temp <- !is.na(indata$Numb)
    indata <- indata[temp,]
    n <- dim(indata)[1]
    y <- indata$Numb
    type <- as.numeric(ordered(indata$type))
    type_position <- as.numeric(ordered(paste(indata$type, indata$Position)))
    ntp <- max(type_position)
    nt <- max(type)
    Plot <- as.numeric(ordered(indata$Plot))
    n.plot <- max(Plot)
    W <- matrix(0, n.plot, n.plot)
    for (i in 1:(n.plot-1))
        W[i, i+1] <- 1
    for (j in 2:n.plot)
        W[j, j-1] <- 1

    data_in <- list(N=n.plot, n=n, nt=nt, ntp=ntp, type=type,
                    typ_pos=type_position, W=W, y=y, plot=Plot,
                    soilC=indata$Total.C, gap=logit(indata$Gap))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(beta0=rnorm(1), beta1=rnorm(ntp),
                           beta2=rnorm(1),
                           delta2=rnorm(nt), beta3=rnorm(1),
                           delta3=rnorm(nt), sigma=runif(3),
                           phi=rnorm(n.plot, 0, 0.1),
                           tau=runif(1), alpha=runif(1))
    para=c("beta0","beta1","beta2","delta2","beta3",
           "delta3","sigma", "phi", "alpha", "tau")
    return(list(data=data_in, inits=inits, paras=para, nchains=n.chains))
}

input.to.stan3 <- inputSeedCar1()
fit3 <- stan_model(model_code=seedStanCar1)
fit2keep3 <- sampling(fit3, data = input.to.stan3$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan3$inits,
                     pars=input.to.stan3$paras, thin=nthin)#,
##                     control=list(adapt_delta=0.95))
print(fit2keep3)
save(fit2keep3, file="multinormCAR.RData")

### Stan model (using sparse CAR prior)
seedStanCar2 <- "
functions {
  /** from https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  * Return the log probability of a proper conditional autoregressive (CAR) prior
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha,
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;

      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }

      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data{
  int n; //sample size
  int nt; //number of tree types
  int ntp; //number of type-position combinations
  int N; //number of plots
  int type[n];
  int typ_pos[n];
  int plot[n];
  int y[n];
  int W_n;
  real soilC[n];
  real gap[n]; //logit of gap
  matrix<lower = 0, upper = 1>[N, N] W;
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] lambda;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(W[i]);
  {
    vector[N] invsqrtD;
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters{
  real beta0;
  real beta1[ntp];
  real beta2;
  real delta2[nt];
  real beta3;
  real delta3[nt];
  real<lower=0> sigma[3];
  vector[N] phi;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;
}
model{
  real log_lambda[n];
  sigma ~ normal(0, 2);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, sigma[1]);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  delta2 ~ normal(0, sigma[2]);
  delta3 ~ normal(0, sigma[3]);
  phi ~ sparse_car(tau, alpha, W_sparse, D_sparse, lambda, N, W_n);
  tau ~ gamma(2, 2);
  for (i in 1:n){
    log_lambda[i] = beta0 + beta1[typ_pos[i]] +
      (beta2+delta2[type[i]])*soilC[i] + (beta3+delta3[type[i]])*gap[i] +
       phi[plot[i]];
  }
  y ~ poisson_log(log_lambda);
}
"

inputSeedCar2 <- function(indata=seedling, n.chains=nchains){
    temp <- !is.na(indata$Numb)
    indata <- indata[temp,]
    n <- dim(indata)[1]
    y <- indata$Numb
    type <- as.numeric(ordered(indata$type))
    type_position <- as.numeric(ordered(paste(indata$type, indata$Position)))
    ntp <- max(type_position)
    nt <- max(type)
    Plot <- as.numeric(ordered(indata$Plot))
    n.plot <- max(Plot)
    W <- matrix(0, n.plot, n.plot)
    for (i in 1:(n.plot-1))
        W[i, i+1] <- 1
    for (j in 2:n.plot)
        W[j, j-1] <- 1

    data_in <- list(N=n.plot, n=n, nt=nt, ntp=ntp, type=type,
                    typ_pos=type_position, W=W, W_n=sum(W)/2,
                    y=y, plot=Plot, soilC=indata$Total.C,
                    gap=logit(indata$Gap))
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(beta0=rnorm(1), beta1=rnorm(ntp),
                           beta2=rnorm(1),
                           delta2=rnorm(nt), beta3=rnorm(1),
                           delta3=rnorm(nt), sigma=runif(3),
                           phi=rnorm(n.plot, 0, 0.1),
                           tau=runif(1), alpha=runif(1))
    para=c("beta0","beta1","beta2","delta2","beta3",
           "delta3","sigma", "phi", "alpha", "tau")
    return(list(data=data_in, inits=inits, paras=para, nchains=n.chains))
}

input.to.stan4 <- inputSeedCar2()
fit4 <- stan_model(model_code=seedStanCar2)
fit2keep4 <- sampling(fit4, data = input.to.stan4$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan4$inits,
                     pars=input.to.stan4$paras, thin=nthin,
                     control=list(adapt_delta=0.95))
print(fit2keep4)

## using rv to calculate variance components
## (multivariate normal CAR with type-position interaction)
seedStan3beta0 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="beta0", permuted=T))))
seedStan3beta1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="beta1", permuted=T))))
seedStan3beta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="beta2", permuted=T))))
seedStan3beta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="beta3", permuted=T))))
seedStan3delta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="delta2", permuted=T))))
seedStan3delta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="delta3", permuted=T))))
seedStan3sigma <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="sigma", permuted=T))))
seedStan3phi <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep3, par="phi", permuted=T))))
alpha0 <- seedStan3beta0 + mean(seedStan3beta1)
alpha1 <- seedStan3beta1 - mean(seedStan3beta1)
alpha2 <- seedStan3beta2 + mean(seedStan3delta2)
alpha3 <- seedStan3beta3 + mean(seedStan3delta3)
del2 <- seedStan3delta2 - mean(seedStan3delta2)
del3 <- seedStan3delta3 - mean(seedStan3delta3)

input.to.stan3$data$gap
input.to.stan3$data$soilC
Int2 <- del2[input.to.stan3$data$type]*input.to.stan3$data$soilC;
Int3 <- del3[input.to.stan3$data$type]*input.to.stan3$data$gap;

sType_Pos <- sd.rv(alpha1)
sSoilC <- abs(alpha2)*sd(input.to.stan4$data$soilC)
sGap <- abs(alpha3)*sd(input.to.stan4$data$gap)
sInt2 <- sd.rv(Int2)
sInt3 <- sd.rv(Int3)
sPhi3 <- sd.rv(seedStan4phi)

vcomp3 <- c(sType_Pos, sSoilC, sGap, sInt2, sInt3, sPhi3)
names(vcomp3) <- c("Type:Position", "Soil C",
                   "Gap", "Type:Soil C", "Type:Gap", "CAR")
tikz(file=paste(plotDIRch5, "seedVC3SparseslowCAR.tex", sep="/"),
     height=2.75, width=3.75, standAlone=F)
par(mar=c(3,3,1,0.25), mgp=c(1.25,0.125,0), las=1,tck=0.01)
mlplot(vcomp3, xlab="standard deviation", xlim=c(0,0.5), cex=0.75)
dev.off()

alp1 <- rvmatrix(alpha1, 5)
row.names(alp1) <- c("Pioneer","Dominant","Companion","Evergreen","Tolerant")
colnames(alp1) <- paste("Pos", 1:5)

par(mfrow=c(1,5), oma=c(3, 3, 1, 1), mar=c(3, 0, 0, 0), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1[1,], xlab=row.names(alp1)[1], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
axis(2, at=1:5, labels=colnames(alp1))
mlplot(alp1[2,], xlab=row.names(alp1)[2], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[3,], xlab=row.names(alp1)[3], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[4,], xlab=row.names(alp1)[4], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[5,], xlab=row.names(alp1)[5], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)

tikz(file=paste(plotDIRch5, "seedInt31slowCARS.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(t(alp1), xlab="", cex=0.5)
dev.off()
tikz(file=paste(plotDIRch5, "seedInt32slowCARS.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1, xlab="", cex=0.5)
dev.off()

## using rv to calculate variance components
## (Sparse CAR with type-position interaction)
seedStan4beta0 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="beta0", permuted=T))))
seedStan4beta1 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="beta1", permuted=T))))
seedStan4beta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="beta2", permuted=T))))
seedStan4beta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="beta3", permuted=T))))
seedStan4delta2 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="delta2", permuted=T))))
seedStan4delta3 <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="delta3", permuted=T))))
seedStan4sigma <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="sigma", permuted=T))))
seedStan4phi <- rvsims(as.matrix(as.data.frame(rstan::extract(
                          fit2keep4, par="phi", permuted=T))))
alpha0 <- seedStan4beta0 + mean(seedStan4beta1)
alpha1 <- seedStan4beta1 - mean(seedStan4beta1)
alpha2 <- seedStan4beta2 + mean(seedStan4delta2)
alpha3 <- seedStan4beta3 + mean(seedStan4delta3)
del2 <- seedStan4delta2 - mean(seedStan4delta2)
del3 <- seedStan4delta3 - mean(seedStan4delta3)

input.to.stan4$data$gap
input.to.stan4$data$soilC
Int2 <- del2[input.to.stan4$data$type]*input.to.stan4$data$soilC;
Int3 <- del3[input.to.stan4$data$type]*input.to.stan4$data$gap;

sType_Pos <- sd.rv(alpha1)
sSoilC <- abs(alpha2)*sd(input.to.stan4$data$soilC)
sGap <- abs(alpha3)*sd(input.to.stan4$data$gap)
sInt2 <- sd.rv(Int2)
sInt3 <- sd.rv(Int3)
sPhi3 <- sd.rv(seedStan4phi)

vcomp4 <- c(sType_Pos, sSoilC, sGap, sInt2, sInt3, sPhi3)
names(vcomp4) <- c("Type:Position", "Soil C",
                   "Gap", "Type:Soil C", "Type:Gap", "CAR")
tikz(file=paste(plotDIRch5, "seedVC3fastCAR.tex", sep="/"),
     height=2.75, width=3.75, standAlone=F)
par(mar=c(3,3,1,0.25), mgp=c(1.25,0.125,0), las=1,tck=0.01)
mlplot(vcomp3, xlab="standard deviation", xlim=c(0,0.5), cex=0.75)
dev.off()

alp1 <- rvmatrix(alpha1, 5)
row.names(alp1) <- c("Pioneer","Dominant","Companion","Evergreen","Tolerant")
colnames(alp1) <- paste("Pos", 1:5)

par(mfrow=c(1,5), oma=c(3, 3, 1, 1), mar=c(3, 0, 0, 0), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1[1,], xlab=row.names(alp1)[1], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
axis(2, at=1:5, labels=colnames(alp1))
mlplot(alp1[2,], xlab=row.names(alp1)[2], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[3,], xlab=row.names(alp1)[3], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[4,], xlab=row.names(alp1)[4], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)
mlplot(alp1[5,], xlab=row.names(alp1)[5], xlim=c(-0.5,1), ylim=c(1,5), axes=F)
axis(1)

tikz(file=paste(plotDIRch5, "seedInt31fastCARS.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(t(alp1), xlab="", cex=0.5)
dev.off()
tikz(file=paste(plotDIRch5, "seedInt32fastCARS.tex", sep="/"),
     height=4, width=2.5, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25,0.125,0),
    las=1,tck=0.01)
mlplot(alp1, xlab="", cex=0.5)
dev.off()
