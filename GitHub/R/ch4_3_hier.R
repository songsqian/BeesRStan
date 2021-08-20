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

neuse <- read.csv(paste(dataDIR, "NeuseEstuary.csv", sep="/"))
neuse[neuse<0] <- NA

neuse$RDate <- as.Date(as.character(neuse$DATE), format="%m/%d/%Y")
neuse$Year <- as.numeric(format(neuse$RDate, "%Y"))

## Sequential updating (Neuse River Estuary example)

sequpdate <- "
  data {
    int<lower=0> N; // sample size (training data, non-zero)
    real y[N]; // response data
    real mu0;
    real<lower=0> n0;
    real<lower=0> alpha;
    real<lower=0> beta;
  }
  parameters {
    real<lower=0> sigma1;
    real<lower=0> lambda;
    real gmu;
    real mu;
  }
  model {
    lambda ~ gamma(alpha, beta);
    mu ~ normal(mu0, 1/sqrt(lambda*n0));
    gmu ~ normal (mu, 1/sqrt(lambda));
    for (i in 1:N){
      y[i] ~ normal(gmu, sigma1);
    }
  }
  generated quantities {
    real<lower=0> sigma2;
    sigma2 = 1/sqrt(lambda);
  }
  "

## Hierarchical model
BHMNormal <- "
  data {
    int<lower=0> N; // sample size (training data, non-zero)
    int<lower=0> G; // number of groups
    real y[N]; // response data
    int<lower=0> gr[N];
  }
  parameters {
    real<lower=0> sigma;
    real<lower=0> hsig;
    real hmu;
    real gmu[G];
  }
  model {
    gmu ~ normal(hmu, hsig);
    for (i in 1:N){
      y[i] ~ normal(gmu[gr[i]], sigma);
    }
  }
  "
## Fitting the hiearchical model with all years
inputHIER <- function(dat=neuse, ycol="SURFCHLA", grp="Year", n.chains=nchains){
    temp <- !is.na(dat[,ycol]) & dat[,ycol]>0
    dat=dat[temp,]
    N <- dim(dat)[1]
    gr <- as.numeric(ordered(dat[,grp]))
    G <- max(gr)
    data <- list(N=N, G=G, y=log(dat[,ycol]), gr=gr)
    inits <- list()
    for (i in 1:n.chains)
        inits[[i]] <- list(sigma=runif(1), hsig=runif(1),
                           hmu=rnorm(1), gmu=rnorm(G))
    para <- c("sigma", "hsig", "hmu", "gmu")
    return(list(data=data, init=inits, paras=para))
}

input.to.stanAll <- inputHIER()
fit2 <- stan_model(model_code=BHMNormal)
fit2keep <- sampling(fit2, data = input.to.stanAll$data,
                     iter = niters, chains = nchains,
                     init=input.to.stanAll$init,
                     pars=input.to.stanAll$para, thin=nthin)

fit2coefAll <- as.data.frame(rstan::extract(fit2keep, permuted=T))
fit2coefAll_rv <- rvsims(as.matrix(fit2coefAll))

BHM_gmu <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep,pars="gmu", permuted=T))))
BHM_hmu <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep,pars="hmu", permuted=T))))
BHM_hsig <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep,pars="hsig", permuted=T))))
BHM_sigma <- rvsims(as.matrix(as.data.frame(
    rstan::extract(fit2keep,pars="sigma", permuted=T))))

## method of moments for estimating prior parameters
prior_pars <- function(mean_mu, v_mu, mean_lambda, v_lambda){
    return(list(alpha=mean_lambda^2/v_lambda,
                beta=mean_lambda/v_lambda,
                mu0=mean_mu,
                n0=(mean_lambda/v_lambda)/(((mean_lambda^2/v_lambda)-1)*v_mu) ))
}

prr_bhm <- prior_pars(mean(BHM_hmu[[1]]), var(BHM_hmu[[1]]),
                      mean(1/BHM_hsig[[1]]^2), sd((1/BHM_hsig[[1]])^2)^2)

## input for sequential updating model
inputSEQ <- function(dat=datfr, ycol, gr, n.chains=nchains,
                     mu0=0,n0=1, alpha=0.1,beta=0.1,
                     sub=dat$gr=="GR"){
      temp <- !is.na(dat[,ycol])
      dat=dat[temp&sub,]
      gr <- as.numeric(ordered(dat[,gr]))
      G <- max(gr)
      N <- dim(dat)[1]
      data <- list(N=N, G=G, y=log(dat[,ycol]), gr=gr,
                   mu0=mu0, n0=n0, alpha=alpha, beta=beta)
      inits <- list()
      for (i in 1:n.chains)
	  inits[[i]] <- list(hmu=rnorm(1), hsig=runif(1),
			     gmu=rnorm(G), sigma=runif(1))
      para <- c("mu", "gmu", "sigma1", "lambda")
      return (list(data=data, init=inits, paras=para))
}

init_prior <- list(mu0=1, n0=10, alpha=1, beta=1)

input.to.stan <- inputSEQ(dat=neuse, ycol="SURFCHLA", gr="Year",
                          mu0=init_prior$mu0,n0=init_prior$n0,
                          alpha=init_prior$alpha, beta=init_prior$beta,
                          sub=neuse$Year==1992)
fit1 <- stan_model(model_code = sequpdate)
fit2keep <- sampling(fit1, data = input.to.stan$data,
                     iter = niters, chains = nchains,
                     init=input.to.stan$init,
                     pars=input.to.stan$para, thin=nthin)
fit2coef <- as.data.frame(rstan::extract(fit2keep, permuted=T))
fit2coef_rv <- rvsims(as.matrix(fit2coef))

updateN <- function(md_code=sequpdate, prior=priors,
                    grp=1990){
    input.to.stanSeq <- inputSEQ(dat=neuse, ycol="SURFCHLA",
                                 gr="Year", sub=neuse$Year==grp,
                                 mu0=prior$mu0, n0=prior$n0,
                                 alpha=prior$alpha, beta=prior$beta)
    fit1 <- stan_model(model_code = md_code)
    fit2 <- sampling(fit1, data = input.to.stanSeq$data,
                 iter = niters, chains = nchains,
                 init=input.to.stanSeq$init,
                 pars=input.to.stanSeq$para, thin=nthin,
                 control=list(adapt_delta=0.99))
    pairs(fit2)
    print(fit2)
    fit2coefSeq <- rstan::extract(fit2, permuted=T)
    fit2coefSeq_rv <- rvsims(as.matrix((as.data.frame(fit2coefSeq))))
    mus <- fit2coefSeq_rv["mu"]
    lmbd <- fit2coefSeq_rv["lambda"]
    priors <- prior_pars(rvmean(mus), rvsd(mus)^2,
                         rvmean(lmbd),rvsd(lmbd)^2)
    return(list(prrs=priors, output_rv=fit2coefSeq_rv))
}

uniqueYr <- as.numeric(levels(ordered(neuse$Year)))
nn <- length(uniqueYr)
tmp <- list()
priors <- init_prior
for (i in 1:nn){
    print(uniqueYr[i])
    tmp[[i]] <- updateN(prior=priors, grp=uniqueYr[i])
    priors <- tmp[[i]]$prrs
}

gmus <- tmp[[1]]$output_rv[2]
seqmu <- tmp[[1]]$output_rv[1]
seqsig1 <- tmp[[1]]$output_rv[3]
seqtau <- 1/sqrt(tmp[[1]]$output_rv[4])

for (i in 2:nn){
    seqmu <- c(seqmu, tmp[[i]]$output_rv[1])
    gmus <- c(gmus, tmp[[i]]$output_rv[2])
    seqsig1 <- c(seqsig1, tmp[[i]]$output_rv[3])
    seqtau <- c(seqtau, tmp[[i]]$output_rv[4])
}

tikz(file=paste(plotDIRch4, "bhm_seqComp1.tex", sep="/"),
     height=2.75, width=5, standAlone=F)
par(mfrow=c(1,2), oma=c(3,3,3,3),mar=c(0,0,0,0),
    mgp=c(1.5, 0.125,0), las=1, tck=0.01)
plot(gmus, axes=F, xlab="years", ylab="log Chla", ylim=c(0.5,3.5))
axis(1, at=1:22, labels=uniqueYr)
axis(2)
axis(3, at=1:22, labels=NA)
box()
text(x=11, y=3.25, labels="Sequential")

plot(BHM_gmu, axes=F, xlab="years", ylab="log Chla", ylim=c(0.5,3.5))
axis(1, at=1:22, labels=NA)
axis(3, at=1:22, labels=uniqueYr)
axis(4)
box()
text(x=11, y=3.25, labels="BHM")
mtext(text="log Chla", side=2, line=1.5, outer=T, las=0)
mtext(text="Year", side=1, line=1.5, outer=T)
dev.off()

tikz(file=paste(plotDIRch4, "bhm_seqComp1Alt.tex", sep="/"),
     height=3, width=3, standAlone=F)
par(oma=c(3,3,3,1),mar=c(0,0,0,0),
    mgp=c(1.5, 0.125,0), las=1, tck=0.01)
mlplot(exp(cbind(gmus, BHM_gmu)), axes=F, xlab="",
       cex=0.75)
axis(1, cex=0.75)
axis(2, at=1:22, labels=uniqueYr, cex=0.75)
box()
mtext(text="Chlorophyll \\textit{a} ($\\mu$g/L)", side=1, line=1.25, outer=T)
dev.off()

seqmu2 <- c(seqmu, BHM_hmu)
seqsig12 <- c(seqsig1, BHM_sigma)
seqtau2 <- c(seqtau, BHM_hsig)

tikz(file=paste(plotDIRch4, "bhm_seqComp2.tex", sep="/"),
     height=3.5, width=4.00, standAlone=F)
par(mfrow=c(3,1), oma=c(3, 3, 3, 3), mar=c(0,3,0,0),
    mgp=c(1.25,0.125,0.0), las=1, tck=0.01)
plot(seqmu2, ylab="$\\mu$", ylim=c(-2,3), axes=F)
axis(2)
axis(4, labels=NA)
box()
axis(3, at=23, labels="BHM")
abline(v=22.5)
plot(seqsig12, ylab="$\\sigma$", axes=F)
axis(2, labels=NA)
axis(4)
box()
abline(v=22.5)
plot(seqtau2, ylab="$\\tau$", ylim=c(0,8), axes=F)
axis(1, at=1:23, labels=c(uniqueYr, "BHM"))
axis(2)
axis(4, labels=NA)
box()
abline(v=22.5)
dev.off()
