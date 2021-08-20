##########################################################################
############# Apportionment - P(W|TL) - in Stan for BEESRStan  ###########
##########################################################################
#use this when generating figures for the book
source("FrontMatter.R")
plotDIRch6 <- paste(plotDIR, "chapter6", "figures", sep="/")


## load packages
packages(rstan)
packages(rstantools)
packages(withr)
packages(rv)
##packages(rstudioapi)
packages(arm)
packages(loo)
packages(bayesplot)
packages(tikzDevice)
packages(MCMCvis)


### the set up
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))
nchains <- min(c(parallel::detectCores(), 8))
niters <- 5000 #50000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)


##############################################
## The model
##############################################
pW_TL <- "
data {
    int<lower=0> N;              // num individuals
    int<lower=1> K;              // num ind predictors
    int<lower=1> J;              // num groups
    int<lower=1> L;              // num group predictors
    int<lower=1,upper=J> jj[N];  // group for individual
    matrix[N,K] x;               // individual predictors
    int<lower=0, upper=1> y[N];                 // outcomes
    matrix[J,L] u;           // group predictors
}

parameters {
    matrix[J,K] z_B;
    corr_matrix[K] Omega;        // prior correlation
    vector<lower=0>[K] tau;      // prior scale
    matrix[L,K] mu_B;            // group coeffs
}

transformed parameters {
    matrix[J, K] B;
    vector[N] theta;

    B = u * mu_B + (z_B * quad_form_diag(Omega, tau));
    theta = rows_dot_product(B[jj] , x);
}

model {
    to_vector(z_B) ~ std_normal();
    tau ~ cauchy(0, 2);
    Omega ~ lkj_corr(2);
    to_vector(mu_B) ~ normal(0, 5);

    y ~ bernoulli_logit(theta);
}"

##############################################



##############################################
## The data and initial values
##############################################
pW_TL_in <- function(y,x,jj,u, n.chains=nchains){
  N <- length(y)
  K <- dim(x)[2]
  J <- max(jj)
  L <- dim(u)[2]
  jj <- jj
  x <- x
  y <- y
  u=u

  data <- list(N=N, K=K, J=J, L=L, jj=jj, x=x, y=y, u=u)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(z_B=(matrix(rnorm(J*L,0,1),ncol=K)),
                       mu_B=(matrix(c(rnorm((L*K)/2,-30,5),
                                      rnorm((L*K)/2,0,5)),nrow=L)),
                       Omega=(diag(1,K)),
                       tau=(runif(2,0,pi/2)))
  paras <- c("B",
             "mu_B",
             "Omega","tau")
  return(list(data=data, init=inits,
              para=paras, nchains=n.chains))
}

## bring in data and subset
data <- read.csv(paste(dataDIR, "GN_length.csv", sep="/")) ## survey data
data <- data[complete.cases(data),]
data_g <- read.csv(paste(dataDIR, "ysi.csv", sep="/")) ## grid level ysi data
data <- data[order(data$Grid),]

## create new grid to fit the model
MyLookupTable <- data.frame(unique(data$Grid),rank(unique(data$Grid)))
colnames(MyLookupTable) <- c("Grid","newGrid")
data$newGrid = MyLookupTable[match(data$Grid,
                                   MyLookupTable$Grid), "newGrid"]
data_g$newGrid = MyLookupTable[match(data_g$Grid,
                                     MyLookupTable$Grid), "newGrid"]
data_g = data_g[complete.cases(data_g),]
head(data); tail(data)
head(data_g); tail(data_g)

## name data
y  <- data$walleye
x  <- matrix(c(rep(1,length(data$LENGTH)),(data$LENGTH/10)),
             ncol=2, byrow=F) #-mean(data$LENGTH/10)
jj <- data$newGrid
u  <- matrix(c(rep(1,max(data_g$newGrid)),data_g$turb),
             ncol=2, byrow=F)
##############################################



##############################################
## Compiling, input data, and running the model
##############################################
fit <- stan_model(model_code = pW_TL)

input.to.stan <- pW_TL_in(y,x,jj,u)

keep <- sampling(fit, data=input.to.stan$data,
                 init=input.to.stan$init,
                 pars=input.to.stan$para,
                 iter=niters,thin=nthin,
                 chains=input.to.stan$nchains,
                 show_messages=F, verbose=F)
                 #,control = list(adapt_delta = 0.99, max_treedepth = 25))

saveRDS(keep,"pW_TL_cor.rds")


# #pairs(keep, pars="Omega") ## divergent transitions below diagonal...
#
#
## check into divergent transitions...
# params <- as.data.frame(extract(keep, permuted=FALSE))
# names(params) <- gsub("chain:1.", "", names(params), fixed = TRUE)
# names(params) <- gsub("[", ".", names(params), fixed = TRUE)
# names(params) <- gsub("]", "", names(params), fixed = TRUE)
# params$iter <- 1:length(params$tau.1)
#
# divergent <- get_sampler_params(keep, inc_warmup=FALSE)[[1]][,'divergent__']
# params$divergent <- divergent
# div_params <- params[params$divergent == 1,]
# nondiv_params <- params[params$divergent == 0,]
#
#
# ## it appears that Neal's funnels has been properly explored
# par(mfrow=c(1,2), mar = c(4, 4, 0.5, 0.5))
# plot(nondiv_params$'mu_B.1,1', log(nondiv_params$tau.1),
#      col="black", pch=16, cex=0.8, xlab="theta.1", ylab="log(tau)")
# points(div_params$'mu_B.1,1', log(div_params$tau.1),
#        col="green", pch=16, cex=0.8)
#
# plot(nondiv_params$'mu_B.1,2', log(nondiv_params$tau.2),
#      col="black", pch=16, cex=0.8, xlab="theta.1", ylab="log(tau)")
# points(div_params$'mu_B.1,2', log(div_params$tau.2),
#        col="green", pch=16, cex=0.8)

#############################################



##############################################
## Processing Stan output and summaries - for book
##############################################
# use this when you get above code working...
#keep <- readRDS("model_pW_TL_cor.rds")
fitcoef <- (rvsims(as.matrix(as.data.frame(
  extract(keep, permute=T)))))

options(digits=2)
summary(fitcoef)


## coef plot
##png("Fig_6.4.5.2_1.png",
##    width=6,height=6, units = 'in', res = 450, pointsize=8, bg="white")
tikz(file=paste(plotDIRch6, "GNapport_hierCoef.tex", sep="/"),
     height=5, width=5, standAlone=F)
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,2,3,1),
    mgp=c(1.25,0.125,0), las=1, tck=0.01)
mlplot(fitcoef[1:21], xlab="Grid-specific intercepts", cex=0.5)
mlplot(fitcoef[22:42], xlab="Grid-specific slopes", cex=0.5)
mtext("Grid-specific parameter estimates", 3, 0, outer=T, cex=1.5)
dev.off()

table(data$newGrid,data$Basin)
## plot all estimated probablity curves, black and white for the paper
x <- seq(0,100, by=1)
clrs <-  (u[,2]/max(u[,2]))
head(data)

##png("Fig_6.4.5.2_2.png",
##    width=4,height=6, units = 'in', res = 450, pointsize=8, bg="white")
tikz(file=paste(plotDIRch6, "GNapport_logistic.tex", sep="/"),
     height=5.5, width=3, standAlone=F)
par(mfrow=c(3,1),mar=c(3,3,1,1),mgp=c(1.25,0.125,0),
    oma=c(4,4,1,1),tck=0.01)
plot(data$LENGTH[data$Grid<8]/10,data$walleye[data$Grid<8],
     ylab="",col="dark gray", xlab="",
     col.axis="black",col.lab="black",cex=1, cex.axis=1.5,
     cex.lab=1.5, bty="n", fg="black",xlim=c(10,70), las=1)
text(60,0.4,"Region 1",cex=1.5)

for(i in 1:7){
    curve(invlogit(rvmean(fitcoef[i]) + rvmean(fitcoef[i+21])*x),
          add=TRUE,col=gray(level=1-clrs[i],alpha=1), lwd=2)
}

plot(data$LENGTH[data$Grid>7 & data$Grid<15]/10,
     data$walleye[data$Grid>7 & data$Grid<15], ylab="",
     col="dark gray", xlab="",
     col.axis="black",col.lab="black",cex=1, cex.axis=1.5,
     cex.lab=1.5, bty="n", fg="black",xlim=c(10,70), las=1)
text(60,0.4,"Region 2",cex=1.5)

for(i in 7:15){
    curve(invlogit(rvmean(fitcoef[i]) + rvmean(fitcoef[i+21])*x),
          add=TRUE,col=gray(level=1-clrs[i],alpha=1), lwd=2)
}

plot(data$LENGTH[data$Grid>14]/10,data$walleye[data$Grid>14],
     ylab="",col="dark gray", xlab="",
     col.axis="black",col.lab="black",cex=1, cex.axis=1.5,
     cex.lab=1.5, bty="n", fg="black",xlim=c(10,70), las=1)
text(60,0.4,"Region 3", cex=1.5)

for(i in 15:21){
    curve(invlogit(rvmean(fitcoef[i]) + rvmean(fitcoef[i+21])*x),
          add=TRUE,col=gray( level=1-clrs[i],alpha=1), lwd=2)
}

mtext("$\\Pr(W \\mid TL,\\hat{\\beta})$",2,2,cex=1., outer=T)
mtext("Total length (cm)",1,2,cex=1., outer=T)
dev.off()

##############################################


# ##############################################
# ## Posterior predictive checks for paper
# ##############################################
# # not exactly sure how to interpret this...
# ## still mis-calibrated...
# yrep <- posterior_predict(zinb2)
# color_scheme_set("gray")
#
#
# png("sturg_ppc.png", units="in",pointsize=12,width=5,height=4,res=600)
# #par(mfrow=c(2,2))
#
# fig1 <- ppc_stat_2d(y = data$CATCH,
#                     yrep = yrep,
#                     stat=c("mean","sd"))
#
# fig2 <- ppc_loo_pit_overlay(
#   y = data$CATCH,
#   yrep = yrep,
#   lw = weights(zinb2_loo$psis_object))
#
# fig3 <- ppc_dens_overlay( y = log(data$CATCH+1),
#                           yrep = log(yrep[1:100,]+1))
#
# fig4 <- ppc_stat_grouped(y = data$CATCH,
#                          yrep = yrep, group = data$YEAR,
#                          stat="mean", binwidth=0.1)
#
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# fig1
# fig2
# fig3
# fig4)
# dev.off()
#
#
# ##############################################

