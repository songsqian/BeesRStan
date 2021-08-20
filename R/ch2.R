source("FrontMatter.R")

## Chapter 2
plotDIRch2 <- paste(plotDIR, "chapter2", "figures", sep="/")

## BMC of the SFD example
n_sims <- 50000

## from chapter 1

post_impft <- function(x=5, n=20, fp=0.07, fn=0.05, k=100){
    theta <- seq(0, 1,, k)
    fpst <- theta*(1-fn) + (1-theta)*fp
    post <- x*log(fpst) + (n-x)*log(1-fpst)
    return(list(pdf=exp(post)/(theta[2]*sum(exp(post))), cdf=cumsum(exp(post))/sum(exp(post))))
}

bmc_sfd_post <- function(prior = sort(runif(n_sims)), x=5, n=20,
                         fp=0.07, fn=0.05, pdfgrid=seq(0,1,0.01)){
    log_likelihood <- x * log(prior*(1-fn) + (1-prior)*fp) +
        (n-x) * log(1-prior*(1-fn)-(1-prior)*fp)
    post <- exp(log_likelihood)/sum(exp(log_likelihood))
    post_cdf <- cumsum(post)
    pdf_interval <- findInterval(prior, pdfgrid)
##    post_pdf <- unlist(by(post_cdf, pdf_interval, function(x)diff(range(x))))
    post_pdf <- unlist(by(post, pdf_interval, sum))
    return(list(theta=prior, cdf=post_cdf, pdf=post_pdf, grid=pdfgrid[-1]))
}
sfdbmc <- bmc_sfd_post()
plot(sfdbmc$cdf ~ sfdbmc$theta, type="l", ylab="cdf", xlab="p")
plot(sfdbmc$pdf ~ sfdbmc$grid, type="l", ylab="pdf", xlab="p")

## Inverse-CDF
set.seed(101)
mu <- 2
sigma <- 1.25
u <- runif(10000)
x <- qnorm(u, mu, sigma)
summary(x)

y <- rnorm(10000, mu, sigma)
summary(y)

hist(x, col=rgb(0.1,0.1,0.1,0.5), xlim=c(-3,6), ylim=c(0, 1800),
     nclass=20, main="")
hist(y, col=rgb(0.8,0.8,0.8,0.5), nclass=20, add=T)
box()

## Inverse-CDF, the SFD example ##
### using BMC
n <- 1000
sfdbmc <- bmc_sfd_post(prior = sort(runif(10000)))
sfd_post <- data.frame(theta=sfdbmc$theta, cdf=sfdbmc$cdf)
u <- runif(n)
tmp <- apply(sfd_post, 1, function(x, unf)
    return(x[2]-unf), unf=u)

## Using evenly spaced pdf from Chapter 1
post_impft_cdf <- function(x=5, n=20, fp=0.07, fn=0.05, k=100){
    theta <- seq(0, 1,, k)
    fpst <- theta*(1-fn) + (1-theta)*fp
    post <- x*log(fpst) + (n-x)*log(1-fpst)
    return(cumsum(exp(post)/sum(exp(post))))
}

post_cdf <- data.frame(theta=seq(0,1,,1000), cdf=post_impft_cdf(k=1000))
u <- runif(n)
tmp <- apply(post_cdf, 1, function(x, unf)
     return(x[2]-unf), unf=u)

theta <- apply(tmp, 1, function(x, theta)
     return(theta[abs(x)==min(abs(x))]),
              theta=post_cdf$theta)
hist(theta)

## Acceptance-Rejection Method Examples
## Standard normal from exponential

n <- 10000
x <- numeric()
for (i in 1:n){
    u1 <- runif(1)
    x[i] <- -log(u1)
    r <- exp(-(1-x[i])^2/2)
    u2 <- runif(1)
    while(u2 > r){
        u1 <- runif(1)
        x[i] <- -log(u1)
        r <- exp(-(1-x[i])^2/2)
        u2 <- runif(1)
    }
}
u <- ifelse(runif(n) < 0.5, -1, 1)
x <- u*x

hist(x)
summary(x)
sd(x)

my_rnorm <- function(n, mu=0, sigma=1){
    x <- -log(runif(n))
    r <- exp(-0.5*(1-x)^2)
    x <- x[r>=runif(n)]
    m <- length(x)
    if (m < n)
        x <- c(x, my_rnorm(n-m, mu, sigma))
    return (mu + sigma * x *  ifelse(runif(n) < 0.5, -1, 1))
}

## relationshis method

my.rbin <- function(m, p, n){

    return(unlist(apply(matrix(as.numeric(runif(n*m)<p), nrow=n), 2, sum)))
}

## numerical integration
nsims <- 10000
## simulated data
### predictive distribution with known sigma
set.seed(2)
sigma <- 2
n <- 25
y <- rnorm(n, 3, sigma)
hat_mu <- mean(y)
rhat_mu <- rnorm(nsims, hat_mu, sigma/sqrt(n))
y_grid <- matrix(seq(1,6,length=500), ncol=1)

pi_yj <- apply(y_grid, 1, function(x)dnorm(x, rhat_mu, sigma/sqrt(25)))
pi_y1 <- apply(pi_yj,2, mean)

plot(pi_y1~y_grid[,1], type="l")

### predictive distribution with unknown mu and sigma
hat_s <- sd(y)
rhat_sig <- sqrt((n-1)*hat_s^2/rchisq(nsims, n-1))
rhat_mu <- rnorm(nsims, hat_mu, rhat_sig/sqrt(n))

pi_yj <- apply(y_grid, 1,
            function(x)dnorm(x, rhat_mu, rhat_sig/sqrt(25)))
pi_y2 <- apply(pi_yj,2, mean)

lines(pi_y2 ~ y_grid[,1], lty=2, type="l")

### 0.75 quantile
q75 <- rhat_mu+rhat_sig * qnorm(0.75)
hist(q75)


## Metropolis-Hasting
## 1. beta distribution
set.seed(10)
n_sims<-4000
x <- numeric()
x[1] <- runif(1) # initial value
a <- 2
b <- 5
for (j in 1:n_sims){
  y <- runif(1)
  alpha <- dbeta(y, a, b)/dbeta(x[j], a, b)
  if (runif(1) < alpha) x[j+1] <- y
  else x[j+1] <- x[j]
}

x <- x[round(1+n_sims/2):n_sims]
tikz(file=paste(plotDIRch2, "met_hast_beta.tex", sep="/"),
     height=2.5, width=5,standAlone=T)
par(mfrow=c(1,2), mar=c(3,3,1,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
plot(x, type="l", ylab="$x$")
hist(x, prob=T, main="", xlab="$x$")
curve(dbeta(x, a, b), add=T)
dev.off()

## Metropolis-Hasting
## 2. Snake fungal disease -- data n=20, x=3
## using unif(0,1) as q()
set.seed(10)
n_sims<-50000
theta <- numeric()
theta[1] <- runif(1) # initial value

n <- 20
x <- 5
log_lk <- function(theta, x=5, n=20, fp=0.07, fn=0.05){
    pp <- theta*(1-fn)+(1-theta)*fp
    llk <- x*log(pp)+(n-x)*log1p(-pp)
    return(llk)
}

for (j in 1:n_sims){
  y <- runif(1)
  alpha <- exp(log_lk(y)-log_lk(theta[j]))
  if (runif(1) < alpha) theta[j+1] <- y
  else theta[j+1] <- theta[j]
}

theta <- theta[round(1+n_sims/2):n_sims]
## numerical program in chapter 1
post_impft <- function(x=5, n=20, fp=0.07, fn=0.05, k=100){
    theta <- seq(0,1,,k)
    fpst <- theta*(1-fn)+(1-theta)*fp
    post <- x*log(fpst)+(n-x)*log1p(-fpst)
    return(exp(post)/(theta[2]*sum(exp(post))))
}

tikz(paste(plotDIRch2, "met_hast_sfd.tex", sep="/"),
     height=2, width=3, standAlone=T)
par(mar=c(3,3,1,0.25), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
hist(theta, prob=T, n=40, xlab="$\\theta$", main="")
lines(seq(0,1,,100), (post2 <- post_impft()))
dev.off()

summary(theta)


## Gibbs sampler examples
## 1. binomial

alpha <- 2
beta <- 2
n <- 20
x <- numeric()
p <- numeric()
p[1] <- rbeta(1, 1, 1)
x[1] <- rbinom(1, n, p[1])
for (j in 2:n_sims){
  p[j] <- rbeta(1, alpha+x[j-1], beta+n-x[j-1])
  x[j] <- rbinom(1, n, p[j])
}

plot(p, x)

## Gibbs sample examples
## 2. PCB in fish
packages(reshape2)
packages(ggplot2)
packages(latex2exp)
packages(rv)
laketrout <- read.csv(paste(dataDIR, "laketrout2.csv", sep="/"))
laketrout$size_cm <- round(laketrout$length * 2.54, 1)
laketrout$gt100 <- laketrout$n * (laketrout$pcb > 1.00)
## # of fish with pcb > 1.00
laketrout$gt190 <- laketrout$n * (laketrout$pcb > 1.90)
## # of fish with pcb > 1.90
laketrout <- laketrout[laketrout$size_cm >0 & !is.na(laketrout$length),]

laketroutM <- melt(laketrout, id.vars="size_cm",
                   measure.vars=c("gt100", "gt190","n"))
laketroutCst <- dcast(laketroutM, size_cm ~ variable, sum)

binmono <- function(yobs, nobs, x, D0, a=10, n_sims=10000, burnin = NULL){
    m <- length(nobs) + 2
    d <- diff(D0)
    f <- matrix(0, ncol=m, nrow=n_sims)
    f[,m] <- 1
    ## initials
    f[1,] <- D0

    for (i in 1:n_sims){
        if (i %% round(0.1*n_sims) == 0)
            print(paste("iter", i, "of", n_sims, "(", i*100/n_sims,"%)"))
        if (i > 1) f[i,] <- f[i-1,]
        for (j in 2:(m-1)){
            ## print(paste("i =", i, "j =",j))
            f[i,j] <- f[i,j-1]+(f[i,j+1]-f[i,j-1])*rbeta(1, a*d[j-1],a*d[j])
            u <- runif(1)
            fstar <- yobs[j-1]/nobs[j-1]
            r <- f[i,j]^yobs[j-1] * (1-f[i,j])^(nobs[j-1]-yobs[j-1]) /
                fstar^yobs[j-1] * (1-fstar)^(nobs[j-1]-yobs[j-1])
            k <- 0
            ## print(c(fstar, f[i,j], r))
            while (u > r){
                k <- k+1
                f[i,j] <- f[i,j-1] + (f[i,j+1]-f[i,j-1])*rbeta(1, a*d[j-1],a*d[j])
                u <- runif(1)
                fstar <- yobs[j-1]/nobs[j-1]
                r <- f[i,j]^yobs[j-1] * (1-f[i,j])^(nobs[j-1]-yobs[j-1]) /
                    fstar^yobs[j-1] * (1-fstar)^(nobs[j-1]-yobs[j-1])
                if (k > 500000) stop(paste("stuck at i = ", i, "j = ", j))
            }
        }
    }
    if (is.null(burnin)) burnin <- n_sims%/%2
    return(rvsims(f[(burnin + 1):n_sims,]))
}

m <- dim(laketroutCst)[1]+2
## very slow ##
frv <- summary(binmono(yobs=laketroutCst$gt100,
                       nobs= laketroutCst$n,
                       D0= seq(0, 1,,m),
                       n_sims=100000))
save(frv, file="ch2_PCBGibbs.RData")
m <- dim(frv)[1]
plotdata <- frv[-c(1,m),]
plotdata$size_cm <- laketroutCst$size_cm
names(plotdata) <- c("mean", "sd", "one", "two.five", "fstQ", "median",
                     "thrdQ", "ninty7.5", "ninty9", "sims","size_cm")
p <- ggplot() + xlab("fish size (cm)") +
    ylab(TeX("$\\Pr(PCB > 1.00 \\mu g/kg)$"))
print(
  p + geom_polygon(data=data.frame(x=c(plotdata$size_cm, rev(plotdata$size_cm)),
                                y=c(plotdata$two.five, rev(plotdata$ninty7.5))),
                 aes(x, y, fill="95%")) +
    geom_polygon(data=data.frame(x=c(plotdata$size_cm, rev(plotdata$size_cm)),
                                 y=c(plotdata$fstQ, rev(plotdata$thrdQ))),
                 aes(x, y, fill="50%")) +
    geom_line(data=plotdata, aes(x=size_cm, y=median))
)

## Snake fungal disease example using Metropolis-Hasting in a Gibbs sampler
alpha <- 1; beta <- 1; ap <- 5; bp <- 20; an <- 2;
bn <- 25; n <- 20; x <- 5

post <- function(theta, fp, fn, x, n)
return ((theta^(alpha-1))*((1-theta)^(beta-1)) *
        (fp^(ap-1))*((1-fp)^(bp-1)) *
        (fn^(an-1))*((1-fn)^(bn-1)) *
        ((theta*(1-fn)+(1-theta)*fp)^x) *
       ((1-theta*(1-fn)-(1-theta)*fp)^(n-x)))

set.seed(101)

met_hast_Gibbs <- function(n_sims=100000, seed=101){
    theta <- fp <- fn <- numeric()
    theta[1] <- rbeta(1, alpha, beta) # initial value
    fp[1] <- rbeta(1, ap, bp)
    fn[1] <- rbeta(1, an, bn)
    for (j in 2:n_sims){
        ## theta
        th <- runif(1)
        r <- post(th, fp[j-1], fn[j-1], x, n) /
            post(theta[j-1], fp[j-1], fn[j-1], x, n)
        theta[j] <- ifelse (runif(1) < r, th, theta[j-1])
        ## fp
        f_p <- runif(1)
        r <- post(theta[j], f_p, fn[j-1], x, n) /
            post(theta[j], fp[j-1], fn[j-1], x, n)
        fp[j] <- ifelse (runif(1) < r, f_p, fp[j-1])
        ## fn
        f_n <- runif(1)
        r <- post(theta[j], fp[j], f_n, x, n) /
            post(theta[j], fp[j], fn[j-1], x, n)
        fn[j] <- ifelse (runif(1) < r, f_n, fn[j-1])
    }

    theta <- theta[round(1+n_sims/2):n_sims]
    fp <- fp[round(1+n_sims/2):n_sims]
    fn <- fn[round(1+n_sims/2):n_sims]
    return(list(theta=theta, fp=fp, fn=fn))
}

sim1 <- met_hast_Gibbs()

tikz(paste(plotDIRch2, "met_hast_Gibbs_sfd.tex", sep="/"),
     height=2.5, width=5.5, standAlone=T)
par(mfrow=c(1,3), mar=c(3, 2, 1, 1), mgp=c(1.25,0.125,0),
    las=1, tck=0.01)
hist(sim1$theta, prob=T, main="", ylab="", xlab="$\\theta$")
curve(dbeta(x, alpha, beta), add=T)
curve(dbeta(x, alpha+5, beta+n-5), add=T, lty=2)
hist(sim1$fp, prob=T, ylab="", main="", xlab="$f_p$")
curve(dbeta(x, ap, bp), add=T)
hist(sim1$fn, prob=T, ylab="", main="", xlab="$f_n$")
curve(dbeta(x, an, bn), add=T)
dev.off()

nn <- length(sim1[[1]])
tmp <- sample(1:nn, size=1000)
par(mfrow=c(1,3))
plot(sim1$theta[tmp], sim1$fp[tmp])
plot(sim1$theta[tmp], sim1$fn[tmp])
plot(sim1$fn[tmp], sim1$fp[tmp])
