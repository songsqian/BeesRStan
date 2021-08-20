source("FrontMatter.R")

## Chapter 3
plotDIRch3 <- paste(plotDIR, "chapter3", "figures", sep="/")

## Neuse River Estuary

neuse <- read.csv(paste(dataDIR, "NeuseEstuary.csv", sep="/"))
neuse[neuse<0] <- NA

neuse$RDate <- as.Date(as.character(neuse$DATE), format="%m/%d/%Y")
neuse$Year <- as.numeric(format(neuse$RDate, "%Y"))
## Acceptance-Rejection Method Examples
## Normal from exponential

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

## Neuse River Estuary Example ##

## prior parameters
## Model 1--TN
pst <- list(alpha =32,
            beta = 28,
            nn = 45,
            mu = 7.6)

NGpost <- function(x, alpha, beta, n0, mu0){
    x_bar <- mean(x)
    n <- length(x)
    s2 <- sd(x)^2
    return(list(nn=n+n0, mu=(n*x_bar+n0*mu0)/(n+n0),
           alpha = alpha+n/2,
           beta = beta+0.5*(n*s2 + (n*n0)*(x_bar-mu0)^2/(n+n0))))
}

tmp <- neuse$Year>=1992 & neuse$Year <= 2000 & ##neuse$SECTION=="UPPER")|
        neuse$SECTION=="MIDDLE"
neuse2 <- neuse[tmp,]

tikz(file=paste(plotDIRch3, "neuseTN.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$TOTN), prob=T, nclass=20,
     axes=F, xlab="Total Nitrogen ($\\mu$g/L)",
     ylim=c(0,1.5), main="", ylab="")
axis(1, at = log(c(10, 100, 250, 500, 1000, 2500, 5000)),
                 labels=c(10, 100,250, 500, 1000, 2500, 5000))
pst <- list(alpha =32,
            beta = 28,
            nn = 45,
            mu = 7.6)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2
curve(dt((x-pr_mut)/pr_sigmat, df=pr_df)/pr_sigmat, add=T, lwd=3)

for (i in 1992:2000){
    tmp <- neuse2$Year==i & !is.na(neuse2$TOTN)
    pst <- NGpost(log(neuse2$TOTN[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
    curve(dt((x-mu_t)/sigma_t, df=df_t)/sigma_t, lty=i-1990, add=T)
}
dev.off()

tikz(file=paste(plotDIRch3, "neuseTN2.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$TOTN), prob=T, nclass=20,
     axes=F, xlab="Total Nitrogen ($\\mu$g/L)",
     ylim=c(0,1), main="", ylab="")
axis(1, at = log(c(10, 100, 250, 500, 1000, 2500, 5000)),
                 labels=c(10, 100,250, 500, 1000, 2500, 5000))
pst <- list(alpha =32,
            beta = 28,
            nn = 45,
            mu = 7.6)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2
segments(x0=pr_mut+pr_sigmat*qt(0.025,df=pr_df),
        x1=pr_mut+pr_sigmat*qt(0.975,df=pr_df),
        y0=0.1, y1=0.1)
segments(x0=pr_mut+pr_sigmat*qt(0.25,df=pr_df),
        x1=pr_mut+pr_sigmat*qt(0.75,df=pr_df),
        y0=0.1, y1=0.1, lwd=3)
points(x=pr_mut+pr_sigmat*qt(0.5,df=pr_df), y=0.1)
j <- 0.1
##curve(dt((x-pr_mut)/pr_sigmat, df=pr_df)/pr_sigmat, add=T, lwd=3)

for (i in 1992:2000){
    j <- j+0.1
    tmp <- neuse2$Year==i & !is.na(neuse2$TOTN)
    pst <- NGpost(log(neuse2$TOTN[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
segments(x0=mu_t+sigma_t*qt(0.025,df=pr_df),
        x1=mu_t+sigma_t*qt(0.975,df=pr_df),
        y0=j, y1=j)
segments(x0=mu_t+sigma_t*qt(0.25,df=pr_df),
        x1=mu_t+sigma_t*qt(0.75,df=pr_df),
        y0=j, y1=j, lwd=3)
points(x=mu_t+sigma_t*qt(0.5,df=pr_df), y=j)
}
axis(2, at=seq(0.1,1,0.1), labels=c("Prior", 1992:2000))
dev.off()

tikz(file=paste(plotDIRch3, "neuseCHLA.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$SURFCHLA), prob=T, nclass=20,
     axes=F, xlab="Chlorophyll a ($\\mu$g/L)",
     ylim=c(0,1), main="", ylab="")
axis(1, at = log(c(0.1,1, 5, 10, 100, 250, 500)),
                 labels=c(0.1,1, 5, 10, 100,250, 500))
pst <- list(alpha =100,
            beta = 40,
            nn = 40,
            mu = 3.4)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2
curve(dt((x-pr_mut)/pr_sigmat, df=pr_df)/pr_sigmat, add=T, lwd=3)

for (i in 1992:2000){
    tmp <- neuse2$Year==i & !is.na(neuse2$SURFCHLA)
    pst <- NGpost(log(neuse2$SURFCHLA[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
    curve(dt((x-mu_t)/sigma_t, df=df_t)/sigma_t, lty=i-1990, add=T)
}
dev.off()

tikz(file=paste(plotDIRch3, "neuseCHLA2.tex", sep="/"),
     width=3.5, height=3, standAlone=F)
par(mar=c(3,3,1,1), mgp=c(1.25, 0.125,0),las=1, tck=0.01)
hist(log(neuse2$SURFCHLA), prob=T, nclass=20,
     axes=F, xlab="Chlorophyll a ($\\mu$g/L)",
     ylim=c(0,0.5), main="", ylab="", border="white",col="gray")
axis(1, at = log(c(0.1,1, 5, 10, 100, 250, 500)),
                 labels=c(0.1,1, 5, 10, 100,250, 500))
pst <- list(alpha =100,
            beta = 40,
            nn = 40,
            mu = 3.4)

pr_mut <- pst$mu
pr_sigmat <- (pst$beta/pst$alpha)*(1+1/pst$nn)
pr_df <- pst$alpha*2

segments(x0=pr_mut+pr_sigmat*qt(0.025,df=pr_df),
        x1=pr_mut+pr_sigmat*qt(0.975,df=pr_df),
        y0=0.05, y1=0.05)
segments(x0=pr_mut+pr_sigmat*qt(0.25,df=pr_df),
        x1=pr_mut+pr_sigmat*qt(0.75,df=pr_df),
        y0=0.05, y1=0.05, lwd=3)
points(x=pr_mut+pr_sigmat*qt(0.5,df=pr_df), y=0.05)
j <- 0.05
for (i in 1992:2000){
    j <- j+0.05
    tmp <- neuse2$Year==i & !is.na(neuse2$SURFCHLA)
    pst <- NGpost(log(neuse2$SURFCHLA[tmp]),
                  pst$alpha, pst$beta, pst$nn, pst$mu)
    mu_t <- pst$mu
    sigma_t <- (pst$beta/pst$alpha)*(1+1/pst$nn)
    df_t <- pst$alpha*2
    segments(x0=mu_t+sigma_t*qt(0.025,df=pr_df),
             x1=mu_t+sigma_t*qt(0.975,df=pr_df),
             y0=j, y1=j)
    segments(x0=mu_t+sigma_t*qt(0.25,df=pr_df),
             x1=mu_t+sigma_t*qt(0.75,df=pr_df),
             y0=j, y1=j, lwd=3)
    points(x=mu_t+sigma_t*qt(0.5,df=pr_df), y=j)
}
axis(2, at=seq(0.05,0.5,0.05), labels=c("Prior", 1992:2000))
dev.off()
