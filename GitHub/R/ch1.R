source("FrontMatter.R")

## Chapter 1
plotDIRch1 <- paste(plotDIR, "chapter1", "figures", sep="/")

x <- 3
n <- 10
pri <- (1:5)/10
prior <- rep(0.2,5)
likelihood <- dbinom(x,n, pri)
post <- likelihood*prior
pi_gvn_y <- post/sum(post)

bayes_binom <- function(x, n, p_i=(1:5)/10, prior=NULL){
  k <- length(p_i)
  if (is.null(prior)) prior <- rep(1/k, k)
  likelihood <- dbinom(x, n, p_i)
  return(likelihood * prior / (sum(likelihood * prior)))
}

pi <- seq(0.1,0.5,0.05)
post <- bayes_binom(x=5, n=22, p_i=pi)

print(cbind(pi, post))
print(pi_gvn_y <- pi[post==max(post)])


tikz(file=paste(plotDIRch1, "betas.tex", sep="/"),
     height=3, width=4, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(c(0, 1), c(0, 3), type="n", xlab="", ylab="")
lines(seq(0,1,,100), dbeta(seq(0,1,,100), 0.5,0.5))
lines(seq(0,1,,100), dbeta(seq(0,1,,100), 1,1), lty=2)
lines(seq(0,1,,100), dbeta(seq(0,1,,100), 1,5), lty=3)
lines(seq(0,1,,100), dbeta(seq(0,1,,100), 2,4), lty=4)
lines(seq(0,1,,100), dbeta(seq(0,1,,100), 5,1), lty=5)

legend(x=0.4, y=3, lty=1:5, legend=c("$\\alpha=\\beta=0.5$",
                                        "$\\alpha=\\beta=1$",
                                        "$\\alpha=1, \\beta=5$",
                                        "$\\alpha=2, \\beta=4$",
                                     "$\\alpha=5, \\beta=1$"),
       cex=0.7, bty="n" )
dev.off()

tikz(file=paste(plotDIRch1, "post_beta.tex", sep="/"),
     height=2.5, width=3, standAlone=F)
par(mar=c(3,3,1,0.5), mgp=c(1.25, 0.125,0), las=1, tck=0.01)
plot(seq(0, 1, 0.01), dbeta(seq(0,1,0.01), 6, 18), type="l",
     xlab="$p$", ylab="$\\pi(p|y)$")
dev.off()

post_impft <- function(x=5, n=20, fp=0.07, fn=0.05, k=100){
    theta <- seq(0, 1,, k)
    fpst <- theta*(1-fn) + (1-theta)*fp
    post <- x*log(fpst) + (n-x)*log(1-fpst)
    return(list(pdf=exp(post)/(theta[2]*sum(exp(post))), cdf=cumsum(exp(post))/sum(exp(post))))
}

tikz(file=paste(plotDIRch1, "postimpft.tex", sep="/"),
                height=2.5, width=3, standAlone=F)
par(mar=c(3, 3, 1, 0.5), mgp=c(1.25, 0.125, 0), las=1, tck=0.01)
plot(seq(0, 1,, 100), post_impft()$pdf, type="l", xlab="$\\theta$",
     ylab="$\\pi(\\theta|y)$")
dev.off()

