##########################################################################
############# Apportionment - P(W|TL) - in Stan for BEESRStan  ###########
##########################################################################
#use this when generating figures for the book
source("FrontMatter.R")
plotDIRch6 <- paste(plotDIR, "chapter6", "figures", sep="/")
packages(arm)
packages(rstan)
packages(rv)

## data and glmer
data <- read.csv(paste(dataDIR, "GN_length.csv", sep="/")) ## survey data
data_g <- read.csv(paste(dataDIR, "ysi.csv", sep="/"))## grid level ysi data
data <- data[order(data$Grid),]
data_g <- data_g[order(data_g$Grid),]

head(data)
head(data_g)

data$len.c <- (data$LENGTH-mean(data$LENGTH))/10
wly.glmer <- glmer(walleye ~ len.c+(1+len.c|Grid),
                   data=data, family="binomial")
summary(wly.glmer)
data$Grid <- as.numeric(ordered(data$Grid))
data$Turb <- data_g$turb[data$Grid]
wly.glmer <- glmer(walleye ~ len.c+Turb+Turb:len.c+
                       (1+len.c|Grid),
                   data=data, family="binomial")
summary(wly.glmer)
fixef(wly.glmer)
ranef(wly.glmer)

temp <- table(data$Grid)>1
M2.coef <- coef (wly.glmer)
a.hat.M2 <- M2.coef[[1]][,1] + M2.coef[[1]][,3]*data_g$turb[temp]
b.hat.M2 <- M2.coef[[1]][,2] + M2.coef[[1]][,4]*data_g$turb[temp]
a.se.M2 <- se.ranef(wly.glmer)[[1]][,1]
b.se.M2 <- se.ranef(wly.glmer)[[1]][,2]

# plot estimated intercepts and slopes

## Figure 10.13
tikz(file=paste(plotDIRch6, "wlyVlength.tex", sep="/"),
     width=4.75, height=3, standAlone=F)
par (mfrow=c(1,2), mar=c(3,3,3,0.25), mgp=c(1.25,0.125,0), las=1, tck=0.01)
lower <- a.hat.M2 - a.se.M2
upper <- a.hat.M2 + a.se.M2
plot (data_g$turb[temp], a.hat.M2,
      ylim=range(lower,upper), cex.lab=0.75, cex.axis=0.75,
      xlab="average turbidity", ylab="regression intercept", pch=20)
curve (fixef(wly.glmer)["(Intercept)"] + fixef(wly.glmer)["Turb"]*x, lwd=1, col="black", add=TRUE)
segments (data_g$turb[temp], lower, data_g$turb[temp], upper,
          lwd=.5, col="gray10")
##text(AveTemp, lower, levels(rtol2$city), adj=c(.5,1),cex=0.5)

lower <- b.hat.M2 - b.se.M2
upper <- b.hat.M2 + b.se.M2

plot (data_g$turb[temp], b.hat.M2, ylim=range(lower,upper),
      cex.lab=0.75, cex.axis=0.75,
      xlab="average turbidity", ylab="regression slope", pch=20)
curve (fixef(wly.glmer)["len.c"] + fixef(wly.glmer)["len.c:Turb"]*x,
       lwd=1, col="black", add=TRUE)
segments (data_g$turb[temp], lower, data_g$turb[temp], upper,
          lwd=.5, col="gray10")
dev.off()

