## SFS 2024 Short Course 
## -- Bayesian Applications in Environmental and 
##    Ecological Studies with R and Stan
## Frontmatters 

knitr::opts_chunk$set(echo = TRUE)
rootDIR <- "https://raw.githubusercontent.com/songsqian/BeesRStan/main/R"
dataDIR <- paste(rootDIR, "Data", sep="/")

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}
source(paste(rootDIR, "BeesRStan.R", sep="/"))

require(rstan)
packages(rv)
packages(car)
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)
