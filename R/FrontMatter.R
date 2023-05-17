### Applications of Bayesian Methods in Environmental Science With R and Stan
### Copyright Song S. Qian
### R scripts
### December 15, 2018
###

## Loading and installing (if not already installed)
##  packages

packages<-function(x, repos="http://cran.r-project.org", ...){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos=repos, ...)
    require(x,character.only=TRUE)
  }
}

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin')
    os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

os <- get_os()
if (os=="linux"){
    base <- "~/google-drive/BEESwithRStan/Book"
} else if (os=="osx") {
    base <- "~/Google\ Drive/BEESwithRStan/Book"
} else {
    base <- "C:/Users/songq/Google\ Drive/BEESwithRStan/Book"
}
RHome <- paste(base, "R", sep="/")
dataDIR <- paste(base, "R","Data", sep="/")
plotDIR <- paste(base, "chapters", sep="/")

#setwd(RHome)

packages(arm)
packages(lattice)
packages(tikzDevice)
source("https://raw.githubusercontent.com/songsqian/BeesRStan/main/R/BeesRStan.R")


