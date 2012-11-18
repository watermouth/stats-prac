# understanding what importance sampling means
# 13:10 - 15:10

# using uniform
# test(n=10000,params=list(normDistMu=0, normDistSigma=1))
# using normal
# test(n=10000,params=list(normDistMu=0, normDistSigma=1, rdistName="normal", mean=0, sd=1))

# test: plot and calculate mean and standard deviation of the approximated distribution
test <- function(n=100, params=NULL){
  temp <- importanceSampling(c(list(n = n),params))
  x <- temp$x
  wt <- temp$wt
  w <- temp$w
  print(paste("n=", n))
  cat("mean(dist)=",sum(wt*x), "\n")
  cat("sd(dist)=", sqrt((sum(wt*x^2) - sum(wt*x)^2)/(n-1)), "\n")
  plot(x,wt)
  browser()
  # extract wt > 1 / n
  idx <- which(wt > 1/n)
  if(length(idx)==0) idx <- seq(along.with=wt)
  x <- x[idx]
  wt <- wt[idx]
  print(paste("n=", length(idx)))
  cat("mean(dist)=",sum(wt*x), "\n")
  cat("sd(dist)=", sqrt((sum(wt*x^2) - sum(wt*x)^2)/(n-1)), "\n")
  plot(x,wt)
}

# approximate normal distribution using particles sampled from uniform distribution
importanceSampling <- function(params){
  # target distribution: normal
  dNormDistMaker <- function(mu,sigma){
    f <- function(x) dnorm(x,mu,sigma)
    return(f)
  }
  normDistMu <- params$normDistMu
  if(is.null(params$normDistMu)){
    normDistMu <- 0
  } 
  normDistSigma <- params$normDistSigma
  if(is.null(params$normDistSigma)){
    normDistSigma <- 1
  }
  dTarget <- dNormDist(normDistMu,normDistSigma)
  # particles sampled from the proposal distributions, and probality weights of the distribution
  temp <- getProposalDists(params)
  x <- temp$x
  w <- temp$w
  # weights for target distribution
  wt <- dTarget(x)
  wt <- wt / w
  wt <- wt / sum(wt)
  return(list(x=x, wt=wt,w=w, normDistMu=normDistMu, normDistSimga=normDistSigma))
}

getProposalDists <- function(params){
  # proposal function
  if(is.null(params$rdistName)) {
    params$rdistName <- "uniform"
    params$min <- -10
    params$max <-  10
  }
  switch(params$rdistName,
         "uniform" = {
           x <- rdist(params$n,min=params$min,max=params$max);
           w <- dunif(x,min=params$min, max=params$max);
         },
         "normal" = {
           x <- rnorm(params$n, mean=params$mean, sd=params$sd)
           w <- dnorm(x,mean=params$mean, sd=params$sd)
         }
  )
  list(x=x,w=w)
}


