#' getMean
#'
#' Takes betas, x.matrix, and link and returns the estimated mean
#' @param betas  p vector of regression parameters
#' @param x.matrix nxp matrix of covariates
#' @param link link between xbeta and the mean, identify for normal, logit for logistic, Gammalog for gamma with log link
#' @return estimated Mean

getMean<-function(betas, x.matrix, link){
  #returns the mean given the betas and x.matrix
  if(link =="identity")return (x.matrix %*% betas)
  if(link =="logit") {
    logit<-x.matrix %*% betas
    pi<-exp(logit)/(1+exp(logit))
    return(pi)
  }
  if(link%in% c("Gammalog", "Poissonlog"))return(exp(x.matrix%*%betas))
}
