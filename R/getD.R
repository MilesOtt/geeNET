#' getD
#'
#' Takes betas and x.matrix and returns matrix of partial derivatives dmu/dbeta
#' @param betas  p vector of regression parameters
#' @param x.matrix nxp matrix of covariates
#' @param link link between xbeta and the mean, identify for normal, logit for logistic, Gammalog for gamma (should work)
#' @return  matrix of partial derivatives dmu/dbeta
getD<-function(betas, x.matrix, link){


  if(link=="identity") return(x.matrix)

  if(link=="logit"){

    p=length(betas)
    D<-x.matrix
    logit<-x.matrix %*% betas
    pi<-exp(logit)/((1+exp(logit))^2)
    for(i in 1:p){
      D[,i]<-x.matrix[,i]*pi
    }
}
    if(link %in% c("Gammalog", "Poissonlog")){
      p=length(betas)
      D<-x.matrix
      mu_i<-exp(x.matrix%*%betas)
      for(i in 1:p){
        D[,i]<-x.matrix[,i]*mu_i
      }
    }
    return(D)

}

