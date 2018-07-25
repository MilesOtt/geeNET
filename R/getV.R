#' getV
#'
#' Takes betas, correlation vector rho, x.matrix, y, adj.matrix, link, W, toep and returns the variance covariance matrix
#' @param betas  p vector of regression parameters
#' @param rho the vector of estimated correlations between those connected in the network
#' @param x.matrix nxp matrix of covariates
#' @param y the vector of dependent variable
#' @param adj.matrix adjacency matrix of the network
#' @param missing.vector TRUE FALSE vector of missing values
#' @param link link between xbeta and the mean, identify for normal, logit for logistic, Poissonlog for Poisson with log link, Gammalog for gamma with log link
#' @param W vector of inverse probability weights
#' @param toep number of rhos to estimate 1-4 currently are the options for this
#' @return estimated Variance Covariance matrix

getV<-function(betas, rho, x.matrix, y, adj.matrix, missing.vector, link, W, toep){
  #returns working covariance matrix
  N=dim(x.matrix)[1]
  p=dim(x.matrix)[2]
  A<-matrix(rep(0, N*N), nrow=N)
  CorrY<- A

  cors<-adj.matrix
  cors[adj.matrix==Inf]<-0
  diag(cors)<-0
  cors.1<-cors
  cors.2<-cors%*%cors
  cors.2[cors.1>0]<-0
  cors.2[cors.2>1]<-1
  cors.3<-cors%*%cors%*%cors
  cors.3[cors.1>0 |cors.2>0]<-0
  cors.3[cors.3>0]<-1
  cors.4<-cors%*%cors%*%cors%*%cors
  cors.4[cors.1>0 |cors.2>0 |cors.3>0]<-0
  cors.4[cors.4>0]<-1

  for(tv in 1:toep){
    if(tv==1) cors=cors.1
    if(tv==2) cors=cors.2
    if(tv==3) cors=cors.3
    if(tv==4) cors=cors.4
    diag(cors)<-0
    cors<-cors[missing.vector==FALSE, missing.vector==FALSE]
    CorrY[cors==1]<-rho[tv]
  }
  diag(CorrY)<-1


  if(link=="identity") {
    residual<-y-getMean(betas, x.matrix, link)
    sigma<-sqrt(t(residual)%*%residual/length(residual))
    diag(A)<-sigma
    return(A%*%CorrY%*%A)

  }

  if(link=="logit"){
    pi<-getMean(betas, x.matrix, link)
    diag(A)<-sqrt(pi*(1-pi)) #need to fix
    return(A%*%CorrY%*%A)

  }

  if(link=="Poissonlog"){
    diag(A)= getMean(betas, x.matrix, link)
    return(A%*%CorrY%*%A)
  }

  if(link=="Gammalog"){
    mui<-getMean(betas, x.matrix, link)
    dispersion<-1/(N-p) * sum((y/mui -1)^2)
    diag(A)<-sqrt(mui^2 * dispersion)
    return(A%*%CorrY%*%A)
  #  #to get estimate of dispersion parameter look at page 23 of
  #  #https://www.statistics.ma.tum.de/fileadmin/w00bdb/www/czado/lec8.pdf
  # general case on page 24 of this:
  #  https://www.stat.ncsu.edu/people/davidian/courses/st732/notes/chap11.pdf
  }

}
