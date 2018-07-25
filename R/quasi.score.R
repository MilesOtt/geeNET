#' quasi.score
#'
#' Fits GEE and returns betas and rho vector
#' @param to.optim a p+1 vector with initial values for betas and rho
#' @param adj.matrix adjacency matrix of the network
#' @param x.matrix nxp matrix of covariates
#' @param missing.vector TRUE FALSE vector of missing values
#' @param link link between xbeta and the mean, identify for normal, logit for logistic, Poissonlog for Poisson with log link, Gammalog for gamma with log link
#' @param y the vector of dependent variable
#' @param W vector of inverse probability weights
#' @param toep number of rhos to estimate 1-4 currently are the options for this
#' @return betas and rho as a vector

quasi.score<-function(to.optim, adj.matrix, x.matrix, missing.vector, link, y, W, toep){
  #returns betas and rho


  p<-length(to.optim)-toep
  betas.old<-to.optim[-((p+1):(p+toep))]
  rho.old<-to.optim[((p+1): (p+toep))]
  rho.new<-rho.old
  betas.new<-betas.old
  start=TRUE
  iterations=1
  while((sum(abs(betas.new-betas.old))>.0001 | start==TRUE ) & iterations<900){
   # print(c(iterations, betas.new))
    iterations<-iterations+1
    start=FALSE
    rho.old<-as.numeric(rho.new)
    betas.old<-betas.new

    #step 1 estimate beta
    V<-getV(betas=betas.old, rho=rho.old, x.matrix=x.matrix, y, adj.matrix=adj.matrix, missing.vector=missing.vector, link=link, W,toep)
    D<-getD(betas=betas.old, x.matrix=x.matrix, link=link)
    Mean<-getMean(betas=betas.old, x.matrix=x.matrix, link=link)


    betas.new<-betas.old+solve(t(D)%*%W%*%solve(V)%*%D)%*%(t(D)%*% solve(V)%*%W%*%(y-Mean))




    #step 2 given new estimate of beta, estimate rho
    Mean<-getMean(betas=betas.new, x.matrix=x.matrix, link=link)
    V<-getV(betas=betas.new, rho=rho.old, x.matrix=x.matrix, y, adj.matrix=adj.matrix, missing.vector, link=link, W, toep)
    residual<-(y-Mean)
    if(link=="identity"){ #****************************
      scale<-V[1,1]
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
      rho.new<-NULL
      for (t in 1: toep){
        if(t==1 )cors=cors.1
        if(t==2) cors=cors.2
        if(t==3) cors=cors.3
        if(t==4) cors=cors.4
        diag(cors)<-0
        cors<-cors[missing.vector==FALSE, missing.vector==FALSE]
        rho.new=c(rho.new, 1/scale * t(residual) %*%cors %*%(residual)/sum(cors))
      }
    }
    if(link %in% c("logit", "Gammalog", "Poissonlog")){ #********************************
      Vars<-sqrt(diag(V))
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
      rho.new<-NULL
      for (t in 1: toep){
        if(t==1 )cors=cors.1
        if(t==2) cors=cors.2
        if(t==3) cors=cors.3
        if(t==4) cors=cors.4
        diag(cors)<-0
        cors<-cors[missing.vector==FALSE, missing.vector==FALSE]
        rho.new<- c(rho.new, (t(residual)/t(Vars))%*%  cors %*% (residual/Vars)/sum(cors))
      }
    }

  }
  return(c(betas.new,rho.new))
}
