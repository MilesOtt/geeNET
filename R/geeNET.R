####

#' geeNET
#'
#' Fits GEE for data with Network Correlation structures
#'
#' @param form the formula for the model to be fit, like in the lm function
#' @param link link between xbeta and the mean, identity for normal, logit for logistic, Poissonlog for Poisson with log link,  Gammalog for gamma with log link
#' @param adj.matrix adjacency matrix of the network
#' @param dataset name of the dataset used for the model
#' @param SE choose either the sandwich or model.based SE estimator
#' @param wt is the vector of inverse probability weights, defaults to NULL which will give all observations equal weights
#' @param toep number of rhos to estimate 1-4 currently are the options for this
#' @return standard regression table output
#' @export
#'


geeNET<-function(form, link, adj.matrix, dataset, SE="sandwich", wt=NULL, toep=1){
  try(if(toep > 4) stop("toep must be less than 5"))
  try(if(toep<1) stop("toep must be at least 1"))
  require(sna)
  #returns coef, SE, Z, CI, p-values



  #------------------ trying to create indicator variables for factors-----#

  #-----------end trying to create indicator variables for factors --------#


  #setting the family

  fam="gaussian"
  if (link=="logit") fam="binomial"
  if (link=="Gammalog") fam=Gamma(link='log')
  if (link=="Poissonlog") fam=poisson(link='log')


  #getting initial beta values

  glm.fit<-glm(form, data=dataset, family=fam)

  initialbetas<-glm.fit$coef
  initialvalues=c(initialbetas, rep(.001, toep))


  #makes identity matrix of dimension of the data set

  W<-diag(dim(dataset)[1])

  if(!is.null(wt))
    diag(W)<-wt

  n<-dim(dataset)[1]
  p<-length(names(glm.fit$coefficients))

  #making xmatrix be in the order of the function call

  x.matrix<-NULL
  for(pp in 1:(p-1)){
    x.matrix<-cbind(x.matrix, as.matrix(dataset[,names(dataset)%in% (names(glm.fit$coefficients))[1+pp]]))
  }

  x.matrix<-cbind(rep(1,n), x.matrix)

  #dealing with missing data

  y.full<-dataset[,names(dataset)==names(glm.fit$model)[1]]
  full.matrix<-cbind(y.full, x.matrix)
  missing.vector<-(is.na(rowSums(cbind(full.matrix, wt)))) #dealing with missing weights
  W<-W[missing.vector==FALSE, missing.vector==FALSE]
  x.matrix<-x.matrix[missing.vector==FALSE,]
  y<-y.full[missing.vector==FALSE]

  #get the coefficients and correlations with the quasi score function
  coefs.cors<-quasi.score(initialvalues, adj.matrix, x.matrix, missing.vector, link, y, W, toep) #get betas, rho

  p=length(coefs.cors)-toep
  coef<-coefs.cors[-((p+1): (p+toep))]
  rho.hat<-coefs.cors[((p+1): (p+toep))]

  # calculate the D and V matrix with final correlation and coefficient estimates
  D<-getD(coef, x.matrix, link)
  V<-getV(coef, rho.hat, x.matrix, y, adj.matrix, missing.vector, link, W, toep)
  Mean<-getMean(coef, x.matrix, link)

  E<-((y-Mean)%*%t(y-Mean))

  #------make sure this is puts 0 values where there is no correlation
  cors<-adj.matrix
  cors[adj.matrix==Inf]<-0
  diag(cors)<-0
  cors.1<-cors
  cors.2<-cors%*%cors
  cors.2[cors.2>1]<-1
  cors.3<-cors%*%cors%*%cors
  cors.3[cors.3>0]<-1
  cors.4<-cors%*%cors%*%cors%*%cors
  cors.4[cors.4>0]<-1
  if(toep==1 )cors=cors.1
  if(toep==2) cors=cors.2+cors.1
  if(toep==3) cors=cors.3+cors.1+cors.2
  if(toep==4) cors=cors.4+cors.1+cors.2+cors.3
  diag(cors)<-1
  cors<-cors[missing.vector==FALSE, missing.vector==FALSE]
  E[cors<1]<-0
  #------------------------------------------------------------#
  Vi<-solve(V)
  Dt<-t(D)

  A<-solve(Dt%*%W%*%Vi%*%D)

  #cacluate the SEs
  ModelBased.SE<-sqrt(diag(A))

  sandwich<-A%*%(Dt%*%Vi%*%W%*%E%*%W%*%Vi%*%D)%*%A

  Sandwich.SE<-sqrt(diag(sandwich))


  #calculate the Z score and CI
  if(SE=="sandwich"){
    Z.Stat<-coef/Sandwich.SE
    CI.LB<-coef-1.959964*Sandwich.SE
    CI.UB<-coef+1.959964*Sandwich.SE
    Standard.Error<-Sandwich.SE# *(n)/(n-p-3)
  }
  if (SE=="model.based"){
    Z.Stat<-coef/ModelBased.SE
    CI.LB<-coef-1.959964*ModelBased.SE
    CI.UB<-coef+1.959964*ModelBased.SE
    Standard.Error<-ModelBased.SE
  }

  Estimate<-coef #renaming this so that it is labeled better in the output



  #putting things together for the function output
  Coefficients<-names(glm.fit$coef)
  p.value<-2*(1-pnorm(abs(Z.Stat)))
  missing.value.statement<-"There were no missing values"
  if (sum(missing.vector)>0)
    missing.value.statement<-paste("There were", sum(missing.vector), "cases with missing values that were deleted." )
  summary<-as.data.frame(cbind(Coefficients, round(cbind(Estimate, Standard.Error, CI.LB, CI.UB, Z.Stat),5), p.value))
  return(list(summary, rho.hat, SE, missing.value.statement))

}
