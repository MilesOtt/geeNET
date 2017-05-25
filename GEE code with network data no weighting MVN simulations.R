
#-----------------------------------------#
#                                         #
# GEE code for networks                   #
#                                         #
# With some simulations                   #
#                                         #
# Work in Progress                        #
#                                         #
# Miles Ott miles_ott@alumni.brown.edu    #
#                                         #
#-----------------------------------------#

#betas = p vector of regression parameters
#x.matrix= nxp matrix of covariates
#link = link between xbeta and the mean

getD<-function(betas, x.matrix, link){ 
  #returns matrix of partial derivatives dmu/dbeta
  
  if(link=="identity") return(x.matrix)
  
  if(link=="logit"){
    
     p=length(betas)	
     D<-x.matrix
     logit<-x.matrix %*% betas
     pi<-exp(logit)/((1+exp(logit))^2)
     for(i in 1:p){
       D[,i]<-x.matrix[,i]*pi
     }	
    return(D)
  }
}

#-------------------------------------------------------------------

#betas = p vector of regression parameters
#rho=estimated correlation between those connected in the network
#x.matrix= nxp matrix of covariates
#dist.matrix is some function of the geodesic distances in the network
#link = link between xbeta and the mean


getV<-function(betas, rho, x.matrix, y, dist.matrix, link){
  #returns working covariance matrix
  N=dim(x.matrix)[1]
  A<-matrix(rep(0, N*N), nrow=N)
  CorrY<-rho^dist.matrix
  CorrY[dist.matrix==Inf]<-0
  
  if(link=="identity") {
    residual<-y-getMean(betas, x.matrix, link)
    sigma<-sqrt(t(residual)%*%residual/length(residual))
    diag(A)<-sigma
    return(A%*%CorrY%*%A)
    
  }
 
  if(link=="logit"){
    pi<-getMean(betas, x.matrix, link)
    diag(A)<-sqrt(pi*(1-pi))
    return(A%*%CorrY%*%A)
    
  }
  
}  
  

#-------------------------------------------------------------------

#betas = p vector of regression parameters
#x.matrix= nxp matrix of covariates
#link = link between xbeta and the mean


getMean<-function(betas, x.matrix, link){
  #returns the mean given the betas and x.matrix
  if(link =="identity")return (x.matrix %*% betas)
  if(link =="logit") {
    logit<-x.matrix %*% betas
    pi<-exp(logit)/(1+exp(logit))
    return(pi)  
  } 
}


#-------------------------------------------------------------------

#to.optim a p+1 vector with initial values for betas and rho
#dist.matrix is some function of the geodesic distances in the network
#x.matrix= nxp matrix of covariates
#link = link between xbeta and the mean
#y = outcome variable

quasi.score<-function(to.optim, dist.matrix, x.matrix, link, y){
  #returns betas and rho
  
  N=dim(dist.matrix)[1]
  p<-length(to.optim)-1
  betas.old<-to.optim[-(p+1)]
  rho.old<-to.optim[(p+1)]
  rho.new<-rho.old
  betas.new<-betas.old+.00001
  
  numit=0
  while(sum(abs(betas.new-betas.old))>.000001 & numit<25){
    numit=numit+1
    rho.old<-as.numeric(rho.new)
    betas.old<-betas.new
    
    #step 1 estimate beta
    V<-getV(betas=betas.old, rho=rho.old, x.matrix=x.matrix, y, dist.matrix=dist.matrix, link=link)
    D<-getD(betas=betas.old, x.matrix=x.matrix, link=link)
    Mean<-getMean(betas=betas.old, x.matrix=x.matrix, link=link)
   
    betas.new<-betas.old+solve(t(D)%*%solve(V)%*%D)%*%t(D)%*%solve(V)%*%(y-Mean) 
   
    #step 2 given new estimate of beta, estimate rho
    Mean<-getMean(betas=betas.new, x.matrix=x.matrix, link=link)
    V<-getV(betas=betas.new, rho=rho.old, x.matrix=x.matrix, y, dist.matrix=dist.matrix, link=link)
    residual<-(y-Mean)
    cors<-dist.matrix

    
    if(link=="identity"){
      cors[dist.matrix==Inf]<-0
      cors[lower.tri(cors)]<-0
      scale<-V[1,1]
      rho.new<-1/scale * t(residual) %*%  cors %*% residual/sum(cors)
    }
    if(link=="logit"){
      #------Using Pearson's Correlation which is constrained by bounds different from -1,1
      Vars<-sqrt(diag(V)) 
      cors<-dist.matrix
      cors[dist.matrix==Inf]<-0
      diag(cors)<-0
      rho.new<- (t(residual)/t(Vars))%*%  cors %*% (residual/Vars)/sum(cors)
      }
  }  
  

  return(c(betas.new,rho.new))
  
  
}


#-------------------------------------------------------
#formula
#link
#dist.matrix
#dataset

gee.net<-function(form, link, dist.matrix, dataset, SE="sandwich"){
  #returns coef, SE, Z, CI, p-values
  
  fam="gaussian"
  if (link=="logit") fam="binomial"
  glm.fit<-glm(form, data=dataset, family=fam)
  
  initialbetas<-glm.fit$coef
  initialvalues=c(initialbetas, .001) 
  
  y<-glm.fit$y
  
  n<-length(y)
  p<-length(names(glm.fit$coefficients))
  x.matrix<-as.matrix(dataset[, names(dataset)%in% (names(glm.fit$coefficients))[2:p]])
  x.matrix<-cbind(rep(1,n), x.matrix)
    

  a<-quasi.score(initialvalues, dist.matrix, x.matrix, link, y) #get betas, rho
  
  p=length(a)-1
  coef<-a[-(p+1)]
  rho.hat<-a[(p+1)]
  
  
  D<-getD(coef, x.matrix, link)
  V<-getV(coef, rho.hat, x.matrix, y, dist.matrix, link)
  Mean<-getMean(coef, x.matrix, link)
  

  E<-((y-Mean)%*%t(y-Mean)) 
  E[dist.matrix==Inf]<-0 #make sure this is puts 0 values where there is no correlation
  Vi<-solve(V)
  Dt<-t(D)
  A<-solve(Dt%*%Vi%*%D)
  
  sandwich<-A%*%(Dt%*%Vi%*%E%*%Vi%*%D)%*%A
  
  Sandwich.SE<-sqrt(diag(sandwich))
  ModelBased.SE<-sqrt(diag(A))


  if(SE=="sandwich"){ 
    Z.Stat<-coef/Sandwich.SE
    CI.LB<-coef-1.96*Sandwich.SE
    CI.UB<-coef+1.96*Sandwich.SE
    Standard.Error<-Sandwich.SE
  }
  if (SE=="model.based"){
    Z.Stat<-coef/ModelBased.SE
    CI.LB<-coef-1.96*ModelBased.SE
    CI.UB<-coef+1.96*ModelBased.SE
    Standard.Error<-ModelBased.SE
  }
  p.value<-2*(1-pnorm(abs(Z.Stat)))
  print(rho.hat)
  summary<-cbind(coef, Standard.Error, CI.LB, CI.UB, Z.Stat, p.value)
  return(summary)
  
}

#--------------------  Simulations  -----------------------------------------
require(sna)
require(MASS)


#-----Initializing  
n.sims<-1000


beta1.gee<-rep(0,100)
beta1.se.gee<-rep(0,100)

beta1.glm<-rep(0,100)
beta1.se.glm<-rep(0,100)


#------ starting the simulation


for ( i in 1:n.sims){

  #-----Generating a Random network

  n=200
  r.net<-rgnm(1, nv=n, m=10*n, mode = "graph")

  distance.matrix<-r.net
  distance.matrix[r.net==0]<-Inf
  diag(distance.matrix)<-0

  #-----Assigning Covariate Values
  Intercept<-rep(1,n)
  Age<-runif(n, 20, 70)
  Female<-rbinom(n, 1, .5)

  x.matrix<-cbind(Intercept, Age, Female)
  #-----Parameter values to be estimated
  beta0= 1
  beta1=.5
  beta2= -2

  betas<-as.vector(c(beta0, beta1, beta2))

  #-----Assigning Outcome variable values
  rho=.1
  cor.mat=rho^distance.matrix
  sigma<-sqrt(400)
  vc.mat<-(sigma^2*diag(n))%*%cor.mat


  Alcohol<- mvrnorm(n = 1, mu=(x.matrix)%*%(betas), 
                  Sigma=vc.mat, tol = .01)
  print(i)
  
  Alch.data<-as.data.frame(cbind(Alcohol, x.matrix))  

  fitglm<-summary(glm(Alcohol~Age+Female, data=Alch.data))$coef
  fitgee<-gee.net(Alcohol~Age+Female, "identity", distance.matrix, data=Alch.data, SE="model.based")
  

    beta1.glm[i]<-fitglm[2,1]

    beta1.gee[i]<-fitgee[2,1]
      
    beta1.se.glm[i]<-fitglm[2,2]
    
    beta1.se.gee[i]<-fitgee[2,2]
}

mean(beta1< beta1.glm+1.96*beta1.se.glm & beta1> beta1.glm-1.96*beta1.se.glm)
mean(beta1< beta1.gee+1.96*beta1.se.gee & beta1> beta1.gee-1.96*beta1.se.gee, na.rm=TRUE)




