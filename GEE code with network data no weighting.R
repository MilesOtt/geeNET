#-----------------------------------------#
#                                         #
# GEE code for networks                   #
#                                         #
# With some example code                  #
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

  while(sum(abs(betas.new-betas.old))>.000001){
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
    
    if(link=="identity"){
      scale<-V[1,1]
      residual<-(y-Mean)
      cors<-dist.matrix
      cors[dist.matrix==Inf]<-0
      rho.new<-1/scale * t(residual) %*%  cors %*% residual/sum(cors)
    }
    if(link=="logit"){
      OR<-NULL
      for(i in 1:(N)){ #change this to be matrix multiplication and not for loops
        for(j in (1):N){
          if(dist.matrix[i,j]>0 & abs(dist.matrix[i,j])<Inf ){
            x.OR<-x.matrix[i,]-x.matrix[j,]
            OR<-c(OR, exp(x.OR%*%t(betas.new)))}
        }
      }
    
      OR.mean<-mean(OR)
      rho.new=.5*(OR.mean-1)/(OR.mean+1) #Still not right.  I'm going to be looking into this more
      }
  }  
  return(c(betas.new,rho.new))
}


#-------------------------------------------------------
#formula
#link
#dist.matrix
#dataset

gee.net<-function(form, link, dist.matrix, dataset){
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

  
  Z.score<-coef/Sandwich.SE
  CI.LB<-coef-1.96*Sandwich.SE
  CI.UB<-coef+1.96*Sandwich.SE
  p.value<-2*(1-pnorm(abs(Z.score)))
  
  summary<-cbind(coef, ModelBased.SE, Sandwich.SE, CI.LB, CI.UB, Z.score, p.value)
  return(list(summary, rho.hat))
  
}

#----------------------------------- Identity Link clustered data try out---------------------------------
require(geepack)
data(dietox)
                 
attach(dietox)

fit.exch <- geeglm (Weight ~ Time+Cu, id =Pig, data = dietox,
                 family=gaussian,corstr="ex", std.err="san.se")


formula<-"Weight ~ Time+Cu"
link.here="identity"
N<-861

distance.matrix<-matrix(rep(Inf, N*N), nrow=N)
for(i in 1:N){
  for(j in 1:N)
    if (Pig[i]==Pig[j])distance.matrix[i,j]<-1
}
diag(distance.matrix)<-0

gee.net(formula, link.here, distance.matrix, data=dietox)


summary(fit.exch)$coef



require(gee)
fit.exch.2 <- gee(formula, family=gaussian,
                   data=dietox, id=Pig, corstr = "exchangeable")

summary(fit.exch.2)$coef
summary(fit.exch.2)$working[1,2] #estimated correlation



#------- Logit Link clustered simulated data try out----------------------------


Pig<-seq(0, 99.75, by=.25)
Pig<-floor(Pig)+1

Feed<-rnorm(400, 0, 3)
b.0<--2
b.1<-1

Pig.logit<-b.0+b.1*Feed +rnorm(400, 0, .5)
Pig.logit[floor(Pig/2)==Pig/2]<-Pig.logit[floor(Pig/2)==Pig/2]+rnorm(200)

p.pig<-exp(Pig.logit)/(1+exp(Pig.logit))
Big<-NULL
for(i in 1:length(Pig)){
  Big<-c(Big, rbinom(1, 1, p.pig[i]))
}

PigDiet<-as.data.frame(cbind(Pig, Big, Feed))

fit.exch <- geeglm (Big ~ Feed, id =Pig, data = PigDiet,
                    family=binomial,corstr="ex", std.err="san.se")


formula<-"Big ~ Feed"
link.here="logit"
N<-dim(PigDiet)[1]

distance.matrix<-matrix(rep(Inf, N*N), nrow=N)
for(i in 1:N){
  for(j in 1:N)
    if (Pig[i]==Pig[j])distance.matrix[i,j]<-1
}
diag(distance.matrix)<-0




summary(fit.exch)$coef



require(gee)
fit.exch.2 <- gee(formula, family=binomial,
                  data=PigDiet, id=Pig, corstr = "exchangeable")

summary(fit.exch)$coef
summary(fit.exch.2)$coef
summary(fit.exch.2)$working[1,2] #estimated correlation

gee.net(formula, link.here, distance.matrix, data=PigDiet)

#------------------------- Logit link clustered data try out----------------

require(geepack)
data(ohio)
ohio.1<-ohio[c(1:100,2001:2100),]
attach(ohio.1)

fit.exch <- geeglm(resp~age, family=binomial(link="logit"),
                   data=ohio.1, id=id, zcor = "exchangeable", std.err="san.se")


formula<-"resp~age"
link="logit"
N<-200

dist.matrix<-matrix(rep(Inf, N*N), nrow=N)
for(i in 1:N){
  for(j in 1:N)
    if (id[i]==id[j])dist.matrix[i,j]<-1
}
diag(dist.matrix)<-0


gee.net(formula, link, dist.matrix, dataset=ohio.1)


summary(fit.exch)$coef



require(gee)
fit.exch.2 <- gee(resp~age, family=binomial(link="logit"),
                  data=ohio.1, id=id, corstr = "exchangeable")

summary(fit.exch.2)$coef
summary(fit.exch.2)$working[1,2] #estimated correlation


#need to fix rho, though it is OK for now


#---------------------- Network data try with UrWeb Data  -----------------

require(foreign)
setwd("C:/Users/Miles/Dropbox/Miles Backup/crystal networks/updated dta files")

pnom<-as.matrix(read.dta("nomination_matrix_129.dta"), nrow=129, ncol=130)
subj<-read.dta("subjectdata.dta")


# setting up network #
w<-matrix(rep(0,129*129), nrow=129, ncol=129)
y<-pnom[1:129,2:130] 

for(i in 1:129){
  for(j in 1:129){
    w[i,j]<-as.numeric(y[i,j])
  }
}
w.un<-symmetrize(w)
Marijuana.ever<-subj$mjever!="No"
Male<-subj$sex=="Male"
frat<-subj$planfrat!="No"
heavy.drinks<-subj$hvydrkdys_smstr
drink.days<-subj$drkdysmo
setwd("C:/Users/Miles/Desktop")

ur.data<-as.data.frame(cbind(Marijuana.ever, Male, frat, heavy.drinks, drink.days))
dist.matrix<-w.un
dist.matrix[w.un==0]<-Inf
diag(dist.matrix)<-0

link="logit"
gee.net(Marijuana.ever~Male+frat+heavy.drinks, link, dist.matrix, dataset=ur.data)

summary(glm(Marijuana.ever~Male+frat+heavy.drinks, 
            family=binomial, data=ur.data))$coef


link="identity"
gee.net(heavy.drinks~Male+frat+drink.days, link, dist.matrix, dataset=ur.data)

summary(glm(heavy.drinks~Male+frat+drink.days, 
            family=gaussian, data=ur.data))$coef

