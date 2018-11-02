#generates multivariate normal random variate with mean 0 and covariance matrix sigma
rmvnorm1=function(n,sigma){
  ev <- eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = T) %*% R
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  which(z < a)
}
#----------------------------------------------------------------------------------------------
#this function calculates the factor to be added in the MH acceptance ratio to correct
#for truncated normal proposal
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#------------------------------------
#generates truncated normal variates based on cumulative normal distribution
tnorm <- function(n,lo,hi,mu,sig){   
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#------------------------------------
#this function generates omega, which is then used to generate the theta matrix
get.omega.theta=function(y,sig2,ncomm,nloc,omega,xmat,betas,phi,jump){
  k=-(1/(2*sig2))
  omega.new=omega.old=omega
  #use tnorm (random walk is bad for very large or very small omega)
  proposed=matrix(tnorm(nloc*ncomm,lo=-6,hi=6,mu=omega.old,sig=jump),nloc,ncomm) 
  correct.tnorm=fix.MH(lo=-6,hi=6,old1=omega.old,new1=proposed,jump=jump)
  proposed[,ncomm]=0
  
  media=xmat%*%betas
  prior.old=k*((omega.old[,-ncomm]-media)^2)
  prior.new=k*((omega.new[,-ncomm]-media)^2)
  
  for (i in 1:(ncomm-1)){
    omega.new[,i]=proposed[,i]
    
    #get corresponding thetas
    eo=exp(omega.old); theta.old=eo/rowSums(eo)
    eo=exp(omega.new); theta.new=eo/rowSums(eo)
    
    #get likel
    prob.old=theta.old%*%phi; 
    prob.new=theta.new%*%phi; 
    llk.old=rowSums(y*log(prob.old))
    llk.new=rowSums(y*log(prob.new))
    
    #accept/reject
    ind=acceptMH(p0=llk.old+prior.old[,i],p1=llk.new+prior.new[,i]+correct.tnorm[,i],
                 x0=omega.old[,i],x1=omega.new[,i])
    omega.old[ind,i]=omega.new[ind,i]
  }
  eo=exp(omega.old); theta=eo/rowSums(eo)
  list(omega=omega.old,theta=theta,accept=omega.old!=omega)
}
#------------------------------------------------------------------------------
#this function samples regression coefficients
get.betas=function(xtx,sig2,omega,nparam,ncomm,txmat){
  prec=(1/sig2)*xtx+diag(1,nparam)
  var1=solve(prec)
  pmedia=(1/sig2)*txmat%*%omega[,-ncomm]
  media=var1%*%pmedia
  t(rmvnorm1(ncomm-1,var1))+media
}
#-----------------------------------------------------------------
#this function generates dirichlet random variables 
#(one random variable for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#------------------------------------------------------------------------------
#this function generates sig2
get.sig2=function(nloc,ncomm,xmat,betas,omega,a.sig2,b.sig2){
  a=(nloc*(ncomm-1)+2*a.sig2)/2
  media=xmat%*%betas
  b=b.sig2+(sum((omega[,-ncomm]-media)^2)/2)
  # 1/tgamma(a=a,b=b,lo=1,hi=10000)
  1/rgamma(1,a,b)
}
#-----------------------------------------------------------------------------------------------
#this function changes jump based on acceptance frequency
#optimal acceptance probability is between 0.2 and 0.4
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.0001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
