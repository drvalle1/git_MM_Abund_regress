# rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_MM_Abund_regress')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data4.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
xmat=data.matrix(read.csv('fake data xmat4.csv',as.is=T))

#useful stuff
ncomm=4
hi=0.999999
lo=0.000001
nspp=ncol(y)
nloc=nrow(y)
npar=ncol(xmat)
txmat=t(xmat)
xtx=txmat%*%xmat

#initial values of parameters
betas=matrix(0,npar,ncomm-1)
omega=matrix(0,nloc,ncomm)
theta=matrix(1/ncomm,nloc,ncomm)
phi=matrix(1/nspp,ncomm,nspp)
sig2=1 

#priors
a.sig2=0.1#nloc*(ncomm-1)/2
b.sig2=0.1#nloc*(ncomm-1)/200

#gibbs details
ngibbs=1000

#objects to store gibbs output
theta.out=matrix(NA,ngibbs,ncomm*nloc)
phi.out=matrix(NA,ngibbs,ncomm*nspp)
betas.out=matrix(NA,ngibbs,npar*(ncomm-1))
sig2.out=matrix(NA,ngibbs,1)
llk=rep(NA,ngibbs)

#MH stuff
accept.output=50
jump1=list(omega=matrix(1,nloc,ncomm))
accept1=list(omega=matrix(0,nloc,ncomm))

#gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #get omega and theta (after integrating out the z's)
  tmp=get.omega.theta(y=y,sig2=sig2,ncomm=ncomm,nloc=nloc,
                      omega=omega,xmat=xmat,betas=betas,phi=phi,jump=jump1$omega)
  omega=tmp$omega
  theta=tmp$theta
  accept1$omega=accept1$omega+tmp$accept
  theta[theta>hi]=hi; theta[theta<lo]=lo
  # omega=omega.true
  # theta=theta.true
  
  #sample z
  tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
  nlk=tmp$nlk
  # nlk=nlk.true
  nks=tmp$nks
  # nks=nks.true

  #sample phi
  phi=rdirichlet1(alpha=nks+1,ncomm=ncomm,nspp=nspp)
  phi[phi>hi]=hi; phi[phi<lo]=lo
  phi[,4]=c(0,0,0,0.0852488) #this ensures that estimated=true last community
  # phi=phi.true
  
  #sample regression coefficients betas
  betas=get.betas(xtx=xtx,sig2=sig2,omega=omega,nparam=npar,
                  ncomm=ncomm,txmat=txmat)

  #sample sig2
  sig2=get.sig2(nloc=nloc,ncomm=ncomm,xmat=xmat,betas=betas,
                omega=omega,a.sig2=a.sig2,b.sig2=b.sig2)

  #adaptation of MH algorithm
  if (i%%accept.output==0 & i<1000){
    k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
    accept1=k$accept1
    jump1=k$jump1
  }
  
  #calculate loglikelihood
  prob=theta%*%phi
  prob[prob>hi]=hi; prob[prob<lo]=lo

  #store results  
  llk[i]=sum(y*log(prob))
  theta.out[i,]=theta
  phi.out[i,]=phi
  betas.out[i,]=betas
  sig2.out[i]=sig2
}

plot(llk,type='l',ylim=range(llk,na.rm=T))
plot(sig2.out,type='l')
