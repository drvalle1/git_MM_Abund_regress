rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_MM_Abund_regress')
source('MM.abund.regression function.R')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data4.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
xmat=data.matrix(read.csv('fake data xmat4.csv',as.is=T))

#run Gibbs sampler
nloc=nrow(y)
ncomm=4
ngibbs=1000
nburn=ngibbs/2
a.sig2=10#nloc*(ncomm-1)/2
b.sig2=1#nloc*(ncomm-1)/20
a.sig2/b.sig2 #precision of 10

res=MM.abund.regression(y=y,xmat=xmat,
                        ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,
                        a.sig2=a.sig2,b.sig2=b.sig2)
  
plot(res$llk,type='l')
plot(res$sig2,type='l')

tmp=colMeans(res$theta)
nloc=nrow(y)
theta=matrix(tmp,nloc,ncomm)
plot(NA,NA,xlim=c(1,nloc),ylim=c(0,1))
for (i in 1:ncomm){
  lines(theta[,i],col=i)
}


