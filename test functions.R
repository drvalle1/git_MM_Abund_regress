rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
sourceCpp('aux1.cpp')

#test multinomial function
n=1000
ncomm=10
z=runif(n)
rmultinom1(z,rep(1/ncomm,ncomm),ncomm)
rmultinom1(z,c(0.9,0.1,rep(0,ncomm-2)),ncomm)
sum(z<0.9); sum(z<1 & z>0.9)
rmultinom1(z,c(0.1,0.9,rep(0,ncomm-2)),ncomm)
sum(z<0.1); sum(z<1 & z>0.1)

rmultinom1(z,c(0.1,0.5,0.4,rep(0,ncomm-3)),ncomm)
sum(z<0.1); sum(z<0.6 & z>0.1); sum(z<1 & z>0.6)

#test ngreater
nloc=1000; ncommun=10
nlk=matrix(sample(1:1000,nloc*ncommun,replace=T),nloc,ncommun)
ngreater1=ngreater(nlk,nloc,ncommun)

tmp=apply(nlk[,ncommun:1],1,cumsum)
tmp1=t(tmp)
tmp2=tmp1[,ncommun:1]

sum(ngreater1!=tmp2)