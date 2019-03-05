rm(list=ls(all=TRUE))
set.seed(4)

nloc=1000
nspp=100
ncommun=4

#generate thetas (each community should dominate in at least some place)
tmp1=c(seq(from=2,to=0,length.out=nloc/2),rep(0,nloc*1/2))
tmp2=c(rep(0,nloc/8),seq(from=0,to=2,length.out=nloc/8),seq(from=2,to=0,length.out=nloc/8),rep(0,nloc*5/8))
tmp3=c(rep(0,nloc/4),seq(from=0,to=2,length.out=nloc/8),seq(from=2,to=0,length.out=nloc/8),rep(0,nloc*4/8))
xmat=cbind(1,tmp1,tmp2,tmp3)
colnames(xmat)=c('interc',paste('cov',1:3,sep=''))

betas=matrix(0,ncol(xmat),ncommun-1)
betas[,1]=c(-3,4,0,0)
betas[,2]=c(-3,0,4,0)
betas[,3]=c(-3,0,0,4)
betas.true=betas
sig2=0.3

media=xmat%*%betas; range(media)
eres=cbind(exp(media),1)
theta.pred=eres/rowSums(eres)
plot(NA,NA,xlim=c(0,nloc),ylim=range(theta.pred))
for (i in 1:ncommun) lines(1:nloc,theta.pred[,i],col=i)

image(media)
omega=matrix(NA,nloc,ncommun)
omega[,ncommun]=0
for (i in 1:(ncommun-1)){
  omega[,i]=rnorm(nloc,mean=media[,i],sd=sqrt(sig2))
}
omega.true=omega
range(omega)
e.omega=exp(omega)
soma=rowSums(e.omega)
theta.true=theta=e.omega/matrix(soma,nloc,ncommun)
apply(theta,1,sum)

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate phi  
tmp=matrix(rnorm(ncommun*nspp,mean=0,sd=2),ncommun,nspp)
tmp[tmp<0.1]=0.1
tmp[,1:(2*ncommun)]=cbind(diag(8,ncommun),diag(8,ncommun))
phi=tmp/matrix(rowSums(tmp),ncommun,nspp)
round(phi[,1:20],2)
table(round(phi,2))

unique(rowSums(phi))
phi.true=phi

#generate actual observations y
nl=floor(runif(nloc,min=100,max=200))
nlk=matrix(NA,nloc,ncommun)
nks=matrix(0,ncommun,nspp)
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  nlk[i,]=rmultinom(1,size=nl[i],prob=theta[i,])
  tmp1=rep(0,ncommun)
  for (k in 1:ncommun){
    tmp=rmultinom(1,size=nlk[i,k],prob=phi[k,])
    nks[k,]=nks[k,]+tmp
    tmp1=tmp1+tmp
  }
  y[i,]=tmp1
}
image(y)

#look at stuff to make sure it makes sense
theta.estim=nlk/matrix(nl,nloc,ncommun)
plot(NA,NA,xlim=c(1,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta.estim[,i],col=i)
nlk.true=nlk

phi.estim=nks/matrix(rowSums(nks),ncommun,nspp,)
plot(phi.true,phi.estim)
nks.true=nks

#export results
setwd('U:\\GIT_models\\git_MM_Abund_regress')
nome=paste(c('fake data','fake data xmat'),ncommun,'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome[1])    
write.csv(xmat,nome[2],row.names=F)