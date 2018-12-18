betas

boxplot(theta)

ind=1:4
ind=c(3,1,2,4)
theta1=theta[,ind]
par(mfrow=c(2,2))
for (i in 1:4) plot(theta.true[,i],theta1[,i])

plot(phi.true,phi[ind,])

betas1=betas[,ind[1:3]]
plot(betas.true,betas1)

par(mfrow=c(1,1))
plot(NA,NA,xlim=c(1,nloc),ylim=c(0,1))
for (i in 1:4) {
  lines(1:nloc,theta.true[,i],col=i)
  lines(1:nloc,theta1[,i],col=i,lty=3)
}