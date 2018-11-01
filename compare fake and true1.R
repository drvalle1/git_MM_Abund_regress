table(gamma.out)

boxplot(theta)

ind=1:4
ind=c(3,2,1,4)
theta1=theta[,ind]
par(mfrow=c(2,2))
for (i in 1:4) plot(theta.true[,i],theta1[,i])

plot(phi.true,phi[ind,])

betas1=betas[,3:1]
plot(betas.true,betas1)