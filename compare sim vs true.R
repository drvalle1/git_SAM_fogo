tmp=data.frame(estim=param$z,true=z.true)
table(tmp)
ordem=c(2,5,1,3,4)

plot(betas.true,param$betas[,ordem])
plot(alpha.true,param$alpha)
plot(param$theta,type='h')