tmp=data.frame(estim=res$z[nrow(res$z),],true=z.true)
table(tmp)
ordem=c(17,20,25,28,20,17,23,17,3,20)

ngroup=30
betas=matrix(res$betas[nrow(res$betas),],ncol(xmat1),ngroup)
rango=range(c(betas.true,betas[,ordem]))
plot(betas.true,betas[,ordem],xlim=rango,ylim=rango)
lines(rango,rango)

alpha=matrix(res$alpha[nrow(res$alpha),],length(unique(plot.id)),ncol(y1)-2)
rango=range(c(alpha,alpha.true))
plot(alpha.true,alpha,xlim=rango,ylim=rango)
lines(rango,rango)

plot(res$theta[nrow(res$theta),],type='h')