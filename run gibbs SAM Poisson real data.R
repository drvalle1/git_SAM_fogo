rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(3)

setwd('U:\\bony\\block experiment\\dados orig')
xmat=read.csv('design matrix.csv',as.is=T)
dat=read.csv('edited data1.csv',as.is=T)
ind=which(colnames(dat)%in%c('plot','ano'))
y1=data.matrix(dat[,-ind])

#basic checking
sum(xmat$plot!=dat$plot)
sum(xmat$plot.id!=dat$plot.id)
sum(xmat$ano!=dat$ano)

#subset just the relevant pieces
plot.id=xmat$plot.id
xmat1=data.matrix(xmat[,c('queima1','queima2','fuel')])

setwd('U:\\GIT_models\\git_SAM_fogo')
source('gibbs SAM Poisson functions.R')
source('gibbs SAM Poisson.R')
sourceCpp('aux1.cpp')

ngibbs=10000
ngroup=30

res=gibbs.SAM.Poisson(y1=y1,xmat1=xmat1,plot.id=plot.id,
                      ngroup=ngroup,ngibbs=ngibbs,nburn=ngibbs/2,
                      alpha.a=0.1,alpha.b=0.1,gamma=0.1,tau2=0.01)

plot(res$llk,type='l')
k=res$z[nrow(res$z),]; table(k)
plot(res$theta[nrow(res$theta),],type='h')

setwd('U:\\bony\\block experiment\\results')
nomes=paste0(c('betas','z','alpha','theta','llk'),'.csv')
write.csv(res$betas,nomes[1],row.names=F)
write.csv(res$z,nomes[2],row.names=F)
write.csv(res$alpha,nomes[3],row.names=F)
write.csv(res$theta,nomes[4],row.names=F)
write.csv(res$llk,nomes[5],row.names=F)
