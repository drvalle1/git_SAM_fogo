rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(3)

setwd('U:\\bony\\block experiment\\dados orig')
xmat=read.csv('design matrix.csv',as.is=T)

#subset just the relevant pieces
plot.id=xmat$plot.id
xmat1=data.matrix(xmat[,c('queima1','queima2','fuel')])

setwd('U:\\GIT_models\\git_SAM_fogo')
source('gibbs SAM Poisson functions.R')
source('gibbs SAM Poisson.R')
sourceCpp('aux1.cpp')

y=read.csv('simulated data.csv',as.is=T)
ind=which(colnames(y)%in%c('plot','ano'))
y1=data.matrix(y[,-ind])

ngibbs=1000
ngroup=30

res=gibbs.SAM.Poisson(y1=y1,xmat1=xmat1,plot.id=plot.id,
                      ngroup=ngroup,ngibbs=ngibbs,nburn=ngibbs/2,
                      alpha.a=0.1,alpha.b=0.1,gamma=0.1,tau2=0.01)

