rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

setwd('U:\\bony\\block experiment\\dados orig')
xmat=read.csv('design matrix.csv',as.is=T)

#subset just the relevant pieces
plot.id=xmat$plot.id
xmat1=data.matrix(xmat[,c('n.queima','fuel')])

setwd('U:\\GIT_models\\git_SAM_fogo')
source('gibbs SAM Poisson functions.R')
sourceCpp('aux1.cpp')

y=read.csv('simulated data.csv',as.is=T)
ind=which(colnames(y)%in%c('plot','ano'))
y1=data.matrix(y[,-ind])

#basic settings
nspp=ncol(y1)
nplot=length(unique(plot.id))
nparam=ncol(xmat1)

#priors
alpha.b=alpha.a=0.1
gamma=0.1

#useful summary statistics
y.ss=matrix(NA,nplot,nspp)
for (i in 1:nplot){
  cond=plot.id==i
  y.ss[i,]=apply(y1[cond,],2,sum)
}

#starting values for parameters
alpha=matrix(mean(y1),nplot,nspp)
ngroup=30
betas=matrix(0,nparam,ngroup)
theta=rep(1/ngroup,ngroup)

#group assignment
tmp=rmultinom(nspp,size=1,prob=rep(1/ngroup,ngroup))
z=apply(tmp==1,2,which)

#stuff for gibbs sampler
ngibbs=10000
nburn=ngibbs/2
vec.betas=matrix(NA,ngibbs,nparam*ngroup)
vec.z=matrix(NA,ngibbs,nspp)
vec.alpha=matrix(NA,ngibbs,nplot*nspp)
vec.theta=matrix(NA,ngibbs,ngroup)
jump1=list(betas=matrix(1,nparam,ngroup))
accept1=list(betas=matrix(0,nparam,ngroup))
  
#start MCMC
accept.output=50
for (i in 1:ngibbs){
  print(i)
  tmp=update.betas(betas=betas,alpha=alpha,z=z,
                     plot.id=plot.id,xmat1=xmat1,y1=y1,
                     nparam=nparam,ngroup=ngroup,jump=jump1$betas)
  betas=tmp$betas
  accept1$betas=accept1$betas+tmp$accept
  
  z=update.z(xmat1=xmat1,y1=y1,
             betas=betas,alpha=alpha,log.theta=log(theta),
             plot.id=plot.id,ngroup=ngroup)
  
  alpha=update.alpha(y.ss=y.ss,xmat=xmat1,
                     betas=betas,z=z,
                     nplot=nplot,plot.id=plot.id,nspp=nspp,
                     alpha.a=alpha.a,alpha.b=alpha.b)
  
  theta=update.theta(z=z,ngroup=ngroup,gamma=gamma)
  
  #adaption of MH algorithm
  if (i%%accept.output==0 & i<nburn){
    tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
    jump1=tmp$jump1
    accept1=tmp$accept1
  }

  #store results
  vec.betas[i,]=betas
  vec.z[i,]=z
  vec.alpha[i,]=alpha
  vec.theta[i,]=theta
}
param=list(z=z,betas=betas,alpha=alpha,theta=theta)
