rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\bony\\block experiment\\dados orig')
dat=read.csv('edited data.csv',as.is=T)
xmat=read.csv('design matrix.csv',as.is=T)

#just to be sure
sum(dat$plot!=xmat$plot)
sum(dat$ano!=xmat$ano)

#subset just the relevant pieces
ind=which(colnames(dat)%in%c('plot','ano'))
dat1=dat[,-ind]
plot.id=xmat$plot.id
xmat1=data.matrix(xmat[,c('n.queima','fuel')])

#basic settings
nspp=ncol(dat1)
nplot=length(unique(plot.id))
nparam=ncol(xmat1)
  
#parameters
alpha.true=alpha=matrix(rnorm(nspp*nplot,mean=1,sd=1)^2,nplot,nspp)
range(alpha)
ngroup=5
betas.true=betas=matrix(rnorm(ngroup*nparam),nparam,ngroup)

#group assignment
tmp=rmultinom(nspp,size=1,prob=rep(1/ngroup,ngroup))
z.true=z=apply(tmp==1,2,which)

#get mean
pmedia=exp(xmat1%*%betas[,z])
alpha1=alpha[plot.id,]
media=alpha1*pmedia

#generate data
nobs=nrow(dat)
tmp=rpois(nobs*nspp,lambda=media)
y=matrix(tmp,nobs,nspp)
colnames(y)=paste0('spp',1:nspp)
ind=which(colnames(dat)%in%c('plot','ano'))
y1=cbind(dat[,ind],y)

setwd('U:\\GIT_models\\git_SAM_fogo')
write.csv(y1,'simulated data.csv',row.names=F)
