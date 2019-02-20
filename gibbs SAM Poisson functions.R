print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#---------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#---------------------------------------------------------------------
update.betas=function(betas,nparam,ngroup,jump,alpha,plot.id,z,xmat1,y1,tau2){
  betas.old=betas.orig=betas
  betas.new=matrix(rnorm(nparam*ngroup,mean=betas.old,sd=jump),nparam,ngroup)  
  alpha1=alpha[plot.id,]
  tmp=data.frame(z=1:ngroup)
  for (i in 1:nparam){
    betas.new=betas.old
    betas.new[i,]=rnorm(ngroup,mean=betas.old[i,],sd=jump[i,])
    
    mult.old=xmat1%*%betas.old
    mult.new=xmat1%*%betas.new
    e.mult.old=exp(mult.old)
    e.mult.new=exp(mult.new)
    
    lprob.old=y1*mult.old[,z]-alpha1*e.mult.old[,z]
    tmp.old=data.frame(z=z,lprob.old=colSums(lprob.old))
    lprob.old1=aggregate(lprob.old~z,data=tmp.old,sum)

    lprob.new=y1*mult.new[,z]-alpha1*e.mult.new[,z]
    tmp.new=data.frame(z=z,lprob.new=colSums(lprob.new))
    lprob.new1=aggregate(lprob.new~z,data=tmp.new,sum)

    #be careful with missing groups    
    fim=merge(lprob.old1,lprob.new1,all=T)
    fim1=merge(fim,tmp,all=T)
    cond=is.na(fim1$lprob.old)
    fim1[cond,c('lprob.old','lprob.new')]=0
    
    #need to add prior
    prior.old=-(1/(2*tau2))*colSums(betas.old^2)
    prior.new=-(1/(2*tau2))*colSums(betas.new^2)
    k=acceptMH(p0=fim1$lprob.old+prior.old,
               p1=fim1$lprob.new+prior.new,
               x0=betas.old[i,],
               x1=betas.new[i,],BLOCK=F)
    betas.old[i,]=k$x
  }
  list(betas=betas.old,accept=betas.old!=betas.orig)
}
#---------------------------------------------------------------------
update.z=function(ngroup,xmat1,betas,plot.id,alpha,y1,log.theta,nspp){
  tmp=exp(xmat1%*%betas)
  alpha1=alpha[plot.id,]
  
  #get log-probabilities
  prob=matrix(NA,nspp,ngroup)
  for (i in 1:ngroup){
    prob[,i]=colSums(dpois(y1,lambda=alpha1*tmp[,i],log=T))+log.theta[i]
  }
  
  #get probabilities
  max1=matrix(apply(prob,1,max),nspp,ngroup)
  prob1=prob-max1
  prob2=exp(prob1)
  soma=matrix(rowSums(prob2),nspp,ngroup)
  prob3=prob2/soma
  
  #sample from multinomial
  z=rmultinom1(prob=prob3,randu=runif(nspp))+1
  z
}
#---------------------------------------------------------------------
update.alpha=function(y.ss,xmat,
                      betas,z,
                      nplot,plot.id,nspp,alpha.a,alpha.b){
  a=y.ss+alpha.a
  tmp=exp(xmat%*%betas)[,z]
  b=matrix(NA,nplot,nspp)
  for (i in 1:nplot){
    cond=plot.id==i
    b[i,]=alpha.b+apply(tmp[cond,],2,sum)
  }
  res=rgamma(nplot*nspp,a,b)
  matrix(res,nplot,nspp)
}
#---------------------------------------------------------------------
update.theta=function(z,ngroup,gamma){
  nk=rep(0,ngroup)  
  tmp=table(z)
  nk[as.numeric(names(tmp))]=tmp
  cumsoma=cumsum(nk[ngroup:1])[ngroup:1]
  vk=c(rbeta(ngroup-1,nk[-ngroup]+1,cumsoma[-1]+gamma),1)
  convertSBtoNormal(v=vk)
}