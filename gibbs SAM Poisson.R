gibbs.SAM.Poisson=function(y1,plot.id,xmat1,alpha.a,alpha.b,gamma,tau2,ngroup,ngibbs,nburn){
  #basic settings
  nspp=ncol(y1)
  nplot=length(unique(plot.id))
  nparam=ncol(xmat1)
  
  #useful summary statistics
  y.ss=matrix(NA,nplot,nspp)
  for (i in 1:nplot){
    cond=plot.id==i
    y.ss[i,]=apply(y1[cond,],2,sum)
  }
  
  #starting values for parameters
  alpha=matrix(mean(y1),nplot,nspp)
  betas=matrix(0,nparam,ngroup)
  theta=rep(1/ngroup,ngroup)
  
  #group assignment
  tmp=rmultinom(nspp,size=1,prob=rep(1/ngroup,ngroup))
  z=apply(tmp==1,2,which)
  
  #stuff for gibbs sampler
  vec.betas=matrix(NA,ngibbs,nparam*ngroup)
  vec.z=matrix(NA,ngibbs,nspp)
  vec.alpha=matrix(NA,ngibbs,nplot*nspp)
  vec.theta=matrix(NA,ngibbs,ngroup)
  vec.llk=matrix(NA,ngibbs,1)
  jump1=list(betas=matrix(1,nparam,ngroup))
  accept1=list(betas=matrix(0,nparam,ngroup))
  
  #start MCMC
  accept.output=50
  for (i in 1:ngibbs){
    print(i)
    tmp=update.betas(betas=betas,alpha=alpha,z=z,
                     plot.id=plot.id,xmat1=xmat1,y1=y1,
                     nparam=nparam,ngroup=ngroup,jump=jump1$betas,tau2=tau2)
    betas=tmp$betas
    accept1$betas=accept1$betas+tmp$accept
    
    z=update.z(xmat1=xmat1,y1=y1,
               betas=betas,alpha=alpha,log.theta=log(theta),
               plot.id=plot.id,ngroup=ngroup,nspp=nspp)
    
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
    
    #calculate llk
    alpha1=alpha[plot.id,]
    emult=exp(xmat1%*%betas)[,z]
    llk=sum(dpois(y1,lambda=alpha1*emult,log=T))
    
    #store results
    vec.betas[i,]=betas
    vec.z[i,]=z
    vec.alpha[i,]=alpha
    vec.theta[i,]=theta
    vec.llk[i,]=llk
  }
  seq1=nburn:ngibbs
  list(betas=vec.betas[seq1,],
       z=vec.z[seq1,],
       alpha=vec.alpha[seq1,],
       theta=vec.theta[seq1,],
       llk=vec.llk[seq1])
}

