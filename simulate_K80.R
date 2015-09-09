library(phyclust)
mutate_K80=function(curbase, kappa ){
  potential_site_v=c(0,1,2,3)[-(curbase+1)]
  totalrate=2+kappa
  if(curbase==0){
    pick=sample(potential_site_v,1,prob=c(kappa/totalrate,1/totalrate,1/totalrate))
  }
  if(curbase==1){
    pick=sample(potential_site_v,1,prob=c(kappa/totalrate,1/totalrate,1/totalrate))
  }
  if(curbase==2){
    pick=sample(potential_site_v,1,prob=c(1/totalrate,1/totalrate,kappa/totalrate))
  }
  if(curbase==3){
    pick=sample(potential_site_v,1,prob=c(1/totalrate,1/totalrate,kappa/totalrate))
  }
  return(pick)
}
evolve_K80=function(base,time,kappa){
  Tt=time
  curbase=base
  cur=rexp(1,rate=2+kappa)
  count=0
  Tt=Tt-cur
  while(Tt>0){
    curbase=mutate_K80(curbase=curbase,kappa=kappa)
    cur=rexp(1,rate=2+kappa)
    Tt=Tt-cur
  }
  return(curbase)
}
evolve_sequence_K80=function(seq,time,kappa){
  evolved_seq=sapply(seq,FUN= function(y) evolve_K80(base=y,time=time,kappa=kappa))
  return(evolved_seq)
}
#generate founders and evolve sequences
evolve_sequence_diff_K80=function(numK,Kprop,timeK,kappaK,foundersT,numSample,nSite){
  seq=gen.seq.HKY(foundersT,pi=c(0.25,0.25,0.25,0.25),kappa=mean(kappaK),L=nSite,rate.scale=10)
  data=read.seqgen(seq)
  X=matrix(nrow=numSample,ncol=nSite)
  
  for(j in 1:numSample){
    pick=sample(1:numK,size=1,prob=Kprop)
    X[j,]=evolve_sequence_K80(seq=data$org[pick,],time=timeK[pick],kappa=kappaK[pick])
  }
  
  return(X)
}

truenumK=5
founders=gen.star.tree(N=truenumK,total.height=5)
kpr=c(0.1,0.2,0.3,0.3,0.1)
ktime=c(0.8,0.15,0.4,0.3,1)
kkappa=c(0.5,1.3,1.2,1.1,1.5)
X=evolve_sequence_diff_K80(numK=truenumK,Kprop=kpr,timeK=ktime,kappaK=kkappa,foundersT=founders,numSample=200,nSite=100)
expinfo=list(proportions=kpr,times=ktime,kappas=kkappa,nsam=nrow(X),nsite=ncol(X))
EMC.2=.EMC
EMC.2$init.procedure <- .init.procedure[2]
EMC.2$init.method=.init.method[1]
EMC.2$substitution.model="K80"
EMC.2$edist.model="D_JC69"
EMC.2$identifier="EV"
maxnumK=10
bic=rep(0,maxnumK)

library(doParallel)
cl=makeCluster(maxnumK)
registerDoParallel(cl)
clock=proc.time()

out=foreach(i=1:maxnumK) %dopar%{
  library(phyclust)
  o=phyclust(X,i,EMC=EMC.2)
  o
}
clock=proc.time()-clock
save(out,clock,expinfo,file="exp313.RData")
