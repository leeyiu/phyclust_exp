
library(phyclust)




evolve_sequence_direct_JC69=function(seq,time,mu){
  evolved_seq=sapply(seq,FUN= function(y) evolve_direct_JC69(base=y,time=time,mu=mu))
  return(evolved_seq)
}

evolve_direct_JC69=function(base,time,mu){
  potential_sitev= potential_site_v=c(0,1,2,3)[-(base+1)]
  Tt=time
  alpha=mu
  pick=runif(n=1,min=0,max=1)
  if(pick<0.25+0.75*exp(-Tt*alpha)) {pick=base
  }else {pick=as.integer(runif(n=1,min=1,max=4))
  if(pick==4) pick=3
  pick=potential_site_v[pick]
  }
  return(pick)
  
}



#generate founders
evolve_sequence_direct_diff_JC69=function(numK,Kprop,timeK,muK,foundersT,numSample,nSite){
  seq=gen.seq.HKY(foundersT,pi=c(0.25,0.25,0.25,0.25),kappa=0.5,L=nSite,rate.scale=5)
  data=read.seqgen(seq)
  X=matrix(nrow=numSample,ncol=nSite)
  
  for(j in 1:numSample){
    pick=sample(1:numK,size=1,prob=Kprop)
    X[j,]=evolve_sequence_direct_JC69(seq=data$org[pick,],time=timeK[pick],mu=muK[pick])
  }
  
  return(X)
}


truenumK=2

founders=gen.star.tree(N=truenumK,total.height=5)
kpr=c(0.2,0.8)
ktime=c(2,1)
kmu=c(1.5,2.5)
X=evolve_sequence_direct_diff_JC69(numK=truenumK,Kprop=kpr,timeK=ktime,muK=kmu,foundersT=founders,numSample=400,nSite=500)
expinfo=list(proportions=kpr,times=ktime,mus=kmu,nsam=nrow(X),nsite=ncol(X))
EMC.2=.EMC
EMC.2$init.procedure <- .init.procedure[2]
EMC.2$substitution.model="JC69"
EMC.2$identifier="EV"
maxnumK=10
bic=rep(0,maxnumK)
library(doParallel)
#parallel
cl=makeCluster(maxnumK)
registerDoParallel(cl)
clock=proc.time()

out=foreach(i=1:maxnumK) %dopar%{
  library(phyclust)
  o=phyclust(X,i,EMC=EMC.2)
  o
}
clock=proc.time()-clock
save(out,clock,expinfo,file="exp314.RData")
