
# calculate Cohen's d 
cd = function(m1,sd1,n1,m2,sd2,n2){
  Swithin<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  eff<-(m1-m2)/Swithin
  vd<-((n1+n2)/(n1*n2))+(eff^2/(2*(n1+n2)))
  return(list(eff=eff,vd=vd))
}

# combine multiple effect sizes 
combine = function(eff=c(),vi=c(),r=.7){
  mu = mean(eff)
  if(length(vi)==1) v = vi
  if(length(vi)>1){
    Svi = vi
    for(i in 1:length(vi)){
      srvv = 0
      for(j in 1:length(vi)) if(i!=j) srvv = srvv + r*sqrt(vi[i])*sqrt(vi[j])
      Svi[i] = Svi[i] + srvv
    }
    v = sum(Svi) * ((1/length(Svi))^2)
  }
  return(list(mu=mu,v=v))
}

# create matrix of effects 
creatematrix = function(levels=c(),nlevels=NA){
  if(!is.na(nlevels)&!is.null(levels)) if(nlevels!=length(levels)) stop("argument 'nlevels' is inconsistent with argument 'levels'!")
  if(is.na(nlevels)) nlevels=length(levels)
  return(rbind(rep(0,nlevels-1),diag(nlevels-1)))
}

