
# a convenience function for installing (if necessary) and loading packages 
install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# calculate Cohen's d 
cd = function(m1,sd1,n1,m2,sd2,n2){
  Swithin<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  eff<-(m1-m2)/Swithin
  vd<-((n1+n2)/(n1*n2))+(eff^2/(2*(n1+n2)))
  vd.0<-((n1+n2)/(n1*n2))
  return(list(eff=eff,vd=vd,vd.0=vd.0))
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
creatematrix = function(levels=c(),nlevels=NULL){
  if(!is.null(nlevels)&!is.null(levels)) if(nlevels!=length(levels)) stop("argument 'nlevels' is inconsistent with argument 'levels'!")
  if(is.null(nlevels)) nlevels=length(levels)
  return(rbind(rep(0,nlevels-1),diag(nlevels-1)))
}


# R2 for metaanalytic models
R2 = function(fit=NULL, fit0=NULL){
  if(is.null(fit$formula.mods)) stop("No moderators found in the model of interest!")
  if(paste(fit$formula.mods,collapse="")=="~1") stop("No moderators found in the model of interest!")
  if(class(fit)[1] == "rma.uni"){
    if(is.null(fit0)) return(fit$R2)
    if(!is.null(fit0)) if(is.null(fit0$formula.mods)) return(fit$R2)
    if(!is.null(fit0)) if(!is.null(fit0$formula.mods)) return(fit$R2-fit0$R2)
  }
  if(class(fit)[1] == "rma.mv"){
    if(is.null(fit0)) fit0 = update(fit,mods=~1)
    return(1 - sum(fit$sigma2)/sum(fit0$sigma2))
  }
}

