
##############################################################

# empty workspace
rm(list=ls())
# load packages
library(readxl)
library(metafor)
library(ggplot2)
library(clubSandwich)

# load custom functions
source("custom functions.R")

##############################################################

# Import and prepare data; calculate Cohen's d

d = data.frame(read_excel("MetaMathGenderData_20240503.xlsx"))
d = d[d$Filter=="1",]
d$GGI = as.numeric(d$GGI)
d$PublicationYear = as.numeric(d$PublicationYear)
d$PublicationYear2010 = as.numeric(d$PublicationYear)-2010
d$MeanAge_months = round(as.numeric(d$MeanAge_months),3)
d$MeanAge_years = d$MeanAge_months/12
d$id = 1:nrow(d)

vars = c("male_mean","male_sd","fem_mean","fem_sd")
for(v in vars) d[,v] = as.numeric(d[,v])
x = cd(d$male_mean,d$male_sd,d$male_sample,d$fem_mean,d$fem_sd,d$fem_sample)
d$eff<-x$eff
d$vi<-x$vd
d = d[abs(d$eff)<2&!is.na(d$eff),]
d = d[!is.na(d$ID_SAMPLE),]
d$eff[d$Sign_score==(-1)] = d$eff[d$Sign_score==(-1)]*(-1)
d$N = d$fem_sample+d$male_sample
d$lb = d$eff + sqrt(d$vi)*qt(.025,d$N)
d$ub = d$eff + sqrt(d$vi)*qt(.975,d$N)

# stima N tot
tb = aggregate(d$N,by=list(d$ID_SAMPLE),FUN=mean)
tb$x
median(tb$x)
sum(tb$x)

# cor anno ed età
tb = aggregate(d[,c("PublicationYear","MeanAge_years")],by=list(d$ID_SAMPLE),FUN=mean)
cor.test(tb[,2],tb[,3],use="complete.obs")

########

colSums(table(d$ID_SAMPLE,d$Math_content2)!=0)

colSums(table(d$ID_SAMPLE,d$CountryCode)!=0)

##############################################################

# FIT ALL MODELS (ESPECIALLY FOR MODERATION ANALYSIS)
# jump below to "load workspace" if already done

# setup parallelized computation
library(foreach)
library(doParallel)
numCores = parallel::detectCores()
cl = makeCluster(numCores)
registerDoParallel(cl)

# define all names and moderator formulas of models
names = c("fit0","fitModsFull","fitMods_PYxMCxCA","fitMods_PY_MCxCA","fitMods_MC_PYxCA","fitMods_CA_PYxMC",
          "fitMods_MCxCA","fitMods_PYxCA","fitMods_PYxMC",
          "fitMods_MC_CA","fitMods_PY_CA","fitMods_PY_MC",
          "fitMods_MC","fitMods_PY","fitMods_CA")
formulas = list("~1","~PublicationYear2010+Math_content2+CountryCode","~PublicationYear2010*Math_content2*CountryCode","~PublicationYear2010+Math_content2*CountryCode","~Math_content2+PublicationYear2010*CountryCode","~CountryCode+PublicationYear2010*Math_content2",
                "~Math_content2*CountryCode","~PublicationYear2010*CountryCode","~PublicationYear2010*Math_content2",
                "~Math_content2+CountryCode","~PublicationYear2010+CountryCode","~PublicationYear2010+Math_content2",
                "~Math_content2","~PublicationYear2010","~CountryCode")

V = impute_covariance_matrix(d$vi,cluster=d$ID_SAMPLE,r=0.7)

# parallelized computation of all models with maximum likelihood estimation for model comparisons
res = foreach(i=1:length(formulas),.packages="metafor") %dopar% {
  rma.mv(yi=eff,V=V,random=~1|ID/ID_SAMPLE/id,data=d,mods=as.formula(formulas[[i]]),method="ML")
}
for(i in 1:length(res)) assign(names[i],res[[i]])
rm(res)
save.image("Workspaces/workspaceMetaMathML.RData")
# see AICs and BICs of all models
tab = data.frame(names=names,AIC=NA,BIC=NA)
for(i in 1:nrow(tab)){
  tab$AIC[i] = round(AIC(get(names[i])),2)
  tab$BIC[i] = round(BIC(get(names[i])),2)
}; tab
tab$names[tab$AIC==min(tab$AIC)]
tab$names[tab$BIC==min(tab$BIC)]

# parallelized computation of all models with REML (restricted maximum likelihood) estimation for better coefficients estimates
res = foreach(i=1:length(formulas),.packages="metafor") %dopar% {
  rma.mv(yi=eff,V=V,random=~1|ID/ID_SAMPLE/id,data=d,mods=as.formula(formulas[[i]]))
}
for(i in 1:length(res)) assign(names[i],res[[i]])
rm(res)
save.image("Workspaces/workspaceMetaMathREML.RData")

##############################################################

#### IF ALREADY DONE THE ABOVE: LOAD WORKSPACE
#load("Workspaces/workspaceMetaMathML.RData")
load("Workspaces/workspaceMetaMathREML.RData")

##############################################################

# INTERCEPT ONLY MODEL

summary(fit0)

# forest(fit)
funnel(fit0,xlab="Standardized Mean Difference")
predict(fit0)

#######################

# PLOT MODERATOR ANALYSIS: "PUBLICATION YEAR"

PublicationYear = seq(0,12,.1)
eff = data.frame(predict(fitMods_PY,PublicationYear))
eff$PublicationYear = PublicationYear

# plot
(ggPY = ggplot(d,aes(x=PublicationYear2010+2010,y=eff))+
  theme(text=element_text(size=16))+
  scale_x_continuous(breaks=seq(2010,2022,2))+
  scale_y_continuous(breaks=seq(-2,2,.2))+
  coord_cartesian(ylim=c(-1,1))+
  geom_point(size=0.8,alpha=.5)+
  geom_smooth(formula="y~x",method="loess",alpha=.5,linetype=2)+
  geom_ribbon(data=eff,aes(x=PublicationYear+2010,y=pred,ymin=ci.lb,ymax=ci.ub),size=1,fill="blue",alpha=.2)+
  geom_line(data=eff,aes(x=PublicationYear+2010,y=pred),size=1,color="blue")+
  xlab("Year of publication")+ylab("Estimated SMD")
)
  
#######################

# PLOT MODERATOR ANALYSIS: "MATH CONTENT 2"

Math_content2 = levels(as.factor(d$Math_content2))
eff = data.frame(predict(fitMods_MC,creatematrix(Math_content2)))
eff$Math_content2 = Math_content2

(ggMC = ggplot(eff,aes(x=Math_content2,y=pred))+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=90))+
  geom_hline(yintercept=0,size=1,linetype=2,color="darkgray")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),width=.2,size=0.8)+
  xlab("Math Content")+ylab("Estimated SMD")
)

#######################

# PLOT MODERATOR ANALYSIS: "COUNTRY AREA"

CountryCode = levels(as.factor(d$CountryCode))
eff = data.frame(predict(fitMods_CA,creatematrix(CountryCode)))
eff$CountryCode = CountryCode

(ggCA = ggplot(eff,aes(x=CountryCode,y=pred))+
  theme(text=element_text(size=16),axis.text.x=element_text(angle=90))+
  geom_hline(yintercept=0,size=1,linetype=2,color="darkgray")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),width=.2,size=0.8)+
  xlab("Country area")+ylab("Estimated SMD")
)

##############################################################

# SEPARATE MODERATOR ANALYSIS: "AGE"
# separate analysis is required because age is coded only for part of the dataset

# select only rows where age data is available
dx = d[!is.na(d$MeanAge_months), ]
# fit model with age as moderator
V = impute_covariance_matrix(dx$vi,cluster=dx$ID_SAMPLE,r=0.7)
fitMods_Age = rma.mv(yi=eff,V=V,random=~1|ID/ID_SAMPLE/id,data=dx,mods=~MeanAge_years)
# see model
summary(fitMods_Age)

# estract predictions for plotting
years = seq(4,28,.5)
eff = data.frame(predict(fitMods_Age,years))
eff$years = years

# plot
(ggAge = ggplot(dx,aes(x=MeanAge_years,y=eff))+
    theme(text=element_text(size=16))+
    scale_x_continuous(breaks=seq(0,30,2))+
    scale_y_continuous(breaks=seq(-2,2,.2))+
    coord_cartesian(ylim=c(-1,1))+
    #geom_errorbar(aes(ymin=lb,ymax=ub),alpha=.5)+
    geom_ribbon(data=eff,aes(x=years,y=pred,ymin=ci.lb,ymax=ci.ub),fill="blue",alpha=.2)+
    geom_point(size=2.5,alpha=.5)+
    geom_smooth(formula="y~x",method="loess",alpha=.5,size=1.1,linetype=2)+
    geom_line(data=eff,aes(x=years,y=pred),size=1.1,color="blue")+
    xlab("Age (years)")+ylab("Estimated SMD")
)

##############################################################

# save.image("Workspaces/workspaceMetaMathREML.RData")

##############################################################
##############################################################
##############################################################
##############################################################

