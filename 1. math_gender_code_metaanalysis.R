
#########################################################

# empty workspace
rm(list=ls())

# Load custom functions
source("_custom functions.R")

# Install and load required packages
required_packages = c("metafor", "clubSandwich", "emmeans", "ggplot2", "gridExtra", "cowplot", "parallel")
install_and_load(required_packages)

# Detect the number of available CPU cores
numCores = detectCores()

#########################################################

#### LOAD DATA AND SEE DESCRIPTIVE STATISTICS
try({ load("workspaceMathGenderMetaAnalysis.RData") })

# Import and prepare data; calculate Cohen's d
dFull = read.table("math_gender_full_dataset.csv",header=T,sep=",")
dFull$eff = as.numeric(dFull$eff)
length(unique(dFull$ID_ARTICLE)); length(unique(dFull$ID_SAMPLE)); nrow(dFull)
d = dFull[!is.na(dFull$GradeM), ]; length(unique(d$ID_ARTICLE)); length(unique(d$ID_SAMPLE)); nrow(d)

#########################################################

#### INTERCEPT-ONLY MODEL (GENERAL MEAN)

# Covariance matrix
V = impute_covariance_matrix(dFull$vi, cluster = dFull$ID_SAMPLE, r = 0.5)

# fit intercept only model
(fit0 = rma.mv(yi = eff, V = V, random = ~1|ID_ARTICLE/ID_SAMPLE/id_effect, data = dFull, control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T))

# forest plot (print to pdf)
pdf("forest_plot.pdf",width=4,height=50)
forest(fit0,slab=paste0(Authors," (sample: '",ID_SAMPLE,"')"),
       efac=.05,
       lwd=.3,
       header=T,
       addfit=T,
       addpred=T)
overall_estimate = fit0$b 
pred_interval = predict(fit0)
polygon(c(overall_estimate, pred_interval$ci.lb, overall_estimate, pred_interval$ci.ub), 
        c(-1-.3, -1, -1+.3, -1), 
        col = "black", border = "black", lwd = .8)
dev.off()

# fit model with and without random variability for sample (fitted with maximum likelihood for model comparability)
fit0sampleML = rma.mv(yi = eff, V = V, random = ~1 | ID_ARTICLE/ID_SAMPLE/id_effect, data = dFull, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
fit0ML = rma.mv(yi = eff, V = V, random = ~1 | ID_ARTICLE/id_effect, data = dFull, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
anova(fit0sampleML,fit0ML) # likleihood ratio test favors model without random variability for sample

# refit intercept-only model with REML for better estimate of random variability
(fit0 = rma.mv(yi = eff, V = V, random = ~1|ID_ARTICLE/id_effect, data = dFull, control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T))

# overall funnel plot
funnel(fit0,level=c(0.90,0.95,0.99),shade=c("white","#CCCCCC",  "#AAAAAA"),back="#DDDDDD",
       cex.lab=1.5,cex.axis=1.5,lwd=2,xlab="Standardized Mean Difference (Cohen's d)",legend=T)

# petpeese test
(fit0Pet = rma.mv(yi = eff, V = V, mods=~sqrt(vi_0), random = ~1|ID_ARTICLE/id_effect, data = dFull, control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T))

#########################################################

#### META-REGRESSION (MODERATORS) MODELS

# Covariance matrix
V = impute_covariance_matrix(d$vi, cluster = d$ID_ARTICLE, r = 0.5)

# starting model
fit0 = rma.mv(yi = eff, V = V, mods=~1, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
fit1 = rma.mv(yi = eff, V = V, mods=~MathContent+GradeM+GeoArea+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
R2(fit1)
# alternative models for main effects / round 1
fit1.1 = rma.mv(yi = eff, V = V, mods=~GradeM+GeoArea+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
fit1.2 = rma.mv(yi = eff, V = V, mods=~MathContent+GeoArea+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
fit1.3 = rma.mv(yi = eff, V = V, mods=~MathContent+GradeM+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
fit1.4 = rma.mv(yi = eff, V = V, mods=~MathContent+GradeM+GeoArea, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
anova(fit1,fit1.1); R2(fit1.1,fit0); R2(fit1,fit0)-R2(fit1.1,fit0); AIC(fit1.1)-AIC(fit1)
anova(fit1,fit1.2); R2(fit1.2,fit0); R2(fit1,fit0)-R2(fit1.2,fit0); AIC(fit1.2)-AIC(fit1)
anova(fit1,fit1.3); R2(fit1.3,fit0); R2(fit1,fit0)-R2(fit1.3,fit0); AIC(fit1.3)-AIC(fit1)
anova(fit1,fit1.4); R2(fit1.4,fit0); R2(fit1,fit0)-R2(fit1.4,fit0); AIC(fit1.4)-AIC(fit1)

# content x grade interaction
fit2 = rma.mv(yi = eff, V = V, mods=~MathContent*GradeM+GeoArea+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
anova(fit2,fit1); R2(fit2,fit0); R2(fit1,fit0)-R2(fit2,fit0); AIC(fit2)-AIC(fit1)
# fit2.1 = rma.mv(yi = eff, V = V, mods=~MathContent*GradeM+GeoArea+GGGI*GradeM, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
# anova(fit2,fit2.1)

# content x grade x gGGGI interaction
# fit2 = rma.mv(yi = eff, V = V, mods=~MathContent*GradeM*GGGI-MathContent:GradeM:GGGI +GeoArea, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
# fit3 = rma.mv(yi = eff, V = V, mods=~MathContent*GradeM*GGGI +GeoArea, random = ~1|ID_ARTICLE/id_effect, data = d, method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
# anova(fit3,fit2)
# Vx = impute_covariance_matrix(d$vi[d$MathContent=="BasicNumeracy"], cluster = d$ID_STUDY[d$MathContent=="BasicNumeracy"], r = 0.5)
# fit2x = rma.mv(yi = eff, V = Vx, mods=~GradeM*GGGI-GradeM:GGGI +GeoArea, random = ~1|ID_ARTICLE/id_effect, data = d[d$MathContent=="BasicNumeracy",], method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
# fit3x = rma.mv(yi = eff, V = Vx, mods=~GradeM*GGGI +GeoArea, random = ~1|ID_ARTICLE/id_effect, data = d[d$MathContent=="BasicNumeracy",], method="ML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
# anova(fit3x,fit2x)

# prepare plot of main effects
fit1 = rma.mv(yi = eff, V = V, mods=~MathContent+GradeM+GeoArea+GGGI, random = ~1|ID_ARTICLE/id_effect, data = d, method="REML", control = list(optimizer = "optimParallel", ncpus = numCores), sparse=T)
ylimits = c(-0.18,0.45)
# math content main effect
empfit1_mc = emmprep(fit1)
(eff_mc = data.frame(emmeans(empfit1_mc,specs="MathContent")))
labelsMC = c("Advanced \n Maths","Basic \n Numeracy", "Broad \n Mathematics", "Computation", "Geometry", "National \n Test")
(gg_mc = ggplot(eff_mc,aes(x=MathContent,y=emmean))+
    geom_hline(yintercept=0,size=1,linetype=2,color="#888888") +
    geom_errorbar(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),width=.2) + 
    geom_point(size=5) + 
    scale_x_discrete(labels=labelsMC) +
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=20),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
    coord_cartesian(ylim=ylimits) +
    ggtitle("Type of content") + ylab("Estimated SMD")
)
# grade main effect
empfit1_gr = emmprep(fit1,at=list(GradeM=seq(0,14,length.out=100)))
eff_gr = data.frame(emmeans(empfit1_gr,specs="GradeM"))
(gg_gr = ggplot(eff_gr,aes(x=GradeM,y=emmean))+
    geom_hline(yintercept=0,size=1,linetype=2,color="#888888") +
    geom_ribbon(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),alpha=.15) + 
    geom_line(size=1.5) + 
    scale_x_continuous(breaks=seq(0,14,2))+
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=20),axis.title.x=element_blank()) +
    coord_cartesian(ylim=ylimits) +
    ggtitle("Grade") + ylab("Estimated SMD")
)
# geographical area main effect
empfit1_ga = emmprep(fit1)
(eff_ga = data.frame(emmeans(empfit1_ga,specs="GeoArea")))
labelsGA = c("Africa","Central \n Europe", "East Asia + \n Oceania", "East Europe + \n Russia", "Middle East", "North America",  "North Europe", "South + Central \n America")
(gg_ga = ggplot(eff_ga,aes(x=GeoArea,y=emmean))+
    geom_hline(yintercept=0,size=1,linetype=2,color="#888888") +
    geom_errorbar(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),width=.2) + 
    geom_point(size=5) + 
    scale_x_discrete(labels=labelsGA) +
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=20),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
    coord_cartesian(ylim=ylimits) +
    ggtitle("Geographical area") + ylab("Estimated SMD")
)
# gGGGI main effect
empfit1_gi = emmprep(fit1,at=list(GGGI=seq(min(d$GGGI),max(d$GGGI),length.out=100)))
eff_gi = data.frame(emmeans(empfit1_gi,specs="GGGI"))
(gg_gi = ggplot(eff_gi,aes(x=GGGI,y=emmean))+
    geom_hline(yintercept=0,size=1,linetype=2,color="#888888") +
    geom_ribbon(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),alpha=.15) + 
    geom_line(size=1.5) + 
    scale_x_continuous(breaks=seq(0,1,0.05)) +
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=20),axis.title.x=element_blank()) +
    coord_cartesian(ylim=ylimits) +
    ggtitle("GGGI") + ylab("Estimated SMD")
)
# plot with 4 panels
plot_grid(gg_mc,gg_gr,gg_ga,gg_gi,nrow=2,align="h")


# analysis of Type of content by Grade interaction
lvl = levels(as.factor(d$MathContent))
res = data.frame(MathContent=lvl,beta=NA,lb=NA,ub=NA,pval=NA)
d$MathContentx = as.factor(d$MathContent)
for(i in 1:length(lvl)){
  d$MathContentx = relevel(d$MathContentx,ref=lvl[i])
  fitx = update(fit2, mods=~MathContentx*GradeM+GeoArea+GGGI)
  res$beta[i] = round(fitx$beta[rownames(fitx$beta)=="GradeM"], 3)
  res$lb[i] = round(fitx$ci.lb[rownames(fitx$beta)=="GradeM"] , 3)
  res$ub[i] = round(fitx$ci.ub[rownames(fitx$beta)=="GradeM"] , 3)
  res$pval[i] = round(fitx$pval[rownames(fitx$beta)=="GradeM"] , 4)
  print(res[i,])
}
# plot of interaction
empfit2_int = emmprep(fit2,at=list(GradeM=seq(0,14,.1)))
eff_int = data.frame(emmeans(empfit2_int,specs="GradeM",by="MathContent"))
(gg_gr = ggplot(eff_int,aes(x=GradeM,y=emmean,group=MathContent,color=MathContent,fill=MathContent))+
    geom_ribbon(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),alpha=.15,color=NA) + 
    geom_hline(yintercept=0,size=1.5,linetype=2,color="#888888") +
    geom_line(size=1.5) + 
    scale_x_continuous(breaks=seq(0,14,2))+
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=28),title=element_text(size=24)) +
    coord_cartesian(ylim=c(-0.2,0.45)) +
    ggtitle("Type of content x Grade interaction") + xlab("Grade") + ylab("Estimated SMD")
)
# see estimated differences at later grades
emmeans(emmprep(fit2,at=list(GradeM=14)),specs="GradeM",by="MathContent")

# plot of interaction plus effects distribution
(gg_gr = ggplot(eff_int,aes(x=GradeM,y=emmean,group=MathContent,color=MathContent,fill=MathContent))+
    geom_point(data=d,aes(x=GradeM,y=eff),size=log(1/d$vi),alpha=.4, position=position_jitter(width=.075,height=0)) +
    geom_ribbon(size=1,aes(ymin=asymp.UCL,ymax=asymp.LCL),alpha=.2,color=NA) + 
    geom_hline(yintercept=0,size=1.5,linetype=2,color="#888888") +
    geom_line(size=2) + 
    scale_x_continuous(breaks=seq(0,14,2))+
    scale_y_continuous(breaks=seq(-2,2,.1)) + 
    theme(text=element_text(size=32),title=element_text(size=23)) +
    coord_cartesian(ylim=c(-0.5,0.7)) +
    ggtitle("Type of content x Grade interaction") + xlab("Grade") + ylab("Estimated SMD")
)


#########################################################

save.image("workspaceMathGenderMetaAnalysis.RData")

#########################################################

