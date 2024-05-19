
# PROVE AGGREGATE

x = data.frame(study = rep(c("a","b"),each=4),
               yi = c(0.0, 0.0, 1.0 ,1.0,  1.0, 2.0, 3.0, 4.0),
               vi =  c(0.5, 0.5, 0.5, 0.5,  0.1, 0.1, 0.4, 1.0))

combine(eff=x$yi[x$study=="a"],vi=x$vi[x$study=="a"],r=0.7)
combine(eff=x$yi[x$study=="b"],vi=x$vi[x$study=="b"],r=0.7)

dat = escalc("GEN",yi=yi,vi=vi, data = x)
aggregate(dat, cluster=study, rho=0.7, weighted = FALSE)

aggregate(dat, cluster=study, rho=0.7, weighted = TRUE)

N = 1000
X = rnorm(N)
Y = rnorm(N) + X*.2
y1 = Y+rnorm(N)
y2 = Y+rnorm(N)
y3 = Y+rnorm(N)

dwide = data.frame(y1,y2,y3,X)
fitlm = lm(cbind(y1,y2,y3)~X,data=dwide)
summary(fitlm)
anova(fitlm)

dlong = data.frame(y = c(y1,y2,y3), X = c(X,X,X), id=rep(1:1000,times=3))
fit = lme4::lmer(y~X+(1|id),data=dlong)
summary(fit)

#####################

dat <- dat.berkey1998
dat$vi = dat$vi*2
V <- metafor::vcalc(dat$vi,cluster=dat$trial,rho=0.99,obs=dat$outcome)
res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat)
res
