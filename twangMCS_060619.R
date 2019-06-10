#set-up commands

setwd("T:/projects/Causal Lifecourse Methods S00333/Data/AnonymisedData")
library(haven)
library(tidyverse)
library(forcats)
library(twang)
library(survey)
set.seed(1)
MCS <- read.csv("MCS_longweightkeyvars_030619.csv", header=FALSE,
                    col.names=c("MCSID",
                                "swgt",
                                "PTTYPE2",
                                "SPTN00",
                                "NH2",
                                "country",
                                "ethmin",
                                "paredu",
                                "smkpreg",
                                "mage",
                                "magect",
                                "amatmh",
                                "apov",
                                "bpov",
                                "cpov",
                                "dpov",
                                "epov",
                                "bmatmh",
                                "cmatmh",
                                "dmatmh",
                                "ematmh",
                                "cumpov",
                                "cumpov5",
                                "anypov",
                                "pstpov",
                                "asmk",
                                "bsmk",
                                "csmk",
                                "dsmk",
                                "esmk",
                                "wales",
                                "scot",
                                "nirlnd",
                                "smkout",
                                "obsfwgt"))
MCS[, c('country')] <- lapply(MCS[,c('country'), drop=F], as.factor)
MCS[, c('magect')] <- lapply(MCS[,c('magect'), drop=F], as.factor)
head(MCS)
smpw <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~obsfwgt, data=MCS, fpc=~NH2)

#simplified version
#wave 1-unweighted
unwgta <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgta)

#wave 1-ATE
psate.a <- ps(apov~ethmin+smkpreg+paredu+magect,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psate.a)
summary(psate.a$gbm.obj,
        n.trees=psate.a$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.a.bal<-bal.table(psate.a)
psate.a.bal
summary(psate.a)
plot(psate.a, plots=1)
plot(psate.a, plots=2)
plot(psate.a, plots=3)
plot(psate.a, plots=4)
plot(psate.a, plots=5)
anum <- svyglm(apov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$apovp <- anum$fitted.values
MCS$aprob <- if_else(MCS$apov==1,MCS$apovp,1-MCS$apovp,missing=NULL)
MCS <- bind_cols(MCS,as.vector(psate.a$ps))
MCS <- mutate(MCS,
              aatep=es.mean.ATE)
MCS <- select(MCS, -(es.mean.ATE))
MCS$aprobd <- if_else(MCS$apov==1,MCS$aatep,1-MCS$aatep,missing=NULL)
MCS <- MCS %>%
  mutate(watea=obsfwgt*(aprob/aprobd),
         atea=aprob/aprobd)

ate.a <-dx.wts(x=MCS$watea,
                  data=MCS,
                  vars=c("ethmin","smkpreg","paredu","magect"),
                  treat.var="apov",
                  perm.test.iters=0, estimand="ATE")
ate.a
bal.table(ate.a)

w.atea <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watea, data=MCS, fpc=~NH2)
est_atea <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=w.atea)
summary(est_atea)

#wave 1-ATT
psatt.a <- ps(apov~ethmin+smkpreg+paredu+magect,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psatt.a)
summary(psatt.a$gbm.obj,
        n.trees=psatt.a$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.a.bal<-bal.table(psatt.a)
psatt.a.bal
summary(psatt.a)
plot(psatt.a, plots=1)
plot(psatt.a, plots=2)
plot(psatt.a, plots=3)
plot(psatt.a, plots=4)
plot(psatt.a, plots=5)
MCS$watta <- get.weights(psatt.a, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.a$ps))
MCS <- mutate(MCS,
              aattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$atta <- if_else(MCS$apov==1,1,(MCS$aattp/(1-MCS$aattp)),missing=NULL)

w.atta <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watta, data=MCS, fpc=~NH2)
est_atta <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=w.atta)
summary(est_atta)

#wave 1-cumulative ATE-equivalent to ATE

#wave 2
#wave 2-unweighted
unwgtb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtb)

#wave 2-ATE
psate.b <- ps(bpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psate.b)
summary(psate.b$gbm.obj,
        n.trees=psate.b$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.b.bal<-bal.table(psate.b)
psate.b.bal
summary(psate.b)
plot(psate.b, plots=1)
plot(psate.b, plots=2)
plot(psate.b, plots=3)
plot(psate.b, plots=4)
plot(psate.b, plots=5)
bnum <- svyglm(bpov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$bpovp <- bnum$fitted.values
MCS$bprob <- if_else(MCS$bpov==1,MCS$bpovp,1-MCS$bpovp,missing=NULL)
MCS <- bind_cols(MCS,as.vector(psate.b$ps))
MCS <- MCS %>%
  mutate(batep=es.mean.ATE) %>%
  select(-(es.mean.ATE))
MCS$bprobd <- if_else(MCS$bpov==1,MCS$batep,1-MCS$batep,missing=NULL)
MCS <- MCS %>%
  mutate(wateb=obsfwgt*(bprob/bprobd),
         ateb=bprob/bprobd)

ate.b <-dx.wts(x=MCS$wateb,
               data=MCS,
               vars=c("ethmin","smkpreg","paredu","magect","apov","asmk","amatmh"),
               treat.var="bpov",
               perm.test.iters=0, estimand="ATE")
ate.b
bal.table(ate.b)

w.ateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wateb, data=MCS, fpc=~NH2)
est_ateb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=w.ateb)
summary(est_ateb)

#wave 2-ATT
psatt.b <- ps(bpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psatt.b)
summary(psatt.b$gbm.obj,
        n.trees=psatt.b$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.b.bal<-bal.table(psatt.b)
psatt.b.bal
summary(psatt.b)
plot(psatt.b, plots=1)
plot(psatt.b, plots=2)
plot(psatt.b, plots=3)
plot(psatt.b, plots=4)
plot(psatt.b, plots=5)
MCS$wattb <- get.weights(psatt.b, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.b$ps))
MCS <- mutate(MCS,
              battp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$attb <- if_else(MCS$bpov==1,1,(MCS$battp/(1-MCS$battp)),missing=NULL)

w.attb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wattb, data=MCS, fpc=~NH2)
est_attb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=w.attb)
summary(est_attb)

#wave 2-cumulative ATE
#should be based on obsfwgt*atea*ateblng
blong <- svyglm(bpov~apov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$bpovplong <- blong$fitted.values
MCS$bproblong <- if_else(MCS$bpov==1,MCS$bpovplong,1-MCS$bpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateblng=bproblong/bprobd,
         cumateb=obsfwgt*atea*ateblng,
         cumpovb=apov+bpov)
table(MCS$cumpovb)
#check balance over cumulative exposure count
MCSpov01 <- filter(MCS,cumpovb<=1)
MCSpov12 <- filter(MCS,cumpovb>=1)
MCSpov12$cumpovb <- MCSpov12$cumpovb-1
table(MCSpov12$cumpovb)
MCSpov02 <- filter(MCS,cumpovb==0|cumpovb==2)
MCSpov02$cumpovb <- if_else(MCSpov02$cumpovb==2,1,0,missing=NULL)
table(MCSpov02$cumpovb)

cumateb01 <-dx.wts(x=MCSpov01$cumateb,
                   data=MCSpov01,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                   treat.var="cumpovb",
                   perm.test.iters=0, estimand="ATE")
cumateb01
bal.table(cumateb01)
cumateb02 <-dx.wts(x=MCSpov02$cumateb,
                   data=MCSpov02,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                   treat.var="cumpovb",
                   perm.test.iters=0, estimand="ATE")
cumateb02
bal.table(cumateb02)
cumateb12 <-dx.wts(x=MCSpov12$cumateb,
                   data=MCSpov12,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                   treat.var="cumpovb",
                   perm.test.iters=0, estimand="ATE")
cumateb12
bal.table(cumateb12)

w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_cumateb <- svyglm(smkout~cumpovb, family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_cumateb)

#effect of any pov to b
#effect of any/first pov in b
#effect of any/persist pov to b
#all use cumulative ate weight to b but different code of exposure

#effect of any pov up to b
MCS$anypovb <- if_else(MCS$apov==1|MCS$bpov==1,1,0,missing=NULL)
table(MCS$anypovb)
anypovb <-dx.wts(x=MCS$cumateb,
                 data=MCS,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                 treat.var="anypovb",
                 perm.test.iters=0, estimand="ATE")
anypovb
bal.table(anypovb)
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_anyb <- svyglm(smkout~anypovb, family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_anyb)

#effect of any/first pov in b
MCS$firstpovb <- if_else(MCS$apov==0 & MCS$bpov==1,2,MCS$anypovb,missing=NULL)
table(MCS$firstpovb)
MCSanypovb <- filter(MCS,MCS$firstpovb==0|MCS$firstpovb==2)
MCSanypovb$firstpovb <- if_else(MCSanypovb$firstpovb==2,1,0,missing=NULL)
table(MCSanypovb$firstpovb)

povb1st <-dx.wts(x=MCSanypovb$cumateb,
                 data=MCSanypovb,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                 treat.var="firstpovb",
                 perm.test.iters=0, estimand="ATE")
povb1st 
bal.table(povb1st)
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_1stb <- svyglm(smkout~as.factor(firstpovb), family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_1stb)

#effect of persistent pov in b
MCS$pstpovb <- if_else(MCS$apov==1 & MCS$bpov==1,2,MCS$anypovb,missing=NULL)
table(MCS$pstpovb)
MCSpstpovb <- filter(MCS,pstpovb==0|pstpovb==2)
MCSpstpovb$pstpovb <- if_else(MCSpstpovb$pstpovb==2,1,0,missing=NULL)
table(MCSpstpovb$pstpovb)
pstpovb <-dx.wts(x=MCSpstpovb$cumateb,
                 data=MCSpstpovb,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh"),
                 treat.var="pstpovb",
                 perm.test.iters=0, estimand="ATE")
pstpovb
bal.table(pstpovb)
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_pstb <- svyglm(smkout~as.factor(pstpovb), family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_pstb)

#wave 3
#wave 3-unweighted
unwgtc <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtc)

#wave 3-ATE
psate.c <- ps(cpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psate.c)
summary(psate.c$gbm.obj,
        n.trees=psate.c$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.c.bal<-bal.table(psate.c)
psate.c.bal
summary(psate.c)
plot(psate.c, plots=1)
plot(psate.c, plots=2)
plot(psate.c, plots=3)
plot(psate.c, plots=4)
plot(psate.c, plots=5)
cnum <- svyglm(cpov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$cpovp <- cnum$fitted.values
MCS$cprob <- if_else(MCS$cpov==1,MCS$cpovp,1-MCS$cpovp,missing=NULL)
MCS <- bind_cols(MCS,as.vector(psate.c$ps))
MCS <- MCS %>%
  mutate(catep=es.mean.ATE) %>%
  select(-(es.mean.ATE))
MCS$cprobd <- if_else(MCS$cpov==1,MCS$catep,1-MCS$catep,missing=NULL)
MCS <- MCS %>%
  mutate(watec=obsfwgt*(cprob/cprobd),
         atec=cprob/cprobd)

ate.c <-dx.wts(x=MCS$watec,
               data=MCS,
               vars=c("ethmin","smkpreg","paredu","magect","apov","asmk","amatmh","bpov","bsmk","bmatmh"),
               treat.var="cpov",
               perm.test.iters=0, estimand="ATE")
ate.c
bal.table(ate.c)

w.atec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watec, data=MCS, fpc=~NH2)
est_atec <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=w.atec)
summary(est_atec)

#wave 3-ATT
psatt.c <- ps(cpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psatt.c)
summary(psatt.c$gbm.obj,
        n.trees=psatt.c$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.c.bal<-bal.table(psatt.c)
psatt.c.bal
summary(psatt.c)
plot(psatt.c, plots=1)
plot(psatt.c, plots=2)
plot(psatt.c, plots=3)
plot(psatt.c, plots=4)
plot(psatt.c, plots=5)
MCS$wattc <- get.weights(psatt.c, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.c$ps))
MCS <- mutate(MCS,
              cattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$attc <- if_else(MCS$cpov==1,1,(MCS$cattp/(1-MCS$cattp)),missing=NULL)

w.attc <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wattc, data=MCS, fpc=~NH2)
est_attc <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=w.attc)
summary(est_attc)

#wave 3 -cumulative ATE
#should be based on obsfwgt*atea*ateblng*ateclng
clong <- svyglm(cpov~apov+bpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$cpovplong <- clong$fitted.values
MCS$cproblong <- if_else(MCS$cpov==1,MCS$cpovplong,1-MCS$cpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateclng=cproblong/cprobd,
         cumatec=obsfwgt*atea*ateblng*ateclng,
         cumpovc=apov+bpov+cpov)
table(MCS$cumpovc)
#check balance over cumulative exposure count
MCSpov01 <- filter(MCS,cumpovc<=1)
MCSpov12 <- filter(MCS,cumpovc>=1 & cumpovc<=2)
MCSpov12$cumpovc <- MCSpov12$cumpovc-1
table(MCSpov12$cumpovc)
MCSpov23 <- filter(MCS,cumpovc>=2)
MCSpov23$cumpovc <- MCSpov23$cumpovc-2
table(MCSpov23$cumpovc)
MCSpov03 <- filter(MCS,cumpovc==0|cumpovc==3)
MCSpov03$cumpovc <- if_else(MCSpov03$cumpovc==3,1,0,missing=NULL)
table(MCSpov03$cumpovc)

cumatec01 <-dx.wts(x=MCSpov01$cumatec,
                   data=MCSpov01,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                   treat.var="cumpovc",
                   perm.test.iters=0, estimand="ATE")
cumatec01
bal.table(cumatec01)
cumatec12 <-dx.wts(x=MCSpov12$cumatec,
                   data=MCSpov12,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                   treat.var="cumpovc",
                   perm.test.iters=0, estimand="ATE")
cumatec12
bal.table(cumatec12)
cumatec23 <-dx.wts(x=MCSpov23$cumatec,
                   data=MCSpov23,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                   treat.var="cumpovc",
                   perm.test.iters=0, estimand="ATE")
cumatec23
bal.table(cumatec23)
cumatec03 <-dx.wts(x=MCSpov03$cumatec,
                   data=MCSpov03,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                   treat.var="cumpovc",
                   perm.test.iters=0, estimand="ATE")
cumatec03
bal.table(cumatec03)

w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_cumatec <- svyglm(smkout~cumpovc, family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_cumatec)

#effect of any pov to c
#effect of any/first pov in c
#effect of any/persist pov to c
#all use cumulative ate weight to c but different code of exposure

#effect of any pov up to c
MCS$anypovc <- if_else(MCS$apov==1|MCS$bpov==1|MCS$cpov==1,1,0,missing=NULL)
table(MCS$anypovc)
anypovc <-dx.wts(x=MCS$cumatec,
                 data=MCS,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                 treat.var="anypovc",
                 perm.test.iters=0, estimand="ATE")
anypovc
bal.table(anypovc)
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_anyc <- svyglm(smkout~anypovc, family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_anyc)

#effect of any/first pov in c
MCS$firstpovc <- if_else(MCS$apov==0 & MCS$bpov==0 & MCS$cpov==1,2,MCS$anypovc,missing=NULL)
table(MCS$firstpovc)
MCS1stpovc <- filter(MCS,MCS$firstpovc==0|MCS$firstpovc==2)
MCS1stpovc$firstpovc <- if_else(MCS1stpovc$firstpovc==2,1,0,missing=NULL) 
table(MCS1stpovc$firstpovc)
povc1st <-dx.wts(x=MCS1stpovc$cumatec,
                 data=MCS1stpovc,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                 treat.var="firstpovc",
                 perm.test.iters=0, estimand="ATE")
povc1st 
bal.table(povc1st)
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_1stc <- svyglm(smkout~as.factor(firstpovc), family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_1stc)

#effect of any/persistent pov in c
MCS$pstpovc <- if_else(MCS$apov==1 & MCS$bpov==1 & MCS$cpov==1,2,MCS$anypovc,missing=NULL)
table(MCS$pstpovc)
MCSpstpovc <- filter(MCS,MCS$pstpovc==0|MCS$pstpovc==2)
MCSpstpovc$pstpovc <- if_else(MCSpstpovc$pstpovc==2,1,0,missing=NULL) 
table(MCSpstpovc$pstpovc)
pstpovc <-dx.wts(x=MCSpstpovc$cumatec,
                 data=MCSpstpovc,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh"),
                 treat.var="pstpovc",
                 perm.test.iters=0, estimand="ATE")
pstpovc
bal.table(pstpovc)
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_pstc <- svyglm(smkout~as.factor(pstpovc), family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_pstc)

#wave 4
#wave 4-unweighted
unwgtd <- svyglm(smkout~dpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtd)

#wave 4-ATE
psate.d <- ps(dpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psate.d)
summary(psate.d$gbm.obj,
        n.trees=psate.d$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.d.bal<-bal.table(psate.d)
psate.d.bal
summary(psate.d)
plot(psate.d, plots=1)
plot(psate.d, plots=2)
plot(psate.d, plots=3)
plot(psate.d, plots=4)
plot(psate.d, plots=5)
dnum <- svyglm(dpov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$dpovp <- dnum$fitted.values
MCS$dprob <- if_else(MCS$dpov==1,MCS$dpovp,1-MCS$dpovp,missing=NULL)
MCS <- bind_cols(MCS,as.vector(psate.d$ps))
MCS <- MCS %>%
  mutate(datep=es.mean.ATE) %>%
  select(-(es.mean.ATE))
MCS$dprobd <- if_else(MCS$dpov==1,MCS$datep,1-MCS$datep,missing=NULL)
MCS <- MCS %>%
  mutate(wated=obsfwgt*(dprob/dprobd),
         ated=dprob/dprobd)

ate.d <-dx.wts(x=MCS$wated,
               data=MCS,
               vars=c("ethmin","smkpreg","paredu","magect","apov","asmk","amatmh","bpov","bsmk","bmatmh","cpov","csmk","cmatmh"),
               treat.var="dpov",
               perm.test.iters=0, estimand="ATE")
ate.d
bal.table(ate.d)

w.ated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wated, data=MCS, fpc=~NH2)
est_ated <- svyglm(smkout~dpov, family=quasibinomial(link = 'logit'), design=w.ated)
summary(est_ated)

#wave 4-ATT
psatt.d <- ps(dpov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psatt.d)
summary(psatt.d$gbm.obj,
        n.trees=psatt.d$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.d.bal<-bal.table(psatt.d)
psatt.d.bal
summary(psatt.d)
plot(psatt.d, plots=1)
plot(psatt.d, plots=2)
plot(psatt.d, plots=3)
plot(psatt.d, plots=4)
plot(psatt.d, plots=5)
MCS$wattd <- get.weights(psatt.d, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.d$ps))
MCS <- mutate(MCS,
              dattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$attd <- if_else(MCS$dpov==1,1,(MCS$dattp/(1-MCS$dattp)),missing=NULL)

w.attd <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wattd, data=MCS, fpc=~NH2)
est_attd <- svyglm(smkout~dpov, family=quasibinomial(link = 'logit'), design=w.attd)
summary(est_attd)

#wave 4-cumulative ATE
#should be based on obsfwgt*atea*ateblng*ateclng*atedlng
dlong <- svyglm(dpov~apov+bpov+cpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$dpovplong <- dlong$fitted.values
MCS$dproblong <- if_else(MCS$dpov==1,MCS$dpovplong,1-MCS$dpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(atedlng=dproblong/dprobd,
         cumated=obsfwgt*atea*ateblng*ateclng*atedlng,
         cumpovd=apov+bpov+cpov+dpov)
table(MCS$cumpovd)
#check balance over cumulative exposure count
MCSpov01 <- filter(MCS,cumpovd<=1)
MCSpov12 <- filter(MCS,cumpovd>=1 & cumpovd<=2)
MCSpov12$cumpovd <- MCSpov12$cumpovd-1
table(MCSpov12$cumpovd)
MCSpov23 <- filter(MCS,cumpovd>=2 & cumpovd<=3)
MCSpov23$cumpovd <- MCSpov23$cumpovd-2
table(MCSpov23$cumpovd)
MCSpov34 <- filter(MCS,cumpovd>=3 & cumpovd<=4)
MCSpov34$cumpovd <- MCSpov34$cumpovd-3
table(MCSpov34$cumpovd)
MCSpov04 <- filter(MCS,cumpovd==0|cumpovd==4)
MCSpov04$cumpovd <- if_else(MCSpov04$cumpovd==4,1,0,missing=NULL)
table(MCSpov04$cumpovd)

cumated01 <-dx.wts(x=MCSpov01$cumated,
                   data=MCSpov01,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                   treat.var="cumpovd",
                   perm.test.iters=0, estimand="ATE")
cumated01
bal.table(cumated01)
cumated12 <-dx.wts(x=MCSpov12$cumated,
                   data=MCSpov12,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                   treat.var="cumpovd",
                   perm.test.iters=0, estimand="ATE")
cumated12
bal.table(cumated12)
cumated23 <-dx.wts(x=MCSpov23$cumated,
                   data=MCSpov23,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                   treat.var="cumpovd",
                   perm.test.iters=0, estimand="ATE")
cumated23
bal.table(cumated23)
cumated34 <-dx.wts(x=MCSpov34$cumated,
                   data=MCSpov34,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                   treat.var="cumpovd",
                   perm.test.iters=0, estimand="ATE")
cumated34
bal.table(cumated34)
cumated04 <-dx.wts(x=MCSpov04$cumated,
                   data=MCSpov04,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                   treat.var="cumpovd",
                   perm.test.iters=0, estimand="ATE")
cumated04
bal.table(cumated04)

w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_cumated <- svyglm(smkout~cumpovd, family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_cumated)

#effect of any pov to d
#effect of any/first pov in d
#effect of any/persist pov to d
#all use cumulative ate weight to d but different code of exposure

#effect of any pov up to d
MCS$anypovd <- if_else(MCS$apov==1|MCS$bpov==1|MCS$cpov==1|MCS$dpov==1,1,0,missing=NULL)
table(MCS$anypovd)
anypovd <-dx.wts(x=MCS$cumated,
                 data=MCS,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                 treat.var="anypovd",
                 perm.test.iters=0, estimand="ATE")
anypovd
bal.table(anypovd)
w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_anyd <- svyglm(smkout~anypovd, family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_anyd)

#effect of any/first pov in d
MCS$firstpovd <- if_else(MCS$apov==0 & MCS$bpov==0 & MCS$cpov==0 & MCS$dpov==1,2,MCS$anypovd,missing=NULL)
table(MCS$firstpovd)
MCS1stpovd <- filter(MCS,MCS$firstpovd==0|MCS$firstpovd==2)
MCS1stpovd$firstpovd <- if_else(MCS1stpovd$firstpovd==2,1,0,missing=NULL) 
table(MCS1stpovd$firstpovd)
povd1st <-dx.wts(x=MCS1stpovd$cumated,
                 data=MCS1stpovd,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                 treat.var="firstpovd",
                 perm.test.iters=0, estimand="ATE")
povd1st 
bal.table(povd1st)
w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_1std <- svyglm(smkout~as.factor(firstpovd), family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_1std)

#effect of any/persistent pov in d
MCS$pstpovd <- if_else(MCS$apov==1 & MCS$bpov==1 & MCS$cpov==1 & MCS$dpov==1,2,MCS$anypovd,missing=NULL)
table(MCS$pstpovd)
MCSpstpovd <- filter(MCS,MCS$pstpovd==0|MCS$pstpovd==2)
MCSpstpovd$pstpovd <- if_else(MCSpstpovd$pstpovd==2,1,0,missing=NULL) 
table(MCSpstpovd$pstpovd)
pstpovd <-dx.wts(x=MCSpstpovd$cumated,
                 data=MCSpstpovd,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh"),
                 treat.var="pstpovd",
                 perm.test.iters=0, estimand="ATE")
pstpovd
bal.table(pstpovd)
w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_pstd <- svyglm(smkout~as.factor(pstpovd), family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_pstd)

#wave 5
#wave 5-unweighted
unwgte <- svyglm(smkout~epov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgte)

#wave 5-ATE
psate.e <- ps(epov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh+dpov+dsmk+dmatmh,
              data=MCS,
              n.trees=30000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psate.e)
summary(psate.e$gbm.obj,
        n.trees=psate.e$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.e.bal<-bal.table(psate.e)
psate.e.bal
summary(psate.e)
plot(psate.e, plots=1)
plot(psate.e, plots=2)
plot(psate.e, plots=3)
plot(psate.e, plots=4)
plot(psate.e, plots=5)
enum <- svyglm(epov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$epovp <- enum$fitted.values
MCS$eprob <- if_else(MCS$epov==1,MCS$epovp,1-MCS$epovp,missing=NULL)
MCS <- bind_cols(MCS,as.vector(psate.e$ps))
MCS <- MCS %>%
  mutate(eatep=es.mean.ATE) %>%
  select(-(es.mean.ATE))
MCS$eprobd <- if_else(MCS$epov==1,MCS$eatep,1-MCS$eatep,missing=NULL)
MCS <- MCS %>%
  mutate(watee=obsfwgt*(eprob/eprobd),
         atee=eprob/eprobd)

ate.e <-dx.wts(x=MCS$watee,
               data=MCS,
               vars=c("ethmin","smkpreg","paredu","magect","apov","asmk","amatmh","bpov","bsmk","bmatmh","cpov","csmk","cmatmh","dpov","dsmk","dmatmh"),
               treat.var="epov",
               perm.test.iters=0, estimand="ATE")
ate.e
bal.table(ate.e)

w.atee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watee, data=MCS, fpc=~NH2)
est_atee <- svyglm(smkout~epov, family=quasibinomial(link = 'logit'), design=w.atee)
summary(est_atee)

#wave 4-ATT
psatt.e <- ps(epov~ethmin+smkpreg+paredu+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk++cmatmh+dpov+dsmk+dmatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
plot(psatt.e)
summary(psatt.e$gbm.obj,
        n.trees=psatt.e$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.e.bal<-bal.table(psatt.e)
psatt.e.bal
summary(psatt.e)
plot(psatt.e, plots=1)
plot(psatt.e, plots=2)
plot(psatt.e, plots=3)
plot(psatt.e, plots=4)
plot(psatt.e, plots=5)
MCS$watte <- get.weights(psatt.e, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.e$ps))
MCS <- mutate(MCS,
              eattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$atte <- if_else(MCS$epov==1,1,(MCS$eattp/(1-MCS$eattp)),missing=NULL)

w.atte <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watte, data=MCS, fpc=~NH2)
est_atte <- svyglm(smkout~epov, family=quasibinomial(link = 'logit'), design=w.atte)
summary(est_atte)

#wave 5-cumulative ATE
#should be based on obsfwgt*atea*ateblng*ateclng*atedlng*ateelng
elong <- svyglm(epov~apov+bpov+cpov+dpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$epovplong <- elong$fitted.values
MCS$eproblong <- if_else(MCS$epov==1,MCS$epovplong,1-MCS$epovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateelng=eproblong/eprobd,
         cumatee=obsfwgt*atea*ateblng*ateclng*atedlng*ateelng,
         cumpove=apov+bpov+cpov+dpov+epov)
table(MCS$cumpove)
#check balance over cumulative exposure count
MCSpov01 <- filter(MCS,cumpove<=1)
MCSpov12 <- filter(MCS,cumpove>=1 & cumpove<=2)
MCSpov12$cumpove <- MCSpov12$cumpove-1
table(MCSpov12$cumpove)
MCSpov23 <- filter(MCS,cumpove>=2 & cumpove<=3)
MCSpov23$cumpove <- MCSpov23$cumpove-2
table(MCSpov23$cumpove)
MCSpov34 <- filter(MCS,cumpove>=3 & cumpove<=4)
MCSpov34$cumpove <- MCSpov34$cumpove-3
table(MCSpov34$cumpove)
MCSpov45 <- filter(MCS,cumpove>=4 & cumpove<=5)
MCSpov45$cumpove <- MCSpov45$cumpove-4
table(MCSpov45$cumpove)
MCSpov05 <- filter(MCS,cumpove==0|cumpove==5)
MCSpov05$cumpove <- if_else(MCSpov05$cumpove==5,1,0,missing=NULL)
table(MCSpov05$cumpove)

cumatee01 <-dx.wts(x=MCSpov01$cumatee,
                   data=MCSpov01,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee01
bal.table(cumatee01)
cumatee12 <-dx.wts(x=MCSpov12$cumatee,
                   data=MCSpov12,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee12
bal.table(cumatee12)
cumatee23 <-dx.wts(x=MCSpov23$cumatee,
                   data=MCSpov23,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee23
bal.table(cumatee23)
cumatee34 <-dx.wts(x=MCSpov34$cumatee,
                   data=MCSpov34,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee34
bal.table(cumatee34)
cumatee45 <-dx.wts(x=MCSpov45$cumatee,
                   data=MCSpov45,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee45
bal.table(cumatee45)
cumatee05 <-dx.wts(x=MCSpov05$cumatee,
                   data=MCSpov05,
                   vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                   treat.var="cumpove",
                   perm.test.iters=0, estimand="ATE")
cumatee05
bal.table(cumatee05)

w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_cumatee <- svyglm(smkout~cumpove, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cumatee)

#effect of any pov to e
#effect of any/first pov in e
#effect of any/persist pov to e
#all use cumulative ate weight to e but different code of exposure

#effect of any pov up to e
MCS$anypove <- if_else(MCS$cumpove>=1,1,0,missing=NULL)
table(MCS$anypove)
anypove <-dx.wts(x=MCS$cumatee,
                 data=MCS,
                 vars=c("ethmin","smkpreg","paredu","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                 treat.var="anypove",
                 perm.test.iters=0, estimand="ATE")
anypove
bal.table(anypove)
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_anye <- svyglm(smkout~anypove, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_anye)

#effect of any/first pov in e
MCS$firstpove <- if_else(MCS$cumpovd==0 & MCS$epov==1,2,MCS$anypove,missing=NULL)
table(MCS$firstpove)
MCS1stpove <- filter(MCS,MCS$firstpove==0|MCS$firstpove==2)
MCS1stpove$firstpove <- if_else(MCS1stpove$firstpove==2,1,0,missing=NULL) 
table(MCS1stpove$firstpove)
pove1st <-dx.wts(x=MCS1stpove$cumatee,
                 data=MCS1stpove,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                 treat.var="firstpove",
                 perm.test.iters=0, estimand="ATE")
pove1st 
bal.table(pove1st)
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_1ste <- svyglm(smkout~as.factor(firstpove), family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_1ste)

#effect of any/persistent pov in e
MCS$pstpove <- if_else(MCS$cumpove==5,2,MCS$anypove,missing=NULL)
table(MCS$pstpove)
MCSpstpove <- filter(MCS,MCS$pstpove==0|MCS$pstpove==2)
MCSpstpove$pstpove <- if_else(MCSpstpove$pstpove==2,1,0,missing=NULL) 
table(MCSpstpove$pstpove)
pstpove <-dx.wts(x=MCSpstpove$cumatee,
                 data=MCSpstpove,
                 vars=c("ethmin","smkpreg","paredu","magect","asmk","amatmh","bsmk","bmatmh","csmk","cmatmh","dsmk","dmatmh"),
                 treat.var="pstpove",
                 perm.test.iters=0, estimand="ATE")
pstpove
bal.table(pstpove)
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_pste <- svyglm(smkout~as.factor(pstpove), family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_pste)

#examine cumated weight for CDE of A not through B-E
#looking for balance in covariates over b-e within strata of a
table(MCS$apov,MCS$cumpove)
MCS <-mutate(MCS, bepov=bpov+cpov+dpov+epov) 
MCSapov0 <- filter(MCS, apov==0)
MCSapov0$cumpove <- if_else(MCSapov0$cumpove>=1,1,0,missing=NULL)
table(MCSapov0$cumpove)
MCSapov1 <- filter(MCS, apov==1)
MCSapov1$cumpove <- if_else(MCSapov1$cumpove<=3,0,1,missing=NULL)
table(MCSapov1$cumpove)
cdea_a0 <-dx.wts(x=MCSapov0$cumatee,
                 data=MCSapov0,
                 vars=c("ethmin","smkpreg","asmk","bsmk","csmk","dsmk"),
                 treat.var="cumpove",
                 perm.test.iters=0, estimand="ATE")
cdea_a0
bal.table(cdea_a0)
cdea_a1 <-dx.wts(x=MCSapov1$cumatee,
                 data=MCSapov1,
                 vars=c("ethmin","smkpreg","asmk","bsmk","csmk","dsmk"),
                 treat.var="cumpove",
                 perm.test.iters=0, estimand="ATE")
cdea_a1
bal.table(cdea_a1)
#estimate CDE of A not through B-E
MCS <- mutate(MCS,abepov=apov*bepov)
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_cdea_be <- svyglm(smkout~apov+bepov+abepov, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cdea_be)
est_cdea_be <- svyglm(smkout~apov+bepov, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cdea_be)



#CDE of B not through C-E
#want a weight that is: obsfwgt*ateb*atecb*atedb*ateeb
cblong <- svyglm(cpov~bpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$cbpovplong <- cblong$fitted.values
MCS$cbproblong <- if_else(MCS$cpov==1,MCS$cbpovplong,1-MCS$cbpovplong,missing=NULL)
dblong <- svyglm(dpov~bpov+cpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$dbpovplong <- dblong$fitted.values
MCS$dbproblong <- if_else(MCS$dpov==1,MCS$dbpovplong,1-MCS$dbpovplong,missing=NULL)
eblong <- svyglm(epov~bpov+cpov+dpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$ebpovplong <- eblong$fitted.values
MCS$ebproblong <- if_else(MCS$epov==1,MCS$ebpovplong,1-MCS$ebpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(atecb=cbproblong/cprobd,
         atedb=dbproblong/dprobd,
         ateeb=ebproblong/eprobd,
         cdeb=obsfwgt*ateb*atecb*atedb*ateeb,
         cepov=cpov+dpov+epov,
         bcepov=bpov*cepov)
#looking for balance in covariates over c-e within strata of b
table(MCS$bpov,MCS$cepov)
MCSbpov0 <- filter(MCS, bpov==0)
MCSbpov0$cepov <- if_else(MCSbpov0$cepov>=2,1,0,missing=NULL)
table(MCSbpov0$cepov)
MCSbpov1 <- filter(MCS, bpov==1)
MCSbpov1$cepov <- if_else(MCSbpov1$cepov<=2,0,1,missing=NULL)
table(MCSbpov1$cepov)
cdeb_b0 <-dx.wts(x=MCSbpov0$cdeb,
                 data=MCSbpov0,
                 vars=c("ethmin","smkpreg","apov","asmk","bsmk","csmk","dsmk"),
                 treat.var="cepov",
                 perm.test.iters=0, estimand="ATE")
cdeb_b0
bal.table(cdeb_b0)
cdeb_b1 <-dx.wts(x=MCSbpov1$cdeb,
                 data=MCSbpov1,
                 vars=c("ethmin","smkpreg","apov","asmk","bsmk","csmk","dsmk"),
                 treat.var="cepov",
                 perm.test.iters=0, estimand="ATE")
cdeb_b1
bal.table(cdeb_b1)
#estimate CDE of B not through C-E
w.cdeb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cdeb, data=MCS, fpc=~NH2)
est_cdeb <- svyglm(smkout~bpov+cepov+bcepov, family=quasibinomial(link = 'logit'), design=w.cdeb)
summary(est_cdeb)
est_cdeb <- svyglm(smkout~bpov+cepov, family=quasibinomial(link = 'logit'), design=w.cdeb)
summary(est_cdeb)

#CDE of C not through D-E
#want a weight that is: obsfwgt*atec*atedc*ateec
dclong <- svyglm(dpov~cpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$dcpovplong <- dclong$fitted.values
MCS$dcproblong <- if_else(MCS$dpov==1,MCS$dcpovplong,1-MCS$dcpovplong,missing=NULL)
eclong <- svyglm(epov~cpov+dpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$ecpovplong <- eclong$fitted.values
MCS$ecproblong <- if_else(MCS$epov==1,MCS$ecpovplong,1-MCS$ecpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(atedc=dcproblong/dprobd,
         ateec=ecproblong/eprobd,
         cdec=obsfwgt*atec*atedc*ateec,
         depov=dpov+epov,
         cdepov=cpov*depov)
#looking for balance in covariates over c-e within strata of b
table(MCS$cpov,MCS$depov)
MCScpov0 <- filter(MCS, cpov==0)
MCScpov0$depov <- if_else(MCScpov0$depov>=1,1,0,missing=NULL)
table(MCScpov0$depov)
MCScpov1 <- filter(MCS, cpov==1)
MCScpov1$depov <- if_else(MCScpov1$depov<=0,0,1,missing=NULL)
table(MCScpov1$depov)
cdec_c0 <-dx.wts(x=MCScpov0$cdec,
                 data=MCScpov0,
                 vars=c("ethmin","smkpreg","apov","asmk","bpov","bsmk","csmk","dsmk"),
                 treat.var="depov",
                 perm.test.iters=0, estimand="ATE")
cdec_c0
bal.table(cdec_c0)
cdec_c1 <-dx.wts(x=MCScpov1$cdec,
                 data=MCScpov1,
                 vars=c("ethmin","smkpreg","apov","asmk","bpov","bsmk","csmk","dsmk"),
                 treat.var="depov",
                 perm.test.iters=0, estimand="ATE")
cdec_c1
bal.table(cdec_c1)
#estimate CDE of C not through D-E
w.cdec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cdec, data=MCS, fpc=~NH2)
est_cdec <- svyglm(smkout~cpov+depov+cdepov, family=quasibinomial(link = 'logit'), design=w.cdec)
summary(est_cdec)
est_cdec <- svyglm(smkout~cpov+depov, family=quasibinomial(link = 'logit'), design=w.cdec)
summary(est_cdec)

#CDE of D not through E
#want a weight that is: obsfwgt*ated*ateed
edlong <- svyglm(epov~dpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$edpovplong <- edlong$fitted.values
MCS$edproblong <- if_else(MCS$epov==1,MCS$edpovplong,1-MCS$edpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateed=edproblong/eprobd,
         cded=obsfwgt*ated*ateed,
         depovi=dpov*epov)
#looking for balance in covariates over c-e within strata of b
table(MCS$dpov,MCS$epov)
MCSdpov0 <- filter(MCS, dpov==0)
MCSdpov0$epov <- if_else(MCSdpov0$epov==1,1,0,missing=NULL)
table(MCSdpov0$epov)
MCSdpov1 <- filter(MCS, dpov==1)
MCSdpov1$epov <- if_else(MCSdpov1$epov==1,1,0,missing=NULL)
table(MCSdpov1$epov)
cded_d0 <-dx.wts(x=MCSdpov0$cded,
                 data=MCSdpov0,
                 vars=c("ethmin","smkpreg","apov","asmk","bpov","bsmk","cpov","csmk","dsmk"),
                 treat.var="epov",
                 perm.test.iters=0, estimand="ATE")
cded_d0
bal.table(cded_d0)
cded_d1 <-dx.wts(x=MCSdpov1$cded,
                 data=MCSdpov1,
                 vars=c("ethmin","smkpreg","apov","asmk","bpov","bsmk","cpov","csmk","dsmk"),
                 treat.var="epov",
                 perm.test.iters=0, estimand="ATE")
cded_d1
bal.table(cded_d1)
#estimate CDE of D not through E
w.cded <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cded, data=MCS, fpc=~NH2)
est_cded <- svyglm(smkout~dpov+epov+depovi, family=quasibinomial(link = 'logit'), design=w.cded)
summary(est_cded)
est_cded <- svyglm(smkout~dpov+epov, family=quasibinomial(link = 'logit'), design=w.cded)
summary(est_cded)

#factorised cumulative exposure?
MCS$cumpovf <- as.factor(MCS$cumpov)
table(MCS$cumpovf)
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_cumateef <- svyglm(smkout~as.factor(cumpovf), family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cumateef)
#early sensitive period?
est_es <- svyglm(smkout~apov+cumpov, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_es)


