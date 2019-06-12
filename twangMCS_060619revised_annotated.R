####intro section: set-up commands
#these commands set up the working environment in R and load the data so that all the subsequent commands can run

#this sets the working directory where my data file is stored
setwd("T:/projects/Causal Lifecourse Methods S00333/Data/AnonymisedData")

#these commands load required R packages for use
#if you do not have these packages installed you will also need to run the following command for each one:
#install.packages("package.name")
library(haven)
library(tidyverse)
library(forcats)
library(twang)
library(survey)
library(cobalt)

#there is a stochastic element to the machine learning methods used to optimise confounder balance
#this command sets R's random number generator seed to a fixed value
#it ensures that the results will be reproducible
set.seed(1)

#this commands loads my CSV data file into R
#my file has no header, hence 'header=FALSE'
#the list of terms enclosed in double-quotes are the variable names
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
#there are two variables that I'm wanting to treat as unordered categorical variables
#these commands identify these two variables to R as 'factor' variables
MCS[, c('country')] <- lapply(MCS[,c('country'), drop=F], as.factor)
MCS[, c('magect')] <- lapply(MCS[,c('magect'), drop=F], as.factor)
#this is just a check showing the first few lines of data to see if it has been read correctly
head(MCS)

#this specifies the sampling design and the weight to adjust for this
smpw <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~obsfwgt, data=MCS, fpc=~NH2)
#it specifies the data file "data=MCS" and the variables that encode the sample design
#SPTN00 represents cluster ids
#PTTYPE2 is the sampling strata
#NH2 is a finite population correction
#swgt is the sampling weight and is not included in above command because:
#obsfwgt combines swgt with a weight for attrition
#this attrition weight was created in another package and combined with the sample weight before loading the data here
#'design=smpw' is used later to denote models using only the sample/attrition weighting

####end of intro section.

####analysis section
#the following is split into a section for each wave and a final cross-wave section
#the section for each wave produces a range of different estimates as follows
#unweighted -ie using only the sampling/attrition weights
#ATE (Average Treatment Effect)
#-i.e. weighting both the exposure/treatment & control group to resemble the total population
#ATT (Average Treatment Effect among the Treated)
#-i.e. weighting the control group to resemble the exposure/treatment group
#cumulative ATE -uses a Marginal Structural Model (MSM) to estimate the cumulative effect of exposure up to that time-point
#this can be interpreted as the 'dose-response' average effect of each additional instance of exposure
#the cumulative ATE assumes the exposure has a linear, cumulative impact on the outcome
#I have also given a few examples of different non-linear codings for the cumulative ATE effects

#in the final cross-wave section weights from across the study period are used to estimate:
#Controlled Direct Effects (CDEs) of exposure, not mediated through later exposure

####wave 1

####wave 1-unweighted
#this is the observational association between poverty at time 1 and the outcome (adolescent smoking)
#'smkout~apov' specifies that smkout is dependent variable and apov is an independent variable
#'family=quasibinomial(link='logit') specifies a logistic regression model
#'design=smpw' calls the sampling/attrition weights that were set up in section 1
#the first command runs the model and stores the output as 'unwgta'
#the second command gives a summary of that output including the regression coefficients and SEs
unwgta <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgta)

####wave 1-ATE
#this command runs a generalised boosted regression model
psate.a <- ps(apov~ethmin+smkpreg+paredu+country+magect,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
#this is a machine learning procedure aimed at producing a propensity score/weight
#it aims to optimise the balance that the weight provides on a specified set of covariates
#the output is stored as a new object called 'psate.a'
#ps is the name of the optimisation procedure and has several arguments
#'apov~ethmin+smkpreg+paredu+magect'
#defines the exposure of interest 'apov' 
#and a set of confounders that you wish to balance across levels of this exposure
#here I've used a set of five baseline confounders:
#ethmin=ethnic minority status (0-1)
#smkpreg=mother smoked during pregnancy (0-1)
#paredu=low parental education (0-1)
#country=UK Country, categorical (1=England,2=Wales,3=Scotland,4=N.Ireland)
#magect=categorisation of mother's age at birth (0=<20,1=20-24,2=25+)
#there is no need to specify interactions (see below)
#'data=MCS' tells it which data file to use
#'n.trees' has to do with how long the optimisation process runs for
#'interaction.depth' tells it what level of interactions to optimise balance over
#i.e. I've said 2, which means it will optimise balance over
#the list of covariates provided and all their second-order interactions
#'stop.method=c("es.mean")' means it aims to minimise the mean effect size (difference)
#between the listed covariates (and their 2nd-order interactions) across levels of the binary exposure
#'estimand="ATE"' means I'm looking for a propensity score/weight 
#that optimises balance for an ATE comparison
#'sampw=MCS$obsfwgt' tells it to use my sampling/attrition weight (even for the 'unweighted' comparison results) 
#the other arguments are just default settings

#you can get several useful pieces of output as standard
#this plots the mean effect size against the n.trees
plot(psate.a)
#basically, if the lowest point on the graph is too close to the right hand side of the plot
#then this indicates that you could get a more optimal model if you use a higher number for n.trees

#this will give you an indication of which variables were most influential in predicting your exposure
#higher 'relative influence'=more influential
summary(psate.a$gbm.obj,
        n.trees=psate.a$desc$es.mean.ATE$n.trees,
        plot=TRUE)

#these next commands produce and store a balance table showing differences between
#your exposure/treatment group and your control group
psate.a.bal<-bal.table(psate.a)
psate.a.bal
#tx=treatment group, ct=control group
#mn=mean, sd=standard deviation
#std.eff.sz=standardised mean difference between the groups
#p=p-value for that difference
#the first part of the table shows the 'unweighted' results, ie using only the sample/attrition weights
#the second part of the table shows what happens when the ATE weight is applied

#this gives you a summary of the balance table
summary(psate.a)
#points of interest are:
#n.treat, n.ctrl=the actual size of each group
#ess.treat, ess.ctrl=the effective sample size of each group after weighting
#the effective sample size takes into account the additional variability introduced by the weights
#max.es=the maximum standardised effect size (difference) across all covariates
#mean.es=the mean standardised effect size (difference) across all covariates
#ks=Kolmogorov-Smirnov statistics
#iter=tells you which iteration of the model produced the lowest/optimal value of the test stat (es.mean)
#this is equivalent to what the plot above was showing us
#in practice I tend to look at this rather than the plot
#as it can be hard to see where the lowest point on the plot is
#again if this value is close to the number specified for n.trees consider re-running with more n.trees

#a couple of other useful standard plots are:
#the distribution of the propensity scores in the treatment vs control group
plot(psate.a, plots=2)
#absolute standard differences between treatment/control for all covars pre/post weighting
plot(psate.a, plots=3)
#the t-test p-value plotted against the rank position of the p-values
#you should see the hollow (weighted) points are close to, or ideally, above the diagonal line
plot(psate.a, plots=4)

#the weight produced above is an Inverse Probability Weight
#that is it=1/P(x|c)
#where P(x|c) is the probability of the individual's observed exposure level conditional on covariates
#however these can be susceptible to problems with excess variability
#stabilised weights are generally preferable
#these are calculated instead as =P(x)/P(x|c)
#where P(x) is just the probability of the individual's observed exposure level (not conditional on covariates)
#therefore we need to... 
#estimate this unconditional probability
anum <- svyglm(apov~1, family=quasibinomial(link = 'logit'), design=smpw)
MCS$apovp <- anum$fitted.values
MCS$aprob <- if_else(MCS$apov==1,MCS$apovp,1-MCS$apovp,missing=NULL)
#extract the propensity score from the covariate model
MCS <- bind_cols(MCS,as.vector(psate.a$ps))
MCS <- mutate(MCS,
              aatep=es.mean.ATE)
MCS <- select(MCS, -(es.mean.ATE))
MCS$aprobd <- if_else(MCS$apov==1,MCS$aatep,1-MCS$aatep,missing=NULL)
#and finally calcualte the stabilised weight
MCS <- MCS %>%
  mutate(watea=obsfwgt*(aprob/aprobd),
         atea=aprob/aprobd)
#atea is the stabilised weight on it's own
#watea combines it with the sample/attrition weight

#we can check that the stabilised weight still produces the desired covariate balance with this command
ate.a <- bal.tab(apov~ethmin+smkpreg+paredu+country+magect,
                 data=MCS,
                 weights=MCS$watea,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
#the first line specifies the exposure and covariates to check balance for
#data=MCS tells it which data file to use
#weights=specifies the new stabilised weight
#m.threshold=specifies a threshold for the standardised difference to be considered 'balanced'
#s.weights=specifies the sample/attrition weights for the 'unweighted' comparison
#estimand="ATE" tells it we want an ATE comparison
print(ate.a,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#this prints the output
#it looks slightly different to what you get with bal.table
#but the key info needed is still there
#m=Mean,SD=standard deviation
#0/1 correspond to values of the exposure variable
#diff=standardised effect sizes
#un='unweighted' ie sample/attrition weights only, while Adj=weighted with the ATE weights

#now we can estimate the ATE of apov on smkout
#we specify use of the new weight as follows
w.atea <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watea, data=MCS, fpc=~NH2)
#estimate the model using that weight
est_atea <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=w.atea)
#and print the output
summary(est_atea)


#wave 1-ATT
#to get an ATT weight we can use essentially the same command as above for the ATE
#we simply change to 
#estimand="ATT"
psatt.a <- ps(apov~ethmin+smkpreg+paredu+country+magect,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
#and we have all the same output options
plot(psatt.a)
summary(psatt.a$gbm.obj,
        n.trees=psatt.a$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.a.bal<-bal.table(psatt.a)
psatt.a.bal
summary(psatt.a)
plot(psatt.a, plots=2)
plot(psatt.a, plots=3)
plot(psatt.a, plots=4)

#for an ATT weight combined with the sample/attrition weight we can extra directly from the output like so
MCS$watta <- get.weights(psatt.a, stop.method = "es.mean")
#if for any reason however we want the ATT weight separated from the sample weight we can extract it as follows
MCS <- bind_cols(MCS,as.vector(psatt.a$ps))
MCS <- mutate(MCS,
              aattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$atta <- if_else(MCS$apov==1,1,(MCS$aattp/(1-MCS$aattp)),missing=NULL)

#now we can use the ATT weight to generate an estimate of the ATT of apov on smkout
w.atta <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watta, data=MCS, fpc=~NH2)
est_atta <- svyglm(smkout~apov, family=quasibinomial(link = 'logit'), design=w.atta)
summary(est_atta)

#wave 1-cumulative ATE
#at wave 1 the ATE weighted estimated is the same as our cumulative ATE estimates
#because we only have 1 wave of exposure to sum over

#wave 2
#wave 2-unweighted
#first we estimate the observational association with wave 2 exposure
unwgtb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtb)

#wave 2-ATE
#key differences from ATE at wave 1 are:
#we have changed the first line to indicate the wave b exposure, bpov
#we have added to the list of covariates:
#apov, the earlier exposure
#asmk and amatmh, two time-varying covariates which we believe to be
#caused by apov but
#possible causes of bpov
#asmk=maternal smoking at wave a
#amatmh=maternal mental health at wave a
#(if we thought instead that asmk & amatmh were causes of apov we would have included them at wave 1)
psate.b <- ps(bpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATE",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
#same output options again
#we want to check that the model has run with enough n.trees to reach the optimal model
#we want to check how well the model produces balance in our newly expanded confounder set
plot(psate.b)
summary(psate.b$gbm.obj,
        n.trees=psate.b$desc$es.mean.ATE$n.trees,
        plot=TRUE)
psate.b.bal<-bal.table(psate.b)
psate.b.bal
summary(psate.b)
plot(psate.b, plots=2)
plot(psate.b, plots=3)
plot(psate.b, plots=4)
#and we want to get a stabilised weight as before
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

#check balance for the stabilised weight
#changing the focal exposure, the weight, and adding the same additional covariates
ate.b <- bal.tab(bpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh,
                 data=MCS,
                 weights=MCS$wateb,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
print(ate.b,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

#then we use it to estimate the ATE of bpov on smkout
w.ateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wateb, data=MCS, fpc=~NH2)
est_ateb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=w.ateb)
summary(est_ateb)

#wave 2-ATT
#again main change from wave 1 here is just to the variable list
#we're putting bpov as the exposure, and adding apov, asmk & amatmh to the confounder list
psatt.b <- ps(bpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh,
              data=MCS,
              n.trees=20000,
              interaction.depth = 2,
              shrinkage=.01,
              perm.test.iters = 0,
              stop.method=c("es.mean"),
              estimand="ATT",
              sampw=MCS$obsfwgt,
              verbose=FALSE)
#output
plot(psatt.b)
summary(psatt.b$gbm.obj,
        n.trees=psatt.b$desc$es.mean.ATT$n.trees,
        plot=TRUE)
psatt.b.bal<-bal.table(psatt.b)
psatt.b.bal
summary(psatt.b)
plot(psatt.b, plots=2)
plot(psatt.b, plots=3)
plot(psatt.b, plots=4)
#get weights/probabilities
MCS$wattb <- get.weights(psatt.b, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.b$ps))
MCS <- mutate(MCS,
              battp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$attb <- if_else(MCS$bpov==1,1,(MCS$battp/(1-MCS$battp)),missing=NULL)

#estimate the ATT effect of bpov on smkout
w.attb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wattb, data=MCS, fpc=~NH2)
est_attb <- svyglm(smkout~bpov, family=quasibinomial(link = 'logit'), design=w.attb)
summary(est_attb)

#wave 2-cumulative ATE
#now we want to estimate the cumulative effect of poverty up to wave b
#ie the 'dose-response' effect of each instance of poverty, either apov or bpov
#for this we use a marginal structural model
#the weight we want should be=obsfwgt*atea*ateblng
#where ateblng is =P(xb|xa)/P(xb|xa,c)
#that is, it is a stabilised weight but the numerator includes prior exposure
blong <- svyglm(bpov~apov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$bpovplong <- blong$fitted.values
MCS$bproblong <- if_else(MCS$bpov==1,MCS$bpovplong,1-MCS$bpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateblng=bproblong/bprobd,
         cumateb=obsfwgt*atea*ateblng,
         cumpovb=apov+bpov)

#we want to check balance over our cumulative exposure count (cumpovb)
#notice we have changed the exposure to cumpovb

#we have removed apov from the list of covariates
#this is because apov is now part of our exposure variable
#we have changed the weight to use the new cumateb
#we have also added pairwise=FALSE which specifies 
#that we compare each level of exposure against all others
cumateb <- bal.tab(as.factor(cumpovb)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh,
                 data=MCS,
                 weights=MCS$cumateb,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE,
                 pairwise=FALSE)
#for the output we add some commands to specify that we want to display all treatment levels
print(cumateb,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat=c("0","1","2"))
#note that in this case we do not necessarily expect all the covariates to be balanced
#the cumateb weight includes both the atea and ateb weights
#thus the variables included in the atea weight should be balanced
#the new variables added for the ateb weight (ie asmk & amatmh) however,
#they will be somewhat balanced, but only to the extent that differences are not attributable to apov
#eg if some of the differences in asmk across levels of cumpovb are caused by apov then
#these will remain imbalanced to some extent (which is what we want)
#differences not caused by apov should be removed though (so you should see some improvement in balance)
#for variables included in both weights eg baseline confounders such as paredu,
#these should be balanced across levels of exposure

#then having checked appropriate balance we can estimate the cumulative ATE up to wave b
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_cumateb <- svyglm(smkout~cumpovb, family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_cumateb)
#we can relax our assumption of linearity for the cumulative effect of poverty as follows
est_cumatefb <- svyglm(smkout~as.factor(cumpovb), family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_cumatefb)
#note that the highest category here gives the effect of persistent poverty through b, relative to no poverty

#we can also use the cumulative ate weight to b to test different codings of exposure
#eg
#effect of any pov up to b
MCS$anypovb <- if_else(MCS$apov==1|MCS$bpov==1,1,0,missing=NULL)

anypovb <- bal.tab(anypovb~ethmin+smkpreg+paredu+country+magect+asmk+amatmh,
                 data=MCS,
                 weights=MCS$cumateb,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
print(anypovb,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#similarly, we expect ethmin-magect to be balanced but may expect some imbalance in asmk/amatmh

#then we can estimate the effect of any poverty up to wave b on smkout
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_anyb <- svyglm(smkout~anypovb, family=quasibinomial(link = 'logit'), design=w.cumateb)
summary(est_anyb)

#effect of experiencing poverty for the 1st time at b
#note here I use the anypovb variable as a catch-all 'other poverty' category
#this ensures the effect of 1st poverty at b is in reference to no poverty (rather than any other pattern)
MCS$firstpovb <- if_else(MCS$apov==0 & MCS$bpov==1,2,MCS$anypovb,missing=NULL)
#as far as checking balance goes we're mainly interested in checking balance between 
#0=no poverty
#2=1st poverty at b
povb1st <- bal.tab(as.factor(firstpovb)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh,
                 data=MCS,
                 weights=MCS$cumateb,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE,
                 focal="2")
#for the output we add some commands to specify that we want to display all treatment levels
print(povb1st,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat="0")

#then after checking balance, we can estimate the effect of 1st poverty experience at b
w.cumateb <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumateb, data=MCS, fpc=~NH2)
est_1stb <- svyglm(smkout~as.factor(firstpovb), family=quasibinomial(link = 'logit'), design=w.cumateb)
#note 1st poverty at b is indicated by category 2 of firstpovb
summary(est_1stb)

#wave 3
#you'll be familiar with the commands now so I'll focus on pointing out differences against earlier waves

#wave 3-'unweighted' observational association with just the sample/atrrition weights
unwgtc <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtc)

#wave 3-ATE
#changed exposure to cpov
#added bpov, bsmk and bmatmh to the covariate list
psate.c <- ps(cpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
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
plot(psate.c, plots=2)
plot(psate.c, plots=3)
plot(psate.c, plots=4)
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

#check balance of stabilised weight
#changed exposure to cpov
#added bpov, bsmk and bmatmh to the covariate list
ate.c <- bal.tab(cpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
                 data=MCS,
                 weights=MCS$watec,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
print(ate.c,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

#estimate ATE of cpov on smkout
w.atec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watec, data=MCS, fpc=~NH2)
est_atec <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=w.atec)
summary(est_atec)

#wave 3-ATT
#changed exposure to cpov
#added bpov, bsmk and bmatmh to the covariate list
psatt.c <- ps(cpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
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
plot(psatt.c, plots=2)
plot(psatt.c, plots=3)
plot(psatt.c, plots=4)
MCS$wattc <- get.weights(psatt.c, stop.method = "es.mean")
MCS <- bind_cols(MCS,as.vector(psatt.c$ps))
MCS <- mutate(MCS,
              cattp=es.mean.ATT)
MCS <- select(MCS, -(es.mean.ATT))
MCS$attc <- if_else(MCS$cpov==1,1,(MCS$cattp/(1-MCS$cattp)),missing=NULL)

#estimate ATT of cpov on smkout
w.attc <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wattc, data=MCS, fpc=~NH2)
est_attc <- svyglm(smkout~cpov, family=quasibinomial(link = 'logit'), design=w.attc)
summary(est_attc)

#wave 3 -cumulative ATE
#following from above weight should be = obsfwgt*atea*ateblng*ateclng
clong <- svyglm(cpov~apov+bpov, family=quasibinomial(link = 'logit'), design=smpw)
MCS$cpovplong <- clong$fitted.values
MCS$cproblong <- if_else(MCS$cpov==1,MCS$cpovplong,1-MCS$cpovplong,missing=NULL)
MCS <- MCS %>%
  mutate(ateclng=cproblong/cprobd,
         cumatec=obsfwgt*atea*ateblng*ateclng,
         cumpovc=apov+bpov+cpov)
#check balance over cumulative exposure count
#note we don't include apov or bpov as confounders as they're part of the exposure
cumatec <- bal.tab(as.factor(cumpovc)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh,
                   data=MCS,
                   weights=MCS$cumatec,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   pairwise=FALSE)
#for the output we add some commands to specify that we want to display all treatment levels
print(cumatec,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat=c("0","1","2","3"))
#remember we expect ethmin-magect to be balanced but time-varying confounders may still be imbalanced

#this is getting to be a lot to go through
#so if you prefer you can get a summary of greatest imbalance across all levels of the exposure like so
print(cumatec,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

#estimate cumulative ATE of poverty up to wave c
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_cumatec <- svyglm(smkout~cumpovc, family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_cumatec)
#we can relax our assumption of linearity for the cumulative effect of poverty as follows
est_cumatefc <- svyglm(smkout~as.factor(cumpovc), family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_cumatefc)
#note that the highest category here gives the effect of persistent poverty through c, relative to no poverty

#we can also use the cumulative ate weight to c to test different codings of exposure
#eg
#effect of any pov up to c
MCS$anypovc <- if_else(MCS$apov==1|MCS$bpov==1|MCS$cpov==1,1,0,missing=NULL)

anypovc <- bal.tab(anypovc~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh,
                   data=MCS,
                   weights=MCS$cumatec,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE)
print(anypovc,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

#then we can estimate the effect of any poverty up to wave c on smkout
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_anyc <- svyglm(smkout~anypovc, family=quasibinomial(link = 'logit'), design=w.cumatec)
summary(est_anyc)

#effect of experiencing poverty for the 1st time at c
#note here I use the anypovb variable as a catch-all 'other poverty' category
#this ensures the effect of 1st poverty at c is in reference to no poverty (rather than any other pattern)
MCS$firstpovc <- if_else(MCS$apov==0 & MCS$bpov==0 & MCS$cpov==1,2,MCS$anypovc,missing=NULL)
#as far as checking balance goes we're mainly interested in checking balance between 
#0=no poverty
#2=1st poverty at b
povc1st <- bal.tab(as.factor(firstpovc)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh,
                   data=MCS,
                   weights=MCS$cumatec,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   focal="2")
#for the output we add some commands to specify that we want to display all treatment levels
print(povc1st,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat="0")

#then after checking balance, we can estimate the effect of 1st poverty experience at c
w.cumatec <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatec, data=MCS, fpc=~NH2)
est_1stc <- svyglm(smkout~as.factor(firstpovc), family=quasibinomial(link = 'logit'), design=w.cumatec)
#note 1st poverty at c is indicated by category 2 of firstpovc
summary(est_1stc)

#wave 4
#wave 4-unweighted
unwgtd <- svyglm(smkout~dpov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgtd)

#wave 4-ATE
psate.d <- ps(dpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
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
plot(psate.d, plots=2)
plot(psate.d, plots=3)
plot(psate.d, plots=4)
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

ate.d <- bal.tab(dpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
                 data=MCS,
                 weights=MCS$wated,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
print(ate.d,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

w.ated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~wated, data=MCS, fpc=~NH2)
est_ated <- svyglm(smkout~dpov, family=quasibinomial(link = 'logit'), design=w.ated)
summary(est_ated)

#wave 4-ATT
psatt.d <- ps(dpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
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
plot(psatt.d, plots=2)
plot(psatt.d, plots=3)
plot(psatt.d, plots=4)
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
#check balance
cumated <- bal.tab(as.factor(cumpovd)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh,
                   data=MCS,
                   weights=MCS$cumated,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   pairwise=FALSE)
#for the output we add some commands to specify that we want to display all treatment levels
print(cumated,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat=c("0","1","2","3","4"))
#remember we expect ethmin-magect to be balanced but time-varying confounders may still be imbalanced

#this is getting to be a lot to go through
#so if you prefer you can get a summary of greatest imbalance across all levels of the exposure like so
print(cumated,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_cumated <- svyglm(smkout~cumpovd, family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_cumated)
est_cumatefd <- svyglm(smkout~as.factor(cumpovd), family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_cumatefd)
#note that the highest category here gives the effect of persistent poverty through d, relative to no poverty

#we can also use the cumulative ate weight to d to test different codings of exposure
#eg
#effect of any pov up to d
MCS$anypovd <- if_else(MCS$apov==1|MCS$bpov==1|MCS$cpov==1|MCS$dpov==1,1,0,missing=NULL)

anypovd <- bal.tab(anypovd~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh,
                   data=MCS,
                   weights=MCS$cumated,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE)
print(anypovd,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

#then we can estimate the effect of any poverty up to wave d on smkout
w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_anyd <- svyglm(smkout~anypovd, family=quasibinomial(link = 'logit'), design=w.cumated)
summary(est_anyd)

#effect of experiencing poverty for the 1st time at d
#note here I use the anypovd variable as a catch-all 'other poverty' category
#this ensures the effect of 1st poverty at d is in reference to no poverty (rather than any other pattern)
MCS$firstpovd <- if_else(MCS$apov==0 & MCS$bpov==0 & MCS$cpov==0 &MCS$dpov==1,2,MCS$anypovd,missing=NULL)
#as far as checking balance goes we're mainly interested in checking balance between 
#0=no poverty
#2=1st poverty at d
povd1st <- bal.tab(as.factor(firstpovd)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh,
                   data=MCS,
                   weights=MCS$cumated,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   focal="2")
#for the output we add some commands to specify that we want to display all treatment levels
print(povd1st,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat="0")

#then after checking balance, we can estimate the effect of 1st poverty experience at d
w.cumated <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumated, data=MCS, fpc=~NH2)
est_1std <- svyglm(smkout~as.factor(firstpovd), family=quasibinomial(link = 'logit'), design=w.cumated)
#note 1st poverty at d is indicated by category 2 of firstpovd
summary(est_1std)

#wave 5
#wave 5-unweighted
unwgte <- svyglm(smkout~epov, family=quasibinomial(link = 'logit'), design=smpw)
summary(unwgte)

#wave 5-ATE
psate.e <- ps(epov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh+dpov+dsmk+dmatmh,
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
plot(psate.e, plots=2)
plot(psate.e, plots=3)
plot(psate.e, plots=4)
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

ate.e <- bal.tab(epov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh+dpov+dsmk+dmatmh,
                 data=MCS,
                 weights=MCS$watee,
                 method="weighting",
                 binary="std",
                 m.threshold=0.1,
                 s.weights=MCS$obsfwgt,
                 estimand="ATE",
                 abs=FALSE,
                 quick=FALSE)
print(ate.e,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

w.atee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~watee, data=MCS, fpc=~NH2)
est_atee <- svyglm(smkout~epov, family=quasibinomial(link = 'logit'), design=w.atee)
summary(est_atee)

#wave 4-ATT
psatt.e <- ps(epov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk++cmatmh+dpov+dsmk+dmatmh,
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
plot(psatt.e, plots=2)
plot(psatt.e, plots=3)
plot(psatt.e, plots=4)
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

#check balance over cumulative exposure count
cumatee <- bal.tab(as.factor(cumpove)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                   data=MCS,
                   weights=MCS$cumatee,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   pairwise=FALSE)
print(cumatee,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat=c("0","1","2","3","4","5"))
#remember we expect ethmin-magect to be balanced but time-varying confounders may still be imbalanced

#this is getting to be a lot to go through
#so if you prefer you can get a summary of greatest imbalance across all levels of the exposure like so
print(cumatee,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_cumatee <- svyglm(smkout~cumpove, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cumatee)
est_cumatefe <- svyglm(smkout~as.factor(cumpove), family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_cumatefe)
#note that the highest category here gives the effect of persistent poverty through e, relative to no poverty

#we can also use the cumulative ate weight to e to test different codings of exposure
#eg
#effect of any pov up to e
MCS$anypove <- if_else(MCS$cumpove>=1,1,0,missing=NULL)

anypove <- bal.tab(anypove~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                   data=MCS,
                   weights=MCS$cumatee,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE)
print(anypove,disp.means=TRUE,disp.sds=TRUE, un=TRUE)

#then we can estimate the effect of any poverty up to wave e on smkout
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_anye <- svyglm(smkout~anypove, family=quasibinomial(link = 'logit'), design=w.cumatee)
summary(est_anye)

#effect of experiencing poverty for the 1st time at e
#note here I use the anypove variable as a catch-all 'other poverty' category
#this ensures the effect of 1st poverty at e is in reference to no poverty (rather than any other pattern)
MCS$firstpove <- if_else(MCS$apov==0 & MCS$bpov==0 & MCS$cpov==0 & MCS$dpov==0 & MCS$epov==1,2,MCS$anypove,missing=NULL)
#as far as checking balance goes we're mainly interested in checking balance between 
#0=no poverty
#2=1st poverty at d
pove1st <- bal.tab(as.factor(firstpove)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                   data=MCS,
                   weights=MCS$cumatee,
                   method="weighting",
                   binary="std",
                   m.threshold=0.1,
                   s.weights=MCS$obsfwgt,
                   estimand="ATE",
                   abs=FALSE,
                   quick=FALSE,
                   focal="2")
#for the output we add some commands to specify that we want to display all treatment levels
print(pove1st,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=FALSE,
      which.treat="0")

#then after checking balance, we can estimate the effect of 1st poverty experience at e
w.cumatee <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cumatee, data=MCS, fpc=~NH2)
est_1ste <- svyglm(smkout~as.factor(firstpove), family=quasibinomial(link = 'logit'), design=w.cumatee)
#note 1st poverty at e is indicated by category 2 of firstpove
summary(est_1ste)

###final section
#here we use weights from across the study period to estimate
#the controlled direct effects (CDEs) of exposure, not mediated via later exposure
#what we want in terms of balance for a CDE is that
#all covariates should be balanced, except when differences are caused by the exposure
#thus, considering apov:
#all pre-treatment covariates should be balanced across apov
#post-treatment covariates should be balanced across later exposure within apov
#similar to how cumalative covariates should be balanced

#examine cumatee weight for CDE of A not through B-E
#check pre-treatment balance
cdea_pre <- bal.tab(apov~ethmin+smkpreg+paredu+country+magect,
                    data=MCS,
                    weights=MCS$cumatee,
                    method="weighting",
                    binary="std",
                    m.threshold=0.1,
                    s.weights=MCS$obsfwgt,
                    estimand="ATE",
                    abs=FALSE,
                    quick=FALSE)
print(cdea_pre,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#looking for balance in covariates over b-e within strata of a
MCS <-mutate(MCS, bepov=bpov+cpov+dpov+epov)
MCS$apovr <- if_else(MCS$apov==1,0,1,missing=NULL)
cdea_post1 <- bal.tab(as.factor(bepov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                    data=MCS,
                    weights=MCS$cumatee,
                    method="weighting",
                    binary="std",
                    m.threshold=0.1,
                    s.weights=MCS$obsfwgt,
                    estimand="ATE",
                    abs=FALSE,
                    subset=as.logical(MCS$apov),
                    quick=FALSE)
print(cdea_post1,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)
cdea_post0 <- bal.tab(as.factor(bepov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cumatee,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$apovr),
                      quick=FALSE)
print(cdea_post0,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

#estimate CDE of A not through B-E
#we include an interaction between apov and bepov
#this is because definition of CDE means it can vary across levels of later exposure
#thus we first test for presence of interaction, adn then re-estimate without if non-sig.
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
#check pre-treatment balance
cdeb_pre <- bal.tab(bpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh,
                    data=MCS,
                    weights=MCS$cdeb,
                    method="weighting",
                    binary="std",
                    m.threshold=0.1,
                    s.weights=MCS$obsfwgt,
                    estimand="ATE",
                    abs=FALSE,
                    quick=FALSE)
print(cdeb_pre,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#looking for balance in covariates over c-e within strata of b
MCS$bpovr <- if_else(MCS$bpov==1,0,1,missing=NULL)
cdeb_post1 <- bal.tab(as.factor(cepov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cdeb,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$bpov),
                      quick=FALSE)
print(cdeb_post1,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)
cdeb_post0 <- bal.tab(as.factor(cepov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cdeb,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$bpovr),
                      quick=FALSE)
print(cdeb_post0,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

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
#check pre-treatment balance
cdec_pre <- bal.tab(cpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh,
                    data=MCS,
                    weights=MCS$cdec,
                    method="weighting",
                    binary="std",
                    m.threshold=0.1,
                    s.weights=MCS$obsfwgt,
                    estimand="ATE",
                    abs=FALSE,
                    quick=FALSE)
print(cdec_pre,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#looking for balance in covariates over d-e within strata of c
MCS$cpovr <- if_else(MCS$cpov==1,0,1,missing=NULL)
cdec_post1 <- bal.tab(as.factor(depov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cdec,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$cpov),
                      quick=FALSE)
print(cdec_post1,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)
cdec_post0 <- bal.tab(as.factor(depov)~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cdec,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$cpovr),
                      quick=FALSE)
print(cdec_post0,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)
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
#check pre-treatment balance
cded_pre <- bal.tab(dpov~ethmin+smkpreg+paredu+country+magect+apov+asmk+amatmh+bpov+bsmk+bmatmh+cpov+csmk+cmatmh,
                    data=MCS,
                    weights=MCS$cded,
                    method="weighting",
                    binary="std",
                    m.threshold=0.1,
                    s.weights=MCS$obsfwgt,
                    estimand="ATE",
                    abs=FALSE,
                    quick=FALSE)
print(cded_pre,disp.means=TRUE,disp.sds=TRUE, un=TRUE)
#looking for balance in covariates over e within strata of d
MCS$dpovr <- if_else(MCS$dpov==1,0,1,missing=NULL)
cded_post1 <- bal.tab(epov~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cded,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$dpov),
                      quick=FALSE)
print(cded_post1,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)
cded_post0 <- bal.tab(epov~ethmin+smkpreg+paredu+country+magect+asmk+amatmh+bsmk+bmatmh+csmk+cmatmh+dsmk+dmatmh,
                      data=MCS,
                      weights=MCS$cded,
                      method="weighting",
                      binary="std",
                      m.threshold=0.1,
                      s.weights=MCS$obsfwgt,
                      estimand="ATE",
                      abs=FALSE,
                      subset=as.logical(MCS$dpovr),
                      quick=FALSE)
print(cded_post0,disp.means=TRUE,disp.sds=TRUE, un=TRUE, 
      multi.summary=TRUE)

#estimate CDE of D not through E
w.cded <- svydesign(ids=~SPTN00, strata=~PTTYPE2, weights=~cded, data=MCS, fpc=~NH2)
est_cded <- svyglm(smkout~dpov+epov+depovi, family=quasibinomial(link = 'logit'), design=w.cded)
summary(est_cded)
est_cded <- svyglm(smkout~dpov+epov, family=quasibinomial(link = 'logit'), design=w.cded)
summary(est_cded)


