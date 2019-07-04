
# This code is an attempt to "redo" the 2017 data cleaning, and subsequent analysis,
# in an attempt to resolve whatever issue was causing all estimates of fitness across
# all levels of a factors to be identical (and standard error)

# data to be used is the second data set sent by Zeb, that now includes measure 
# of the distance from seed source to common garden in North Dakota


setwd("C:/Users/Mason Kulbaba/Dropbox/git/geum-aster")

#load data
#dat<- read.csv("NV_CG_Experiment2.csv")

dat<- read.csv("cleaned_data_for_aster.csv")

#subset data for 2017 analysis
dat2<- dat[c("Family.Unique",   "Block.ID", "HabitatType", "Region", "Population",
             "Dist.from.cg.km","Germination.Y.N","Survival.Y.N","Survival.Y.N.2017", 
             "Flower.Y.N.2016","Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017",
             "Fruit.Y.N.2016","Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017",
             "sm", "sm.2", "Surv2017", "sm2017")]



######################
vars<- c("Germination.Y.N", "Survival.Y.N","Surv2017", "Flower.Y.N.2016","Flower.Y.N.2017",
         "No.Flowers.2016","Total.Flowers.2017","No.Fruit.2016", "No.Fruit.2017","sm", "sm.2", "sm2017")


#look into distributions for nodes

#recall nodes of graphical model
vars

#Isolate non bernoulli varibles, and prepare to test for 

flwno<- dat2$No.Flowers.2016


flwno2<- dat2$Total.Flowers.2017

frtno<- dat2$No.Fruit.2016

frtno2<- dat2$No.Fruit.2017

seeds<-dat2$sm2017

sm<-dat$sm

sm2<- dat$sm.2

sm2017<- dat$sm2017

library(MASS)

#library(fitdistrplus)
#library(gamlss) #to include ZIP etc.

#2016 flw no
fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.087611463
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2

#2017 flw no
fl2.1<- fitdistr(flwno2, "normal")
fl2.2<- fitdistr(flwno2, "negative binomial")#size: 0.30514117
fl2.3<- fitdistr(flwno2, "poisson")

AIC(fl2.1, fl2.2, fl2.3)
fl2.2

#2016 fruit number
frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.027821465
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)
frt.2

#2017 fruit number
frt2.1<- fitdistr(frtno2, "normal")
frt2.2<- fitdistr(frtno2, "negative binomial")#size: 0.23720330
frt2.3<- fitdistr(frtno2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2

#seeds set (2016 + 2017)
seed.1<- fitdistr(seeds, "normal")
seed.2<- fitdistr(seeds, "negative binomial")#size: 8.544092e-02
seed.3<- fitdistr(seeds, "poisson")

AIC(seed.1, seed.2, seed.3)
seed.2

#sm
sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.0058024283
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2

#sm.2
sm.2.1<- fitdistr(sm2, "normal")
sm.2.2<- fitdistr(sm2, "negative binomial")#size: 0.08457787
sm.2.3<- fitdistr(sm2, "poisson")

AIC(sm.2.1, sm.2.2, sm.2.3)
sm.2.2


#sm2017
sm2017.2.1<- fitdistr(sm2017, "normal")
sm2017.2.2<- fitdistr(sm2017, "negative binomial")#size: 8.544092e-02
sm2017.2.3<- fitdistr(sm2017, "poisson")

AIC(sm2017.2.1, sm2017.2.2, sm2017.2.3)
sm2017.2.2

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2017 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of fitness variable for 2016 data
fit <- grepl("sm2017", as.character(redata2017$varb))
fit<- as.numeric(fit)

redata2017$fit <- fit

#check
with(redata2017, sort(unique(as.character(varb)[fit == 0])))
with(redata2017, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2017<- data.frame(redata2017, root=1)

#check classes of redata2017
sapply(redata2017, class)

#make block.id a factor

redata2017$Block.ID<- as.factor(redata2017$Block.ID)


#load aster package
library(aster)


#set up custom family list

famlist <- list(fam.bernoulli(),
                fam.negative.binomial(0.087611463),
                fam.negative.binomial(0.30514117),
                fam.negative.binomial(0.027821465),
                fam.negative.binomial(0.23720330), 
                fam.negative.binomial(0.0058024283),
                fam.negative.binomial(0.08457787),
                fam.negative.binomial(8.544092e-02))





pred<- c(0,1,2,2,3,4,5,6,7,8,9,1)

fam<- c(1,1,1,1,1,2,3,4,5,6,7,8)
#sapply(fam.default(), as.character)[fam]

#fixed effect model for 2017 with only fitness: note the use of 'famlist'
aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata2017,famlist = famlist)

summary(aouta, show.graph=T, info.tol = 1e-10)

#include HabitatType in model
aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2017, famlist=famlist)


summary(aout, show.graph = TRUE, info.tol=1e-10)

anova(aouta, aout)#Region is significant

aoutb<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2017, famlist=famlist)

summary(aoutb, show.graph=T, info.tol=1e-11)

anova(aouta, aoutb)#block on it's own is not significant



###############################################

# Split redata2017 into region-specific data, add block.id, and estiamte
# mean fitness to habitat-specific aster models, and generate fitness estiamtes

levels(redata2017$Region)


redata.gla<- subset(redata2017, Region=="GL_alvar")
redata.gla<- droplevels(redata.gla)

redata.mba<- subset(redata2017, Region=="MB_alvar")
redata.mba<- droplevels(redata.mba)

redata.p<- subset(redata2017, Region=="Prairie")
redata.p<- droplevels(redata.p)

#make sure block.id still a factor
redata.gla$Block.ID<- as.factor(redata.gla$Block.ID)
redata.mba$Block.ID<- as.factor(redata.mba$Block.ID)
redata.p$Block.ID<- as.factor(redata.p$Block.ID)

## split origional data file (dat2) into alvar and prairie data to
# estiamte distribution parameters for alvar and prairie subsets

dat2.gla<- subset(dat2, Region=="GL_alvar")
dat2.gla<- droplevels(dat2.gla)

dat2.mba<- subset(dat2, Region=="MB_alvar")
dat2.mba<- droplevels(dat2.mba)

dat2.p<- subset(dat2, Region=="Prairie")
dat2.p<- droplevels(dat2.p)


#begin estimating 'alvar' distributions

flwno1<- dat2.gla$No.Flowers.2016

flwno2<- dat2.gla$Total.Flowers.2017

frt1<- dat2.gla$No.Fruit.2016

frt2<- dat2.gla$No.Fruit.2017

sm<- dat2.gla$sm

sm2<- dat2.gla$sm.2

sm2017<- dat2.gla$sm2017


#flwno1
flwno1.1 <- fitdistr(flwno1, "normal")
flwno1.2 <- fitdistr(flwno1, "negative binomial")#size: 0.150132301
flwno1.3 <- fitdistr(flwno1, "poisson")

AIC(flwno1.1, flwno1.2, flwno1.3)
flwno1.2

#flwno2
flwno2.1 <- fitdistr(flwno2, "normal")
flwno2.2 <- fitdistr(flwno2, "negative binomial")#size: 0.49295158
flwno2.3 <- fitdistr(flwno2, "poisson")

AIC(flwno2.1, flwno2.2, flwno2.3)
flwno2.2

#frt1
frt1.1 <- fitdistr(frt1, "normal")
frt1.2 <- fitdistr(frt1, "negative binomial")#size: 0.048512293
frt1.3 <- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2

#frt2
frt2.1 <- fitdistr(frt2, "normal")
frt2.2 <- fitdistr(frt2, "negative binomial")#size: 0.3951908
frt2.3 <- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2

#sm
sm.1 <- fitdistr(sm, "normal")
sm.2 <- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3 <- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2

#sm.2
sm2.1 <- fitdistr(sm2, "normal")
sm2.2 <- fitdistr(sm2, "negative binomial")#size: 1.420398e-01
sm2.3 <- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2


#sm2017
sm2017.1 <- fitdistr(sm2017, "normal")
sm2017.2 <- fitdistr(sm2017, "negative binomial")#size:  1.433965e-01
sm2017.3 <- fitdistr(sm2017, "poisson")

AIC(sm2017.1, sm2017.2, sm2017.3)
sm2017.2

#make new famlist for alvar data
famlist.gla <- list(fam.bernoulli(),
                fam.negative.binomial(0.150132301),
                fam.negative.binomial(0.49295158),
                fam.negative.binomial(0.048512293),
                fam.negative.binomial(0.3951908), 
                fam.negative.binomial(0.009840565),
                fam.negative.binomial(1.420398e-01),
                fam.negative.binomial(1.433965e-01))


#alvar aster analysis with only fitness data
aout.a1<- aster(resp~varb, pred, fam, varb, id, root, data=redata.gla,famlist = famlist.gla)

summary(aout.a1, show.graph=T, info.tol=1e-10)



aout<- aster(resp~varb +fit:(Block.ID), pred, fam, varb, id, root, data=redata.gla,famlist = famlist.gla)

summary(aout, show.graph=T, info.tol = 1e-10)

anova(aout.a1, aout)#block not significant for GL_alvar, but that's ok



# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-10)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Block.ID=levels(redata.gla$Block.ID),
                   Germination.Y.N=1, Survival.Y.N=1,Surv2017=1, Flower.Y.N.2016=1,
                   No.Flowers.2016=1, Flower.Y.N.2017=1, Total.Flowers.2017=1, 
                   No.Fruit.2016=1, No.Fruit.2017=1, sm=1, sm.2=1, sm2017=1, root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm2017")

#charlie add
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm2017")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)


#Generate fintess estimates and standard errors for each block
nReg<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nReg, nnode, nReg))
dim(amat)# makes an 12 x 12 x 12 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm2017", vars)
for(k in 1:nReg)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol=1e-10)

#combine estimates with standard error, and then round
#to three decimal places
gl.a<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(gl.a)<- as.character(fred$Block.ID)


colnames(gl.a)<- c("Expected Fitness", "SE")

#Expected fitness for blocks in GL_alvar region

#as block does not explain a sifnificant amount of variation, can omit and used following estimates
round(gl.a, 3) 

summary(gl.a)#median= 615.8 -> corresponds with block 9: 576.436 (125.247)


#########################################################################
# MB_alvar analysis

#begin estimating 'alvar' distributions

flwno1<- dat2.mba$No.Flowers.2016

flwno2<- dat2.mba$Total.Flowers.2017

frt1<- dat2.mba$No.Fruit.2016

frt2<- dat2.mba$No.Fruit.2017

sm<- dat2.mba$sm

sm2<- dat2.mba$sm.2

sm2017<- dat2.mba$sm2017


#flwno1
flwno1.1 <- fitdistr(flwno1, "normal")
flwno1.2 <- fitdistr(flwno1, "negative binomial")#size: 0.06758479
flwno1.3 <- fitdistr(flwno1, "poisson")

AIC(flwno1.1, flwno1.2, flwno1.3)
flwno1.2

#flwno2
flwno2.1 <- fitdistr(flwno2, "normal")
flwno2.2 <- fitdistr(flwno2, "negative binomial")#size: 0.56971309
flwno2.3 <- fitdistr(flwno2, "poisson")

AIC(flwno2.1, flwno2.2, flwno2.3)
flwno2.2

#frt1
frt1.1 <- fitdistr(frt1, "normal")
frt1.2 <- fitdistr(frt1, "negative binomial")#size: 0.006770819
frt1.3 <- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2

#frt2
frt2.1 <- fitdistr(frt2, "normal")
frt2.2 <- fitdistr(frt2, "negative binomial")#size: 0.30670838
frt2.3 <- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2

#sm
sm.1 <- fitdistr(sm, "normal")
sm.2 <- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3 <- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2

#sm.2
sm2.1 <- fitdistr(sm2, "normal")
sm2.2 <- fitdistr(sm2, "negative binomial")#size: 1.420398e-01
sm2.3 <- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2

#sm2017
sm2017.1 <- fitdistr(sm2017, "normal")
sm2017.2 <- fitdistr(sm2017, "negative binomial")#size: 8.883848e-02
sm2017.3 <- fitdistr(sm2017, "poisson")

AIC(sm2017.1, sm2017.2, sm2017.3)
sm2017.2

#make new famlist for alvar data
famlist.mba <- list(fam.bernoulli(),
                    fam.negative.binomial(0.06758479),
                    fam.negative.binomial(0.56971309),
                    fam.negative.binomial(0.006770819),
                    fam.negative.binomial(0.30670838), 
                    fam.negative.binomial(0.009840565),
                    fam.negative.binomial(1.420398e-01),
                    fam.negative.binomial(8.883848e-02))


#alvar aster analysis with only fitness data
aout.a1<- aster(resp~varb, pred, fam, varb, id, root, data=redata.mba,famlist = famlist.mba)

summary(aout.a1, show.graph=T, info.tol=1e-10)



aout<- aster(resp~varb +fit:(Block.ID), pred, fam, varb, id, root, data=redata.mba,famlist = famlist.mba)

summary(aout, show.graph=T, info.tol = 1e-10)

anova(aout.a1, aout)#block not sig. for mb_alvar, but that's ok



# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-10)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Block.ID=levels(redata.mba$Block.ID),
                   Germination.Y.N=1, Survival.Y.N=1,Surv2017=1, Flower.Y.N.2016=1,
                   No.Flowers.2016=1, Flower.Y.N.2017=1, Total.Flowers.2017=1, 
                   No.Fruit.2016=1, No.Fruit.2017=1, sm=1, sm.2=1, sm2017=1,root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm2017")

#charlie add
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm2017")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)


#Generate fintess estimates and standard errors for each block
nReg<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nReg, nnode, nReg))
dim(amat)# makes an 12 x 12 x 12 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm2017", vars)
for(k in 1:nReg)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol=1e-10)

#combine estimates with standard error, and then round
#to three decimal places
mb.a<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(mb.a)<- as.character(fred$Block.ID)


colnames(mb.a)<- c("Expected Fitness", "SE")

#Expected fitness for blocks in GL_alvar region

#as block does not explain a sifnificant amount of variation, can omit and used following estimates
round(mb.a, 3) 

summary(mb.a)#median: 94.55 ->corresponds to block 4: 93.6 (64.409)

#######################################################################################
## Now do Prairie analysis

#begin estimating 'alvar' distributions

flwno1<- dat2.p$No.Flowers.2016

flwno2<- dat2.p$Total.Flowers.2017

frt1<- dat2.p$No.Fruit.2016

frt2<- dat2.p$No.Fruit.2017

sm<- dat2.p$sm

sm2<- dat2.p$sm.2

sm2017<- dat2.p$sm2017


#flwno1
flwno1.1 <- fitdistr(flwno1, "normal")
flwno1.2 <- fitdistr(flwno1, "negative binomial")#size: 0.028979806
flwno1.3 <- fitdistr(flwno1, "poisson")

AIC(flwno1.1, flwno1.2, flwno1.3)
flwno1.2

#flwno2
flwno2.1 <- fitdistr(flwno2, "normal")
flwno2.2 <- fitdistr(flwno2, "negative binomial")#size: 0.13857497
flwno2.3 <- fitdistr(flwno2, "poisson")

AIC(flwno2.1, flwno2.2, flwno2.3)
flwno2.2

#frt1
frt1.1 <- fitdistr(frt1, "normal")
frt1.2 <- fitdistr(frt1, "negative binomial")#size: 0.003323081
frt1.3 <- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2

#frt2
frt2.1 <- fitdistr(frt2, "normal")
frt2.2 <- fitdistr(frt2, "negative binomial")#size: 0.11502653
frt2.3 <- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2

#sm
sm.1 <- fitdistr(sm, "normal")
sm.2 <- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3 <- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2

#sm.2
sm2.1 <- fitdistr(sm2, "normal")
sm2.2 <- fitdistr(sm2, "negative binomial")#size: 0.034831540
sm2.3 <- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2

#sm2017
sm2017.1 <- fitdistr(sm2017, "normal")
sm2017.2 <- fitdistr(sm2017, "negative binomial")#size: 0.03511662
sm2017.3 <- fitdistr(sm2017, "poisson")

AIC(sm2017.1, sm2017.2, sm2017.3)
sm2017.2

#make new famlist for alvar data
famlist.p <- list(fam.bernoulli(),
                    fam.negative.binomial(0.028979806),
                    fam.negative.binomial(0.13857497),
                    fam.negative.binomial(0.003323081),
                    fam.negative.binomial(0.11502653), 
                    fam.negative.binomial(0.009840565),
                    fam.negative.binomial(0.034831540),
                    fam.negative.binomial(0.03511662))


#alvar aster analysis with only fitness data
aout.a1<- aster(resp~varb, pred, fam, varb, id, root, data=redata.p,famlist = famlist.p)

summary(aout.a1, show.graph=T, info.tol=1e-10)



aout<- aster(resp~varb +fit:(Block.ID), pred, fam, varb, id, root, data=redata.p,famlist = famlist.p)

summary(aout, show.graph=T, info.tol = 1e-10)

anova(aout.a1, aout)#block not sig. for prairie, but that's ok



# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-10)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Block.ID=levels(redata.p$Block.ID),
                   Germination.Y.N=1, Survival.Y.N=1,Surv2017=1, Flower.Y.N.2016=1,
                   No.Flowers.2016=1, Flower.Y.N.2017=1, Total.Flowers.2017=1, 
                   No.Fruit.2016=1, No.Fruit.2017=1, sm=1, sm.2=1,sm2017=1, root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm2017")

#charlie add
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm2017")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)


#Generate fintess estimates and standard errors for each block
nReg<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nReg, nnode, nReg))
dim(amat)# makes an 12 x 12 x 12 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm2017", vars)
for(k in 1:nReg)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol=1e-10)

#combine estimates with standard error, and then round
#to three decimal places
pr<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(pr)<- as.character(fred$Block.ID)


colnames(pr)<- c("Expected Fitness", "SE")

#Expected fitness for blocks in GL_alvar region

#as block does not explain a sifnificant amount of variation, can omit and used following estimates
round(pr, 3) 

summary(pr)#median: 68.02 -> corresponds to block 8: 70.379 (75.872)



