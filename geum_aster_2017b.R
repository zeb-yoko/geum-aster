
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

#make variable for "any fruit" in 2016 or 2017

dat2$any.rep[dat2$sm2017 > 0 | dat2$No.Fruit.2016 >0 | dat2$No.Fruit.2017>0]=1

#replace NAs with 0 is "any.seeds" variable

dat2$any.rep[is.na(dat2$any.rep)] <- 0


######################
vars<- c("Germination.Y.N", "Survival.Y.N","Surv2017", "Flower.Y.N.2016","Flower.Y.N.2017",
         "No.Flowers.2016","Total.Flowers.2017","No.Fruit.2016", "No.Fruit.2017","sm", "sm.2")


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


#load aster package
library(aster)


#set up custom family list

famlist <- list(fam.bernoulli(),
                fam.negative.binomial(0.087611463),
                fam.negative.binomial(0.30514117),
                fam.negative.binomial(0.027821465),
                fam.negative.binomial(0.23720330), 
                fam.negative.binomial(0.0058024283),
                fam.negative.binomial(0.08457787))





pred<- c(0,1,2,2,3,4,5,6,7,8,9)

fam<- c(1,1,1,1,1,2,3,4,5,6,7)
#sapply(fam.default(), as.character)[fam]

#fixed effect model for 2017 with only fitness: note the use of 'famlist'
aouta<- aster(resp~varb + 0, pred, fam, varb, id, root, data=redata2017,famlist = famlist)

summary(aouta, show.graph=T, info.tol = 1e-10)

#include HabitatType in model
aout<- aster(resp~varb + fit:(HabitatType), pred, fam, varb, id, root, data=redata2017, famlist=famlist)


summary(aout, show.graph = TRUE, info.tol=1e-10)

anova(aouta, aout)#HabitatType is significant

aoutb<- aster(resp~varb + fit:(HabitatType + Block.ID), pred, fam, varb, id, root, data=redata2017, famlist=famlist)

summary(aoutb, show.graph=T, info.tol=1e-11)

anova(aouta, aout, aoutb)#block does not add anything.


# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-10)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(HabitatType=levels(redata2017$HabitatType),
                   Germination.Y.N=1, Survival.Y.N=1,Surv2017=1, No.Flowers.2016=1, Total.Flowers.2017=1, 
                   No.Fruit.2016=1, No.Fruit.2017=1, any.rep=1, sm2017=1, root = 1)

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

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm2017")


#Generate fintess estimates and standard errors for each block
nHab<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nHab, nnode, nHab))
dim(amat)# makes an 3 x 9 x 3 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm2017", vars)
for(k in 1:nHab)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol=1e-10)

#combine estimates with standard error, and then round
#to three decimal places
Hab<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(Hab)<- as.character(fred$HabitatType)


colnames(Hab)<- c("Expected Fitness", "SE")

#Expected fitness for each region

#as block does not explain a sifnificant amount of variation, can omit and used following estimates
round(Hab, 3) 





