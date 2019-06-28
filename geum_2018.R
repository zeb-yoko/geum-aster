#Load data and run code that cleans 2016 & 2017 data.

setwd("C:/Users/Mason Kulbaba/Dropbox/git/geum-aster")

#load data
dat<- read.csv("cleaned_data_for_aster.csv")

#sum seed mass across all three years
dat$sm.tot<- dat$sm + dat$sm.2 + dat$sm.3


#subset data for 2017 analysis
dat2<- dat[c("Family.Unique",   "Block.ID", "HabitatType", "Region", "Population",
                          "Dist.from.cg.km","Germination.Y.N","Survival.Y.N","Surv2017", "Surv2018", 
                          "Flower.Y.N.2016","Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017",
                          "Fruit.Y.N.2016","Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017",
                          "sm", "sm.2", "sm.3", "sm.tot", "Flowering.Y.N.2018",
             "Total.Flowers.2018", "Fruit.Y.N.2018", 
             "No.Fruit.2018" )]



#make variable for "any reproduction" in 2016 through 2018
dat2$any.rep[dat2$sm.tot >0]=1

#replace NAs with 0 is "any.seeds" variable
dat2$any.rep[is.na(dat2$any.rep)] <- 0




#set response variables -> these represent variables in graphical model
vars<- c( "Germination.Y.N","Survival.Y.N","Surv2017", "Surv2018","Flower.Y.N.2016", 
          "Flower.Y.N.2017", "Flowering.Y.N.2018","No.Flowers.2016","Total.Flowers.2017",
          "Total.Flowers.2018","No.Fruit.2016","No.Fruit.2017", "No.Fruit.2018","sm", "sm.2",
          "sm.3")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#Designation of fitness variable for 2016 data
fit <- grepl("sm.tot", as.character(redata$varb))
fit<- as.numeric(fit)

redata$fit <- fit

#check
with(redata, sort(unique(as.character(varb)[fit == 0])))
with(redata, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata<- data.frame(redata, root=1)

#check class of each variable
sapply(redata, class)

#make sure Block.ID is a factor
redata$Block.ID <- as.factor(redata$Block.ID)

##########################################

# Estimate distribtuions from 2018 data (reuse 2016 and 2017 distributions)

flwno1<- dat2$No.Flowers.2016

flwno2<- dat2$Total.Flowers.2017

frt1<- dat2$No.Fruit.2016

frt2<- dat2$No.Fruit.2017

flwno<- dat2$Total.Flowers.2018

frtno<- dat2$No.Fruit.2018

sm<- dat2$sm

sm2<- dat2$sm.2

sm3<- dat2$sm.3

seeds<-dat2$sm.tot

sm3<- dat2$sm.3

library(MASS)

#2016 flower number
fl1.1<- fitdistr(flwno1, "normal")
fl1.2<- fitdistr(flwno1, "negative binomial")#size: 0.087611463
fl1.3<- fitdistr(flwno1, "poisson")

AIC(fl1.1, fl1.2, fl1.3)
fl1.2


#2017 flower number
fl2.1<- fitdistr(flwno2, "normal")
fl2.2<- fitdistr(flwno2, "negative binomial")#size: 0.30514117
fl2.3<- fitdistr(flwno2, "poisson")

AIC(fl2.1, fl2.2, fl2.3)
fl2.2


#flower number
fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.32113855
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


#2016 fruit number
frt1.1<- fitdistr(frt1, "normal")
frt1.2<- fitdistr(frt1, "negative binomial")#size: 0.027821465
frt1.3<- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2


#2017 fruit number
frt2.1<- fitdistr(frt2, "normal")
frt2.2<- fitdistr(frt2, "negative binomial")#size: 0.23720330
frt2.3<- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2


#2018 fruit number
frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.22983979
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)
frt.2

#total seed mass
seeds.1<- fitdistr(seeds, "normal")
seeds.2<- fitdistr(seeds, "negative binomial")#size: 9.216314e-02
seeds.3<- fitdistr(seeds, "poisson")

AIC(seeds.1, seeds.2, seeds.3)
seeds.2

#2016 seed mass
sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.0058024283
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2


#2017 seed mass
sm2.1<- fitdistr(sm2, "normal")
sm2.2<- fitdistr(sm2, "negative binomial")#size: 0.08457787
sm2.3<- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2


#2018 seed mass
sm3.1<- fitdistr(sm3, "normal")
sm3.2<- fitdistr(sm3, "negative binomial")#size: 6.633839e-02
sm3.3<- fitdistr(sm3, "poisson")

AIC(sm3.1, sm3.2, sm3.3)
sm3.2

#NOTE: sm, sm.2 size estimates for neg.binomial dist from code: geum_aster_2017b.R

#load aster package
library(aster)


#set custom famlist
                #1
famlist <- list(fam.bernoulli(),
                #2
                fam.negative.binomial(0.087611463),#2016 flower number
                #3
                fam.negative.binomial(0.30514117),#2017 flower number
                #4
                fam.negative.binomial(0.32113855),#2018 flower number
                #5
                fam.negative.binomial(0.027821465),#2016 fruit number
                #6
                fam.negative.binomial(0.23720330), #2017 fruit number
                #7
                fam.negative.binomial(0.22983979), #2018 fruit number
                #8
                fam.negative.binomial(0.0058024283),#sm
                #9
                fam.negative.binomial(0.08457787),#sm.2
                #10
                fam.negative.binomial(6.633839e-02))#sm.3 

#Node:   1  2  3  4  5  6  7  8  9  10  11  12  13   14   15   16
pred<- c(0, 1, 2, 3, 2, 3, 4, 5, 6, 7,  8,  9,  10,  11,  12,  13)
fam<-  c(1, 1, 1, 1, 1, 1, 1, 2, 3, 4,  5,  6,  7,   8,   9,   10)


#The above graphical model produces a direction of recession/constancy error

# try slightly reducing the model: remove all seed masses, except for 2018

#Node:   1  2  3  4  5  6  7  8  9  10  11  12  13   14   15   16
pred<- c(0, 1, 2, 3, 2, 3, 4, 5, 6, 7,  8,  9,  10,  11,  12,  13)
fam<-  c(1, 1, 1, 1, 1, 1, 1, 2, 3, 4,  5,  6,  7,   8,   9,   10)


aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aouta, show.graph=TRUE,info.tol = 1e-20)















aout<- aster(resp~varb + fit:(HabitatType), pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aout, show.graph=TRUE, info.tol = 1e-16)

anova(aouta, aout)

aoutb<- aster(resp~ varb + fit:(Region +  Block.ID), pred, fam, varb, id, root, data = redata)

summary(aoutb, show.graph=T, info.tol = 1e-10)

aoutc<- aster(resp~ varb + fit:(Region +  Block.ID+ Block.ID*Region), pred, fam, varb, id, root, data = redata)
summary(aoutc, show.graph=T)

anova(aout, aoutb, aoutc)

#############################

#To incorporate block effects into fitness estimates for each Retion (or any other factor), split 
# data into three regions

glalvar<- subset(redata, Region=="GL_alvar")
glalvar<- droplevels(glalvar)

mbalvar<- subset(redata, Region=="MB_alvar")
mbalvar<- droplevels(mbalvar)

pra<- subset(redata, Region=="Prairie")
pra<- droplevels(pra)

####aster analyses with block for each region-specific data set

aout.gla<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=glalvar,method='nlm', maxiter=5000)

aout.mba<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=mbalvar,method='nlm', maxiter=5000)

aout.pra<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=pra,method='nlm', maxiter=5000)

summary(aout.gla, show.graph = T, info.tol = 1e-11)
summary(aout.mba, show.graph = T, info.tol = 1e-11)
summary(aout.pra, show.graph = T, info.tol = 1e-11)

#Make design matrix
fred <- data.frame( Block.ID=levels(redata$Block.ID),
                   Germination.Y.N=1, Survival.Y.N=1, Survival.Y.N.2018=1,
                   Flowering.Y.N.2018=1,Total.Flowers.2018=1, 
                   Fruit.Y.N.2018=1, No.Fruit.2018=1, sm3=1, root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm3")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm3")


#Generate fintess estimates and standard errors for each block
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 8x 12 matrix (12 blocs and 8 nodes of graphicla model)


#only want means for k'th individual that contribute to expected
#fitness, and want to add only sm3 entries

foo<- grepl("sm3", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1


#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout.gla, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-10)

pout.amat.mba<- predict(aout.mba, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-10)

pout.amat.pra<- predict(aout.pra, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-10)
#combine estimates with standard error, and then round
#to three decimal places for each region
block.gla<- cbind(pout.amat.gla$fit, pout.amat.gla$se.fit)
rownames(block.gla)<- as.character(fred$Block.ID)
colnames(block.gla)<- c("Expected Fitness", "SE")

block.gla<- round(block.gla, 3) 
block.gla

#use block with median fitness value to represent region
summary(block.gla)#median = 311.98 (calculated due to even number of blocks), so use:

# Block 4: 303.627 (49.436)

block.mba<- cbind(pout.amat.mba$fit, pout.amat.mba$se.fit)
rownames(block.mba)<- as.character(fred$Block.ID)
colnames(block.mba)<- c("Expected Fitness", "SE")

block.mba<- round(block.mba, 3) 
block.mba

#use block with median fitness value to represent region
summary(block.mba)#median =117.05, corresponds to block 6 93.55 (49.185)


block.pra<- cbind(pout.amat.pra$fit, pout.amat.pra$se.fit)
rownames(block.pra)<- as.character(fred$Block.ID)
colnames(block.pra)<- c("Expected Fitness", "SE")

block.pra<- round(block.pra, 3) 
block.pra

summary(block.pra)# median = 70.04; corresponds to block 7: 77.076 (29.768)



