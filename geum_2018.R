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
                          "sm", "sm.2", "sm.3",  "Flowering.Y.N.2018",
             "Total.Flowers.2018", "Fruit.Y.N.2018", 
             "No.Fruit.2018","sm.tot" )]




#set response variables -> these represent variables in graphical model
#vars<- c( "Germination.Y.N","Survival.Y.N","Surv2017", "Surv2018","Flower.Y.N.2016", 
 #         "Flower.Y.N.2017", "Flowering.Y.N.2018","No.Flowers.2016","Total.Flowers.2017",
  #        "Total.Flowers.2018","No.Fruit.2016","No.Fruit.2017", "No.Fruit.2018","sm", "sm.2",
   #       "sm.3", "sm.tot")



vars<- c( "Germination.Y.N","Survival.Y.N","Surv2017", "Surv2018",
          "No.Flowers.2016","Total.Flowers.2017",
          "Total.Flowers.2018","No.Fruit.2016","No.Fruit.2017", 
          "No.Fruit.2018", "sm.tot")


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

# Estimate distribtuions for data

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

#total seed mass
seeds.1<- fitdistr(seeds, "normal")
seeds.2<- fitdistr(seeds, "negative binomial")#size: 9.216314e-02
seeds.3<- fitdistr(seeds, "poisson")

AIC(seeds.1, seeds.2, seeds.3)
seeds.2

#NOTE: sm, sm.2 size estimates for neg.binomial dist from code: geum_aster_2017b.R

#load aster package
library(aster)


#Note: right below is "full" graphical model, that we could not avoid direction of
#       recession/constancy. Therefore, Flowering.Y.N (2016 - 2018) removed, and 
#       analysis works correctly.

#Node:   1  2  3  4  5  6  7  8  9  10  11  12  13   14   15   16   17
#pred<- c(0, 1, 2, 3, 2, 3, 4, 5, 6, 7,  8,  9,  10,  11,  12,  13,   1)
#fam<-  c(1, 1, 1, 1, 1, 1, 1, 2, 3, 4,  5,  6,  7,   8,   9,   10,  11)


############################################################
#Below graphical model excludes Flowering.Y.N variables


vars<- c( "Germination.Y.N","Survival.Y.N","Surv2017", "Surv2018",
          "No.Flowers.2016","Total.Flowers.2017",
          "Total.Flowers.2018","No.Fruit.2016","No.Fruit.2017", 
          "No.Fruit.2018","sm", "sm.2", "sm.3", "sm.tot")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#Designation of fitness variable
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
                fam.negative.binomial(6.633839e-02),#sm.3 
                #11
                fam.negative.binomial(9.216314e-02))#sm.tot

pred<- c(0,1,2,3,2,3,4,5,6,7,8,9,10,1)
fam<- c(1,1,1,1,2,3,4,5,6,7,8,9,10,11)


aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aouta, show.graph=TRUE,info.tol = 1e-11)



aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata, famlist=famlist)

summary(aout, show.graph=TRUE, info.tol = 1e-11)

anova(aouta, aout) #significant effect of Region


#############################

#Estimate region-specific fitness as median value of fitness across all blocks (within each reagion)

#subset redata into region-specific data sets
redata.gla<- subset(redata, Region=="GL_alvar")
redata.gla<- droplevels(redata.gla)

redata.mba<- subset(redata, Region=="MB_alvar")
redata.mba<- droplevels(redata.mba)

redata.p<- subset(redata, Region=="Prairie")
redata.p<- droplevels(redata.p)

#################################################################################
# Estiamte GL_alvar fitness


dat.gla<- subset(dat2, Region=="GL_alvar")

# Estimate distribtuions for data

flwno1<- dat.gla$No.Flowers.2016

flwno2<- dat.gla$Total.Flowers.2017

frt1<- dat.gla$No.Fruit.2016

frt2<- dat.gla$No.Fruit.2017

flwno<- dat.gla$Total.Flowers.2018

frtno<- dat.gla$No.Fruit.2018

sm<- dat.gla$sm

sm2<- dat.gla$sm.2

sm3<- dat.gla$sm.3

seeds<-dat.gla$sm.tot



library(MASS)

#2016 flower number
fl1.1<- fitdistr(flwno1, "normal")
fl1.2<- fitdistr(flwno1, "negative binomial")#size: 0.150132301
fl1.3<- fitdistr(flwno1, "poisson")

AIC(fl1.1, fl1.2, fl1.3)
fl1.2


#2017 flower number
fl2.1<- fitdistr(flwno2, "normal")
fl2.2<- fitdistr(flwno2, "negative binomial")#size: 0.49295158
fl2.3<- fitdistr(flwno2, "poisson")

AIC(fl2.1, fl2.2, fl2.3)
fl2.2


#flower number
fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.52971311
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


#2016 fruit number
frt1.1<- fitdistr(frt1, "normal")
frt1.2<- fitdistr(frt1, "negative binomial")#size: 0.048512293
frt1.3<- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2


#2017 fruit number
frt2.1<- fitdistr(frt2, "normal")
frt2.2<- fitdistr(frt2, "negative binomial")#size: 0.3951908
frt2.3<- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2


#2018 fruit number
frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.36867383
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)
frt.2


#2016 seed mass
sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2


#2017 seed mass
sm2.1<- fitdistr(sm2, "normal")
sm2.2<- fitdistr(sm2, "negative binomial")#size:  1.420398e-01
sm2.3<- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2


#2018 seed mass
sm3.1<- fitdistr(sm3, "normal")
sm3.2<- fitdistr(sm3, "negative binomial")#size:  1.078905e-01
sm3.3<- fitdistr(sm3, "poisson")

AIC(sm3.1, sm3.2, sm3.3)
sm3.2

#total seed mass
seeds.1<- fitdistr(seeds, "normal")
seeds.2<- fitdistr(seeds, "negative binomial")#size: 1.540548e-01
seeds.3<- fitdistr(seeds, "poisson")

AIC(seeds.1, seeds.2, seeds.3)
seeds.2


#set famlist for GL_alvar

famlist.gla <- list(fam.bernoulli(),
                #2
                fam.negative.binomial(0.150132301),#2016 flower number
                #3
                fam.negative.binomial(0.49295158),#2017 flower number
                #4
                fam.negative.binomial(0.52971311),#2018 flower number
                #5
                fam.negative.binomial(0.048512293),#2016 fruit number
                #6
                fam.negative.binomial(0.3951908), #2017 fruit number
                #7
                fam.negative.binomial(0.36867383), #2018 fruit number
                #8
                fam.negative.binomial(0.009840565),#sm
                #9
                fam.negative.binomial(1.420398e-01),#sm.2
                #10
                fam.negative.binomial(1.078905e-01),#sm.3 
                #11
                fam.negative.binomial(1.540548e-01))#sm.tot




aout.gla<- aster(resp~varb, pred, fam, varb, id, root, data=redata.gla, famlist=famlist.gla)

summary(aout.gla, show.graph = T, info.tol = 1e-11)


aout<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata.gla, famlist=famlist.gla)

summary(aout, show.graph = T, info.tol = 1e-11)


anova(aout.gla, aout)


#Make design matrix
fred <- data.frame( Block.ID=levels(redata$Block.ID),
                    Germination.Y.N=1, Survival.Y.N=1, Surv2017=1,Surv2018=1,
                    No.Flowers.2016=1, Total.Flowers.2017=1, Flowering.Y.N.2018=1,
                    Total.Flowers.2018=1, No.Fruit.2016=1, No.Fruit.2017=1,
                    No.Fruit.2018=1, sm=1, sm.2=1, sm.3=1, sm.tot=1,root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm.tot")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm.tot")


#Generate fintess estimates and standard errors for each block
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 14 x 12 matrix (12 blocs and 8 nodes of graphicla model)


#only want means for k'th individual that contribute to expected
#fitness, and want to add only sm3 entries

foo<- grepl("sm.tot", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1


#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-11)


gla<- cbind(pout.amat.gla$fit, pout.amat.gla$se.fit)
rownames(gla)<- as.character(fred$Block.ID)
colnames(gla)<- c("Expected Fitness", "SE")

gla<- round(gla, 3) 
gla

#use block with median fitness value to represent region

summary(gla)# median= 912.3 -> corresponds to block10: 926.324 (268.496)

#######################################################################################
# Estimate MB_alvar fitness


dat.mba<- subset(dat2, Region=="MB_alvar")

# Estimate distribtuions for data

flwno1<- dat.mba$No.Flowers.2016

flwno2<- dat.mba$Total.Flowers.2017

frt1<- dat.mba$No.Fruit.2016

frt2<- dat.mba$No.Fruit.2017

flwno<- dat.mba$Total.Flowers.2018

frtno<- dat.mba$No.Fruit.2018

sm<- dat.mba$sm

sm2<- dat.mba$sm.2

sm3<- dat.mba$sm.3

seeds<-dat.mba$sm.tot



library(MASS)

#2016 flower number
fl1.1<- fitdistr(flwno1, "normal")
fl1.2<- fitdistr(flwno1, "negative binomial")#size: 0.06758479
fl1.3<- fitdistr(flwno1, "poisson")

AIC(fl1.1, fl1.2, fl1.3)
fl1.2


#2017 flower number
fl2.1<- fitdistr(flwno2, "normal")
fl2.2<- fitdistr(flwno2, "negative binomial")#size: 0.56971309
fl2.3<- fitdistr(flwno2, "poisson")

AIC(fl2.1, fl2.2, fl2.3)
fl2.2


#flower number
fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.55131225
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


#2016 fruit number
frt1.1<- fitdistr(frt1, "normal")
frt1.2<- fitdistr(frt1, "negative binomial")#size: 0.006770819
frt1.3<- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2


#2017 fruit number
frt2.1<- fitdistr(frt2, "normal")
frt2.2<- fitdistr(frt2, "negative binomial")#size: 0.30670838
frt2.3<- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2


#2018 fruit number
frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.26732056
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)
frt.2


#2016 seed mass
sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2


#2017 seed mass
sm2.1<- fitdistr(sm2, "normal")
sm2.2<- fitdistr(sm2, "negative binomial")#size: 8.597147e-02
sm2.3<- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2


#2018 seed mass
sm3.1<- fitdistr(sm3, "normal")
sm3.2<- fitdistr(sm3, "negative binomial")#size:  6.134588e-02
sm3.3<- fitdistr(sm3, "poisson")

AIC(sm3.1, sm3.2, sm3.3)
sm3.2

#total seed mass
seeds.1<- fitdistr(seeds, "normal")
seeds.2<- fitdistr(seeds, "negative binomial")#size: 0.11070444
seeds.3<- fitdistr(seeds, "poisson")

AIC(seeds.1, seeds.2, seeds.3)
seeds.2


#set famlist for GL_alvar

famlist.mba <- list(fam.bernoulli(),
                    #2
                    fam.negative.binomial(0.06758479),#2016 flower number
                    #3
                    fam.negative.binomial(0.56971309),#2017 flower number
                    #4
                    fam.negative.binomial(0.55131225),#2018 flower number
                    #5
                    fam.negative.binomial(0.006770819),#2016 fruit number
                    #6
                    fam.negative.binomial(0.30670838), #2017 fruit number
                    #7
                    fam.negative.binomial(0.26732056), #2018 fruit number
                    #8
                    fam.negative.binomial(0.009840565),#sm
                    #9
                    fam.negative.binomial(8.597147e-02),#sm.2
                    #10
                    fam.negative.binomial(6.134588e-02),#sm.3 
                    #11
                    fam.negative.binomial(0.11070444))#sm.tot




aout.mba<- aster(resp~varb, pred, fam, varb, id, root, data=redata.mba, famlist=famlist.mba)

summary(aout.mba, show.graph = T, info.tol = 1e-11)


aout<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata.mba, famlist=famlist.mba)

summary(aout, show.graph = T, info.tol = 1e-11)


anova(aout.mba, aout)#block not significant, but that's ok


#Make design matrix
fred <- data.frame( Block.ID=levels(redata$Block.ID),
                    Germination.Y.N=1, Survival.Y.N=1, Surv2017=1,Surv2018=1,
                    No.Flowers.2016=1, Total.Flowers.2017=1, Flowering.Y.N.2018=1,
                    Total.Flowers.2018=1, No.Fruit.2016=1, No.Fruit.2017=1,
                    No.Fruit.2018=1, sm=1, sm.2=1, sm.3=1, sm.tot=1,root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm.tot")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm.tot")


#Generate fintess estimates and standard errors for each block
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 14 x 12 matrix (12 blocs and 8 nodes of graphicla model)


#only want means for k'th individual that contribute to expected
#fitness, and want to add only sm3 entries

foo<- grepl("sm.tot", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1


#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-11)


mba<- cbind(pout.amat.gla$fit, pout.amat.gla$se.fit)
rownames(mba)<- as.character(fred$Block.ID)
colnames(mba)<- c("Expected Fitness", "SE")

mba<- round(mba, 3) 
mba

#use block with median fitness value to represent region

summary(mba)# median= 242.30 -> corresponds to block6: 248.6 (194.65)

#######################################################################################
# Estimate Prairie fitness


dat.p<- subset(dat2, Region=="Prairie")

# Estimate distribtuions for data

flwno1<- dat.p$No.Flowers.2016

flwno2<- dat.p$Total.Flowers.2017

frt1<- dat.p$No.Fruit.2016

frt2<- dat.p$No.Fruit.2017

flwno<- dat.p$Total.Flowers.2018

frtno<- dat.p$No.Fruit.2018

sm<- dat.p$sm

sm2<- dat.p$sm.2

sm3<- dat.p$sm.3

seeds<-dat.p$sm.tot



library(MASS)

#2016 flower number
fl1.1<- fitdistr(flwno1, "normal")
fl1.2<- fitdistr(flwno1, "negative binomial")#size: 0.028979806
fl1.3<- fitdistr(flwno1, "poisson")

AIC(fl1.1, fl1.2, fl1.3)
fl1.2


#2017 flower number
fl2.1<- fitdistr(flwno2, "normal")
fl2.2<- fitdistr(flwno2, "negative binomial")#size: 0.13857497
fl2.3<- fitdistr(flwno2, "poisson")

AIC(fl2.1, fl2.2, fl2.3)
fl2.2


#flower number
fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.11756156
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


#2016 fruit number
frt1.1<- fitdistr(frt1, "normal")
frt1.2<- fitdistr(frt1, "negative binomial")#size: 0.003323081
frt1.3<- fitdistr(frt1, "poisson")

AIC(frt1.1, frt1.2, frt1.3)
frt1.2


#2017 fruit number
frt2.1<- fitdistr(frt2, "normal")
frt2.2<- fitdistr(frt2, "negative binomial")#size: 0.11502653
frt2.3<- fitdistr(frt2, "poisson")

AIC(frt2.1, frt2.2, frt2.3)
frt2.2


#2018 fruit number
frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.09239198
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)
frt.2


#2016 seed mass
sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.009840565
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)
sm.2


#2017 seed mass
sm2.1<- fitdistr(sm2, "normal")
sm2.2<- fitdistr(sm2, "negative binomial")#size: 0.034831540
sm2.3<- fitdistr(sm2, "poisson")

AIC(sm2.1, sm2.2, sm2.3)
sm2.2


#2018 seed mass
sm3.1<- fitdistr(sm3, "normal")
sm3.2<- fitdistr(sm3, "negative binomial")#size:  0.024156866
sm3.3<- fitdistr(sm3, "poisson")

AIC(sm3.1, sm3.2, sm3.3)
sm3.2

#total seed mass
seeds.1<- fitdistr(seeds, "normal")
seeds.2<- fitdistr(seeds, "negative binomial")#size: 3.642185e-02
seeds.3<- fitdistr(seeds, "poisson")

AIC(seeds.1, seeds.2, seeds.3)
seeds.2


#set famlist for Prairie analysis

famlist.p <- list(fam.bernoulli(),
                    #2
                    fam.negative.binomial(0.028979806),#2016 flower number
                    #3
                    fam.negative.binomial(0.13857497),#2017 flower number
                    #4
                    fam.negative.binomial(0.11756156),#2018 flower number
                    #5
                    fam.negative.binomial(0.003323081),#2016 fruit number
                    #6
                    fam.negative.binomial(0.11502653), #2017 fruit number
                    #7
                    fam.negative.binomial(0.09239198), #2018 fruit number
                    #8
                    fam.negative.binomial(0.009840565),#sm
                    #9
                    fam.negative.binomial(0.034831540),#sm.2
                    #10
                    fam.negative.binomial(0.024156866),#sm.3 
                    #11
                    fam.negative.binomial(3.642185e-02))#sm.tot




aout.p<- aster(resp~varb, pred, fam, varb, id, root, data=redata.p, famlist=famlist.p)

summary(aout.p, show.graph = T, info.tol = 1e-11)


aout<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata.p, famlist=famlist.p)

summary(aout, show.graph = T, info.tol = 1e-11)


anova(aout.p, aout)#block not significant, but that's ok


#Make design matrix
fred <- data.frame( Block.ID=levels(redata$Block.ID),
                    Germination.Y.N=1, Survival.Y.N=1, Surv2017=1,Surv2018=1,
                    No.Flowers.2016=1, Total.Flowers.2017=1, Flowering.Y.N.2018=1,
                    Total.Flowers.2018=1, No.Fruit.2016=1, No.Fruit.2017=1,
                    No.Fruit.2018=1, sm=1, sm.2=1, sm.3=1, sm.tot=1,root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm.tot")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm.tot")


#Generate fintess estimates and standard errors for each block
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 14 x 12 matrix (12 blocs and 8 nodes of graphicla model)


#only want means for k'th individual that contribute to expected
#fitness, and want to add only sm3 entries

foo<- grepl("sm.tot", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1


#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-11)


pra<- cbind(pout.amat.gla$fit, pout.amat.gla$se.fit)
rownames(pra)<- as.character(fred$Block.ID)
colnames(pra)<- c("Expected Fitness", "SE")

pra<- round(pra, 3) 
pra

#use block with median fitness value to represent region

summary(pra)# median= 132.26 -> corresponds to block4: 94.445 (100.454)


