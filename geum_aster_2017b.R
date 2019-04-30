
# This code is an attempt to "redo" the 2017 data cleaning, and subsequent analysis,
# in an attempt to resolve whatever issue was causing all estimates of fitness across
# all levels of a factors to be identical (and standard error)

# data to be used is the second data set sent by Zeb, that now includes measure 
# of the distance from seed source to common garden in North Dakota


setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/aster-analysis/")

#load data
dat<- read.csv("NV_CG_Experiment2wdist2.csv")


#subset data for 2017 analysis
dat2<- dat[c("Family.Unique",   "Block.ID", "HabitatType", "Region", "Population",
             "Dist.from.cg.km","Germination.Y.N","Survival.Y.N","Survival.Y.N.2017", 
             "Flower.Y.N.2016","Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017",
             "Fruit.Y.N.2016","Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017",
             "sm", "sm.2")]


#replace NAs with 0 (zero)
dat2[is.na(dat2)] <- 0

#make sum of seedmass in 2017 and 2017 a new variable
#dat2$sm.2017<-dat2$sm + dat2$sm.2


########################################################
########################################################
## Begin cleaning data-> start with 2016 cleaning code##
########################################################
########################################################


#survival but no germination
subset(dat2, Germination.Y.N==0 & Survival.Y.N==1)# no errors

#flowered but no survival
subset(dat2, Flower.Y.N.2016==1 & Survival.Y.N==0)# one error

#fix error: set survival to 1 for this plant
dat2$Survival.Y.N[dat2$Flower.Y.N.2016==1 & dat2$Survival.Y.N==0]=1

#flw no >0 but flw.y.n ==0
subset(dat2, Flower.Y.N.2016 ==0 & No.Flowers.2016 >0)# 8 errors

#fix errors: set flw.y.n to 1 for these plants
dat2$Flower.Y.N.2016[dat2$Flower.Y.N.2016 ==0 & dat2$No.Flowers.2016 >0]=1

#fruit.y.n=1 but flow.no ==0
subset(dat2, Fruit.Y.N.2016==1 & No.Flowers.2016==0)# 1 error

#fix error: set No.Flowers.2016 to 1
dat2$No.Flowers.2016[dat2$Fruit.Y.N.2016==1 & dat2$No.Flowers.2016==0]=1

#frt. no >0 but fruit.y.n =0
subset(dat2, No.Fruit.2016 >0 & Fruit.Y.N.2016==0)# no errors

#sm >0 but fruit no. =0
subset(dat2, sm >0 & No.Fruit.2016==0)#no errors

#####################################################################################
#####################################################################################
## The above cleaned data was tested and gives appropriate resutls for estimates of##
## mean fitness across different levels of factors. Now move to 2017 data cleaning###
#####################################################################################
#####################################################################################

#suvive to 2017 but surv.2016=0
subset(dat2, Survival.Y.N==0 & Survival.Y.N.2017==1)# no errors

##################### NEW
dat2$Survival.Y.N[dat2$Survival.Y.N==0 & dat2$Survival.Y.N.2017==1]=1
####################

#flw.no >0, but flw.y.n=0
subset(dat2, Total.Flowers.2017 >0 & Flower.Y.N.2017==0)# 10 errors

#set flw.y.n=1, but this introduces 3 new erorrs to survival, so recorrect
dat2$Flower.Y.N.2017[dat2$Total.Flowers.2017 >0 & dat2$Flower.Y.N.2017==0]=1








#flowered but no survival
subset(dat2, Survival.Y.N.2017==0 & Flower.Y.N.2017==1)# 17 errors

#fix errors: set survival.y.n to 1
dat2$Survival.Y.N.2017[dat2$Survival.Y.N.2017==0 & dat2$Flower.Y.N.2017==1]=1





#recheck for surv.2017=0 & flw.y.n.2017 =1
subset(dat2, Survival.Y.N.2017==0 & Flower.Y.N.2017==1)# 3 new errors

#fix errors: set survival.y.n to 1 to correct these new errors
dat2$Survival.Y.N.2017[dat2$Survival.Y.N.2017==0 & dat2$Flower.Y.N.2017==1]=1

#fruit.y.n =1, flw.no=0
subset(dat2, Fruit.Y.N.2017==1 & Total.Flowers.2017 ==0)# 7 errors

#fix errors: set total.flw =1: this should be checked by Zeb
dat2$Total.Flowers.2017[dat2$Fruit.Y.N.2017==1 & dat2$Total.Flowers.2017 ==0]=1

#frt.no >0, but frt.y.n=0
subset(dat2, Fruit.Y.N.2017==0 & No.Fruit.2017 > 0)# many errors

#fix errors: set frt.y.n=1 for these erors -> this should be checked by Zeb
dat2$Fruit.Y.N.2017[dat2$Fruit.Y.N.2017==0 & dat2$No.Fruit.2017 > 0]=1

#sm.2 >0 but frt.no =0
subset(dat2, sm.2 > 0 & No.Fruit.2017 ==0)# 5 error

#fix errors: set frut.No=1, but will introduce 1 new error with frut.y.n, so recheck
dat2$No.Fruit.2017[dat2$sm.2 > 0 & dat2$No.Fruit.2017 ==0]=1

#recheck frt.no >0, but frt.y.n=0
subset(dat2, Fruit.Y.N.2017==0 & No.Fruit.2017 > 0)# 1 new introduced error

#fix error: set frt.y.n=1
dat2$Fruit.Y.N.2017[dat2$Fruit.Y.N.2017==0 & dat2$No.Fruit.2017 > 0]=1

###########################################
#Begin 2017 analysis

#preliminaries: make 2017 seed mass variable (2016 seedmass + 2017 seedmass)

dat2$sm2017<- dat2$sm + dat2$sm.2


vars<- c("Germination.Y.N", "Survival.Y.N","Flower.Y.N.2016", "No.Flowers.2016", 
         "Fruit.Y.N.2016","No.Fruit.2016", "Survival.Y.N.2017","Flower.Y.N.2017",
         "Total.Flowers.2017","Fruit.Y.N.2017", "No.Fruit.2017","sm2017")

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

pred<- c(0,1,2,3,4,5,6,7)
fam<- c(1,1,1,1,2,1,2,2)
#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model for 2015 with only fitness
aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata2017)

summary(aouta, show.graph=T, info.tol=1e-9)

aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2017, method='nlm')


summary(aout, show.graph = TRUE, info.tol=1e-10)


anova(aouta, aout)

# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-10)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Region=levels(redata2017$Region),
                   Germination.Y.N=1, Survival.Y.N=1,Survival.Y.N.2017=1, Flower.Y.N.2017=1,
                   Total.Flowers.2017=1, Fruit.Y.N.2017=1, No.Fruit.2017=1, sm.2=1, root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm.2")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)

renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm.2")


#Generate fintess estimates and standard errors for each block
nRegion<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nRegion, nnode, nRegion))
dim(amat)# makes an 3 x 6 x 3 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm.2", vars)
for(k in 1:nRegion)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol=1e-10)

#combine estimates with standard error, and then round
#to three decimal places
RegionType<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(RegionType)<- as.character(fred$Region)


colnames(RegionType)<- c("Expected Fitness", "SE")

#Oof! This produces identical estiamtes for all region types...that doesn't make sense!
round(RegionType, 3) 

#Charlie Geyer to the rescue!

#Something goofy going on with renewdata file
renewdata$Region
aout$formula
aout$coefficients
renewdata$fit #here is the problem

#should be
as.numeric(as.character(renewdata$varb) == "sm.2")

# when that is fixed, it works
renewdata$fit <- as.numeric(as.character(renewdata$varb) == "sm.2")
pout.amat<- predict(aout, newdata = renewdata, varvar= varb,
                    idvar = id, root = root, se.fit = TRUE, amat = amat, info.tol = 1e-10)
RegionType<- cbind(pout.amat$fit, pout.amat$se.fit)
rownames(RegionType)<- as.character(fred$Region)
colnames(RegionType)<- c("Expected Fitness", "SE")
round(RegionType, 3) 

#Now that this is fixed, can replace "region" with any other factor to estimate W
