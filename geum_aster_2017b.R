
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
             "Germination.Y.N","Survival.Y.N","Survival.Y.N.2017", "Flower.Y.N.2016",
             "Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017", "Fruit.Y.N.2016",
             "Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017", "sm", "sm.2")]

#make function to look for NAs across all columns with sapply()

fun<- function(x) table(is.na(x))

sapply(dat2, fun) # good, not factors have NA values (note: not using Family.NonUnique)
sapply(dat2, class) #Block is class integer; make it a factor

#make block.id a factor
dat2$Block.ID<- as.factor(dat2$Block.ID)

#replace NAs with 0 (zero)
dat2[is.na(dat2)] <- 0

#make sum of seedmass in 2017 and 2017 a new variable
#dat2$sm.2017<-dat2$sm + dat2$sm.2



##############################################
# Set up simplified graph. models that can be recalled
# quickly, for the purpose of code cleaning


#this is the 2016 graph. model:
gm2016<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016", "Fruit.Y.N.2016",
 "No.Fruit.2016","sm")

#this is the simplified graph. model for 2017, for cleaning purposes. The later graph. model
# with contain a "branch".
gm2017<- c("Germination.Y.N","Survival.Y.N","Survival.Y.N.2017","Flower.Y.N.2017","Total.Flowers.2017",
           "Fruit.Y.N.2017", "No.Fruit.2017", "sm.2")


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

#Below is the check for 2016

#set response variables -> these represent variables in graphical model
vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016","No.Fruit.2016", "sm")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2016 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

write.csv(redata2016, file="redata2016.csv", row.names = FALSE, quote = FALSE)

#Designation of fitness variable for 2016 data
fit <- grepl("sm", as.character(redata2016$varb))
fit<- as.numeric(fit)

redata2016$fit <- fit

#check
with(redata2016, sort(unique(as.character(varb)[fit == 0])))
with(redata2016, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2016<- data.frame(redata2016, root=1)

#make sure Block.ID, Region, Habitat Type, and Population a factor
#redata2016$Block.ID <- as.factor(redata2016$Block.ID)
#redata2016$Region <- as.factor(redata2016$Region)
#redata2016$Population <- as.factor(redata2016$Population)
#redata2016$HabitatType<- as.factor(redata2016$HabitatType)

#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes
pred<- c(0,1,2,3,4,5)
fam<- c(1,1,1,2,2,2) #might want to play with these distributions, especially seedmass

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model for 2015 with only fitness variable: this is the "basic" model
# that we compare with later models, to determine significance of facotrs (e.g. Habitat Type, etc.)
aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph=TRUE)

# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Region=levels(redata2016$Region),
                   Germination.Y.N=1, Survival.Y.N=1, Flower.Y.N.2016=1, No.Flowers.2016=1,
                   No.Fruit.2016=1, sm=1, root = 1)

# reshape the "made up data" just as the actual data
renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")

# make character string from "varb" of renewdata, without actual values (e.g., the layers of varb in renewdata)
layer<- gsub("[0-9]", "", as.character(renewdata$varb))

# add layer to renewdata
renewdata<- data.frame(renewdata, layer= layer)

# add Seedmass.2016 in new layer col of renewdata as numeric, called fit
fit<- as.numeric(layer=="sm")

# add fit to renewdata
renewdata<- data.frame(renewdata, fit = fit)


#Generate fintess estimates and standard errors for each block
nRegion<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nRegion, nnode, nRegion))
dim(amat)# makes an 3 x 6 x 3 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nRegion)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat)

#combine estimates with standard error, and then round
#to three decimal places
RegionType<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(RegionType)<- as.character(fred$Region)


colnames(RegionType)<- c("Expected Fitness", "SE")

round(RegionType, 3) 

##############################################################################
##############################################################################
## The above code cleaned and tested the 2016 data, now move on to 2017 data##
##############################################################################
##############################################################################

#recall 2017 graph. model: includes a germ and survival after transplant (2016),
# so start at survival to 2017 after survial after transplant 
gm2017

#suvive to 2017 but surv.2016=0
subset(dat2, Survival.Y.N==0 & Survival.Y.N.2017==1)# no errors

##################### NEW
dat2$Survival.Y.N[dat2$Survival.Y.N==0 & dat2$Survival.Y.N.2017==1]=1
####################


#flowered but no survival
subset(dat2, Survival.Y.N.2017==0 & Flower.Y.N.2017==1)# 17 errors

#fix errors: set survival.y.n to 1
dat2$Survival.Y.N.2017[dat2$Survival.Y.N.2017==0 & dat2$Flower.Y.N.2017==1]=1

#flw.no >0, but flw.y.n=0
subset(dat2, Total.Flowers.2017 >0 & Flower.Y.N.2017==0)# 10 errors

#set flw.y.n=1, but this introduces 3 new erorrs to survival, so recorrect
dat2$Flower.Y.N.2017[dat2$Total.Flowers.2017 >0 & dat2$Flower.Y.N.2017==0]=1

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


#################################################################################
##  This should be enough cleaning to check 2017 data with simple aster model  ##
#################################################################################


#still not giving reasonable results. Take closer look at factor variables
    #NOTE: predictors are all being dropped! why?!?!?


#isolate factors
fs<- dat2[c("Family.Unique", "Block.ID", "HabitatType", "Region", "Population")]

sapply(fs, class) #all factors

table(is.na(fs)) # no NAs

sapply(fs, table)# all seems fine with factors.


vars<- c("Germination.Y.N", "Survival.Y.N","Survival.Y.N.2017", "Flower.Y.N.2017"   
         ,"Total.Flowers.2017" ,"Fruit.Y.N.2017","No.Fruit.2017","sm.2")

#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2017 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of fitness variable for 2016 data
fit <- grepl("sm.2", as.character(redata2017$varb))
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
