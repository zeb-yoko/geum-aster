#Load data and run code that cleans 2016 & 2017 data.

setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/aster-analysis/")

#load data
dat2<- read.csv("NV_CG_Experiment2wdist2.csv")


#subset data for 2017 analysis
#dat2<- dat[c("Family.Unique",   "Block.ID", "HabitatType", "Region", "Population",
 #            "Germination.Y.N","Survival.Y.N","Survival.Y.N.2017", "Flower.Y.N.2016",
  #           "Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017", "Fruit.Y.N.2016",
   #          "Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017", "sm", "sm.2")]




########################################################
########################################################
## Begin cleaning data-> start with 2016 cleaning code##
########################################################
########################################################

#replace NAs with 0 (zero)
dat2[is.na(dat2)] <- 0

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

######################################################################################
# Now begin cleaning the 2018 data


#survival to 2018 but not 2017
subset(dat2, Survival.Y.N.2018==1 & Survival.Y.N.2017==0)# no errors

#survival to 2018 but not from greenhouse
subset(dat2, Survival.Y.N.2018==1 & Survival.Y.N==0)# 7 errors

#fix above errors
dat2$Survival.Y.N[dat2$Survival.Y.N.2018==1 & dat2$Survival.Y.N==0]=1

#flower in 2018 but no survival
subset(dat2, Flowering.Y.N.2018==1 & Survival.Y.N.2018==0)# 7 errors

#fix above errors
dat2$Survival.Y.N.2018[dat2$Flowering.Y.N.2018==1 & dat2$Survival.Y.N.2018==0]=1

#flowering =0, but total flowers >0
subset(dat2, Flowering.Y.N.2018==0 & Total.Flowers.2018 > 0)# 6 errors

#fix above errors
dat2$Flowering.Y.N.2018[dat2$Flowering.Y.N.2018==0 & dat2$Total.Flowers.2018 > 0]=1

#no flowers produced but fruit set=1
subset(dat2, Fruit.Y.N.2018==1 & Total.Flowers.2018==0)#6 errors

#correct above errors (make fruit y.n. =0)
dat2$Fruit.Y.N.2018[dat2$Fruit.Y.N.2018==1 & dat2$Total.Flowers.2018==0]=0

#Number of fruits > 0, but fruit y.n.=0
subset(dat2, Fruit.Y.N.2018==0 & No.Fruit.2018 >0)#22 errors

#fix above errors
dat2$Fruit.Y.N.2018[dat2$Fruit.Y.N.2018==0 & dat2$No.Fruit.2018 >0]=1

#seed mass >0, but fruit number =0
subset(dat2, seedmass.2018.g. >0 & No.Fruit.2018==0)# 2 errors

#fix above errors, make fruit number =1
dat2$No.Fruit.2018[dat2$seedmass.2018.g. >0 & dat2$No.Fruit.2018==0]=1

#Recheck for flower2018=1 but survival=0 for introduced errors
subset(dat2, Flowering.Y.N.2018==1 & Survival.Y.N.2018==0)# 3 new errors

#fix above errors
dat2$Survival.Y.N.2018[dat2$Flowering.Y.N.2018==1 & dat2$Survival.Y.N.2018==0]=1

#That should conclude the data cleaning.

#subset what's needed for 2018 analysis

dat3<- dat2[c("Family.Unique","Block.ID", "HabitatType", "Region", "Population",
              "Germination.Y.N","Survival.Y.N", "Survival.Y.N.2018", "Flowering.Y.N.2018",
              "Total.Flowers.2018", "Fruit.Y.N.2018", "No.Fruit.2018", "seedmass.2018.g."  )]


#write.table(dat3, file="dat3_2018.csv", sep=",", row.names = F, quote = F)

#conver seedmass.g to seedmass in mg

dat3$sm3<- round((dat3$seedmass.2018.g. * 1000), 0)


#set response variables -> these represent variables in graphical model
vars<- c( "Germination.Y.N","Survival.Y.N", "Survival.Y.N.2018", "Flowering.Y.N.2018",
          "Total.Flowers.2018", "Fruit.Y.N.2018", "No.Fruit.2018", "sm3")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata <- reshape(dat3, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#Designation of fitness variable for 2016 data
fit <- grepl("sm3", as.character(redata$varb))
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


#load aster package
library(aster)

pred<- c(0,1,2,3,4,5,6,7)
fam<- c(1,1,1,1,2,1,2,2) #might want to play with these distributions, especially seedmass

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

layer <- gsub("[0-9]", "", as.character(redata$varb))
unique(layer)

redata <- data.frame(redata, layer = layer)
with(redata, class(layer))

aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata, method = 'nlm')

summary(aouta, show.graph=TRUE, info.tol = 1e-10)

aout<- aster(resp~varb + fit:(Block.ID + Region), pred, fam, varb, id, root, data=redata,method='nlm', maxiter=5000)

summary(aout, show.graph=TRUE, info.tol = 1e-11)

anova(aouta, aout)

aoutb<- aster(resp~ varb + fit:(Block.ID + Region + Block.ID*Region), pred, fam, varb, id, root, data = redata, method='nlm')

summary(aoutb, show.graph=T, info.tol = 1e-10)

anova(aout, aoutb)

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



