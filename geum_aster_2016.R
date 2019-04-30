
# Aster code to estimate mean expected fitness of Geum from unique habitat types, Regions,
# and populaitons over three years


setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/aster-analysis/")

#load data
dat<- read.csv("NV_CG_Experiment2.csv")


#subset data for 2016 analysis
dat2<- dat[c("Family.NonUnique", "Family.Unique",   "Block.ID", "HabitatType", "Region", "Population", "Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016", "Fruit.Y.N.2016", "No.Fruit.2016", "sm")]

#replace NA with zero
dat2[is.na(dat2)] <- 0


#The below conversion of seedmass from g -> mg has been done in the excel file,
# to avoid a formatting issue that fouling later aster analyses. New variable "sm"
# is seed mass for 2016 plants in milligrams

#Aster can't handle decimal data
#Conver seedmass from g to mg by " x 1000"
#dat2$Sm.2016mg<- (dat2$Seedmass.2016 * 1000)


#Note: found error when plant produced 1 fruit, but data has "0" for flower number
subset(dat2, No.Flowers.2016==0 & Fruit.Y.N.2016==1)#individual 890

# Replace with flower number =1 for this case
dat2$No.Flowers.2016[dat2$Flower.Y.N.2016==1&dat2$Fruit.Y.N.2016==1]= 1

#check
subset(dat2, No.Flowers.2016==0 & Fruit.Y.N.2016==1) #error corrected

#Additional Errors

#Flower.Y.N.2016 = 0 but No.Flowers.2016 > 0?
subset(dat2, Flower.Y.N.2016==0 & No.Flowers.2016 >0)# 8 erros

#correct
dat2$Flower.Y.N.2016[dat2$Survival.Y.N==1 & dat2$No.Flowers.2016 >0]=1

#check correction
subset(dat2, Flower.Y.N.2016==0 & No.Flowers.2016 >0)# erros corrected

#Third error: survival =0, but Flower.y.N.2016 ==1
subset(dat2, Survival.Y.N==0 & Flower.Y.N.2016==1 )#10 errors

#correct
dat2$Survival.Y.N[dat2$Germination.Y.N==1 & dat2$Flower.Y.N.2016==1]=1

#check correction
subset(dat2, Survival.Y.N==0 & Flower.Y.N.2016==1 )# errors corrected

#dat2$Seedmass.2016<- as.integer(dat2$Seedmass.2016)

#set response variables -> these represent variables in graphical model
vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016","Fruit.Y.N.2016", "No.Fruit.2016", "sm")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2016 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

#write.csv(redata2016, file="redata2016.csv", row.names = FALSE, quote = FALSE)

#Designation of fitness variable for 2016 data
fit <- grepl("sm", as.character(redata2016$varb))
fit<- as.numeric(fit)

redata2016$fit <- fit

#check
with(redata2016, sort(unique(as.character(varb)[fit == 0])))
with(redata2016, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2016<- data.frame(redata2016, root=1)


#check class of redata columns

sapply(redata2016, class)

#make block.id a factor

redata2016$Block.ID<- as.factor(redata2016$Block.ID)

#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes
pred<- c(0,1,2,3,4,5,6)
fam<- c(1,1,1,2,1,2,2) #might want to play with these distributions, especially seedmass

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model for 2016 with only fitness variable: this is the "basic" model
# that we compare with later models, to determine significance of facotrs (e.g. Habitat Type, etc.)
aout2016a<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016)

summary(aout2016a, show.graph=TRUE)

#Add block ID
aout2016b<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016)

summary(aout2016b, show.graph = TRUE)

#likelihood test for Block.ID
anova(aout2016a, aout2016b)# significant effect of block

#add habitat type
aout2016c<- aster(resp~varb + fit:(HabitatType), pred, fam, varb, id, root, data=redata2016)

summary(aout2016c)

anova(aout2016a, aout2016c)#habitat type significant

#test effects of Region
aout2016d<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016)

summary(aout2016d)

#likelihood test
anova(aout2016a, aout2016d)# significant effect of Region

#add Region
aout2016e<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016)

summary(aout2016e, info.tol = 1e-15)

anova(aout2016a, aout2016e)#significant effect of population




# The above aster analyses and likelihood-ratio test indicate significant effects of:
#   Block, Habitat type, Region, and population


#   Now want to generate mean expected fintess (and standard errors) for these factors


###############################################################################
#     Generate mean fitness (and stand errors) estimates for Habitat Type 
###############################################################################


# Start with Region Type 
aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016)

aout1<- aster(resp~varb + fit:(Region + Block.ID), pred, fam, varb, id, root, data=redata2016)

aout2<- aster(resp~varb + fit:(Region + Block.ID + Region*Block.ID), pred, fam, varb, id, root, data=redata2016)

anova(aout, aout1, aout2)# block + region sig, but not interaction. 

#So split the data into Region Specific data sets

redata2016.gla<- subset(redata2016, Region=="GL_alvar")
redata2016.gla<- droplevels(redata2016.gla)

redata2016.mba<- subset(redata2016, Region=="MB_alvar")
redata2016.mba<- droplevels(redata2016.mba)

redata2016.pra<- subset(redata2016, Region=="Prairie")
redata2016.pra<- droplevels(redata2016.pra)
  

#aster analyses with Block.ID for all three region types

aout.gla<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.gla)

aout.mba<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.mba)

aout.pra<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.pra)

summary(aout.gla, show.graph = T, info.tol = 1e-16)

summary(aout.mba, show.graph = T, info.tol = 1e-16)

summary(aout.pra, show.graph = T, info.tol = 1e-16)

  
#Setting up design matrix that will eventually hold estimates from 'predict'

#Note, we can use the same desing matrix for all three models (3 regions)

fred <- data.frame(Block.ID=levels(redata2016.gla$Block.ID),
                   Germination.Y.N=1, Survival.Y.N=1, Flower.Y.N.2016=1, No.Flowers.2016=1,
                 Fruit.Y.N.2016=1, No.Fruit.2016=1, sm=1, root = 1)

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
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 7 x 12 matrix (12 blocks types and 7 nodes of graphicla model)

#only want prediction for k'th individual that contribute to expected
#fitness, and want to add only seedmass (sm) entries

foo<- grepl("sm", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to "sm"

#generate predicted valuses using region-specific aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout.gla, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-16)

pout.amat.mba<- predict(aout.mba, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-16)

pout.amat.pra<- predict(aout.pra, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-16)

#combine estimates with standard error, and then round
#to three decimal places
Region.gla<- cbind(pout.amat.gla$fit, pout.amat.gla$se.fit)

Region.mba<- cbind(pout.amat.mba$fit, pout.amat.mba$se.fit)

Region.pra<- cbind(pout.amat.pra$fit, pout.amat.pra$se.fit)


rownames(Region.gla)<- as.character(fred$Block.ID)
rownames(Region.mba)<- as.character(fred$Block.ID)
rownames(Region.pra)<- as.character(fred$Block.ID)

colnames(Region.gla)<- c("Expected Fitness", "SE")
colnames(Region.mba)<- c("Expected Fitness", "SE")
colnames(Region.pra)<- c("Expected Fitness", "SE")

round(Region.gla, 3) 
round(Region.mba, 3) 
round(Region.pra, 3) 

summary(Region.gla)# median = 2.345 corresponds to block 6: 3.027 (1.245)
summary(Region.mba)# median = 2.50E-10 which is basically zero, so go with smalles non-zero? Block 9 1.6 (2.670)
summary(Region.pra)# median = 8.055348e-12 (basically zero, as above), so use block 0.061 (0.155)

###############################################################################
#     Generate mean fitness (and stand errors) estimates for Regions 
###############################################################################


# Start with Habitat Type 
aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph = TRUE)

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

###############################################################################
#     Generate mean fitness (and stand errors) estimates for Populations 
###############################################################################


# Start with Habitat Type 
aout<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph = TRUE, info.tol = 1e-16)

# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol = 1e-16)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Population=levels(redata2016$Population),
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
nPop<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nPop, nnode, nPop))
dim(amat)# makes an 2 x 6 x 2 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nPop)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-16)

#combine estimates with standard error, and then round
#to three decimal places
pop<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(pop)<- as.character(fred$Population)


colnames(pop)<- c("Expected Fitness", "SE")

round(pop, 3) 
write.csv(pop, file="2016pop.csv", quote = F)

###############################################################################
#     Generate mean fitness (and stand errors) estimates for Family.NonUnique
###############################################################################

redata2016$Family.NonUnique<- as.factor(redata2016$Family.NonUnique)

# Start with Habitat Type 
aout<- aster(resp~varb + fit:(Family.NonUnique), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph = TRUE, info.tol = 1e-20)#direction of recession issue

# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE, info.tol=1e-2000)

# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Family.NonUnique=levels(redata2016$Family.NonUnique),
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
nfam<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nfam, nnode, nfam))
dim(amat)# makes an 41 x 6 x 41 matrix (41 family IDs types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nfam)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-50)

#The above produces an error: "direction of recession -> cannot compute standard error"
# I will have to look into this.

#For more info on this error, see Charlie's great explaination: http://www.stat.umn.edu/geyer/8931aster/slides/s9.pdf

# and also this posting on the Aster Forum: https://groups.google.com/forum/#!msg/aster-analysis-user-group/oivUr6Q99rE/cxLljfoCfA0J;context-place=forum/aster-analysis-user-group