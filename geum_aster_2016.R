
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

#test habitat type + block.ID
aout2016c2<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016)

summary(aout2016c2)

#likelihood test
anova(aout2016b, aout2016c2)# significant effect of Habitat Type + Block.ID

#add Region
aout2016d<- aster(resp~varb + fit:(Block.ID + HabitatType), pred, fam, varb, id, root, data=redata2016)

summary(aout2016d)

anova(aout2016c, aout2016d)#significant effect of Region + Block.ID

aout2016e<- aster(resp~varb + fit:(Block.ID + Region), pred, fam, varb, id, root, data=redata2016)

#summary(aout2016e, show.graph = T, info.tol = 1e-16)# NOTE: direction of recession error with Pop

anova(aout2016c2, aout2016e)#significant effect of Block + Population



# The above aster analyses and likelihood-ratio test indicate significant effects of:
#   Block, Habitat type, Region, and population


#   Now want to generate mean expected fintess (and standard errors) for these factors


###############################################################################
#     Generate mean fitness (and stand errors) estimates for Habitat Type 
###############################################################################


# Start with Habitat Type 
aout<- aster(resp~varb + fit:(Block.ID + Region + Population), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph = TRUE)

# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE)


#exploring directions of recession
fred <- eigen(aout$fisher, symmetric = TRUE)
dor <- fred$vectors[ , fred$values == min(fred$values)]
names(dor) <- names(aout$coefficients)
dor <- zapsmall(dor / max(dor))
dor

modmat <- aout$modmat
dim(modmat)

modmat <- as.vector(modmat)
modmat <- matrix(modmat, ncol = length(dor))
dor.phi <- modmat %*% dor
dor.phi <- as.vector(dor.phi)


unique(dor.phi)

sum(dor.phi)

foo <- data.frame(Population = as.character(redata2016$Population),
                   id = redata2016$id,
                  varb = as.character(redata2016$varb),
                  resp = redata2016$resp, stringsAsFactors = FALSE)
 foo <- foo[dor.phi == 1, ]
 
 foo

 unique(foo$Population) #population WA-BLK
 
 unique(foo$varb) #seed mass
 
 unique(foo$id) #these individuals
 
 unique(foo$resp)# just zeros
 
 #each individual in Population WA-BLK failed to set any seeds (all zeros)
 
 # Various options exist: merge this pop with anohter, delete this pop from data
 # and just "say" it has zero fitness, change the graphical model: in this case, make 
 # fruit number terminal fitness node?
 
    #looking at data, all failed to set fruit, too. But, some (2 inds) produced
    #a single flower...maybe go with this?
 
 
 #take a closer look at data
 tapply(dat2$sm, dat2$Population, max) # hmmm 9 families set a max of zero seeds
 
 tapply(dat2$No.Fruit.2016, dat2$Population, max)# same nine families failed to set >0 fruit, too
 
 tapply(dat2$No.Flowers.2016, dat2$Population, max)# only three families did not produce >0 number of flowers 
 

 #####################################################################
 #  So, try No.Flowers.2016 as fitness node (i.e. truncate model)
 #####################################################################
 
 #set response variables -> these represent variables in graphical model
 vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016")
 
 
 #reshape data so that all response variables are located in a single vector in a new data
 #set called "redata"
 redata2016 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")
 

 #Designation of fitness variable for 2016 data
 fit <- grepl("No.Flowers.2016", as.character(redata2016$varb))
 fit<- as.numeric(fit)
 
 redata2016b$fit <- fit
 
 #check
 with(redata2016, sort(unique(as.character(varb)[fit == 0])))
 with(redata2016, sort(unique(as.character(varb)[fit == 1])))
 
 
 #add a variable "root" to redata files, where value is 1
 redata2016b<- data.frame(redata2016, root=1)
 

 
 #load aster package
 library(aster)
 
 #set graphical mode and dist. for fitness nodes
 pred<- c(0,1,2,3)
 fam<- c(1,1,1,2) #might want to play with these distributions, especially seedmass
 
 #describe dist. of preds.
 sapply(fam.default(), as.character)[fam]
 
 
 aout<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016)
 
 summary(aout, show.graph = TRUE) #same problem, this time pop "AB-RO" all failed to flower!
 
 
 
 
 # so hard way, is to redisign the data, so that aster consideres the whole data set
 # as a single individual
 
 #But we have to re-specify the graphical model, which is difficult
 
 #kill the graph model
 outies <- dor.phi == 1
 subdata <- redata2016[! outies, ]
 
 #save real ID values
 id <- subdata$id
 
 #call it all one indiv
 subdata$id <- 1
 
 #remake the graph
 idx <- seq(1, nrow(subdata))
  varb <- as.character(subdata$varb)
  pred <- rep(NA, length(idx))
  fam <- rep(NA, length(idx))
  pred[varb == "sm"] <- 0
  fam[varb == "sm"] <- 1
  head(idx[varb == "sm"])
 

  
  
  sum(varb == "Germination.Y.N") == sum(varb == "Survival.Y.N")#true
  
  pred[varb == "Germination.Y.N"] <- idx[varb == "Survival.Y.N"]
  fam[varb == "Germination.Y.N."] <- 1
  
  sum(varb == "Survival.Y.N") == sum(varb == "Flower.Y.N.2016")#true
  
  pred[varb == "Survival.Y.N"] <- idx[varb == "Flower.Y.N.2016"]#true
  fam[varb == "Survival.Y.N"] <- 1
  
  sum(varb == "Flower.Y.N.2016") == sum(varb == "No.Flowers.2016")#true
  
  pred[varb == "Flower.Y.N.2016"] <- idx[varb == "No.Flowers.2016"]
  fam[varb == "Flower.Y.N.2016"] <- 1
  
  
  sum(varb == "Flower.Y.N.2016") == sum(varb == "No.Flowers.2016")#true
  
  pred[varb == "Flower.Y.N.2016"] <- idx[varb == "No.Flowers.2016"]
  fam[varb == "Flower.Y.N.2016"] <- 1
  
  sum(varb == "No.Flowers.2016") == sum(varb == "No.Fruit.2016")# true
  
  pred[varb == "No.Flowers.2016"] <- idx[varb == "No.Fruit.2016"]
  fam[varb == "No.Flowers.2016"] <- 2
  
  sum(varb == "No.Fruit.2016") == sum(varb == "sm")# false
  
  bar <- match(id[varb == "sm"], id[varb == "No.Fruit.2016"])
  pred[varb == "sm"] <- idx[varb == "No.Fruit.2016"][bar]
  fam[varb == "sm"] <- 2
  
#Now remake "varb"
  subvarb <- paste(as.character(subdata$varb), id, sep = "")
  subdata <- data.frame(subdata, subvarb = subvarb)
  
  aout.sub <- aster(resp ~ varb + fit : Population,
                    +pred, fam, subvarb, id, root, data = subdata)
  summary(aout.sub)
  
  
# make up  data for hypothetical individual that meet "typical" criteria:
# Therefore, "make up" covariate data for hypothetical individuals that are comparable and obtain mean values for them

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(HabitatType=levels(redata2016$HabitatType),Family.NonUnique=redata2016$Family.NonUnique,
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
nHabitatType<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nHabitatType, nnode, nHabitatType))
dim(amat)# makes an 2 x 6 x 2 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nHabitatType)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat)

#combine estimates with standard error, and then round
#to three decimal places
Habtype<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(Habtype)<- as.character(fred$HabitatType)


colnames(Habtype)<- c("Expected Fitness", "SE")

round(Habtype, 3) 

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