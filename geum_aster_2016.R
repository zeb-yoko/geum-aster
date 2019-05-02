
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
vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016", "No.Fruit.2016", "sm")


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

###############################################################
#look into distributions for nodes

#recall nodes of graphical model
vars

#Germ, Survival, Flower.Y.N, all bernoulli

flwno<- dat2$No.Flowers.2016
flw.no<- subset(flwno, flwno >0)
frtno<- dat2$No.Fruit.2016
frt.no<-subset(frtno, frtno >0)
sm<- dat2$sm
sm.0<- subset(sm, sm > 0)

library(MASS)

fl.1<- fitdistr(flw.no, "normal")
fl.2<- fitdistr(flw.no, "negative binomial")#size: 1.12571436
fl.3<- fitdistr(flw.no, "poisson")

AIC(fl.1, fl.2, fl.3)

frt.1<- fitdistr(frt.no, "normal")
frt.2<- fitdistr(frt.no, "negative binomial")#size: 3.5074887
frt.3<- fitdistr(frt.no, "poisson")

AIC(frt.1, frt.2, frt.3)

sm.1<- fitdistr(sm.0, "normal")
sm.2<- fitdistr(sm.0, "negative binomial")#size: 1.0341103
sm.3<- fitdistr(sm.0, "poisson")

AIC(sm.1, sm.2, sm.3)

#####################################
#need to check dist for region-specific 

#load aster package
library(aster)

#set family list

famlist<- list(fam.bernoulli(), fam.negative.binomial(1.12571436), fam.negative.binomial(3.5074887),
               fam.negative.binomial(1.0341103))

#set graphical mode and dist. for fitness nodes
pred<- c(0,1,2,3,4,5)
fam<- c(1,1,1,2,3,4) #might want to play with these distributions, especially seedmass

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model for 2016 with only fitness variable: this is the "basic" model
# that we compare with later models, to determine significance of facotrs (e.g. Habitat Type, etc.)
aout2016a<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016a, show.graph=TRUE)

#Add block ID
aout2016b<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016b, show.graph = TRUE)

#likelihood test for Block.ID
anova(aout2016a, aout2016b)# significant effect of block

#add habitat type
aout2016c<- aster(resp~varb + fit:(HabitatType), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016c)

anova(aout2016a, aout2016c)#habitat type significant

#test effects of Region
aout2016d<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016d)

#likelihood test
anova(aout2016a, aout2016d)# significant effect of Region

#add Region
aout2016e<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016e, info.tol = 1e-15)#direction of recession/constancy!

anova(aout2016a, aout2016e)#significant effect of population




# The above aster analyses and likelihood-ratio test indicate significant effects of:
#   Block, Habitat type, Region, and population


#   Now want to generate mean expected fintess (and standard errors) for these factors


###############################################################################
#     Generate mean fitness (and stand errors) estimates for Habitat Type 
###############################################################################


# Start with Region Type 
aout<- aster(resp~varb + fit:(Region), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

aout1<- aster(resp~varb + fit:(Region + Block.ID), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

aout2<- aster(resp~varb + fit:(Region + Block.ID + Region*Block.ID), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout, show.graph = T)
summary(aout1, show.graph= T)
summary(aout2, show.graph = T, info.tol = 1e-15)

anova(aout, aout1, aout2)# block + region sig, but not interaction. 

#So, split the data into Region Specific data sets

redata2016.gla<- subset(redata2016, Region=="GL_alvar")
redata2016.gla<- droplevels(redata2016.gla)

dat.gla<- subset(dat2, Region=="GL_alvar")
gla.flw<- subset(dat.gla, No.Flowers.2016 >0 )
gla.flw<- gla.flw$No.Flowers.2016

gla.flw.0<- fitdistr(gla.flw, "negative binomial")# size: 1.10997286

gla.frt<- subset(dat.gla, No.Fruit.2016 > 0)
gla.frt<- gla.frt$No.Fruit.2016

gla.frt.0<- fitdistr(gla.frt, "negative binomial")# size: 3.3698966

gla.sm<- subset(dat.gla, sm > 0)
gla.sm<- gla.sm$sm

gla.sm.0<- fitdistr(gla.sm, "negative binomial")# size: 1.0217064

redata2016.mba<- subset(redata2016, Region=="MB_alvar")
redata2016.mba<- droplevels(redata2016.mba)

dat.mba<- subset(dat2, Region=="MB_alvar")
mba.flw<- subset(dat.mba, No.Flowers.2016 > 0)
mba.flw<- mba.flw$No.Flowers.2016

mba.flw.0<- fitdistr(mba.flw, "negative binomial")# size: 1.7177893

mba.frt<- subset(dat.mba, No.Fruit.2016 >0)
mba.frt<- mba.frt$No.Fruit.2016

mba.frt.0<- fitdistr(mba.frt, "negative binomial")# size: 100.00010

mba.sm<- subset(dat.mba, sm > 0)
mba.sm<- mba.sm$sm

mba.sm.0<- fitdistr(mba.sm, "negative binomial")# size: 10.57473


redata2016.pra<- subset(redata2016, Region=="Prairie")
redata2016.pra<- droplevels(redata2016.pra)

dat.pra<- subset(dat2, Region=="Prairie")
pra.flw<- subset(dat.pra, No.Flowers.2016 > 0)
pra.flw<- pra.flw$No.Flowers.2016

pra.flw.0<- fitdistr(pra.flw, "negative binomial")# size: 2.2972941

pra.frt<- subset(dat.pra, No.Fruit.2016 > 0)
pra.frt<- pra.frt$No.Fruit.2016

pra.frt.0<- fitdistr(pra.frt, "negative binomial")# size: 100.0001238

pra.sm<- subset(dat.pra, sm > 0)
pra.sm<- pra.sm$sm

pra.sm.0<- fitdistr(pra.sm, "negative binomial")# size: 1.481091

#aster analyses with Block.ID for all three region types

famlist.gla<- list(fam.bernoulli(), fam.negative.binomial(1.10997286), fam.negative.binomial(3.3698966), fam.negative.binomial(1.0217064))

aout.gla<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.gla, famlist = famlist.gla)


famlist.mba<- list(fam.bernoulli(), fam.negative.binomial(1.7177893), fam.negative.binomial(100.00010), fam.negative.binomial(10.57473))

aout.mba<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.mba, famlist = famlist.mba)


famlist.pra<- list(fam.bernoulli(), fam.negative.binomial(2.2972941), fam.negative.binomial(100.0001238), fam.negative.binomial(1.481091))

aout.pra<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.pra, famlist = famlist.pra)

summary(aout.gla, show.graph = T, info.tol = 1e-14)

summary(aout.mba, show.graph = T, info.tol = 1e-15)

summary(aout.pra, show.graph = T, info.tol = 1e-14)

  
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
dim(amat)# makes an 12 x 6 x 12 matrix (12 blocks types and 7 nodes of graphicla model)

#only want prediction for k'th individual that contribute to expected
#fitness, and want to add only seedmass (sm) entries

foo<- grepl("sm", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to "sm"

#generate predicted valuses using region-specific aout object, with renewdata, and amat format
pout.amat.gla<- predict(aout.gla, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-14)

pout.amat.mba<- predict(aout.mba, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-15)

pout.amat.pra<- predict(aout.pra, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-15)

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

Region.gla<-round(Region.gla, 3) 
Region.mba<- round(Region.mba, 3) 
Region.pra<-round(Region.pra, 3) 

Region.gla
Region.mba
Region.pra

summary(Region.gla)# median = 2.345 corresponds to block 6: 3.027 (1.245)
summary(Region.mba)# median = 2.50E-10 which is basically zero
summary(Region.pra)# median = 8.055348e-12 (basically zero, as above)



###############################################################################
#     Generate mean fitness (and stand errors) estimates for 2016 Habitatat Types 
###############################################################################


# Start with Habitat Type 
aout<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aouta)

aouta<- aster(resp~varb + fit:(HabitatType), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout, show.graph = TRUE)


aoutb<- aster(resp~varb + fit:(HabitatType + Block.ID), pred, fam, varb, id, root, data=redata2016, famlist = famlist)
summary(aoutb)

aoutc<- aster(resp~varb + fit:(HabitatType + Block.ID + HabitatType*Block.ID), pred, fam, varb, id, root, data=redata2016, famlist = famlist)
summary(aoutc, show.graph=T, info.tol = 1e-16)

anova(aout, aouta, aoutb, aoutc)

#Significant interaction between Block and Habitat Type, so split data into Alvar 
# and prairie datasets and perform separate analyses

#separate data and estimate size parameter of negative binomaial distributions
# for flower number, fruit number, and seedmass

#Start with Alvars
redata2016.alv<- subset(redata2016, HabitatType=="Alvar")
redata2016.alv<- droplevels(redata2016.alv)

dat.alv<- subset(dat2, HabitatType=="Alvar")
alv.flw<- subset(dat.alv, No.Flowers.2016 > 0)
alv.flw<- alv.flw$No.Flowers.2016

alv.flw.0<-fitdistr(alv.flw, "negative binomial")# size: 1.10997286

alv.frt<- subset(dat.alv, No.Fruit.2016 > 0)
alv.frt<- alv.frt$No.Fruit.2016

alv.frt.0<-fitdistr(alv.frt, "negative binomial")# size: 3.3698966

alv.sm<- subset(dat.alv, sm >0)
alv.sm<- alv.sm$sm

alv.sm.0<- fitdistr(alv.sm, "negative binomial")# size: 1.0217064

famlist.alv<- list(fam.bernoulli(), fam.negative.binomial(1.10997286), 
                   fam.negative.binomial(3.3698966), fam.negative.binomial(1.0217064))

aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016.alv, famlist = famlist.alv)
summary(aouta)

aout.alv<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.alv, famlist = famlist.alv)
summary(aout, info.tol = 1e-14)

anova(aouta, aout.alv)

#Next do Prairie
redata2016.pra<- subset(redata2016, HabitatType=="Prairie")
redata2016.pra<- droplevels(redata2016.pra)

dat.pra<- subset(dat2, HabitatType=="Prairie")
pra.flw<- subset(dat.pra, No.Flowers.2016 > 0)
pra.flw<- pra.flw$No.Flowers.2016

pra.flw.0<-fitdistr(pra.flw, "negative binomial")# size: 1.8744208

pra.frt<- subset(dat.pra, No.Fruit.2016 > 0)
pra.frt<- pra.frt$No.Fruit.2016

pra.frt.0<-fitdistr(pra.frt, "negative binomial")# size: 100.0002113

pra.sm<- subset(dat.pra, sm >0)
pra.sm<- pra.sm$sm

pra.sm.0<- fitdistr(pra.sm, "negative binomial")# size: 1.593404

famlist.pra<- list(fam.bernoulli(), fam.negative.binomial(1.8744208), 
                   fam.negative.binomial(100.0002113), fam.negative.binomial(1.593404))

aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016.pra, famlist = famlist.pra)
summary(aouta)

aout.pra<- aster(resp~varb + fit:(Block.ID), pred, fam, varb, id, root, data=redata2016.pra, famlist = famlist.pra)
summary(aout, info.tol = 1e-14)

anova(aouta, aout.pra)

# make data.frame of indivudals for each habitat type (Alvar and Prairie)

fred <- data.frame(Block.ID=levels(redata2016$Block.ID),
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
nBlock<- nrow(fred)#all data has same number of blocks so any file will do
nnode<- length(vars)
amat<- array(0, c(nBlock, nnode, nBlock))
dim(amat)# makes an 12 x 6 x 12 matrix (12 blocks and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nBlock)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat.alv<- predict(aout.alv, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-14)

pout.amat.pra<- predict(aout.pra, newdata= renewdata, varvar= varb,
                        idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-15)
#combine estimates with standard error, and then round
#to three decimal places
alv.block<- cbind(pout.amat.alv$fit, pout.amat.alv$se.fit)
pra.block<- cbind(pout.amat.pra$fit, pout.amat.pra$se.fit)

rownames(alv.block)<- as.character(fred$Block.ID)
rownames(pra.block)<- as.character(fred$Block.ID)


colnames(alv.block)<- c("Expected Fitness", "SE")
colnames(pra.block)<- c("Expected Fitness", "SE")


alv.block<- round(alv.block, 3)
pra.block<- round(pra.block, 3) 

alv.block
pra.block

summary(alv.block)#median = 2.345, corresponds with block 8 1.664 (0.699)
summary(pra.block)#median = 0.000, so keep at zero

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