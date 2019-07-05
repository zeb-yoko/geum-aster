
# Begin working with Fitness Surface Work. Starting with 2016 Data


setwd("C:/Users/Mason Kulbaba/Dropbox/git/geum-aster")

#load data
#dat<- read.csv("NV_CG_Experiment2.csv")

dat<- read.csv("fitness_landscape_data_cleaned.csv")


#subset data for 2016 analysis
dat2<- dat[c( "Family.Unique",   "Block.ID", "HabitatType","Dist.from.cg.km", "No.Days.to.Germ",
             "Region", "Population", "Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016",
             "No.Flowers.2016", "Fruit.Y.N.2016", "No.Fruit.2016", "sm")]

#remove rows with NAs in "days to germ." variable

dat3<- subset(dat2, No.Days.to.Germ >0)

dat3<- subset(dat3, Germination.Y.N==1 | Germination.Y.N==0)

dat2<-dat3

#set response variables -> these represent variables in graphical model
vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016", "No.Fruit.2016", "sm")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2016 <- reshape(dat3, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")

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

frtno<- dat2$No.Fruit.2016

sm<- dat2$sm


library(MASS)

fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.161672878
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size:  0.045111784
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)

frt.2

sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.009215800 
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)

sm.2

#####################################
# Overall Model work

#load aster package
library(aster)

#set family list

famlist<- list(fam.bernoulli(), 
               fam.negative.binomial(0.161672878), 
               fam.negative.binomial(0.045111784),
               fam.negative.binomial(0.009215800))

#set graphical mode and dist. for fitness nodes
pred<- c(0,1,2,3,4,5)
fam<- c(1,1,1,2,3,4) 

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]

#fixed effect model for 2016 with only fitness variable: this is the "basic" model
# that we compare with later models, to determine significance of facotrs (e.g. Habitat Type, etc.)
aout2016<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout2016, show.graph=TRUE)

#add distance from seed source
aout1<- aster(resp~varb+0 + fit:Region, pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout1, show.graph=T)

anova(aout2016, aout1)#Region not significant here, but we knew that and it's ok.


####################################################################################

#Divide into region-specific data set and fit node distributions

#Start with GL_alvar

dat.gla<- subset(dat2, Region=="GL_alvar")
dat.gla<- droplevels(dat.gla)


flwno<- dat.gla$No.Flowers.2016

frtno<- dat.gla$No.Fruit.2016

sm<- dat.gla$sm


library(MASS)

fl.1<- fitdistr(flwno, "normal")
fl.2<- fitdistr(flwno, "negative binomial")#size: 0.23784984
fl.3<- fitdistr(flwno, "poisson")

AIC(fl.1, fl.2, fl.3)
fl.2


frt.1<- fitdistr(frtno, "normal")
frt.2<- fitdistr(frtno, "negative binomial")#size: 0.06728733
frt.3<- fitdistr(frtno, "poisson")

AIC(frt.1, frt.2, frt.3)

frt.2

sm.1<- fitdistr(sm, "normal")
sm.2<- fitdistr(sm, "negative binomial")#size: 0.013494317
sm.3<- fitdistr(sm, "poisson")

AIC(sm.1, sm.2, sm.3)

sm.2

#set gla famlist

famlist.gla<- list(fam.bernoulli(), 
               fam.negative.binomial(0.23784984), 
               fam.negative.binomial(0.06728733),
               fam.negative.binomial(0.013494317))

redata.gla<- subset(redata2016, Region=="GL_alvar")
redata.gla<- droplevels(redata.gla)



aout1<- aster(resp~varb, pred, fam, varb, id, root, data=redata.gla, famlist = famlist.gla)

summary(aout1, show.graph=TRUE, info.tol=1e-16)



#Add dist to seed source and Number of days to germ
aout3<- aster(resp~varb+0+Dist.from.cg.km + No.Days.to.Germ + I(Dist.from.cg.km^2) + I(No.Days.to.Germ^2) + I(2*Dist.from.cg.km*No.Days.to.Germ), pred, fam, varb, id, root, 
            maxiter=8000,data=redata.gla, famlist = famlist.gla)

summary(aout3, show.graph=T, info.tol=1e-13)


#check for coefficients of above model

aout3$coefficients

aout<- aout3


######################################################################
# Estimate Selection Gradient (distance from source)

pout <- predict(aout)
pout <- matrix(pout, nrow = nrow(aout1$x), ncol = ncol(aout1$x))
colnames(pout) <- colnames(aout1$x)
mufit <- pout[, grep("sm", colnames(pout))]


#only needed when >1 year of data
#mufit <- apply(mufit, 1, "sum")

#calcualte mean fitness
wmu <- mufit/mean(mufit)

#perform linear analysis
wmout <- lm(wmu ~ dat2$Dist.from.cg.km)

pre_w<- predict(wmout)

summary(wmout)


#Now try with two predictors dist to source and days to germination

#see new aout model on line 125

#extract two coeff

a1 <- aout$coefficients["Dist.from.cg.km"]
a2 <- aout$coefficients["No.Days.to.Germ"]
a <- c(a1, a2)

A11 <- aout$coefficients["I(Dist.from.cg.km^2)"]
A22 <- aout$coefficients["I(No.Days.to.Germ^2)"]
A12 <- aout$coefficients["I(2 * Dist.from.cg.km * No.Days.to.Germ)"]
A <- matrix(c(A11, A12, A12, A22), 2, 2)

eigen(A, symmetric = TRUE, only.values = TRUE)$values


max8 <- (-solve(A, a)/2)
print(max8)


#Plot dist from source & days to germ, on fitness contours
plot(dat2$Dist.from.cg.km, dat2$No.Days.to.Germ, xlab = "Dist.", ylab = "Days to Germ")
 ufoo <- par("usr")
 nx <- 101
 ny <- 101
 z <- matrix(NA, nx, ny)
 x <- seq(ufoo[1], ufoo[2], length = nx)
 y <- seq(ufoo[3], ufoo[4], length = ny)
 points(max8[1], max8[2], pch = 19)
 for (i in 1:nx) {
   for (j in 1:ny) {
     b <- c(x[i], y[j])
     z[i, j] <- sum(a * b) + as.numeric(t(b) %*% A %*%
                                           + b)
     }
   }
 b <- as.numeric(max8)
 contour(x, y, z, add = TRUE)
 contour(x, y, z, levels = c(0.325), add = TRUE)
 
 
 
 
 
 #OK, cool. But let's be thorough and compare this with the Lande and Arnold (1984) way:
 
 dat2$relfit <- dat2$sm/mean(dat2$sm)
  lout <- lm(relfit ~ Dist.from.cg.km + No.Days.to.Germ + I(Dist.from.cg.km^2) +
                I(No.Days.to.Germ^2) + I(2*Dist.from.cg.km*No.Days.to.Germ), data = dat2)
  summary(lout)


  a1 <- lout$coefficients["Dist.from.cg.km"]
  a2 <- lout$coefficients["No.Days.to.Germ"]
  a <- c(a1, a2)
  
  A11 <- lout$coefficients["I(Dist.from.cg.km^2)"]
  A22 <- lout$coefficients["I(No.Days.to.Germ^2)"]
  A12 <- lout$coefficients["I(2 * Dist.from.cg.km * No.Days.to.Germ)"]
  A <- matrix(c(A11, A12, A12, A22), 2, 2)
  
  eigen(A, symmetric = TRUE, only.values = TRUE)$values
  
  
  max8 <- (-solve(A, a)/2)
  print(max8)
  
  
  #plot OLS (Lande and Arnold) way
 plot(dat2$Dist.from.cg.km, dat2$No.Days.to.Germ, xlab = "Dist.", ylab = "Days to Germ")
  ufoo <- par("usr")
  nx <- 101
  ny <- 101
  z <- matrix(NA, nx, ny)
  x <- seq(ufoo[1], ufoo[2], length = nx)
  y <- seq(ufoo[3], ufoo[4], length = ny)
  points(max8[1], max8[2], pch = 19)
  for (i in 1:nx) {
    for (j in 1:ny) {
      b <- c(x[i], y[j])
      z[i, j] <- sum(a * b) + as.numeric(t(b) %*% A %*%
                                           + b)
    }
  }
  b <- as.numeric(max8)
  contour(x, y, z, add = TRUE)
  contour(x, y, z, levels = c(0.325), add = TRUE)
  

######################################################################



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


aouta<- aster(resp~varb, pred, fam, varb, id, root, data=redata2016, famlist = famlist)

aout<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016, famlist = famlist)

summary(aout, show.graph = TRUE, info.tol = 1e-15)

anova(aouta, aout)

#Generate fitness estimates for populations

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
dim(amat)# makes an 22 x 6 x 22 matrix (2 habitat types and 6 nodes of graphicla model)

#only want means for k'th individual that contribute to expected
#fitness, and want to add only Seedmass.2016 entries

foo<- grepl("sm", vars)
for(k in 1:nPop)
  amat[k, foo, k]<- 1

#check
foo #yes, only last node is "true"; corresponds to Seedmass.2016

#generate predicted valuses using aout object, with renewdata, and amat format
pout.amat<- predict(aout, newdata= renewdata, varvar= varb,
                    idvar= id, root = root, se.fit=TRUE, amat = amat, info.tol = 1e-15)

#combine estimates with standard error, and then round
#to three decimal places
pop<- cbind(pout.amat$fit, pout.amat$se.fit)


rownames(pop)<- as.character(fred$Population)


colnames(pop)<- c("Expected Fitness", "SE")

pop<-round(pop, 3) 

pop