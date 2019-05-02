#Exploring fitness surfaces with 2016 data, and distance from seed source


setwd("C:/Users/Mason Kulbaba/Dropbox/Rscripts/aster-analysis/")

#load data
dat<- read.csv("NV_CG_Experiment2wdist2.csv")


#subset data for 2017 analysis
#dat2<- dat[c("Family.Unique",   "Block.ID", "HabitatType", "Region", "Population",
 #            "Dist.from.cg.km","Germination.Y.N","Survival.Y.N","Survival.Y.N.2017", 
  #           "Flower.Y.N.2016","Flower.Y.N.2017","No.Flowers.2016","Total.Flowers.2017",
   #          "Fruit.Y.N.2016","Fruit.Y.N.2017", "No.Fruit.2016","No.Fruit.2017",
    #         "sm", "sm.2")]

dat2<- dat[c("Family.NonUnique", "Family.Unique", "Block.ID", "HabitatType", "Region",
             "Population","Dist.from.cg.km", "Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016", "Fruit.Y.N.2016", "No.Fruit.2016", "sm")]

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

aout2016b<- aster(resp~varb+ fit:Dist.from.cg.km, pred, fam, varb, id, root, data=redata2016, famlist = famlist)
summary(aout2016b, info.tol=1e-11)

aout2016<- aster(resp~varb+ fit:(Region + Dist.from.cg.km), pred, fam, varb, id, root, data=redata2016, famlist = famlist)
summary(aout2016, info.tol=1e-11)

aout2016d<- aster(resp~varb+ fit:(Region + Dist.from.cg.km + Region*Dist.from.cg.km), pred, fam, varb, id, root, data=redata2016, famlist = famlist)
summary(aout2016d, info.tol=1e-16)


anova(aout2016a, aout2016b, aout2016, aout2016d)




