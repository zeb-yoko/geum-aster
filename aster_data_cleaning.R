
########################################################
########################################################
## Begin cleaning data-> start with 2016 cleaning code##
########################################################
########################################################

#replace NAs with 0 (zero)
#dat2[is.na(dat2)] <- 0

dat2<- read.csv("NV_CG_Experiment2wdist2.csv")


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

#suvive to 2017 but surv.2016=0
subset(dat2, Survival.Y.N==0 & Survival.Y.N.2017==1)# no errors

##################### NEW
dat2$Survival.Y.N[dat2$Survival.Y.N==0 & dat2$Survival.Y.N.2017==1]=1
####################

#flw.no >0, but flw.y.n=0
subset(dat2, Total.Flowers.2017 >0 & Flower.Y.N.2017==0)# 10 errors

#following changes confirmed by Zeb 04/30/2019
dat2$Total.Flowers.2017[dat2$Family.Unique=="MB-CRN_10" & dat2$Block.ID==1]=0
dat2$Total.Flowers.2017[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==7]=0

#set flw.y.n=1 for these above cases
dat2$Flower.Y.N.2017[dat2$Total.Flowers.2017 >0 & dat2$Flower.Y.N.2017==0]=1


#fruit.y.n =1, flw.no=0
subset(dat2, Fruit.Y.N.2017==1 & Total.Flowers.2017 ==0)# 7 errors

#change total number of flowers (2017) to 1. Confirmed by Zeb 04/30/2019
dat2$Total.Flowers.2017[dat$Fruit.Y.N.2017==1 & dat2$Total.Flowers.2017 ==0]=1


#frt number >0, but fruit.y.n==0
subset(dat2, No.Fruit.2017 > 0 & Fruit.Y.N.2017==0)

#the following changes confirmed by Zeb 04/30/2019
dat2$No.Fruit.2017[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==7]=0

#Change remaining issues (from line98) to Fruit.Y.N.2017=1, confirmed by Zeb 04/30/2019
dat2$Fruit.Y.N.2017[dat2$No.Fruit.2017 > 0 & dat2$Fruit.Y.N.2017==0]=1

#seed mass >0 but fruit number 2017 =0
subset(dat2, sm.2 > 0 & No.Fruit.2017==0)

#the following changes confirmed by Zeb 04/30/2019
dat2$sm.2[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==7]=0

# Removing individuals: CAR-NBA_4, MB-MR_32, NAP-CE_6, SD-MUD_10, SD-PMG_NA

dat3<-dat2[!(dat2$sm.2 >0 & dat2$No.Fruit.2017==0),]

dat2<-dat3

#now make new survival to 2017 after winter of 2016 variable

dat2$Surv2017[dat2$Flower.Y.N.2017==1]=1

dat2$Surv2017[is.na(dat2$Surv2017)] <- 0

###################################################################################
# Now begin cleaning the 2018 data


#survival to 2018 but not 2017
subset(dat2, Survival.Y.N.2018==1 & Surv2017==0)# 29 errors

#correct above errors
dat2$Surv2017[dat2$Survival.Y.N.2018==1 & dat2$Surv2017==0]=1


#survival to 2018 but not from greenhouse
subset(dat2, Survival.Y.N.2018==1 & Survival.Y.N==0)#  errors

#flower in 2018 but no survival
subset(dat2, Flowering.Y.N.2018==1 & Survival.Y.N.2018==0)# 7 errors

#fix above errors
dat2$Survival.Y.N.2018[dat2$Flowering.Y.N.2018==1 & dat2$Survival.Y.N.2018==0]=1

#flowering =0, but total flowers >0
subset(dat2, Flowering.Y.N.2018==0 & Total.Flowers.2018 > 0)# 6 errors

#fix above errors
#individual corrections confirmed by Zeb 05/02/2019
dat2$Flowering.Y.N.2018[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0
dat2$Total.Flowers.2018[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0
dat2$seedmass.2018.g.[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0

dat2$ Flowering.Y.N.2018[dat2$Family.Unique=="MAN-KIP_15" & dat2$Block.ID==7]=0
dat2$Total.Flowers.2018[dat2$Family.Unique=="MAN-KIP_15" & dat2$Block.ID==7]=0
dat2$seedmass.2018.g.[dat2$Family.Unique=="MAN-KIP_15" & dat2$Block.ID==7]=0

dat2$Flowering.Y.N.2018[dat2$Family.Unique=="MB-MR_13" & dat2$Block.ID==1]=0
dat2$Total.Flowers.2018[dat2$Family.Unique=="MB-MR_13" & dat2$Block.ID==1]=0
dat2$seedmass.2018.g.[dat2$Family.Unique=="MB-MR_13" & dat2$Block.ID==1]=0

dat2$Flowering.Y.N.2018[dat2$Family.Unique=="MB-MR_38" & dat2$Block.ID==1]=0
dat2$Total.Flowers.2018[dat2$Family.Unique=="MB-MR_38" & dat2$Block.ID==1]=0

dat2$Flowering.Y.N.2018[dat2$Family.Unique=="NAP-ASS_6" & dat2$Block.ID==9]=1

dat2$Flowering.Y.N.2018[dat2$Family.Unique=="SD-MUD_18" & dat2$Block.ID==6]=1



#no flowers produced but fruit set=1
subset(dat2, Fruit.Y.N.2018==1 & Total.Flowers.2018==0)#6 errors
########################
########################


#correct above errors

#individual corrections confirmed by Zeb 05/02/2019
dat2$Total.Flowers.2018[dat2$Family.Unique=="AB-LL_40" & dat2$Block.ID==1]=2
dat2$No.Fruit.2018[dat2$Family.Unique=="AB-LL_40" & dat2$Block.ID==1]=2
dat2$seedmass.2018.g.[dat2$Family.Unique=="AB-LL_40" & dat2$Block.ID==1]=0.1372

dat2$Total.Flowers.2018[dat2$Family.Unique=="CAR-PSR_14" & dat2$Block.ID==12]=1
dat2$No.Fruit.2018[dat2$Family.Unique=="CAR-PSR_14" & dat2$Block.ID==12]=1
dat2$seedmass.2018.g.[dat2$Family.Unique=="CAR-PSR_14" & dat2$Block.ID==12]=0.0563

dat2$No.Fruit.2018[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==12]=0
dat2$seedmass.2018.g.[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==12]=0

dat2$Total.Flowers.2018[dat2$Family.Unique=="MB-CRN_21" & dat2$Block.ID==9]=2
dat2$No.Fruit.2018[dat2$Family.Unique=="MB-CRN_21" & dat2$Block.ID==9]=2
dat2$seedmass.2018.g.[dat2$Family.Unique=="MB-CRN_21" & dat2$Block.ID==9]=0.1306

dat2$Total.Flowers.2018[dat2$Family.Unique=="MB-CRN_29" & dat2$Block.ID==5]=5
dat2$No.Fruit.2018[dat2$Family.Unique=="MB-CRN_29" & dat2$Block.ID==5]=5
dat2$seedmass.2018.g.[dat2$Family.Unique=="MB-CRN_29" & dat2$Block.ID==5]=0.2158

dat2$ Fruit.Y.N.2018[dat2$Family.Unique=="MB-CRN_2" & dat2$Block.ID==10]=0

dat2$ Fruit.Y.N.2018[dat2$Family.Unique=="MB-CRN_17" & dat2$Block.ID==11]=0

dat2$Total.Flowers.2018[dat2$Family.Unique=="MB-MR_40" & dat2$Block.ID==12]=1
dat2$No.Fruit.2018[dat2$Family.Unique=="MB-MR_40" & dat2$Block.ID==12]=1
dat2$seedmass.2018.g.[dat2$Family.Unique=="MB-MR_40" & dat2$Block.ID==12]=0.0595

dat2$Total.Flowers.2018[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0
dat2$No.Fruit.2018[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0
dat2$seedmass.2018.g.[dat2$Family.Unique=="AB-LL_10" & dat2$Block.ID==12]=0

#Number of fruits > 0, but fruit y.n.=0
subset(dat2, Fruit.Y.N.2018==0 & No.Fruit.2018 >0)#22 errors

#fix above errors

#individual correction confirmed by zeb 05/02/2019
dat2$Fruit.Y.N.2018[dat2$Family.Unique=="MB-MR_13" & dat2$Block.ID==1]=0
dat2$No.Fruit.2018[dat2$Family.Unique=="MB-MR_13" & dat2$Block.ID==1]=0

#fix remaining errors (fruit.y.n=1)
dat2$Fruit.Y.N.2018[dat2$Fruit.Y.N.2018==0 & dat2$No.Fruit.2018 >0]=1


#seed mass >0, but fruit number =0
subset(dat2, seedmass.2018.g. >0 & No.Fruit.2018==0)# 2 errors

###########################
#temp. removal############
##########################

dat3<-dat2[!(dat2$seedmass.2018.g. >0 & dat2$No.Fruit.2018==0),]

dat2<-dat3

#Recheck for flower2018=1 but survival=0 for introduced errors
subset(dat2, Flowering.Y.N.2018==1 & Survival.Y.N.2018==0)# 3 new errors

#fix above errors
dat2$Survival.Y.N.2018[dat2$Flowering.Y.N.2018==1 & dat2$Survival.Y.N.2018==0]=1


#survival from 2017 (over winter) to 2018
dat2$Surv2018[dat2$Flowering.Y.N.2018==1]=1

dat2$Surv2018[is.na(dat2$Surv2018)] <- 0

#check for Surv2018 errors
subset(dat2, Surv2018==1 & Survival.Y.N==0)

subset(dat2, Survival.Y.N==1 & Germination.Y.N==0)

subset(dat2, Flowering.Y.N.2018==1 & Surv2018==0)

subset(dat2, Fruit.Y.N.2018==1 & Total.Flowers.2018 == 0)
#######################################################################################

subset(dat2, Survival.Y.N==1 & Germination.Y.N==0)# 0 errors

subset(dat2, Surv2018==1 & Survival.Y.N==0)# 0 errors

subset(dat2, Flowering.Y.N.2018==1 & Surv2018==0)# 0 errors

subset(dat2, Fruit.Y.N.2018==1 & Total.Flowers.2018==0)# 0 errors



########################################################################
#Isolate "other" covariate data

cov<- dat2

#what's actually in here?
names(cov)

#isolate what we actually want
cov2<-cov[c("Family.Unique", "Block.ID", "Germination.Date","No.Days.to.Germ","True.Leaf", 
            "No.Days.to.TrueLeaf","Germ.to.TrueLeaf")]



