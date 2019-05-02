##Resolved Models##
##load libraries##
library(lme4); library(lmerTest); library(effects)
library(tidyverse); library(rptR);library(emmeans)
#install.packages("QGglmm")
library(QGglmm)
##new package##
#install.packages("fitdistrplus")
library(fitdistrplus)
df <- read.csv('NV_CG_Experiment-mg.csv')

##germination##
##############################
##check distribution first##
hist(df$No.Days.to.Germ)
class(df$No.Days.to.Germ)
range(df$No.Days.to.Germ)
##df w/out NA's for germination##
##for fitdist##
germ <- df[!is.na(df$No.Days.to.Germ),]
descdist(germ$No.Days.to.Germ, boot = 100)
f1g <- fitdist(germ$No.Days.to.Germ, "pois")
##works##
plot(f1g)
summary(f1g)
##germination##
germ.mod <- glmer(No.Days.to.Germ~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = df,
						family=poisson(link=log))
summary(germ.mod)
hist(residuals(germ.mod))
##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(germ.mod))[, c('grp','vcov')]
intercept <- fixef(germ.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, is it worth breaking apart?##
herit2
?QGparams
col.classes = c("character", "character", "numeric")
col.names = c("Trait", "Year", "Heritability")
h2 <-read.table(text = "",colClasses = col.classes, col.names =col.names)
h2[1,1] <- "Germination"
h2[1,2] <- "2015"
h2[1,3] <- herit2$h2.obs
#############################

##trueleaf##
#############################
##check distribution first##
hist(df$No.Days.to.TrueLeaf)
class(df$No.Days.to.TrueLeaf)
range(df$No.Days.to.TrueLeaf)
##df w/out NA's for TrueLeafination##
##for fitdist##
TrueLeaf <- df[!is.na(df$No.Days.to.TrueLeaf),]
descdist(log(TrueLeaf$No.Days.to.TrueLeaf), boot = 100)
f1g <- fitdist(log(TrueLeaf$No.Days.to.TrueLeaf), "lnorm")
##works##
plot(f1g)
summary(f1g)
##TrueLeaf##
TrueLeaf.mod <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = df,
							 family=poisson(link=log))

summary(TrueLeaf.mod)
hist(residuals(TrueLeaf.mod))
##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(TrueLeaf.mod))[, c('grp','vcov')]
intercept <- fixef(TrueLeaf.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, is it worth breaking apart?##
herit2
h2[2,1] <- "TrueLeaf"
h2[2,2] <- "2015"
h2[2,3] <- herit2$h2.obs
h2
#############################

##DTFF 2016##
#############################
##check distribution first##
#make date number?
flr.16 <- filter(df, Flower.Y.N.2016 >= 1)
#trueleaf#
hist(flr.16$no.Planting.to.DTFF)
FLR <- flr.16[!is.na(flr.16$no.Planting.to.DTFF),]
descdist(FLR$no.Planting.to.DTFF, boot = 100)
f1g <- fitdist(FLR$no.Planting.to.DTFF, "pois")
##works##
plot(f1g)
summary(f1g)
dtff.mod<- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr.16,
					  family = poisson(link=log))
hist(residuals(dtff.mod))
summary(dtff.mod)
hist(residuals(dtff.mod))
##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.mod))[, c('grp','vcov')]
intercept <- fixef(dtff.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <-vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, is it worth breaking apart?##
herit2
##NOTE--BASICALLY ZERO BECAUSE MODEL EXPLAINS ESSENTIALLY NOTHING##
library(rptR)
rpt.dtff<-rptPoisson(formula = no.Planting.to.DTFF~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtff
h2[3,1] <- "Date to first flower"
h2[3,2] <- "2016"
h2[3,3] <- 0
h2
############################

##No. flowers 2016##
############################
n.flr.mod <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = df,
						 family=poisson(link=log))
hist(residuals(n.flr.mod))

vars <- as.data.frame(VarCorr(n.flr.mod))[, c('grp','vcov')]
intercept <- fixef(n.flr.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, is it worth breaking apart?##
herit2
rpt.dtff<-rptPoisson(formula = No.Flowers.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtff
##NOTE--BASICALLY ZERO BECAUSE MODEL EXPLAINS ESSENTIALLY NOTHING##
h2[4,1] <- "Number of Flowers"
h2[4,2] <- "2016"
h2[4,3] <- 0
h2
###########################

##Number of fruit 2016##
###########################
hist(df$No.Fruit.2016)
df1 <- filter(df, No.Fruit.2016 >= 0)
##only option is poisson##
f1g <- fitdist(df1$No.Fruit.2016, "pois")
plot(f1g)
n.fruit.mod <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), data = df1,
							family=poisson(link=log))
n.fruit.out <-	summary(n.fruit.mod)
n.fruit.out
#store residuals
n.fruit.resid <- residuals(n.fruit.mod)
rpt.dtff<-rptPoisson(formula = No.Fruit.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtff
vars <- as.data.frame(VarCorr(n.fruit.mod))[, c('grp','vcov')]
intercept <- fixef(n.fruit.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
herit2
##add to table
h2[5,1] <- "Number of Fruit"
h2[5,2] <- "2016"
h2[5,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2016##
###########################
hist(df$sm)
hist(df1$sm)
df1<-df[!is.na(df$sm),]
f1g <- fitdist(df1$sm, "pois")
f1g <- fitdist(df1$sm, "nbinom")
plot(f1g)
n.seed.mod <- glmer.nb(Seedmass16.mg~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = df1)#,
n.seed.out <-	summary(n.seed.mod)
n.seed.out
hist(residuals(n.seed.mod))
rpt.seedmass<-rpt(formula = Seedmass16.mg~Region + (1 | Population) + 
				  	(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = df1, datatype = "nbinom", link = "log", nboot =0, ratio =T, adjusted =F)
rpt.seedmass
vars <- as.data.frame(VarCorr(n.seed.mod))[, c('grp','vcov')]
intercept <- fixef(n.seed.mod)['(Intercept)']
vars
intercept
theta <- 1.5239
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit2
##add to table
h2[6,1] <- "Seedmass"
h2[6,2] <- "2016"
h2[6,3] <- herit2$h2.obs
h2
###########################

########2017 Season#############

##DTFF 2017##
###########################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(df$DTFF.Ordinal.Day.2017)
flr.17$DTFF.Ordinal.Day.2017 <-as.numeric(flr.17$DTFF.Ordinal.Day.2017)
#DTFF17#
#make date number?
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(flr.17$DTFF.Ordinal.Day.2017)
flr17 <- flr.17[!is.na(flr.17$DTFF.Ordinal.Day.2017),]
f1g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "norm")
f2g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "pois")
dtff.mod<- glmer(DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr.17,
					  ##first day = 107
					  family = poisson(link=log))
hist(residuals(dtff.mod))
rpt.dtff<-rptPoisson(formula = DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr.17, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtff
##NOTE--BASICALLY ZERO BECAUSE MODEL EXPLAINS ESSENTIALLY NOTHING##

##add to table
h2[7,1] <- "Date to First Flower"
h2[7,2] <- "2017"
h2[7,3] <- 0
h2
###########################

##Date to bolt 2017##
###########################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(df$DtB.O.Day.2017)
flr.17$DtB.O.Day.2017 <-as.numeric(flr.17$DtB.O.Day.2017)
flr17 <- flr.17[!is.na(flr.17$DtB.O.Day.2017),]
f1g <- fitdist(flr17$DtB.O.Day.2017, "norm")
f2g <- fitdist(flr17$DtB.O.Day.2017, "pois")
plot(f1g)
plot(f2g)
##model##
dtb.mod<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(1 | Family.Unique) + (1 | Block.ID), data = flr.17,
					 ##first day = ~114
					 family = poisson(link=log))
hist(residuals(dtb.mod))
rpt.dtb<-rpt(formula = DtB.O.Day.2017~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr17, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtb
vars <- as.data.frame(VarCorr(dtb.mod))[, c('grp','vcov')]
intercept <- fixef(dtb.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit2
##add to table
h2[8,1] <- "Date to bolt"
h2[8,2] <- "2017"
h2[8,3] <- herit2$h2.obs
h2
###########################

##Date to Fruit##
###########################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(log(df$Fruit.O.Day.2017))
flr17 <- flr.17[!is.na(flr.17$Fruit.O.Day.2017),]
flr.17$Fruit.O.Day.2017 <-as.numeric(flr.17$Fruit.O.Day.2017)
f1g <- fitdist(flr17$Fruit.O.Day.2017, "norm")
f2g <- fitdist(flr17$Fruit.O.Day.2017, "pois")
plot(f1g)
plot(f2g)
dtfr.mod<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = ~121
					  family = poisson(link=log))
rpt.dtfr<-rptPoisson(formula = Fruit.O.Day.2017~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), 
						  grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
						  data = flr17, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtfr
vars <- as.data.frame(VarCorr(dtfr.mod))[, c('grp','vcov')]
intercept <- fixef(dtfr.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp,  model = "Poisson.log")
herit2
##add to table
h2[9,1] <- "Date to fruit"
h2[9,2] <- "2017"
h2[9,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2017##
###########################
flr.17 <- filter(df, Seedmass17.mg >=1)
#View(flr.17)
hist(df$Seedmass17.mg)
summary(df$Seedmass17.mg)
flr17 <- flr.17[!is.na(flr.17$Seedmass17.mg),]
f1g <- fitdist(flr17$Seedmass17.mg, "norm")
f2g <- fitdist(flr17$Seedmass17.mg, "pois")
f3g <- fitdist(flr17$Seedmass17.mg, "nbinom")
plot(f1g)
plot(f2g)
plot(f3g)
##negative binomial rather than poisson, bc density on right tail high##
##when ran glmer.nb: Negative Binomial(1.347)  ( log )##
seed17.mod<- glmer(Seedmass17.mg~1 +Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = flr.17,
						  family = negative.binomial(1.347))
summary(seed17.mod)
hist(residuals(seed17.mod))
rpt.seed17<-rpt(formula = DtB.O.Day.2017~Region + (1 | Population) + 
				 	(1 | Family.Unique) + (1 | Block.ID), 
				 grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
				 data = flr17,  datatype = "nbinom", link = "log", nboot =0, ratio =T, adjusted =F)
rpt.seed17
?rptR
vars <- as.data.frame(VarCorr(seed17.mod))[, c('grp','vcov')]
intercept <- fixef(seed17.mod)['(Intercept)']
vars
intercept
theta <- 1.347
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit2
##add to table
h2[10,1] <- "Seedmass (mg)"
h2[10,2] <- "2017"
h2[10,3] <- herit2$h2.obs
h2
###########################

########2018 Season###########

##DTFF##
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(df$DTFF.18.Oday)
hist(flr.18$DTFF.18.Oday)
flr18 <- flr.18[!is.na(flr.18$DTFF.18.Oday),]
f1g <- fitdist(flr18$DTFF.18.Oday, "norm")
f2g <- fitdist(flr18$DTFF.18.Oday, "pois")
plot(f1g)
plot(f2g)

dtff18.mod<- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr.18,
						 family = poisson(link = log))
summary(dtff18.mod)
hist(residuals(dtff18.mod))
vars <- as.data.frame(VarCorr(dtff18.mod))[, c('grp','vcov')]
intercept <- fixef(dtff18.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp,  model = "Poisson.log")
herit2
##add to table
h2[11,1] <- "Date to first flower"
h2[11,2] <- "2018"
h2[11,3] <- herit2$h2.obs
h2
###########################

##Date to bolt##
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(df$DtB.Oday.2018)
flr18 <- flr.18[!is.na(flr.18$DtB.Oday.2018),]
f1g <- fitdist(flr18$DtB.Oday.2018, "norm")
f2g <- fitdist(flr18$DtB.Oday.2018, "pois")

##Roughly normal-NOT poisson##
dtb18.mod<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						##went with (neg) Gamma dist
						family = poisson(link=log))
dtb18.out <-	summary(dtb18.mod)
dtb18.out
hist(residuals(dtb18.mod))

vars <- as.data.frame(VarCorr(dtb18.mod))[, c('grp','vcov')]
intercept <- fixef(dtb18.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp,  model = "Poisson.log")
herit2
##add to table
h2[12,1] <- "Date to bolt"
h2[12,2] <- "2018"
h2[12,3] <- herit2$h2.obs
h2
###########################

##Date to Fruit 2018##
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$Date.to.Fruit.Oday.2018)
flr18 <- flr.18[!is.na(flr.18$Date.to.Fruit.Oday.2018),]
f1g <- fitdist(flr18$Date.to.Fruit.Oday.2018, "norm")
f2g <- fitdist(flr18$Date.to.Fruit.Oday.2018, "pois")
plot(f1g)
plot(f2g)

dtfr18.mod<- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = poisson(link=log))
dtfr18.mod2<- lmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18)
						 
summary(dtfr18.mod2)
hist(residuals(dtfr18.mod2))
dtfr18.out <-	summary(dtfr18.mod)
dtfr18.out
hist(residuals(dtfr18.mod))
vars <- as.data.frame(VarCorr(dtfr18.mod2))[, c('grp','vcov')]
intercept <- fixef(dtfr18.mod2)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp,  model = "Gaussian")
herit2
##add to table
h2[13,1] <- "Date to Fruit"
h2[13,2] <- "2018"
h2[13,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2018##
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$Seedmass18.mg)
flr18 <- flr.18[!is.na(flr.18$Seedmass18.mg),]
f1g <- fitdist(flr18$Seedmass18.mg, "norm")
f2g <- fitdist(flr18$Seedmass18.mg, "pois")
fg3 <- fitdist(log(flr18$Seedmass18.mg), "norm")
fg4 <- fitdist(log(flr18$Seedmass18.mg), "nbinom")
plot(f1g)
plot(f2g)
plot(fg3)
plot(fg4)
seeds18.mod<- glmer(Seedmass18.mg~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						  family = negative.binomial(1.2657))
seeds18.out <-	summary(seeds18.mod)
seeds18.out
hist(residuals(seeds18.mod))
vars <- as.data.frame(VarCorr(seeds18.mod))[, c('grp','vcov')]
intercept <- fixef(seeds18.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
theta <- 1.2657
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = 1.2657, model = "negbin.log")
herit2
##add to table
h2[14,1] <- "seedmass"
h2[14,2] <- "2018"
h2[14,3] <- herit2$h2.obs
h2
###########################

##Number of fruit  2018##
###########################
f1g <- fitdist(flr18$No.Fruit.2018, "norm")
f2g <- fitdist(flr18$No.Fruit.2018, "pois")
plot(f1g)
plot(f2g)
hist(flr.18$No.Fruit.2018)
nfr18.mod<- glmer(No.Fruit.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						##went with poisson
						family=poisson(link=log))
nfr18.out <-	summary(nfr18.mod)
nfr18.out
hist(residuals(nfr18.mod))
vars <- as.data.frame(VarCorr(nfr18.mod))[, c('grp','vcov')]
intercept <- fixef(nfr18.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp
#latent scale heritability#
lh2 <- va/vp
lh2
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
herit2
##add to table
h2[15,1] <- "Number of Fruit"
h2[15,2] <- "2018"
h2[15,3] <- herit2$h2.obs
h2
###########################


str(h2)
write.csv(h2, "heritabilities from pheno-fit3.csv")
