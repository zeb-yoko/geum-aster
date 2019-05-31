##Resolved Models##
##load libraries##
library(lme4); library(lmerTest); library(effects)
library(tidyverse); library(rptR);library(emmeans)
#install.packages("QGglmm")
library(QGglmm)
##new package##
#install.packages("fitdistrplus")
library(fitdistrplus)
##NOTE: Distributions for variables in ASTER analysis taken from there##
##phenology variables primarily determined in aster-distros script##
#old datasheet
#df <- read.csv('NV_CG_Experiment-mg.csv')
##cleaned data from aster model datasheet in git repo##
df <- read.csv("full_clean_geum_experiment.csv")


###########################################################
##INDIVIDUAL MODELS PER TRAIT TO CALCULATE HERITABILITIES##
###########################################################

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
##model statement##
germ.mod <- glmer(No.Days.to.Germ~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = df,
						family=poisson(link=log))
##View outputs##
summary(germ.mod)
hist(residuals(germ.mod))
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(germ.mod))[, c('grp','vcov')]
intercept <- fixef(germ.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
herit2

##Create a table to compile heritabilities## 
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
f1g <- fitdist(TrueLeaf$No.Days.to.TrueLeaf, "lnorm")
f2g <- fitdist(TrueLeaf$No.Days.to.TrueLeaf, "pois")
##works##
plot(f1g)
summary(f1g)
plot(f2g)
summary(f2g)
##TrueLeaf##
TrueLeaf.mod <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = df,
							 family=poisson(link=log))

summary(TrueLeaf.mod)
hist(residuals(TrueLeaf.mod))
##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(TrueLeaf.mod))[, c('grp','vcov')]
intercept <- fixef(TrueLeaf.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
herit2
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart?##
herit2
h2[2,1] <- "TrueLeaf"
h2[2,2] <- "2015"
h2[2,3] <- herit2$h2.obs
h2
#############################

##DTFF 2016##
##CHECK DISTRO--Gamma working, but not resolved in QGparams##
#############################
##check distribution first##
#make date number?
flr.16 <- filter(df, Flower.Y.N.2016 >= 1)
#trueleaf#
hist(flr.16$no.Planting.to.DTFF)
FLR <- flr.16[!is.na(flr.16$no.Planting.to.DTFF),]
descdist(FLR$no.Planting.to.DTFF, boot = 100)
f1g <- fitdist(FLR$no.Planting.to.DTFF, "pois")
f2g <- fitdist(FLR$no.Planting.to.DTFF, "gamma")
boxplot(FLR$no.Planting.to.DTFF~Region, data = FLR)
f4g <- fitdistr(FLR$no.Planting.to.DTFF, "gamma")
f4g
k <- f4g$estimate[1]
theta <- 1/(f4g$estimate[2])
plot(f2g)
##works##
plot(f1g)
summary(f1g)
f1g <- fitdist(FLR$no.Planting.to.DTFF, "pois")
##works##
plot(f1g)
summary(f1g)
##working poisson model--crappy, explains little variance##
dtff.pos<- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr.16,
					  family = poisson(link=log))
hist(residuals(dtff.pos))
summary(dtff.pos)
dtff.gam <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr.16,
					  family = Gamma(link=log))
hist(residuals(dtff.gam))
summary(dtff.gam)
hist(residuals(dtff.gam))
dtff.gam
##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.gam))[, c('grp','vcov')]
intercept <- fixef(dtff.gam)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
herit2
##IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma
##Pull from fitdistr##
f4g <- fitdistr(FLR$no.Planting.to.DTFF, "gamma")
f4g
k <- f4g$estimate[1]
theta <- 1/(f4g$estimate[2])

e <- exp(1)
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){exp(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
									d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam
str(dtff.gam)							

##run rptPoisson to see how much variance explained by effects##
rpt.dtff<-rptGamma(formula = no.Planting.to.DTFF~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), 
							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
rpt.dtff
##NOTE--BASICALLY ZERO BECAUSE MODEL EXPLAINS ESSENTIALLY NOTHING##

h2[3,1] <- "Date to first flower"
h2[3,2] <- "2016"
h2[3,3] <- herit2$h2.obs
h2
############################

##No. flowers 2016##
############################
n.flr.mod <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = df,
						 family=neg.bin(theta = 1.12571436))
n.flr.mod2 <- glmer.nb(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = df)
n.flr.mod2
n.flr.mod3 <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = df,
						 family=neg.bin(theta = 0.1874))
hist(residuals(n.flr.mod3))
hist(residuals(n.flr.mod2))
hist(residuals(n.flr.mod))

vars <- as.data.frame(VarCorr(n.flr.mod))[, c('grp','vcov')]
intercept <- fixef(n.flr.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  1.12571436
##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
herit2

##Run rptPoisson to see how much variance explained by model effects##
#rpt.nflr<-rptPoisson(formula = No.Flowers.2016~Region + (1 | Population) + 
#								(1 | Family.Unique) + (1 | Block.ID), 
#							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
#							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
#rpt.nflr

##NOTE--With Negative binomial and theta from Mason--residuals Bad##

h2[4,1] <- "Number of Flowers"
h2[4,2] <- "2016"
h2[4,3] <- herit2$h2.obs
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
							family=neg.bin(theta = 3.5074887))
n.fruit.out <-	summary(n.fruit.mod)
n.fruit.out
#store residuals
n.fruit.resid <- residuals(n.fruit.mod)

vars <- as.data.frame(VarCorr(n.fruit.mod))[, c('grp','vcov')]
intercept <- fixef(n.fruit.mod)['(Intercept)']

##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <- 3.5074887
##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
herit2

#rpt.nfruit<-rptPoisson(formula = No.Fruit.2016~Region + (1 | Population) + 
#								(1 | Family.Unique) + (1 | Block.ID), 
#							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
#							data = flr.16, link = "log", nboot =0, ratio =T, adjusted =F)
#rpt.nfruit

##add to table
h2[5,1] <- "Number of Fruit"
h2[5,2] <- "2016"
h2[5,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2016##
###########################
View(df)
hist(df$sm)
hist(df1$Seedmass.2016)
df1<-df[!is.na(df$sm),]
hist(df1$sm)
f1g <- fitdist(df1$sm, "pois")
f1g <- fitdist(df1$sm, "nbinom")
plot(f1g)
n.seed.mod <- glmer(sm~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = df1,
							family = neg.bin(theta = 1.0341103))
n.seed.out <-	summary(n.seed.mod)
n.seed.out
n.seed.mod2 <- glmer.nb(Seedmass16.mg~Region + (1 | Population) + 
							  	(1 | Family.Unique) + (1 | Block.ID), data = df1)#,
n.seed.out2 <-	summary(n.seed.mod2)

n.seed.out
hist(residuals(n.seed.mod))
##negbin is not a datatype rptR works with##
#rpt.seedmass<-rpt(formula = sm~Region + (1 | Population) + 
#				  	(1 | Family.Unique) + (1 | Block.ID), 
#							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
#							data = df1, datatype = "nbinom", link = "log", nboot =0, ratio =T, adjusted =F)
#rpt.seedmass
vars <- as.data.frame(VarCorr(n.seed.mod))[, c('grp','vcov')]
intercept <- fixef(n.seed.mod)['(Intercept)']
vars
intercept
##Negative binomial needs additional (dispersion) parameter, theta##
theta <- 1.0341103
##latent region mean##
mu <- intercept
##Additive Variance--note 4x for half-sib design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##Total variance##
vp <- sum(vars[,"vcov"])
vp
##latent scale narrow-sense heritability##
lh2 <- va/vp
lh2

##put in QGparams--get observed scale h2##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit2
##add to table
h2[6,1] <- "Seedmass"
h2[6,2] <- "2016"
h2[6,3] <- 0  #herit2$h2.obs
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
					  family = Gamma(link=log))
hist(residuals(dtff.mod))
#rpt.dtff<-rptPoisson(formula = DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
#								(1 | Family.Unique) + (1 | Block.ID), 
#							grname = c("Fixed", "Block.ID", "Population", "Family.Unique"),
#							data = flr.17, link = "log", nboot =0, ratio =T, adjusted =F)
#rpt.dtff
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
##model statement##
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
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
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
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
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
##pull data from model##
vars <- as.data.frame(VarCorr(seed17.mod))[, c('grp','vcov')]
intercept <- fixef(seed17.mod)['(Intercept)']
vars
intercept
##additonal negative binomial parameter##
theta <- 1.347

##latent region mean##
mu <- intercept
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
vp <- sum(vars[,"vcov"])
vp

##latent scale heritability##
lh2 <- va/vp
lh2

##put in QGparams--get observed scale h2##
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

##model statement##
dtff18.mod<- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr.18,
						 family = poisson(link = log))
summary(dtff18.mod)
hist(residuals(dtff18.mod))
vars <- as.data.frame(VarCorr(dtff18.mod))[, c('grp','vcov')]
intercept <- fixef(dtff18.mod)['(Intercept)']

##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##

##region mean##
mu <- intercept
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*vars[vars[["grp"]] == "Family.Unique", "vcov"]
va
##total variance in trait##
vp <- sum(vars[,"vcov"])
vp
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2

##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp, model = "Poisson.log")
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as h2 of trait 
#across all regions, could be worth breaking apart##
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
f1g <- fitdistr(flr18$DtB.Oday.2018, "normal")
f2g <- fitdistr(flr18$DtB.Oday.2018, "poisson")
f3g <- fitdistr(flr18$DtB.Oday.2018, "negative binomial")
f1g
AIC(f1g, f2g)
plot(f1g)
plot(f2g)
##Roughly normal? OR poisson??##
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
