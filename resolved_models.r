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
df <- read.csv("cleaned_NV_CG_experiment_data.csv")
View(df)

###########################################################
##INDIVIDUAL MODELS PER TRAIT TO CALCULATE HERITABILITIES##
###########################################################

##germination##1
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

##trueleaf##2
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
f3g <- fitdist(TrueLeaf$No.Days.to.TrueLeaf, "gamma")
plot(f1g)
summary(f1g)
plot(f2g)
summary(f2g)
plot(f3g)
f1 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "lognormal")
f2 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "poisson")
f3 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "Gamma")
AIC(f1,f2,f3)

##Model Statement##
TrueLeaf.mod <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = df,
							 family=poisson(link=log))
##output##
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
##NOTE: with mu being regional(fixed effect) mean, heritability is calculated as 
#h2 of trait across all regions, could be worth breaking apart##
herit2

##Add to Table##
h2[2,1] <- "TrueLeaf"
h2[2,2] <- "2015"
h2[2,3] <- herit2$h2.obs
h2
#############################

##DTFF 2016##3
#############################
##check distribution first##
##remove NAs for descdist##
flr.16 <- filter(df, Flower.Y.N.2016 >= 1)
hist(flr.16$no.Planting.to.DTFF)
FLR <- flr.16[!is.na(flr.16$no.Planting.to.DTFF),]
descdist(FLR$no.Planting.to.DTFF, boot = 100)
f1g <- fitdist(FLR$no.Planting.to.DTFF, "pois")
f2g <- fitdist(FLR$no.Planting.to.DTFF, "gamma")
plot(f1g)
plot(f2g)
##Gamma seems best, will need parameters (from fitdistr)##
g <- fitdistr(FLR$no.Planting.to.DTFF, "gamma")
g
k <- g$estimate[1]
theta <- 1/(g$estimate[2])

##Model Statement##
dtff.gam <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr.16,
					  family = Gamma(link=log))
##model output##
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

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(FLR$no.Planting.to.DTFF, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
									d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam

h2[3,1] <- "Date to first flower"
h2[3,2] <- "2016"
h2[3,3] <- herit.gam$h2.obs
h2
############################

##No. flowers 2016##4
############################
nfl<-df[!is.na(df$No.Flowers.2016),]
f1 <- fitdistr(nfl$No.Flowers.2016, "normal")
f2 <- fitdistr(nfl$No.Flowers.2016, "Poisson")
f3 <- fitdistr(nfl$No.Flowers.2016, "gamma")
f4 <- fitdistr(nfl$No.Flowers.2016, "negative binomial")
hist(log(nfl$No.Flowers.2016))
hist(nfl$No.Flowers.2016)

AIC(f1,f2,f3,f4)
f4 #theta 1.29674244

n.flr.mod <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = nfl,
						 family=neg.bin(theta = 1.29674244))
summary(n.flr.mod)
n.flr.mod2 <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = nfl,
						  family = gaussian(link=log))
n.flr.mod2

hist(residuals(n.flr.mod))
hist(residuals(n.flr.mod2))

vars <- as.data.frame(VarCorr(n.flr.mod))[, c('grp','vcov')]
intercept <- fixef(n.flr.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  1.29674244
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

##Number of fruit 2016##5
###########################
hist(df$No.Fruit.2016)
nfr<-df[!is.na(df$No.Fruit.2016),]
f1 <- fitdistr(nfr$No.Fruit.2016, "normal")
f2 <- fitdistr(nfr$No.Fruit.2016, "Poisson")
f3 <- fitdistr(nfr$No.Fruit.2016, "gamma")
f4 <- fitdistr(nfr$No.Fruit.2016, "negative binomial")
AIC(f1,f2,f3,f4)
##Negative binomial best fit##
f4 #theta 0.16237431

n.fruit.mod <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), data = nfr,
							family=neg.bin(theta = 0.16237431))
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
theta <- 0.16237431
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
h2[5,1] <- "Number of Fruit"
h2[5,2] <- "2016"
h2[5,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2016##6
###########################
View(df)
df1 <-df[!is.na(df$sm),]
f1 <- fitdistr(df1$sm, "normal")
f2 <- fitdistr(df1$sm, "Poisson")
f3 <- fitdistr(df1$sm, "gamma")
f4 <- fitdistr(df1$sm, "negative binomial")
AIC(f1,f2,f3,f4)
hist(df1$sm)
f4 #theta = 0.013268236
n.seed.mod <- glmer(sm~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = df1,
							family = neg.bin(theta = 0.013268236))
n.seed.out <-	summary(n.seed.mod)
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
theta <- 0.013268236
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

##DTFF 2017##7
###########################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(df$DTFF.Ordinal.Day.2017)
flr.17$DTFF.Ordinal.Day.2017 <-as.numeric(flr.17$DTFF.Ordinal.Day.2017)
#DTFF17#
#make date number?
flr17 <- df[!is.na(df$DTFF.Ordinal.Day.2017),]
descdist(flr17$DTFF.Ordinal.Day.2017)
f1 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "normal")
f2 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "poisson")
f3 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "Gamma")
#f4 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "negative binomial")## error-fails
AIC(f1,f2,f3)

dtff.mod<- glmer(DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.mod))

##intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.mod))[, c('grp','vcov')]
intercept <- fixef(dtff.mod)['(Intercept)']
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

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
								 d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam
##add to table
h2[7,1] <- "Date to First Flower"
h2[7,2] <- "2017"
h2[7,3] <- herit.gam$h2.obs
h2
###########################

##Date to bolt 2017##8
###########################
flr17 <- df[!is.na(df$DtB.O.Day.2017),]
f1 <- fitdistr(flr17$DtB.O.Day.2017, "normal")
f2 <- fitdistr(flr17$DtB.O.Day.2017, "poisson")
f3 <- fitdistr(flr17$DtB.O.Day.2017, "Gamma")
#f4 <- fitdistr(flr17$DtB.O.Day.2017, "negative binomial")## error-fails
AIC(f1,f2,f3)
##model statement##
dtb.mod<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.mod))

rpt.dtb<-rptPoisson(formula = DtB.O.Day.2017~Region + (1 | Population) + 
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

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr17$DtB.O.Day.2017, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
								 d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam
##add to table
h2[8,1] <- "Date to bolt"
h2[8,2] <- "2017"
h2[8,3] <- herit.gam$h2.obs
h2
###########################

##Date to Fruit##9
###########################
flr17 <- df[!is.na(df$Fruit.O.Day.2017),]
hist(flr17$Fruit.O.Day.2017)
f1 <- fitdistr(flr17$Fruit.O.Day.2017, "normal")
f2 <- fitdistr(flr17$Fruit.O.Day.2017, "poisson")
f3 <- fitdistr(flr17$Fruit.O.Day.2017, "Gamma")
#f4 <- fitdistr(flr17$Fruit.O.Day.2017, "negative binomial")## error-fails
AIC(f1,f2,f3)
###################Gamma model######
dtfr.mod<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))

hist(residuals(dtfr.mod))
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
#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr17$Fruit.O.Day.2017, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
								 d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam
##add to table
h2[9,1] <- "Date to fruit"
h2[9,2] <- "2017"
h2[9,3] <- herit.gam$h2.obs
h2
###########################

##No. flowers 2017##10
############################
flr.17 <- filter(df, Total.Flowers.2017 >=1)
#View(flr.17)
flr17 <- df[!is.na(df$Total.Flowers.2017),]
hist(flr17$Total.Flowers.2017)
fi <- fitdistr(flr17$Total.Flowers.2017, "normal")
fj <- fitdistr(flr17Total.Flowers.2017, "poisson")
fk <- fitdistr(flr17$Total.Flowers.2017, "Gamma")
fl <- fitdistr(flr17$Total.Flowers.2017, "negative binomial")## error-fails
AIC(fi,fl)
fl #theta = 2.4967562
fitdistr(flr17$Total.Flowers.2017, "negative binomial")

n.flr.mod <- glmer(Total.Flowers.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = 2.4967562))
hist(residuals(n.flr.mod))

vars <- as.data.frame(VarCorr(n.flr.mod))[, c('grp','vcov')]
intercept <- fixef(n.flr.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  2.4967562
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

##NOTE--With Negative binomial and theta from Mason--residuals ok--model failed to converge##

h2[10,1] <- "Number of Flowers"
h2[10,2] <- "2017"
h2[10,3] <- herit2$h2.obs
h2
############################

##No. fruit 2017##11
############################
flr17 <- df[!is.na(df$No.Fruit.2017),]
summary(flr17$No.Fruit.2017)
hist(flr17$No.Fruit.2017)
f1g <- fitdistr(flr17$No.Fruit.2017, "normal")
f2g <- fitdistr(flr17$No.Fruit.2017, "poisson")
f3g <- fitdistr(flr17$No.Fruit.2017, "negative binomial")
AIC(f1g,f2g,f3g)
f3g #theta 2.1804933
n.frt.mod <- glmer(No.Fruit.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = .237))

hist(residuals(n.frt.mod))

vars <- as.data.frame(VarCorr(n.frt.mod))[, c('grp','vcov')]
intercept <- fixef(n.frt.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  0.237
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

##NOTE--With Negative binomial and theta from Mason--residuals ok--model failed to converge##

h2[11,1] <- "Number of Fruit"
h2[11,2] <- "2017"
h2[11,3] <- herit2$h2.obs
h2
###########################

##Seedmass 2017##12
###########################
View(df)
flr17 <- df[!is.na(df$sm.2),]
hist(flr17$sm.2)
f1g <- fitdistr(flr17$sm.2, "normal")
f2g <- fitdistr(flr17$sm.2, "poisson")
f3g <- fitdistr(flr17$sm.2, "Gamma")
f4g <- fitdistr(flr17$sm.2, "negative binomial")
AIC(f1g, f2g, f4g)
f4g #theta = 0.50654046
##negative binomial rather than poisson, bc density on right tail high##
##when ran glmer.nb: Negative Binomial(1.347)  ( log )##
seed17.mod<- glmer(sm.2~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
							 family = negative.binomial(0.50654046))
summary(seed17.mod)
hist(residuals(seed17.mod))

##pull data from model##
vars <- as.data.frame(VarCorr(seed17.mod))[, c('grp','vcov')]
intercept <- fixef(seed17.mod)['(Intercept)']
vars
intercept
##additonal negative binomial parameter##
theta <- 0.50654046

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
h2[12,1] <- "Seedmass (mg)"
h2[12,2] <- "2017"
h2[12,3] <- herit2$h2.obs
h2
##############################

########2018 Season###########

##DTFF##13
###########################
hist(df$DTFF.18.Oday)
flr18 <- df[!is.na(df$DTFF.18.Oday),]
descdist(flr18$DTFF.18.Oday, boot = 100)
f1g <- fitdist(flr18$DTFF.18.Oday, "norm")
f2g <- fitdist(flr18$DTFF.18.Oday, "pois")
f3g <- fitdist(flr18$DTFF.18.Oday, "gamma")
f4g <- fitdist(flr18$DTFF.18.Oday, "lnorm")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)

f1 <- fitdistr(flr18$DTFF.18.Oday, "normal")
f2 <- fitdistr(flr18$DTFF.18.Oday, "Poisson")
f3 <- fitdistr(flr18$DTFF.18.Oday, "gamma")
f4 <- fitdistr(flr18$DTFF.18.Oday, "lognormal")
f5 <- fitdistr(flr18$DTFF.18.Oday, "negative binomial")
AIC(f1,f2,f3,f4)

##model statement##
dtff18.mod<- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link = log))
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

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr18$DTFF.18.Oday, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
								 d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##

##add to table
h2[13,1] <- "Date to first flower"
h2[13,2] <- "2018"
h2[13,3] <- herit.gam$h2.obs
h2
###########################

##Date to bolt##14
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(df$DtB.Oday.2018)
flr18 <- flr.18[!is.na(flr.18$DtB.Oday.2018),]
f1g <- fitdistr(flr18$DtB.Oday.2018, "normal")
f2g <- fitdistr(flr18$DtB.Oday.2018, "poisson")
#f3g <- fitdistr(flr18$DtB.Oday.2018, "negative binomial")
f4g <- fitdistr(flr18$DtB.Oday.2018, "Gamma")
f1g
AIC(f1g, f2g, f4g)
plot(f1g)
plot(f2g)
##Roughly normal? OR poisson??##
dtb18.mod<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = flr18,
							family = gaussian)
dtb18.out <-	summary(dtb18.mod)
dtb18.out
hist(residuals(dtb18.mod))

##look at fitted values vs original data##
plot(fitted(dtb18.mod),flr18$DtB.Oday.2018)
##check correlation of data##
cor(fitted(dtb18.mod),flr18$DtB.Oday.2018)
##other models:##
dtb18.mod2<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						family = poisson(link = "log"))
dtb18.out2 <-	summary(dtb18.mod2)
dtb18.out2
hist(residuals(dtb18.mod2))

dtb18.mod3<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link = "log"))
dtb18.out3 <-	summary(dtb18.mod3)
dtb18.out3
hist(residuals(dtb18.mod3))

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
herit2 <- QGparams(mu = mu, var.a = va, var.p = vp,  model = "Gaussian")
herit2
##add to table
h2[14,1] <- "Date to bolt"
h2[14,2] <- "2018"
h2[14,3] <- herit2$h2.obs
h2
###########################

##Date to Fruit 2018##15
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$Date.to.Fruit.Oday.2018)
flr18 <- flr.18[!is.na(flr.18$Date.to.Fruit.Oday.2018),]
f1g <- fitdist(flr18$Date.to.Fruit.Oday.2018, "norm")
f2g <- fitdist(flr18$Date.to.Fruit.Oday.2018, "pois")
f3g <- fitdist(flr18$Date.to.Fruit.Oday.2018, "gamma")
plot(f1g)
plot(f2g)
plot(f3g)
f1 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "normal")
f2 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "poisson")
f3 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "gamma")
AIC(f1, f2, f3)

dtfr18.mod<- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link=log))

summary(dtfr18.mod)
hist(residuals(dtfr18.mod))
dtfr18.out <-	summary(dtfr18.mod)
dtfr18.out
hist(residuals(dtfr18.mod))

##pull variables from model##
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
#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "gamma")
g
##Shape parameter (k)##
k <- g$estimate[1]
##fitdistr gives rate parameter (beta), which is inverse of scale(theta)##
##using theta because scale term seems more common##
theta <- 1/(g$estimate[2])
e <- exp(1)
##define functions for QGparams##
inv.link <- function(x){exp(x)}
var.func <- function(x){k*theta^2}
d.inv.link <- function(x){e^(x)}
custom.functions <- list(inv.link =inv.link, var.func=var.func,
								 d.inv.link = d.inv.link)
herit.gam <- QGparams(mu = mu, var.a = va, var.p = vp, 
							 custom.model = custom.functions) 
herit.gam

##add to table
h2[15,1] <- "Date to Fruit"
h2[15,2] <- "2018"
h2[15,3] <- herit.gam$h2.obs
h2
###########################

##No. flowers 2018##16
############################
View(df)
flr.18 <- filter(df, Total.Flowers.2018 >=1)
#View(flr.18)
flr18 <- df[!is.na(df$Total.Flowers.2018),]
hist(flr18$Total.Flowers.2018)
f1g <- fitdist(flr18$Total.Flowers.2018, "norm")
f2g <- fitdist(flr18$Total.Flowers.2018, "pois")
f3g <- fitdist(flr18$Total.Flowers.2018, "nbinom")
f4g <- fitdist(flr18$Total.Flowers.2018, "gamma")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
f1 <- fitdistr(flr18$Total.Flowers.2018, "normal")
f2 <- fitdistr(flr18$Total.Flowers.2018, "poisson")
f3 <- fitdistr(flr18$Total.Flowers.2018, "negative binomial")
f3 #theta 0.31204072
##with zeros removed: theta = 1.9207413
AIC(f1, f2, f3)

n.flr.mod <- glmer(Total.Flowers.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr.18,
						 family=neg.bin(theta = 1.9207413))
hist(residuals(n.flr.mod))

vars <- as.data.frame(VarCorr(n.flr.mod))[, c('grp','vcov')]
intercept <- fixef(n.flr.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  1.9207413
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

h2[16,1] <- "Number of Flowers"
h2[16,2] <- "2018"
h2[16,3] <- herit2$h2.obs
h2
############################

##No. fruit 2018##17
############################
flr.18 <- filter(df, No.Fruit.2018 >=1)
flr18 <- df[!is.na(df$No.Fruit.2018),]
summary(flr18$No.Fruit.2018)
hist(flr18$No.Fruit.2018)
##distributions with zeros in##
f1g <- fitdistr(flr18$No.Fruit.2018, "normal")
f2g <- fitdistr(flr18$No.Fruit.2018, "poisson")
f3g <- fitdistr(flr18$No.Fruit.2018, "negative binomial")
##f4g <- fitdistr(flr18$No.Fruit.2018, "gamma") doesn't work##
AIC(f1g,f2g,f3g)
f3g #theta = 0.22335483

##subset out zeros##
f1g <- fitdistr(flr.18$No.Fruit.2018, "normal")
f2g <- fitdistr(flr.18$No.Fruit.2018, "poisson")
f3g <- fitdistr(flr.18$No.Fruit.2018, "negative binomial")
f4g <- fitdistr(flr.18$No.Fruit.2018, "gamma")
AIC(f1g,f2g,f3g, f4g)
f3g #theta = 2.5503853
f4g 

n.frt.mod <- glmer(No.Fruit.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr.18,
						 family=neg.bin(theta = 2.5503853))
hist(residuals(n.frt.mod))

vars <- as.data.frame(VarCorr(n.frt.mod))[, c('grp','vcov')]
intercept <- fixef(n.frt.mod)['(Intercept)']
##verify data loaded in to objects##
vars
##verify data loaded in to objects##
intercept

##View latent-scale values region mean##
##values are for variables from model, not yet converted to observation scale##
theta <-  2.5503853
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

##NOTE--With Negative binomial and theta from Mason--residuals ok--model failed to converge##

h2[17,1] <- "Number of Fruit"
h2[17,2] <- "2018"
h2[17,3] <- herit2$h2.obs
h2
###########################



##Seedmass 2018##18
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$sm.3)
flr18 <- df[!is.na(df$sm.3),]
hist(flr18$sm.3)
##distributions with zeros in##
f1g <- fitdistr(flr18$sm.3, "normal")
f2g <- fitdistr(flr18$sm.3, "poisson")
f3g <- fitdistr(flr18$sm.3, "negative binomial")
f4g <- fitdistr(flr18$sm.3, "gamma") 
AIC(f1g,f2g,f3g, f4g)
f3g #theta = 0.22335483
f4g #theta = 1.67350953

##subset out zeros##
f1g <- fitdistr(flr.18$sm.3, "normal")
f2g <- fitdistr(flr.18$sm.3, "poisson")
f3g <- fitdistr(flr.18$sm.3, "negative binomial")
f4g <- fitdistr(flr.18$sm.3, "gamma")
AIC(f1g,f2g,f3g, f4g)
f3g #theta = 2.5503853
f4g 

f1g <- fitdist(flr18$sm.3, "norm")
f2g <- fitdist(flr18$sm.3, "pois")
f3g <- fitdist(flr18$sm.3, "nbinom")
f4g <- fitdist(flr18$sm.3, "gamma") 
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)

seeds18.mod<- glmer(sm.3~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = flr.18,
						  family = negative.binomial(2.5503853))
seeds18.out <-	summary(seeds18.mod)
seeds18.out
hist(residuals(seeds18.mod))
vars <- as.data.frame(VarCorr(seeds18.mod))[, c('grp','vcov')]
intercept <- fixef(seeds18.mod)['(Intercept)']
vars
intercept
#latent region mean#
mu <- intercept
theta <- 2.5503853
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
h2[18,1] <- "seedmass"
h2[18,2] <- "2018"
h2[18,3] <- herit2$h2.obs
h2
###########################

str(h2)
write.csv(h2, "heritabilities_for_traits.csv")
