##load libraries##
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("effects")
#install.packages("tidyverse")
#install.packages("rptR")
#install.packages("emmeans")
#install.packages("QGglmm")
library(QGglmm)
#install.packages("fitdistrplus")
library(fitdistrplus)
library(lme4); library(lmerTest); library(effects)
library(tidyverse); library(rptR);library(emmeans)
library(boot)
library(simpleboot)
##NOTE: Distributions for variables in ASTER analysis taken from there##
	##shape parameters changed, however##
##cleaned data from aster model datasheet in git repo##
setwd(thelio-zeb)
df <- read.csv("cleaned_NV_CG_heritability_data.csv")
#View(df)

###########################################################
##INDIVIDUAL MODELS PER TRAIT TO CALCULATE HERITABILITIES##
###########################################################

##germination## trait no. 1
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
							(1 | Family.Unique) + (1 | Block.ID), data = germ,
						family=poisson(link=log))
##View outputs##
summary(germ.mod)
hist(residuals(germ.mod))

##model statement##
##separating out family.unique (Va) by region to compare##
germ.mod2 <- glmer(No.Days.to.Germ~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = germ,
						family=poisson(link=log))
##confidence intervals from bootstrap--cluster ONLY##
#C1<- confint.merMod(germ.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")
##View outputs##
germ.out <- summary(germ.mod2)
hist(residuals(germ.mod2))
fixef(germ.mod2)

##model statement##
##Centering population by region as well##
germ.mod3 <- glmer(No.Days.to.Germ~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = germ,
						family=poisson(link=log))
##View outputs##
summary(germ.mod3)
hist(residuals(germ.mod3))
A1 <- AIC(germ.mod, germ.mod2, germ.mod3)
##conclusion: separating populations by region mean is not informative##
A1

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(germ.mod2))
vars
print(VarCorr(germ.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(germ.mod2)['(Intercept)']*(nrow(dplyr::filter(germ, germ$Region =="GL_alvar"))/nrow(germ))
gla
##Prairie region mean##
pra <-fixef(germ.mod2)['RegionPrairie']*(nrow(dplyr::filter(germ, germ$Region == "Prairie"))/nrow(germ))
##MB_alvar region mean##
mba <- fixef(germ.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(germ, germ$Region =="MB_alvar"))/nrow(germ))

##look at model values to make sure mu's makes sense##
fixef(germ.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(germ.mod2)[1] * model.matrix(germ.mod2)[, 1] +  fixef(germ.mod2)[2] * model.matrix(germ.mod2)[, 2] +
	fixef(germ.mod2)[3] * model.matrix(germ.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##
VarF
##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(germ.mod2)[1] * model.matrix(germ.mod2)[, 1]) ##gives 0 for GLA
##Manitoba alvar variance##
varMBA <- var(fixef(germ.mod2)[2] * model.matrix(germ.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(germ.mod2)[3] * model.matrix(germ.mod2)[, 3])
varGLA
varMBA
varPRA

##############Narrow-sense Heritabilities###########
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, model = "Poisson.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, model = "Poisson.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, model = "Poisson.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##Create a table to compile heritabilities## 
col.classes = c("character", "character", "numeric", "numeric", "character", "numeric", "numeric")
col.names = c("Trait", "Year", "Heritability", "Evolvability", "Region", "CVa","Va")
h2 <-read.table(text = "",colClasses = col.classes, col.names =col.names)
h2[1,1] <- "Germination"
h2[1,2] <- "2015"
h2[1,3] <- herit.gla$h2.obs
h2[1,4] <- ev.gla
h2[1,5] <- "GLA"
h2[1,6] <- cv.gla
h2[1,7] <- herit.gla$var.a.obs

h2[2,1] <- "Germination"
h2[2,2] <- "2015"
h2[2,3] <- herit.mba$h2.obs
h2[2,4] <- ev.mba 
h2[2,5] <-	"MBA"
h2[2,6] <- cv.mba
h2[2,7] <- herit.mba$var.a.obs

h2[3,1] <- "Germination"
h2[3,2] <- "2015"
h2[3,3] <- herit.pra$h2.obs
h2[3,4] <- ev.pra
h2[3,5] <- "PRA"
h2[3,6] <- cv.pra
h2[3,7] <- herit.pra$var.a.obs
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
f1 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "lognormal")
f2 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "poisson")
f3 <- fitdistr(TrueLeaf$No.Days.to.TrueLeaf, "Gamma")
AIC(f1,f2,f3)
##poisson visually looks best as description of data##
##Model Statement##
TrueLeaf.mod <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							 	(1 | Family.Unique) + (1 | Block.ID), data = TrueLeaf,
							 family=poisson(link=log))
##output##
summary(TrueLeaf.mod)
hist(residuals(TrueLeaf.mod))

TrueLeaf.mod2 <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							 	(Region| Family.Unique) + (1 | Block.ID), data = TrueLeaf,
							 family=poisson(link=log))
##output##
summary(TrueLeaf.mod2)
hist(residuals(TrueLeaf.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C2 <- confint.merMod(TrueLeaf.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

TrueLeaf.mod3 <- glmer(No.Days.to.TrueLeaf~Region + (Region | Population) + 
							 	(Region | Family.Unique) + (1 | Block.ID), data = TrueLeaf,
							 family=poisson(link=log))
##output##
summary(TrueLeaf.mod3)
hist(residuals(TrueLeaf.mod3))

A2 <-AIC(TrueLeaf.mod, TrueLeaf.mod2, TrueLeaf.mod3)
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(TrueLeaf.mod2))
vars
print(VarCorr(TrueLeaf.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(TrueLeaf.mod2)['(Intercept)']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region =="GL_alvar"))/nrow(TrueLeaf))
gla
##Prairie region mean##
pra <-fixef(TrueLeaf.mod2)['RegionPrairie']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region == "Prairie"))/nrow(TrueLeaf))
##MB_alvar region mean##
mba <- fixef(TrueLeaf.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region =="MB_alvar"))/nrow(TrueLeaf))

##look at model values to make sure mu's makes sense##
fixef(TrueLeaf.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(TrueLeaf.mod2)[1] * model.matrix(TrueLeaf.mod2)[, 1] +  fixef(TrueLeaf.mod2)[2] * model.matrix(TrueLeaf.mod2)[, 2] +
	fixef(TrueLeaf.mod2)[3] * model.matrix(TrueLeaf.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(TrueLeaf.mod2)[1] * model.matrix(TrueLeaf.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(TrueLeaf.mod2)[2] * model.matrix(TrueLeaf.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(TrueLeaf.mod2)[3] * model.matrix(TrueLeaf.mod2)[, 3])

##############Narrow-sense Heritabilities###########
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, model = "Poisson.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, model = "Poisson.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, model = "Poisson.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##Create a table to compile heritabilities## 
h2[4,1] <- "TrueLeaf"
h2[4,2] <- "2015"
h2[4,3] <- herit.gla$h2.obs
h2[4,4] <- ev.gla
h2[4,5] <- "GLA"
h2[4,6] <- cv.gla
h2[4,7] <- herit.gla$var.a.obs

h2[5,1] <- "TrueLeaf"
h2[5,2] <- "2015"
h2[5,3] <- herit.mba$h2.obs
h2[5,4] <- ev.mba 
h2[5,5] <-	"MBA"
h2[5,6] <- cv.mba
h2[5,7] <- herit.mba$var.a.obs

h2[6,1] <- "TrueLeaf"
h2[6,2] <- "2015"
h2[6,3] <- herit.pra$h2.obs
h2[6,4] <- ev.pra
h2[6,5] <- "PRA"
h2[6,6] <- cv.pra
h2[6,7] <- herit.pra$var.a.obs
############################################

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
##NOTE: using .gam instead of .mod for gamma-distributed models##
dtff.gam <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.gam)
hist(residuals(dtff.gam))
dtff.gam

##Model Statement##
##NOTE: using .gam instead of .mod for gamma-distributed models##
dtff.gam2 <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.gam2)
hist(residuals(dtff.gam2))
##confidence intervals from bootstrap--cluster ONLY##
#C3 <- confint.merMod(dtff.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

##Model Statement##
##NOTE: using .gam instead of .mod for gamma-distributed models##
dtff.gam3 <- glmer(no.Planting.to.DTFF~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.gam3)
hist(residuals(dtff.gam3))
A3 <- AIC(dtff.gam, dtff.gam2,dtff.gam3)

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.gam2))
vars
print(VarCorr(dtff.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff.gam2)['(Intercept)']*(nrow(dplyr::filter(FLR, FLR$Region =="GL_alvar"))/nrow(FLR))
gla
##Prairie region mean##
pra <-fixef(dtff.gam2)['RegionPrairie']*(nrow(dplyr::filter(FLR, FLR$Region == "Prairie"))/nrow(FLR))
##MB_alvar region mean##
mba <- fixef(dtff.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(FLR, FLR$Region =="MB_alvar"))/nrow(FLR))

##look at model values to make sure mu's makes sense##
fixef(dtff.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtff.gam2)[1] * model.matrix(dtff.gam2)[, 1] +  fixef(dtff.gam2)[2] * model.matrix(dtff.gam2)[, 2] +
	fixef(dtff.gam2)[3] * model.matrix(dtff.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtff.gam2)[1] * model.matrix(dtff.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtff.gam2)[2] * model.matrix(dtff.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtff.gam2)[3] * model.matrix(dtff.gam2)[, 3])

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
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[7,1] <- "DTFF"
h2[7,2] <- "2016"
h2[7,3] <- herit.gla$h2.obs
h2[7,4] <- ev.gla
h2[7,5] <- "GLA"
h2[7,6] <- cv.gla
h2[7,7] <- herit.gla$var.a.obs

h2[8,1] <- "DTFF"
h2[8,2] <- "2016"
h2[8,3] <- herit.mba$h2.obs
h2[8,4] <- ev.mba 
h2[8,5] <-	"MBA"
h2[8,6] <- cv.mba
h2[8,7] <- herit.mba$var.a.obs

h2[9,1] <- "DTFF"
h2[9,2] <- "2016"
h2[9,3] <- herit.pra$h2.obs
h2[9,4] <- ev.pra
h2[9,5] <- "PRA"
h2[9,6] <- cv.pra
h2[9,7] <- herit.pra$var.a.obs
############################################

##No. flowers 2016##4
##########################################################
nfl<-df[!is.na(df$No.Flowers.2016),]
f1 <- fitdistr(nfl$No.Flowers.2016, "normal")
f2 <- fitdistr(nfl$No.Flowers.2016, "Poisson")
f3 <- fitdistr(nfl$No.Flowers.2016, "gamma")
f4 <- fitdistr(nfl$No.Flowers.2016, "negative binomial")
hist(log(nfl$No.Flowers.2016))
hist(nfl$No.Flowers.2016)

AIC(f1,f2,f3,f4)
f4 #theta 1.27066
theta <- f4$estimate[1]
theta
n.flr.mod <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = nfl,
						 family=neg.bin(theta = theta))

n.flr.mod2 <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = nfl,
						 family=neg.bin(theta = theta))
##confidence intervals from bootstrap--cluster ONLY##
#C4<- confint.merMod(n.flr.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.flr.mod3 <- glmer(No.Flowers.2016~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = nfl,
						 family=neg.bin(theta = theta))
A4 <- AIC(n.flr.mod, n.flr.mod2, n.flr.mod3)

####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.mod2))
vars
print(VarCorr(n.flr.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.mod2)['(Intercept)']*(nrow(dplyr::filter(nfl, nfl$Region =="GL_alvar"))/nrow(nfl))
gla
##Prairie region mean##
pra <-fixef(n.flr.mod2)['RegionPrairie']*(nrow(dplyr::filter(nfl, nfl$Region == "Prairie"))/nrow(nfl))
##MB_alvar region mean##
mba <- fixef(n.flr.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(nfl, nfl$Region =="MB_alvar"))/nrow(nfl))

##look at model values to make sure mu's makes sense##
fixef(n.flr.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1] +  fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2] +
	fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[10,1] <- "Number of Flowers"
h2[10,2] <- "2016"
h2[10,3] <- herit.gla$h2.obs
h2[10,4] <- ev.gla
h2[10,5] <- "GLA"
h2[10,6] <- cv.gla
h2[10,7] <- herit.gla$var.a.obs

h2[11,1] <- "Number of Flowers"
h2[11,2] <- "2016"
h2[11,3] <- herit.mba$h2.obs
h2[11,4] <- ev.mba 
h2[11,5] <-	"MBA"
h2[11,6] <- cv.mba
h2[11,7] <- herit.mba$var.a.obs

h2[12,1] <- "Number of Flowers"
h2[12,2] <- "2016"
h2[12,3] <- herit.pra$h2.obs
h2[12,4] <- ev.pra
h2[12,5] <- "PRA"
h2[12,6] <- cv.pra
h2[12,7] <- herit.pra$var.a.obs
############################################

##Number of fruit 2016##5
########################################
hist(df$No.Fruit.2016)
nfr<-df[!is.na(df$No.Fruit.2016),]
descdist(nfr$No.Fruit.2016)
f1 <- fitdistr(nfr$No.Fruit.2016, "normal")
f2 <- fitdistr(nfr$No.Fruit.2016, "Poisson")
#f3 <- fitdistr(nfr$No.Fruit.2016, "gamma")#doesn't work with this data
f4 <- fitdistr(nfr$No.Fruit.2016, "negative binomial")
AIC(f1,f2,f4)
##Negative binomial best fit##
theta <- f4$estimate[1]
n.fruit.mod <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out <-	summary(n.fruit.mod)
hist(residuals(n.fruit.mod))

n.fruit.mod2 <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(Region | Family.Unique) + (1 | Block.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out2 <-	summary(n.fruit.mod2)
hist(residuals(n.fruit.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C5<- confint.merMod(n.fruit.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.fruit.mod3 <- glmer(No.Fruit.2016~Region + (Region | Population) + 
								(Region | Family.Unique) + (1 | Block.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out3 <-	summary(n.fruit.mod3)
hist(residuals(n.fruit.mod3))

A5 <- AIC(n.fruit.mod, n.fruit.mod2, n.fruit.mod3)

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.fruit.mod2))
vars
print(VarCorr(n.fruit.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.fruit.mod2)['(Intercept)']*(nrow(dplyr::filter(nfr, nfr$Region =="GL_alvar"))/nrow(nfr))
gla
##Prairie region mean##
pra <-fixef(n.fruit.mod2)['RegionPrairie']*(nrow(dplyr::filter(nfr, nfr$Region == "Prairie"))/nrow(nfr))
##MB_alvar region mean##
mba <- fixef(n.fruit.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(nfr, nfr$Region =="MB_alvar"))/nrow(nfr))

##look at model values to make sure mu's makes sense##
fixef(n.fruit.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.fruit.mod2)[1] * model.matrix(n.fruit.mod2)[, 1] +  fixef(n.fruit.mod2)[2] * model.matrix(n.fruit.mod2)[, 2] +
	fixef(n.fruit.mod2)[3] * model.matrix(n.fruit.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.fruit.mod2)[1] * model.matrix(n.fruit.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.fruit.mod2)[2] * model.matrix(n.fruit.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.fruit.mod2)[3] * model.matrix(n.fruit.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[13,1] <- "Number of Fruit"
h2[13,2] <- "2016"
h2[13,3] <- herit.gla$h2.obs
h2[13,4] <- ev.gla
h2[13,5] <- "GLA"
h2[13,6] <- cv.gla
h2[13,7] <- herit.gla$var.a.obs

h2[14,1] <- "Number of Fruit"
h2[14,2] <- "2016"
h2[14,3] <- herit.mba$h2.obs
h2[14,4] <- ev.mba 
h2[14,5] <-	"MBA"
h2[14,6] <- cv.mba
h2[14,7] <- herit.mba$var.a.obs

h2[15,1] <- "Number of Fruit"
h2[15,2] <- "2016"
h2[15,3] <- herit.pra$h2.obs
h2[15,4] <- ev.pra
h2[15,5] <- "PRA"
h2[15,6] <- cv.pra
h2[15,7] <- herit.pra$var.a.obs
############################################

##Seedmass 2016##6
#######################################
df1 <-df[!is.na(df$sm),]
f1 <- fitdistr(df1$sm, "normal")
f2 <- fitdistr(df1$sm, "Poisson")
f3 <- fitdistr(df1$sm, "gamma") #doesn't work with this data
f4 <- fitdistr(df1$sm, "negative binomial")
AIC(f1,f2,f4)
hist(df1$sm)
theta = f4$estimate[1]
n.seed.mod <- glmer(sm~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out <-	summary(n.seed.mod)
n.seed.out
hist(residuals(n.seed.mod))

n.seed.mod2 <- glmer(sm~Region + (1 | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out2 <-	summary(n.seed.mod2)
n.seed.out2
hist(residuals(n.seed.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C6<- confint.merMod(n.seed.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.seed.mod3 <- glmer(sm~Region + (Region | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out3 <-	summary(n.seed.mod3)
n.seed.out3
hist(residuals(n.seed.mod3))
A6 <- AIC(n.seed.mod, n.seed.mod2, n.seed.mod3)

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.seed.mod2))
vars
print(VarCorr(n.seed.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.seed.mod2)['(Intercept)']*(nrow(dplyr::filter(df1, df1$Region =="GL_alvar"))/nrow(df1))
gla
##Prairie region mean##
pra <-fixef(n.seed.mod2)['RegionPrairie']*(nrow(dplyr::filter(df1, df1$Region == "Prairie"))/nrow(df1))
##MB_alvar region mean##
mba <- fixef(n.seed.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(df1, df1$Region =="MB_alvar"))/nrow(df1))

##look at model values to make sure mu's makes sense##
fixef(n.seed.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.seed.mod2)[1] * model.matrix(n.seed.mod2)[, 1] +  fixef(n.seed.mod2)[2] * model.matrix(n.seed.mod2)[, 2] +
	fixef(n.seed.mod2)[3] * model.matrix(n.seed.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.seed.mod2)[1] * model.matrix(n.seed.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.seed.mod2)[2] * model.matrix(n.seed.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.seed.mod2)[3] * model.matrix(n.seed.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[16,1] <- "seedmass (mg)"
h2[16,2] <- "2016"
h2[16,3] <- herit.gla$h2.obs
h2[16,4] <- ev.gla
h2[16,5] <- "GLA"
h2[16,6] <- cv.gla
h2[16,7] <- herit.gla$var.a.obs

h2[17,1] <- "seedmass (mg)"
h2[17,2] <- "2016"
h2[17,3] <- herit.mba$h2.obs
h2[17,4] <- ev.mba 
h2[17,5] <-	"MBA"
h2[17,6] <- cv.mba
h2[17,7] <- herit.mba$var.a.obs

h2[18,1] <- "seedmass (mg)"
h2[18,2] <- "2016"
h2[18,3] <- herit.pra$h2.obs
h2[18,4] <- ev.pra
h2[18,5] <- "PRA"
h2[18,6] <- cv.pra
h2[18,7] <- herit.pra$var.a.obs
############################################

########2017 Season#############

##DTFF 2017##7
###########################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(df$DTFF.Ordinal.Day.2017)
flr.17$DTFF.Ordinal.Day.2017 <-as.numeric(flr.17$DTFF.Ordinal.Day.2017)
#make date number?
flr17 <- df[!is.na(df$DTFF.Ordinal.Day.2017),]
descdist(flr17$DTFF.Ordinal.Day.2017)
f1 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "normal")
f2 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "poisson")
f3 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "Gamma")
#f4 <- fitdistr(flr17$DTFF.Ordinal.Day.2017, "negative binomial")## error-fails
AIC(f1,f2,f3)

dtff.gam<- glmer(DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam))

dtff.gam2<- glmer(DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam2))
##confidence intervals from bootstrap--cluster ONLY##
#C7 <- confint.merMod(dtff.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtff.gam3 <- glmer(DTFF.Ordinal.Day.2017~Region + (Region | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam3))
A7 <- AIC(dtff.gam, dtff.gam2, dtt.gam3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.gam2))
vars
print(VarCorr(dtff.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff.gam2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtff.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtff.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(dtff.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtff.gam2)[1] * model.matrix(dtff.gam2)[, 1] +  fixef(dtff.gam2)[2] * model.matrix(dtff.gam2)[, 2] +
	fixef(dtff.gam2)[3] * model.matrix(dtff.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtff.gam2)[1] * model.matrix(dtff.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtff.gam2)[2] * model.matrix(dtff.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtff.gam2)[3] * model.matrix(dtff.gam2)[, 3])

##############Narrow-sense Heritabilities###########(Gamma)
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
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[19,1] <- "DTFF"
h2[19,2] <- "2017"
h2[19,3] <- herit.gla$h2.obs
h2[19,4] <- ev.gla
h2[19,5] <- "GLA"
h2[19,6] <- cv.gla
h2[19,7] <- herit.gla$var.a.obs

h2[20,1] <- "DTFF"
h2[20,2] <- "2017"
h2[20,3] <- herit.mba$h2.obs
h2[20,4] <- ev.mba 
h2[20,5] <-	"MBA"
h2[20,6] <- cv.mba
h2[20,7] <- herit.mba$var.a.obs

h2[21,1] <- "DTFF"
h2[21,2] <- "2017"
h2[21,3] <- herit.pra$h2.obs
h2[21,4] <- ev.pra
h2[21,5] <- "PRA"
h2[21,6] <- cv.pra
h2[21,7] <- herit.gla$var.a.obs
############################################

##Date to bolt 2017##8
###########################
flr17 <- df[!is.na(df$DtB.O.Day.2017),]
f1 <- fitdistr(flr17$DtB.O.Day.2017, "normal")
f2 <- fitdistr(flr17$DtB.O.Day.2017, "poisson")
f3 <- fitdistr(flr17$DtB.O.Day.2017, "Gamma")
#f4 <- fitdistr(flr17$DtB.O.Day.2017, "negative binomial")## error-fails
AIC(f1,f2,f3)
##model statement##
dtb.gam<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.gam))

dtb.gam2<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.gam2))
##confidence intervals from bootstrap--cluster ONLY##
#C8<- confint.merMod(dtb.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtb.gam3 <- glmer(DtB.O.Day.2017~Region + (Region | Population) + 
					 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.gam3))
A8 <- AIC(dtb.gam,dtb.gam2,dtb.gam3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtb.gam2))
vars
print(VarCorr(dtb.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtb.gam2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtb.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtb.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(dtb.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtb.gam2)[1] * model.matrix(dtb.gam2)[, 1] +  fixef(dtb.gam2)[2] * model.matrix(dtb.gam2)[, 2] +
	fixef(dtb.gam2)[3] * model.matrix(dtb.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtb.gam2)[1] * model.matrix(dtb.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtb.gam2)[2] * model.matrix(dtb.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtb.gam2)[3] * model.matrix(dtb.gam2)[, 3])

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
##############Narrow-sense Heritabilities###########(Gamma)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[22,1] <- "Date to Bolt"
h2[22,2] <- "2017"
h2[22,3] <- herit.gla$h2.obs
h2[22,4] <- ev.gla
h2[22,5] <- "GLA"
h2[22,6] <- cv.gla
h2[22,7] <- herit.gla$var.a.obs

h2[23,1] <- "Date to Bolt"
h2[23,2] <- "2017"
h2[23,3] <- herit.mba$h2.obs
h2[23,4] <- ev.mba 
h2[23,5] <-	"MBA"
h2[23,6] <- cv.mba
h2[23,7] <- herit.mba$var.a.obs

h2[24,1] <- "Date to Bolt"
h2[24,2] <- "2017"
h2[24,3] <- herit.pra$h2.obs
h2[24,4] <- ev.pra
h2[24,5] <- "PRA"
h2[24,6] <- cv.pra
h2[24,7] <- herit.pra$var.a.obs
############################################

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
dtfr.gam<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.gam))

dtfr.gam2<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.gam2))
##confidence intervals from bootstrap--cluster ONLY##
#C9 <- confint.merMod(dtfr.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtfr.gam3<- glmer(Fruit.O.Day.2017~Region + (Region | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.gam3))

A9 <- AIC(dtfr.gam, dtfr.gam2, dtfr.gam3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtfr.gam2))
vars
print(VarCorr(dtfr.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtfr.gam2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtfr.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtfr.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(dtfr.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtfr.gam2)[1] * model.matrix(dtfr.gam2)[, 1] +  fixef(dtfr.gam2)[2] * model.matrix(dtfr.gam2)[, 2] +
	fixef(dtfr.gam2)[3] * model.matrix(dtfr.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtfr.gam2)[1] * model.matrix(dtfr.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtfr.gam2)[2] * model.matrix(dtfr.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtfr.gam2)[3] * model.matrix(dtfr.gam2)[, 3])

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

##############Narrow-sense Heritabilities###########(Gamma)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[25,1] <- "Date to Fruit"
h2[25,2] <- "2017"
h2[25,3] <- herit.gla$h2.obs
h2[25,4] <- ev.gla
h2[25,5] <- "GLA"
h2[25,6] <- cv.gla
h2[25,7] <- herit.gla$var.a.obs

h2[26,1] <- "Date to Fruit"
h2[26,2] <- "2017"
h2[26,3] <- herit.mba$h2.obs
h2[26,4] <- ev.mba 
h2[26,5] <-	"MBA"
h2[26,6] <- cv.mba
h2[26,7] <- herit.mba$var.a.obs

h2[27,1] <- "Date to Fruit"
h2[27,2] <- "2017"
h2[27,3] <- herit.pra$h2.obs
h2[27,4] <- ev.pra
h2[27,5] <- "PRA"
h2[27,6] <- cv.pra
h2[27,7] <- herit.pra$var.a.obs
############################################

##No. flowers 2017##10
############################
flr.17 <- filter(df, Total.Flowers.2017 >=1)
#View(flr.17)
flr17 <- df[!is.na(df$Total.Flowers.2017),]
hist(flr17$Total.Flowers.2017)
fi <- fitdistr(flr17$Total.Flowers.2017, "normal")
fj <- fitdistr(flr17$Total.Flowers.2017, "poisson")
#fk <- fitdistr(flr17$Total.Flowers.2017, "Gamma")# Doesn't work
fl <- fitdistr(flr17$Total.Flowers.2017, "negative binomial")## error-fails
AIC(fi,fj,fl)
theta <- fl$estimate[1] #1.34227663
fitdistr(flr17$Total.Flowers.2017, "negative binomial")

n.flr.mod <- glmer(Total.Flowers.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.mod))

n.flr.mod2 <- glmer(Total.Flowers.2017~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C10 <- confint.merMod(n.flr.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.flr.mod3 <- glmer(Total.Flowers.2017~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.mod3))
A10 <- AIC(n.flr.mod, n.flr.mod2, n.flr.mod3)
######

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.mod2))
vars
print(VarCorr(n.flr.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(n.flr.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(n.flr.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(n.flr.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1] +  fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2] +
	fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[28,1] <- "Number of Flowers"
h2[28,2] <- "2017"
h2[28,3] <- herit.gla$h2.obs
h2[28,4] <- ev.gla
h2[28,5] <- "GLA"
h2[28,6] <- cv.gla
h2[28,7] <- herit.gla$var.a.obs

h2[29,1] <- "Number of Flowers"
h2[29,2] <- "2017"
h2[29,3] <- herit.mba$h2.obs
h2[29,4] <- ev.mba 
h2[29,5] <-	"MBA"
h2[29,6] <- cv.mba
h2[29,7] <- herit.mba$var.a.obs

h2[30,1] <- "Number of Flowers"
h2[30,2] <- "2017"
h2[30,3] <- herit.pra$h2.obs
h2[30,4] <- ev.pra
h2[30,5] <- "PRA"
h2[30,6] <- cv.pra
h2[30,7] <- herit.pra$var.a.obs
############################################

##No. fruit 2017##11
############################
flr17 <- df[!is.na(df$No.Fruit.2017),]
summary(flr17$No.Fruit.2017)
hist(flr17$No.Fruit.2017)
f1g <- fitdistr(flr17$No.Fruit.2017, "normal")
f2g <- fitdistr(flr17$No.Fruit.2017, "poisson")
f3g <- fitdistr(flr17$No.Fruit.2017, "negative binomial")
AIC(f1g,f2g,f3g)
f3g 
theta <- f3g$estimate[1]

n.frt.mod <- glmer(No.Fruit.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod))

n.frt.mod2 <- glmer(No.Fruit.2017~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C11 <- confint.merMod(n.frt.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.frt.mod3 <- glmer(No.Fruit.2017~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod3))

A11 <- AIC(n.frt.mod, n.frt.mod2, n.frt.mod3)
######

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.frt.mod2))
vars
print(VarCorr(n.frt.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.frt.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(n.frt.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(n.frt.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(n.frt.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.frt.mod2)[1] * model.matrix(n.frt.mod2)[, 1] +  fixef(n.frt.mod2)[2] * model.matrix(n.frt.mod2)[, 2] +
	fixef(n.frt.mod2)[3] * model.matrix(n.frt.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.frt.mod2)[1] * model.matrix(n.frt.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.frt.mod2)[2] * model.matrix(n.frt.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.frt.mod2)[3] * model.matrix(n.frt.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[31,1] <- "Number of Fruit"
h2[31,2] <- "2017"
h2[31,3] <- herit.gla$h2.obs
h2[31,4] <- ev.gla
h2[31,5] <- "GLA"
h2[31,6] <- cv.gla
h2[31,7] <- herit.gla$var.a.obs

h2[32,1] <- "Number of Fruit"
h2[32,2] <- "2017"
h2[32,3] <- herit.mba$h2.obs
h2[32,4] <- ev.mba 
h2[32,5] <-	"MBA"
h2[32,6] <- cv.mba
h2[32,7] <- herit.mba$var.a.obs

h2[33,1] <- "Number of Fruit"
h2[33,2] <- "2017"
h2[33,3] <- herit.pra$h2.obs
h2[33,4] <- ev.pra
h2[33,5] <- "PRA"
h2[33,6] <- cv.pra
h2[33,7] <- herit.pra$var.a.obs
############################################

##Seedmass 2017##12
########################################
flr17 <- df[!is.na(df$sm.2),]
df$sm.2
hist(flr17$sm.2)
f1g <- fitdistr(flr17$sm.2, "normal")
f2g <- fitdistr(flr17$sm.2, "poisson")
#f3g <- fitdistr(flr17$sm.2, "Gamma") # doesn't work
f4g <- fitdistr(flr17$sm.2, "negative binomial")
AIC(f1g, f2g, f4g)
f4g #
theta <- f4g$estimate[1]

##negative binomial rather than poisson, bc density on right tail high##
seed17.mod<- glmer(sm.2~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod)
hist(residuals(seed17.mod))

seed17.mod2<- glmer(sm.2~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod2)
hist(residuals(seed17.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C12 <- confint.merMod(seed17.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

seed17.mod3<- glmer(sm.2~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod3)
hist(residuals(seed17.mod))

A12 <- AIC(seed17.mod, seed17.mod2, seed17.mod3)
A12
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(seed17.mod2))
vars
print(VarCorr(seed17.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(seed17.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(seed17.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(seed17.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mu's makes sense##
fixef(seed17.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(seed17.mod2)[1] * model.matrix(seed17.mod2)[, 1] +  fixef(seed17.mod2)[2] * model.matrix(seed17.mod2)[, 2] +
	fixef(seed17.mod2)[3] * model.matrix(seed17.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(seed17.mod2)[1] * model.matrix(seed17.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(seed17.mod2)[2] * model.matrix(seed17.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(seed17.mod2)[3] * model.matrix(seed17.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[34,1] <- "Seedmass (mg)"
h2[34,2] <- "2017"
h2[34,3] <- herit.gla$h2.obs
h2[34,4] <- ev.gla
h2[34,5] <- "GLA"
h2[34,6] <- cv.gla
h2[34,7] <- herit.gla$var.a.obs

h2[35,1] <- "Seedmass (mg)"
h2[35,2] <- "2017"
h2[35,3] <- herit.mba$h2.obs
h2[35,4] <- ev.mba 
h2[35,5] <-	"MBA"
h2[35,6] <- cv.mba
h2[35,7] <- herit.mba$var.a.obs

h2[36,1] <- "Seedmass (mg)"
h2[36,2] <- "2017"
h2[36,3] <- herit.pra$h2.obs
h2[36,4] <- ev.pra
h2[36,5] <- "PRA"
h2[36,6] <- cv.pra
h2[36,7] <- herit.pra$var.a.obs
############################################

########2018 Season###########

##DTFF##13
###########################
hist(df$DTFF.18.Oday)
flr18 <- df[!is.na(df$DTFF.18.Oday),]
descdist(flr18$DTFF.18.Oday, boot = 100)
f1 <- fitdistr(flr18$DTFF.18.Oday, "normal")
f2 <- fitdistr(flr18$DTFF.18.Oday, "Poisson")
f3 <- fitdistr(flr18$DTFF.18.Oday, "gamma")
f4 <- fitdistr(flr18$DTFF.18.Oday, "lognormal")
f5 <- fitdistr(flr18$DTFF.18.Oday, "negative binomial")
AIC(f1,f2,f3,f4)

##model statement##
dtff18.gam<- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.gam)
hist(residuals(dtff18.gam))

dtff18.gam2 <- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.gam2)
hist(residuals(dtff18.gam2))
##confidence intervals from bootstrap--cluster ONLY##
#C13 <- confint.merMod(dtff18.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtff18.gam3<- glmer(DTFF.18.Oday~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.gam3)
hist(residuals(dtff18.gam3))
A13 <- AIC(dtff18.gam, dtff18.gam2, dtff18.gam3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff18.gam2))
vars
print(VarCorr(dtff18.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff18.gam2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtff18.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtff18.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(dtff18.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtff18.gam2)[1] * model.matrix(dtff18.gam2)[, 1] +  fixef(dtff18.gam2)[2] * model.matrix(dtff18.gam2)[, 2] +
	fixef(dtff18.gam2)[3] * model.matrix(dtff18.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtff18.gam2)[1] * model.matrix(dtff18.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtff18.gam2)[2] * model.matrix(dtff18.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtff18.gam2)[3] * model.matrix(dtff18.gam2)[, 3])

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g13 <- fitdistr(flr18$DTFF.18.Oday, "gamma")
g13
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
##############Narrow-sense Heritabilities###########(Gamma)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[37,1] <- "DTFF"
h2[37,2] <- "2018"
h2[37,3] <- herit.gla$h2.obs
h2[37,4] <- ev.gla
h2[37,5] <- "GLA"
h2[37,6] <- cv.gla
h2[37,7] <- herit.gla$var.a.obs

h2[38,1] <- "DTFF"
h2[38,2] <- "2018"
h2[38,3] <- herit.mba$h2.obs
h2[38,4] <- ev.mba 
h2[38,5] <-	"MBA"
h2[38,6] <- cv.mba
h2[38,7] <- herit.mba$var.a.obs

h2[39,1] <- "DTFF"
h2[39,2] <- "2018"
h2[39,3] <- herit.pra$h2.obs
h2[39,4] <- ev.pra
h2[39,5] <- "PRA"
h2[39,6] <- cv.pra
h2[39,7] <- herit.pra$var.a.obs
############################################

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
##Roughly normal? OR poisson??##
dtb18.mod<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						family = gaussian)
dtb18.out <-	summary(dtb18.mod)
hist(residuals(dtb18.mod))

dtb18.mod2<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						family = gaussian)
dtb18.out2 <-	summary(dtb18.mod2)
hist(residuals(dtb18.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C14 <- confint.merMod(dtb18.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtb18.mod3 <- glmer(DtB.Oday.2018~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						family = gaussian)
dtb18.out3 <-	summary(dtb18.mod3)
hist(residuals(dtb18.mod3))

A14 <- AIC(dtb18.mod, dtb18.mod2, dtb18.mod3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtb18.mod2))
vars
print(VarCorr(dtb18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtb18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtb18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtb18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(dtb18.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtb18.mod2)[1] * model.matrix(dtb18.mod2)[, 1] +  fixef(dtb18.mod2)[2] * model.matrix(dtb18.mod2)[, 2] +
	fixef(dtb18.mod2)[3] * model.matrix(dtb18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtb18.mod2)[1] * model.matrix(dtb18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtb18.mod2)[2] * model.matrix(dtb18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtb18.mod2)[3] * model.matrix(dtb18.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Gaussian)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, model = "gaussian")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, model = "gaussian")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, model = "gaussian")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[40,1] <- "Date to Bolt"
h2[40,2] <- "2018"
h2[40,3] <- herit.gla$h2.obs
h2[40,4] <- ev.gla
h2[40,5] <- "GLA"
h2[40,6] <- cv.gla
h2[40,7] <- herit.gla$var.a.obs

h2[41,1] <- "Date to Bolt"
h2[41,2] <- "2018"
h2[41,3] <- herit.mba$h2.obs
h2[41,4] <- ev.mba 
h2[41,5] <-	"MBA"
h2[41,6] <- cv.mba
h2[41,7] <- herit.mba$var.a.obs

h2[42,1] <- "Date to Bolt"
h2[42,2] <- "2018"
h2[42,3] <- herit.pra$h2.obs
h2[42,4] <- ev.pra
h2[42,5] <- "PRA"
h2[42,6] <- cv.pra
h2[42,7] <- herit.pra$var.a.obs
############################################

##Date to Fruit 2018##15
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$Date.to.Fruit.Oday.2018)
flr18 <- flr.18[!is.na(flr.18$Date.to.Fruit.Oday.2018),]
f1 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "normal")
f2 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "poisson")
f3 <- fitdistr(flr18$Date.to.Fruit.Oday.2018, "gamma")
AIC(f1, f2, f3)

dtfr18.gam<- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr18.gam)
hist(residuals(dtfr18.gam))
dtfr18.out <-	summary(dtfr18.gam)

dtfr18.gam2 <- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr18.gam2)
hist(residuals(dtfr18.gam2))
dtfr18.out <-	summary(dtfr18.gam2)
##confidence intervals from bootstrap--cluster ONLY##
#C15 <- confint.merMod(dtfr18.gam2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

dtfr18.gam3 <- glmer(Date.to.Fruit.Oday.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr18.gam3)
hist(residuals(dtfr18.gam3))
dtfr18.out <-	summary(dtfr18.gam3)

A15 <- AIC(dtfr18.gam, dtfr18.gam2, dtfr18.gam3)

#####
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtfr18.gam2))
vars
print(VarCorr(dtfr18.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtfr18.gam2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtfr18.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtfr18.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(dtfr18.gam2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtfr18.gam2)[1] * model.matrix(dtfr18.gam2)[, 1] +  fixef(dtfr18.gam2)[2] * model.matrix(dtfr18.gam2)[, 2] +
	fixef(dtfr18.gam2)[3] * model.matrix(dtfr18.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtfr18.gam2)[1] * model.matrix(dtfr18.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtfr18.gam2)[2] * model.matrix(dtfr18.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtfr18.gam2)[3] * model.matrix(dtfr18.gam2)[, 3])

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr18$Fruit.O.Day.2017, "gamma")
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

##############Narrow-sense Heritabilities###########(Gamma)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##Create a table to compile heritabilities## 
h2[43,1] <- "Date to Fruit"
h2[43,2] <- "2018"
h2[43,3] <- herit.gla$h2.obs
h2[43,4] <- ev.gla
h2[43,5] <- "GLA"
h2[43,6] <- cv.gla
h2[43,7] <- herit.gla$var.a.obs

h2[44,1] <- "Date to Fruit"
h2[44,2] <- "2018"
h2[44,3] <- herit.mba$h2.obs
h2[44,4] <- ev.mba 
h2[44,5] <-	"MBA"
h2[44,6] <- cv.mba
h2[44,7] <- herit.mba$var.a.obs

h2[45,1] <- "Date to Fruit"
h2[45,2] <- "2018"
h2[45,3] <- herit.pra$h2.obs
h2[45,4] <- ev.pra
h2[45,5] <- "PRA"
h2[45,6] <- cv.pra
h2[45,7] <- herit.pra$var.a.obs
############################################

##No. flowers 2018##16
############################
flr.18 <- filter(df, Total.Flowers.2018 >=1)
#View(flr.18)
flr18 <- df[!is.na(df$Total.Flowers.2018),]
hist(flr18$Total.Flowers.2018)
f1 <- fitdistr(flr18$Total.Flowers.2018, "normal")
f2 <- fitdistr(flr18$Total.Flowers.2018, "poisson")
f3 <- fitdistr(flr18$Total.Flowers.2018, "negative binomial")
f3 
theta <- f3$estimate[1]
AIC(f1, f2, f3)

n.flr.mod <- glmer(Total.Flowers.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.mod))
nflr.mod.out <- summary(n.flr.mod)

n.flr.mod2 <- glmer(Total.Flowers.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.mod2))
nflr.mod.out2<- summary(n.flr.mod2)
##confidence intervals from bootstrap--cluster ONLY##
#C15 <- confint.merMod(n.flr.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.flr.mod3 <- glmer(Total.Flowers.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.mod3))
nflr.mod.out3<- summary(n.flr.mod3)

A16 <- AIC(n.flr.mod, n.flr.mod2, n.flr.mod3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.mod2))
vars
print(VarCorr(n.flr.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(n.flr.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(n.flr.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(n.flr.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1] +  fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2] +
	fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.mod2)[1] * model.matrix(n.flr.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.mod2)[2] * model.matrix(n.flr.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.mod2)[3] * model.matrix(n.flr.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[46,1] <- "Number of Flowers"
h2[46,2] <- "2018"
h2[46,3] <- herit.gla$h2.obs
h2[46,4] <- ev.gla
h2[46,5] <- "GLA"
h2[46,6] <- cv.gla
h2[46,7] <- herit.gla$var.a.obs

h2[47,1] <- "Number of Flowers"
h2[47,2] <- "2018"
h2[47,3] <- herit.mba$h2.obs
h2[47,4] <- ev.mba 
h2[47,5] <-	"MBA"
h2[47,6] <- cv.mba
h2[47,7] <- herit.mba$var.a.obs

h2[48,1] <- "Number of Flowers"
h2[48,2] <- "2018"
h2[48,3] <- herit.pra$h2.obs
h2[48,4] <- ev.pra
h2[48,5] <- "PRA"
h2[48,6] <- cv.pra
h2[48,7] <- herit.pra$var.a.obs
############################################

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
f3g 
theta <- f3g$estimate[1]

n.frt.mod <- glmer(No.Fruit.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod))
frt.out <- summary(n.frt.mod)

n.frt.mod2 <- glmer(No.Fruit.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod2))
frt.out2 <- summary(n.frt.mod2)
##confidence intervals from bootstrap--cluster ONLY##
#C17 <- confint.merMod(n.frt.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

n.frt.mod3 <- glmer(No.Fruit.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.frt.mod3))
frt.out3 <- summary(n.frt.mod3)

A17 <- AIC(n.frt.mod, n.frt.mod2, n.frt.mod3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.frt.mod2))
vars
print(VarCorr(n.frt.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.frt.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(n.frt.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(n.frt.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(n.frt.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.frt.mod2)[1] * model.matrix(n.frt.mod2)[, 1] +  fixef(n.frt.mod2)[2] * model.matrix(n.frt.mod2)[, 2] +
	fixef(n.frt.mod2)[3] * model.matrix(n.frt.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.frt.mod2)[1] * model.matrix(n.frt.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.frt.mod2)[2] * model.matrix(n.frt.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.frt.mod2)[3] * model.matrix(n.frt.mod2)[, 3])

##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[49,1] <- "Number of Fruit"
h2[49,2] <- "2018"
h2[49,3] <- herit.gla$h2.obs
h2[49,4] <- ev.gla
h2[49,5] <- "GLA"
h2[49,6] <- cv.gla
h2[49,7] <- herit.gla$var.a.obs

h2[50,1] <- "Number of Fruit"
h2[50,2] <- "2018"
h2[50,3] <- herit.mba$h2.obs
h2[50,4] <- ev.mba 
h2[50,5] <-	"MBA"
h2[50,6] <- cv.mba
h2[50,7] <- herit.mba$var.a.obs

h2[51,1] <- "Number of Fruit"
h2[51,2] <- "2018"
h2[51,3] <- herit.pra$h2.obs
h2[51,4] <- ev.pra
h2[51,5] <- "PRA"
h2[51,6] <- cv.pra
h2[51,7] <- herit.pra$var.a.obs
############################################

##Seedmass 2018##18
###########################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$sm.3)
flr18 <- df[!is.na(df$seedmass.2018.g.),]
flr18$sm.3<-flr18$seedmass.2018.g.*1000
flr18$sm.3<-as.integer(flr18$sm.3)
hist(flr18$sm.3)
##distributions with zeros in##
descdist(flr18$sm.3)
f1g <- fitdistr(flr18$sm.3, "normal")
f2g <- fitdistr(flr18$sm.3, "poisson")
f3g <- fitdistr(flr18$sm.3, "negative binomial")
f4g <- fitdistr(flr18$sm.3, "gamma") 
AIC(f1g,f2g,f3g)
f3g
theta <- f3g$estimate[1]

seeds18.mod<- glmer(sm.3~Region + (1 | Population) + 
						  	(1 | Family.Unique) + (1 | Block.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out <-	summary(seeds18.mod)
hist(residuals(seeds18.mod))

seeds18.mod2 <- glmer(sm.3~Region + (1 | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out2 <-	summary(seeds18.mod2)
seeds18.out2
hist(residuals(seeds18.mod2))
##confidence intervals from bootstrap--cluster ONLY##
#C18<- confint.merMod(seeds18.mod2, level = 0.95, method = 'boot', nsim=100, boot.type = "basic")

seeds18.mod3 <- glmer(sm.3~Region + (Region | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out3 <-	summary(seeds18.mod3)
hist(residuals(seeds18.mod3))
A18 <- AIC(seeds18.mod, seeds18.mod2, seeds18.mod3)
#####

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(seeds18.mod2))
vars
print(VarCorr(seeds18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[1,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[2,4]

##Family.Unique variance in Prairie##
va.pra <- vars[3,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(seeds18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(seeds18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
pra
##MB_alvar region mean##
mba <- fixef(seeds18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mu's makes sense##
fixef(seeds18.mod2)

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(seeds18.mod2)[1] * model.matrix(seeds18.mod2)[, 1] +  fixef(seeds18.mod2)[2] * model.matrix(seeds18.mod2)[, 2] +
	fixef(seeds18.mod2)[3] * model.matrix(seeds18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
View(model.matrix(seeds18.mod2))
varGLA <- var(fixef(seeds18.mod2)[1] * model.matrix(seeds18.mod2)[, 1])
varGLA
##Manitoba alvar variance##
varMBA <- var(fixef(seeds18.mod2)[2] * model.matrix(seeds18.mod2)[, 2])
varMBA
##prairie fixed effect##
varPRA <- var(fixef(seeds18.mod2)[3] * model.matrix(seeds18.mod2)[, 3])
varPRA
##############Narrow-sense Heritabilities###########(Negative binomial)
####GLA####
#values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla
####

####MBA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance##
##		Va/family	Pop	Block		Vfixed##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba
####

####PRA####
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va

##total variance (vp)##
##		Va/family	Pop	Block		Vfixed##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  theta = theta, model = "negbin.log")
herit.pra
####

##################Evolvability####################
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- herit.ma$var.a.obs/(herit.gla$mean.obs^2)
ev.mba <- herit.mba$var.a.obs/(herit.mba$mean.obs^2)
ev.pra <- herit.pra$var.a.obs/(herit.pra$mean.obs^2)
##################

##################Coefficient of additive variance
cv.gla <- 100* (sqrt(herit.gla$var.a.obs)/herit.gla$mean.obs)
cv.mba <- 100* (sqrt(herit.mba$var.a.obs)/herit.mba$mean.obs)
cv.pra <- 100* (sqrt(herit.pra$var.a.obs)/herit.pra$mean.obs)
##################

##table to compile heritabilities## 
h2[52,1] <- "Seedmass (mg)"
h2[52,2] <- "2018"
h2[52,3] <- herit.gla$h2.obs
h2[52,4] <- ev.gla
h2[52,5] <- "GLA"
h2[52,6] <- cv.gla
h2[52,7] <- herit.gla$var.a.obs

h2[53,1] <- "Seedmass (mg)"
h2[53,2] <- "2018"
h2[53,3] <- herit.mba$h2.obs
h2[53,4] <- ev.mba 
h2[53,5] <-	"MBA"
h2[53,6] <- cv.mba
h2[53,7] <- herit.mba$var.a.obs

h2[54,1] <- "Seedmass (mg)"
h2[54,2] <- "2018"
h2[54,3] <- herit.pra$h2.obs
h2[54,4] <- ev.pra
h2[54,5] <- "PRA"
h2[54,6] <- cv.pra
h2[54,7] <- herit.pra$var.a.obs
############################################

##EVALUATE THAT MODELS RAN BY VIEWING COMPARISON##
##IF OBJECT NOT FOUND, MODELS FAILED--ASSIGN TRAIT NA VALUE (FOR NOW)##
tryCatch(A1,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A2,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A3,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A4,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A5,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A6,error =function(e){cat("error", conditionMessage(e),"\n")}) #MODEL BREAK
tryCatch(A7,error =function(e){cat("error", conditionMessage(e),"\n")}) #MODEL BREAK
tryCatch(A8,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A9,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A10,error =function(e){cat("error", conditionMessage(e),"\n")}) #MODEL 3 SEEMS BETTER?
tryCatch(A11,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A12,error =function(e){cat("error", conditionMessage(e),"\n")}) #MODEL 3 SEEMS BETTER?
tryCatch(A13,error =function(e){cat("error", conditionMessage(e),"\n")}) #MODEL BREAK
tryCatch(A14,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A15,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A16,error =function(e){cat("error", conditionMessage(e),"\n")})
tryCatch(A17,error =function(e){cat("error", conditionMessage(e),"\n")})

h2[16,3:7]<-NA
h2[17,3:7]<-NA
h2[18,3:7]<-NA
#h2[19,3:4]<-NA
#h2[20,3:4]<-NA
#h2[21,3:4]<-NA
#h2[37,3:4]<-NA
#h2[38,3:4]<-NA
#h2[39,3:4]<-NA
write.csv(h2, "Region-level_heritabilities_evolvability_no_CI.csv", row.names = F)

#write.csv(C1, "germination_CI.csv", row.names = F)
#write.csv(C2, "trueleaf_CI.csv", row.names = F)
#write.csv(C3, "DTFF16_CI.csv", row.names =F)
#write.csv(C4, "No.Flr16_CI.csv", row.names =F)
#write.csv(C5, "No.Frt16_CI.csv", row.names =F)
#write.csv(C6, "No.Seeds16_CI.csv", row.names =F)
#write.csv(C7, "DTFF17_CI.csv", row.names =F)
#write.csv(C8, "DTB17_CI.csv", row.names =F)
#write.csv(C9, "DTFR17_CI.csv", row.names =F)
#write.csv(C10, "No.Flr17_CI.csv", row.names =F)
#write.csv(C11, "No.Frt17_CI.csv", row.names =F)
#write.csv(C12, "No.Seeds17_CI.csv", row.names =F)
#write.csv(C13, "DTFF18_CI.csv", row.names =F)
#write.csv(C14, "DTB18_CI.csv", row.names =F)
#write.csv(C15, "DTFR18_CI.csv", row.names =F)
#write.csv(C16, "No.Flr18_CI.csv", row.names =F)
#write.csv(C17, "No.Frt18_CI.csv", row.names =F)
#write.csv(C18, "No.Seeds18_CI.csv", row.names =F)

#rm(list=ls())

