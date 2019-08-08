##load libraries##
library(lme4); library(lmerTest); library(effects)
library(tidyverse); library(rptR);library(emmeans)
#install.packages("QGglmm")
library(QGglmm)
#install.packages("fitdistrplus")
library(fitdistrplus)
##NOTE: Distributions for variables in ASTER analysis taken from there##
	##shape parameters changed, however##
##cleaned data from aster model datasheet in git repo##
df <- read.csv("cleaned_NV_CG_experiment_data.csv")
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
							(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = germ,
						family=poisson(link=log))
##View outputs##
summary(germ.mod)
hist(residuals(germ.mod))

##model statement##
##separating out family.unique (Va) by region to compare##
germ.mod2 <- glmer(No.Days.to.Germ~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = germ,
						family=poisson(link=log))

##bootstrap methods##  #NOTE: Either not running, or taking forever##
#bt <- bootMer(germ.mod2,FUN = fixef, nsim = 2 )
#bt
#summary(bt)
#class(bt)
#boot(bt)
library(boot)
library(simpleboot)
#ci <- boot.ci(bt, conf = .95, type = 'basic')
#bt2 <- bootMer(germ.mod2,FUN = predict, nsim = 20 )
#conf <- confint.merMod(germ.mod2, level = 0.95, method = 'boot', nsim=1, boot.type = "basic")
#ci
##View outputs##
summary(germ.mod2)
hist(residuals(germ.mod2))
fixef(germ.mod2)

##model statement##
##Centering population by region as well##
germ.mod3 <- glmer(No.Days.to.Germ~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = germ,
						family=poisson(link=log))
##View outputs##
summary(germ.mod3)
hist(residuals(germ.mod3))
A1 <- AIC(germ.mod, germ.mod2, germ.mod3)
##conclusion: separating populations by region mean is not informative##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(germ.mod2))
vars
print(VarCorr(germ.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(germ.mod2)['(Intercept)']*(nrow(dplyr::filter(germ, germ$Region =="GL_alvar"))/nrow(germ))
gla
##Prairie region mean##
pra <-fixef(germ.mod2)['RegionPrairie']*(nrow(dplyr::filter(germ, germ$Region == "Prairie"))/nrow(germ))
##MB_alvar region mean##
mba <- fixef(germ.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(germ, germ$Region =="MB_alvar"))/nrow(germ))

##look at model values to make sure mus makes sense##
fixef(germ.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

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
##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(germ.mod2)
##Null model (No fixed effect) to get residual error variance##
m0 <- glmer(No.Days.to.Germ~1 + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = germ,
						 family=poisson(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, model = "Poisson.log")
herit.gla
##find Ve?##
h2.obs <- 0.101713

ve <- (va-(h2.obs*vp)/h2.obs)
ve
va
vp
?lme4::family
##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(germ.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, model = "Poisson.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(germ.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, model = "Poisson.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
col.classes = c("character", "character", "numeric", "numeric", "character")
col.names = c("Trait", "Year", "Heritability", "Evolvability", "Region")
h2 <-read.table(text = "",colClasses = col.classes, col.names =col.names)
h2[1,1] <- "Germination"
h2[1,2] <- "2015"
h2[1,3] <- herit.gla$h2.obs
h2[1,4] <- ev.gla
h2[1,5] <- "GLA"

h2[2,1] <- "Germination"
h2[2,2] <- "2015"
h2[2,3] <- herit.mba$h2.obs
h2[2,4] <- ev.mba 
h2[2,5] <-	"MBA"

h2[3,1] <- "Germination"
h2[3,2] <- "2015"
h2[3,3] <- herit.pra$h2.obs
h2[3,4] <- ev.pra
h2[3,5] <- "PRA"
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
							 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = TrueLeaf,
							 family=poisson(link=log))
##output##
summary(TrueLeaf.mod)
hist(residuals(TrueLeaf.mod))

TrueLeaf.mod2 <- glmer(No.Days.to.TrueLeaf~Region + (1 | Population) + 
							 	(Region| Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = TrueLeaf,
							 family=poisson(link=log))
##output##
summary(TrueLeaf.mod2)
hist(residuals(TrueLeaf.mod2))

TrueLeaf.mod3 <- glmer(No.Days.to.TrueLeaf~Region + (Region | Population) + 
							 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = TrueLeaf,
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
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(TrueLeaf.mod2)['(Intercept)']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region =="GL_alvar"))/nrow(TrueLeaf))
gla
##Prairie region mean##
pra <-fixef(TrueLeaf.mod2)['RegionPrairie']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region == "Prairie"))/nrow(TrueLeaf))
##MB_alvar region mean##
mba <- fixef(TrueLeaf.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(TrueLeaf, TrueLeaf$Region =="MB_alvar"))/nrow(TrueLeaf))

##look at model values to make sure mus makes sense##
fixef(TrueLeaf.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

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

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(TrueLeaf.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.TrueLeaf~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = TrueLeaf,
#				family=poisson(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp
ve <-vars[1,4]
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, model = "Poisson.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(TrueLeaf.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, model = "Poisson.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(TrueLeaf.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, model = "Poisson.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[4,1] <- "TrueLeaf"
h2[4,2] <- "2015"
h2[4,3] <- herit.gla$h2.obs
h2[4,4] <- ev.gla
h2[4,5] <- "GLA"

h2[5,1] <- "TrueLeaf"
h2[5,2] <- "2015"
h2[5,3] <- herit.mba$h2.obs
h2[5,4] <- ev.mba 
h2[5,5] <-	"MBA"

h2[6,1] <- "TrueLeaf"
h2[6,2] <- "2015"
h2[6,3] <- herit.pra$h2.obs
h2[6,4] <- ev.pra
h2[6,5] <- "PRA"
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
dtff.16.gam <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.16.gam)
hist(residuals(dtff.16.gam))
dtff.16.gam

##Model Statement##
##NOTE: using .gam instead of .mod for gamma-distributed models##
dtff.16.gam2 <- glmer(no.Planting.to.DTFF~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.16.gam2)
hist(residuals(dtff.16.gam2))

##Model Statement##
##NOTE: using .gam instead of .mod for gamma-distributed models##
dtff.16.gam3 <- glmer(no.Planting.to.DTFF~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = FLR,
						family = Gamma(link=log))
##model output##
summary(dtff.16.gam3)
hist(residuals(dtff.16.gam3))
A3 <- AIC(dtff.16.gam, dtff.16.gam2,dtff.16.gam3)
######
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.mod2))
vars
print(VarCorr(dtff.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff.mod2)['(Intercept)']*(nrow(dplyr::filter(FLR, FLR$Region =="GL_alvar"))/nrow(FLR))
gla
##Prairie region mean##
pra <-fixef(dtff.mod2)['RegionPrairie']*(nrow(dplyr::filter(FLR, FLR$Region == "Prairie"))/nrow(FLR))
##MB_alvar region mean##
mba <- fixef(dtff.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(FLR, FLR$Region =="MB_alvar"))/nrow(FLR))

##look at model values to make sure mus makes sense##
fixef(dtff.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtff.mod2)[1] * model.matrix(dtff.mod2)[, 1] +  fixef(dtff.mod2)[2] * model.matrix(dtff.mod2)[, 2] +
	fixef(dtff.mod2)[3] * model.matrix(dtff.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtff.mod2)[1] * model.matrix(dtff.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtff.mod2)[2] * model.matrix(dtff.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtff.mod2)[3] * model.matrix(dtff.mod2)[, 3])

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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.model = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.model = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.model = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[7,1] <- "DTFF"
h2[7,2] <- "2016"
h2[7,3] <- herit.gla$h2.obs
h2[7,4] <- ev.gla
h2[7,5] <- "GLA"

h2[8,1] <- "DTFF"
h2[8,2] <- "2016"
h2[8,3] <- herit.mba$h2.obs
h2[8,4] <- ev.mba 
h2[8,5] <-	"MBA"

h2[9,1] <- "DTFF"
h2[9,2] <- "2016"
h2[9,3] <- herit.pra$h2.obs
h2[9,4] <- ev.pra
h2[9,5] <- "PRA"
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
n.flr.16.mod <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfl,
						 family=neg.bin(theta = theta))

n.flr.16.mod2 <- glmer(No.Flowers.2016~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfl,
						 family=neg.bin(theta = theta))

n.flr.16.mod3 <- glmer(No.Flowers.2016~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfl,
						 family=neg.bin(theta = theta))
A4 <- AIC(n.flr.16.mod, n.flr.16.mod2, n.flr.16.mod3)

####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.16.mod2))
vars
print(VarCorr(n.flr.16.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.16.mod2)['(Intercept)']*(nrow(dplyr::filter(nfl, nfl$Region =="GL_alvar"))/nrow(nfl))
gla
##Prairie region mean##
pra <-fixef(n.flr.16.mod2)['RegionPrairie']*(nrow(dplyr::filter(nfl, nfl$Region == "Prairie"))/nrow(nfl))
##MB_alvar region mean##
mba <- fixef(n.flr.16.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(nfl, nfl$Region =="MB_alvar"))/nrow(nfl))

##look at model values to make sure mus makes sense##
fixef(n.flr.16.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.16.mod2)[1] * model.matrix(n.flr.16.mod2)[, 1] +  fixef(n.flr.16.mod2)[2] * model.matrix(n.flr.16.mod2)[, 2] +
	fixef(n.flr.16.mod2)[3] * model.matrix(n.flr.16.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.16.mod2)[1] * model.matrix(n.flr.16.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.16.mod2)[2] * model.matrix(n.flr.16.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.16.mod2)[3] * model.matrix(n.flr.16.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.16.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.nfl~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfl,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.16.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.16.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[10,1] <- "Number of Flowers"
h2[10,2] <- "2016"
h2[10,3] <- herit.gla$h2.obs
h2[10,4] <- ev.gla
h2[10,5] <- "GLA"

h2[11,1] <- "Number of Flowers"
h2[11,2] <- "2016"
h2[11,3] <- herit.mba$h2.obs
h2[11,4] <- ev.mba 
h2[11,5] <-	"MBA"

h2[12,1] <- "Number of Flowers"
h2[12,2] <- "2016"
h2[12,3] <- herit.pra$h2.obs
h2[12,4] <- ev.pra
h2[12,5] <- "PRA"
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
n.fruit.16.mod <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out <-	summary(n.fruit.16.mod)
hist(residuals(n.fruit.16.mod))

n.fruit.16.mod2 <- glmer(No.Fruit.2016~Region + (1 | Population) + 
								(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out2 <-	summary(n.fruit.16.mod2)
hist(residuals(n.fruit.16.mod2))

n.fruit.16.mod3 <- glmer(No.Fruit.2016~Region + (Region | Population) + 
								(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfr,
							family=neg.bin(theta = theta))
n.fruit.out3 <-	summary(n.fruit.16.mod3)
hist(residuals(n.fruit.16.mod3))

A5 <- AIC(n.fruit.16.mod, n.fruit.16.mod2, n.fruit.16.mod3)

####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.fruit.16.mod2))
vars
print(VarCorr(n.fruit.16.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.fruit.16.mod2)['(Intercept)']*(nrow(dplyr::filter(nfr, nfr$Region =="GL_alvar"))/nrow(nfr))
gla
##Prairie region mean##
pra <-fixef(n.fruit.16.mod2)['RegionPrairie']*(nrow(dplyr::filter(nfr, nfr$Region == "Prairie"))/nrow(nfr))
##MB_alvar region mean##
mba <- fixef(n.fruit.16.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(nfr, nfr$Region =="MB_alvar"))/nrow(nfr))

##look at model values to make sure mus makes sense##
fixef(n.fruit.16.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.fruit.16.mod2)[1] * model.matrix(n.fruit.16.mod2)[, 1] +  fixef(n.fruit.16.mod2)[2] * model.matrix(n.fruit.16.mod2)[, 2] +
	fixef(n.fruit.16.mod2)[3] * model.matrix(n.fruit.16.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.fruit.16.mod2)[1] * model.matrix(n.fruit.16.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.fruit.16.mod2)[2] * model.matrix(n.fruit.16.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.fruit.16.mod2)[3] * model.matrix(n.fruit.16.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.16.mod2)
##Null model (No fixed effect) to get residual error variance##
m0 <- glmer(No.Days.to.nfr~1 + (1 | Population) + 
					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = nfr,
				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.16.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.16.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[13,1] <- "Number of Fruit"
h2[13,2] <- "2016"
h2[13,3] <- herit.gla$h2.obs
h2[13,4] <- ev.gla
h2[13,5] <- "GLA"

h2[14,1] <- "Number of Fruit"
h2[14,2] <- "2016"
h2[14,3] <- herit.mba$h2.obs
h2[14,4] <- ev.mba 
h2[14,5] <-	"MBA"

h2[15,1] <- "Number of Fruit"
h2[15,2] <- "2016"
h2[15,3] <- herit.pra$h2.obs
h2[15,4] <- ev.pra
h2[15,5] <- "PRA"
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
						  	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out <-	summary(n.seed.mod)
n.seed.out
hist(residuals(n.seed.mod))

n.seed.mod2 <- glmer(sm~Region + (1 | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out2 <-	summary(n.seed.mod2)
n.seed.out2
hist(residuals(n.seed.mod2))

n.seed.mod3 <- glmer(sm~Region + (Region | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = df1,
						  family = neg.bin(theta = theta))
n.seed.out3 <-	summary(n.seed.mod3)
n.seed.out3
hist(residuals(n.seed.mod3))
A6 <- AIC(n.seed.mod, n.seed.mod2, n.seed.mod3)
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.seed.mod2))
vars
print(VarCorr(n.seed.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.seed.mod2)['(Intercept)']*(nrow(dplyr::filter(df1, df1$Region =="GL_alvar"))/nrow(df1))
gla
##Prairie region mean##
pra <-fixef(n.seed.mod2)['RegionPrairie']*(nrow(dplyr::filter(df1, df1$Region == "Prairie"))/nrow(df1))
##MB_alvar region mean##
mba <- fixef(n.seed.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(df1, df1$Region =="MB_alvar"))/nrow(df1))

##look at model values to make sure mus makes sense##
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

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.seed.mod2)
##Null model (No fixed effect) to get residual error variance##
m0 <- glmer(No.Days.to.df1~1 + (1 | Population) + 
					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = df1,
				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.seed.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.seed.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[16,1] <- "seedmass (mg)"
h2[16,2] <- "2016"
h2[16,3] <- herit.gla$h2.obs
h2[16,4] <- ev.gla
h2[16,5] <- "GLA"

h2[17,1] <- "seedmass (mg)"
h2[17,2] <- "2016"
h2[17,3] <- herit.mba$h2.obs
h2[17,4] <- ev.mba 
h2[17,5] <-	"MBA"

h2[18,1] <- "seedmass (mg)"
h2[18,2] <- "2016"
h2[18,3] <- herit.pra$h2.obs
h2[18,4] <- ev.pra
h2[18,5] <- "PRA"
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
					  	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam))

dtff.gam2<- glmer(DTFF.Ordinal.Day.2017~Region + (1 | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam2))

dtff.gam3 <- glmer(DTFF.Ordinal.Day.2017~Region + (Region | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = 107
					  family = Gamma(link=log), control = glmerControl(optCtrl = list(maxfun=10000000)))
##not converging--model fails?##
hist(residuals(dtff.gam3))
A7 <- AIC(dtff.gam, dtff.gam2, dtt.gam3)
#####
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff.gam2))
vars
print(VarCorr(dtff.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff.gam2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtff.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtff.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(dtff.gam2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[19,1] <- "DTFF"
h2[19,2] <- "2017"
h2[19,3] <- herit.gla$h2.obs
h2[19,4] <- ev.gla
h2[19,5] <- "GLA"

h2[20,1] <- "DTFF"
h2[20,2] <- "2017"
h2[20,3] <- herit.mba$h2.obs
h2[20,4] <- ev.mba 
h2[20,5] <-	"MBA"

h2[21,1] <- "DTFF"
h2[21,2] <- "2017"
h2[21,3] <- herit.pra$h2.obs
h2[21,4] <- ev.pra
h2[21,5] <- "PRA"
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
dtb.17.gam<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.17.gam))

dtb.17.gam2<- glmer(DtB.O.Day.2017~Region + (1 | Population) + 
					 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.17.gam2))

dtb.17.gam3 <- glmer(DtB.O.Day.2017~Region + (Region | Population) + 
					 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					 ##first day = ~114
					 family = Gamma(link=log))
hist(residuals(dtb.17.gam3))
A8 <- AIC(dtb.17.gam,dtb.17.gam2,dtb.17.gam3)
#####
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtb.17.gam2))
vars
print(VarCorr(dtb.17.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtb.17.gam2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtb.17.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtb.17.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(dtb.17.gam2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtb.17.gam2)[1] * model.matrix(dtb.17.gam2)[, 1] +  fixef(dtb.17.gam2)[2] * model.matrix(dtb.17.gam2)[, 2] +
	fixef(dtb.17.gam2)[3] * model.matrix(dtb.17.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtb.17.gam2)[1] * model.matrix(dtb.17.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtb.17.gam2)[2] * model.matrix(dtb.17.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtb.17.gam2)[3] * model.matrix(dtb.17.gam2)[, 3])

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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtb.17.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtb.17.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[22,1] <- "Date to Bolt"
h2[22,2] <- "2017"
h2[22,3] <- herit.gla$h2.obs
h2[22,4] <- ev.gla
h2[22,5] <- "GLA"

h2[23,1] <- "Date to Bolt"
h2[23,2] <- "2017"
h2[23,3] <- herit.mba$h2.obs
h2[23,4] <- ev.mba 
h2[23,5] <-	"MBA"

h2[24,1] <- "Date to Bolt"
h2[24,2] <- "2017"
h2[24,3] <- herit.pra$h2.obs
h2[24,4] <- ev.pra
h2[24,5] <- "PRA"
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
dtfr.17.mod<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.17.mod))

dtfr.17.mod2<- glmer(Fruit.O.Day.2017~Region + (1 | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.17.mod2))

dtfr.17.mod3<- glmer(Fruit.O.Day.2017~Region + (Region | Population) + 
					  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
					  ##first day = ~121
					  family = Gamma(link=log))
hist(residuals(dtfr.17.mod3))

A9 <- AIC(dtfr.17.mod, dtfr.17.mod2, dtfr.17.mod3)

#####
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtfr.17.mod2))
vars
print(VarCorr(dtfr.17.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtfr.17.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(dtfr.17.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(dtfr.17.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(dtfr.17.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtfr.17.mod2)[1] * model.matrix(dtfr.17.mod2)[, 1] +  fixef(dtfr.17.mod2)[2] * model.matrix(dtfr.17.mod2)[, 2] +
	fixef(dtfr.17.mod2)[3] * model.matrix(dtfr.17.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtfr.17.mod2)[1] * model.matrix(dtfr.17.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtfr.17.mod2)[2] * model.matrix(dtfr.17.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtfr.17.mod2)[3] * model.matrix(dtfr.17.mod2)[, 3])

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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtfr.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtfr.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[25,1] <- "Date to Fruit"
h2[25,2] <- "2017"
h2[25,3] <- herit.gla$h2.obs
h2[25,4] <- ev.gla
h2[25,5] <- "GLA"

h2[26,1] <- "Date to Fruit"
h2[26,2] <- "2017"
h2[26,3] <- herit.mba$h2.obs
h2[26,4] <- ev.mba 
h2[26,5] <-	"MBA"

h2[27,1] <- "Date to Fruit"
h2[27,2] <- "2017"
h2[27,3] <- herit.pra$h2.obs
h2[27,4] <- ev.pra
h2[27,5] <- "PRA"
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

n.flr.17.mod <- glmer(Total.Flowers.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.17.mod))

n.flr.17.mod2 <- glmer(Total.Flowers.2017~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.17.mod2))

n.flr.17.mod3 <- glmer(Total.Flowers.2017~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = 1.34227663))
hist(residuals(n.flr.17.mod3))
A10 <- AIC(n.flr.17.mod, n.flr.17.mod2, n.flr.17.mod3)
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.17.mod2))
vars
print(VarCorr(n.flr.17.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.17.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(n.flr.17.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(n.flr.17.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(n.flr.17.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.17.mod2)[1] * model.matrix(n.flr.17.mod2)[, 1] +  fixef(n.flr.17.mod2)[2] * model.matrix(n.flr.17.mod2)[, 2] +
	fixef(n.flr.17.mod2)[3] * model.matrix(n.flr.17.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.17.mod2)[1] * model.matrix(n.flr.17.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.17.mod2)[2] * model.matrix(n.flr.17.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.17.mod2)[3] * model.matrix(n.flr.17.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.17.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr17~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[28,1] <- "Number of Flowers"
h2[28,2] <- "2017"
h2[28,3] <- herit.gla$h2.obs
h2[28,4] <- ev.gla
h2[28,5] <- "GLA"

h2[29,1] <- "Number of Flowers"
h2[29,2] <- "2017"
h2[29,3] <- herit.mba$h2.obs
h2[29,4] <- ev.mba 
h2[29,5] <-	"MBA"

h2[30,1] <- "Number of Flowers"
h2[30,2] <- "2017"
h2[30,3] <- herit.pra$h2.obs
h2[30,4] <- ev.pra
h2[30,5] <- "PRA"
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

n.fruit.17.mod <- glmer(No.Fruit.2017~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.17.mod))

n.fruit.17.mod2 <- glmer(No.Fruit.2017~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.17.mod2))

n.fruit.17.mod3 <- glmer(No.Fruit.2017~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.17.mod3))

A11 <- AIC(n.fruit.17.mod, n.fruit.17.mod2, n.fruit.17.mod3)
######
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.fruit.17.mod2))
vars
print(VarCorr(n.fruit.17.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.fruit.17.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(n.fruit.17.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(n.fruit.17.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(n.fruit.17.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.fruit.17.mod2)[1] * model.matrix(n.fruit.17.mod2)[, 1] +  fixef(n.fruit.17.mod2)[2] * model.matrix(n.fruit.17.mod2)[, 2] +
	fixef(n.fruit.17.mod2)[3] * model.matrix(n.fruit.17.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.fruit.17.mod2)[1] * model.matrix(n.fruit.17.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.fruit.17.mod2)[2] * model.matrix(n.fruit.17.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.fruit.17.mod2)[3] * model.matrix(n.fruit.17.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.17.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr17~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[31,1] <- "Number of Fruit"
h2[31,2] <- "2017"
h2[31,3] <- herit.gla$h2.obs
h2[31,4] <- ev.gla
h2[31,5] <- "GLA"

h2[32,1] <- "Number of Fruit"
h2[32,2] <- "2017"
h2[32,3] <- herit.mba$h2.obs
h2[32,4] <- ev.mba 
h2[32,5] <-	"MBA"

h2[33,1] <- "Number of Fruit"
h2[33,2] <- "2017"
h2[33,3] <- herit.pra$h2.obs
h2[33,4] <- ev.pra
h2[33,5] <- "PRA"
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
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod)
hist(residuals(seed17.mod))

seed17.mod2<- glmer(sm.2~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod2)
hist(residuals(seed17.mod2))

seed17.mod3<- glmer(sm.2~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
						 family = negative.binomial(theta = theta))
summary(seed17.mod3)
hist(residuals(seed17.mod))

A12 <- AIC(seed17.mod, seed17.mod2, seed17.mod3)
A12
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(seed17.mod2))
vars
print(VarCorr(seed17.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(seed17.mod2)['(Intercept)']*(nrow(dplyr::filter(flr17, flr17$Region =="GL_alvar"))/nrow(flr17))
gla
##Prairie region mean##
pra <-fixef(seed17.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr17, flr17$Region == "Prairie"))/nrow(flr17))
##MB_alvar region mean##
mba <- fixef(seed17.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr17, flr17$Region =="MB_alvar"))/nrow(flr17))

##look at model values to make sure mus makes sense##
fixef(seed17.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

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

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seed17.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr17~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr17,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seed17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seed17.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[34,1] <- "Seedmass (mg)"
h2[34,2] <- "2017"
h2[34,3] <- herit.gla$h2.obs
h2[34,4] <- ev.gla
h2[34,5] <- "GLA"

h2[35,1] <- "Seedmass (mg)"
h2[35,2] <- "2017"
h2[35,3] <- herit.mba$h2.obs
h2[35,4] <- ev.mba 
h2[35,5] <-	"MBA"

h2[36,1] <- "Seedmass (mg)"
h2[36,2] <- "2017"
h2[36,3] <- herit.pra$h2.obs
h2[36,4] <- ev.pra
h2[36,5] <- "PRA"
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
dtff18.mod<- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.mod)
hist(residuals(dtff18.mod))

dtff18.mod2 <- glmer(DTFF.18.Oday~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.mod2)
hist(residuals(dtff18.mod2))

dtff18.mod3<- glmer(DTFF.18.Oday~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link = log))
summary(dtff18.mod3)
hist(residuals(dtff18.mod3))

#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtff18.mod2))
vars
print(VarCorr(dtff18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtff18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtff18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtff18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(dtff18.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtff18.mod2)[1] * model.matrix(dtff18.mod2)[, 1] +  fixef(dtff18.mod2)[2] * model.matrix(dtff18.mod2)[, 2] +
	fixef(dtff18.mod2)[3] * model.matrix(dtff18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtff18.mod2)[1] * model.matrix(dtff18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtff18.mod2)[2] * model.matrix(dtff18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtff18.mod2)[3] * model.matrix(dtff18.mod2)[, 3])

#################Gamma custom###########
###IF RUNNING GAMMA DISTRIBUTION, NEED 'CUSTOM' MODEL DESGIN IN QGPARAMS##
##per wikipedia gamma consists of two parameters: shape parameter (k) and scale (theta)##
##https://stats.stackexchange.com/questions/96972/how-to-interpret-parameters-in-glm-with-family-gamma

##Pull parameters from fitdistr##
g <- fitdistr(flr18$DTFF18.18.Oday, "gamma")
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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtff18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[37,1] <- "DTFF"
h2[37,2] <- "2018"
h2[37,3] <- herit.gla$h2.obs
h2[37,4] <- ev.gla
h2[37,5] <- "GLA"

h2[38,1] <- "DTFF"
h2[38,2] <- "2018"
h2[38,3] <- herit.mba$h2.obs
h2[38,4] <- ev.mba 
h2[38,5] <-	"MBA"

h2[39,1] <- "DTFF"
h2[39,2] <- "2018"
h2[39,3] <- herit.pra$h2.obs
h2[39,4] <- ev.pra
h2[39,5] <- "PRA"
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
dtb.18.mod<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
							(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						family = gaussian)
dtb.18.out <-	summary(dtb.18.mod)
hist(residuals(dtb.18.mod))

dtb.18.mod2<- lmer(DtB.Oday.2018~Region + (1 | Population) + 
							(Region | Family.Unique) + (1 | Block.ID), data = flr18,
						)
dtb.18.out2 <-	summary(dtb.18.mod2)
hist(residuals(dtb.18.mod2))

dtb.18.gam2<- glmer(DtB.Oday.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link = 'log'))
dtb.18.gout2 <-	summary(dtb.18.gam2)
hist(residuals(dtb.18.gam2))

dtb.18.mod3 <- glmer(DtB.Oday.2018~Region + (Region | Population) + 
							(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						family = gaussian)
dtb.18.out3 <-	summary(dtb.18.mod3)
hist(residuals(dtb.18.mod3))

A14 <- AIC(dtb.18.mod, dtb.18.mod2, dtb.18.mod3)
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtb.18.mod2))
vars
print(VarCorr(dtb.18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtb.18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtb.18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtb.18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(dtb.18.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtb.18.mod2)[1] * model.matrix(dtb.18.mod2)[, 1] +  fixef(dtb.18.mod2)[2] * model.matrix(dtb.18.mod2)[, 2] +
	fixef(dtb.18.mod2)[3] * model.matrix(dtb.18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtb.18.mod2)[1] * model.matrix(dtb.18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtb.18.mod2)[2] * model.matrix(dtb.18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtb.18.mod2)[3] * model.matrix(dtb.18.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtb.18.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr18~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
#				family=gaussian(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp,  model = "gaussian")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtb.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp,  model = "gaussian")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtb.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp,  model = "gaussian")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[40,1] <- "Date to Bolt"
h2[40,2] <- "2018"
h2[40,3] <- herit.gla$h2.obs
h2[40,4] <- ev.gla
h2[40,5] <- "GLA"

h2[41,1] <- "Date to Bolt"
h2[41,2] <- "2018"
h2[41,3] <- herit.mba$h2.obs
h2[41,4] <- ev.mba 
h2[41,5] <-	"MBA"

h2[42,1] <- "Date to Bolt"
h2[42,2] <- "2018"
h2[42,3] <- herit.pra$h2.obs
h2[42,4] <- ev.pra
h2[42,5] <- "PRA"
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

dtfr.18.gam<- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr.18.gam)
hist(residuals(dtfr.18.gam))
dtfr.18.out <-	summary(dtfr.18.gam)

dtfr.18.gam2 <- glmer(Date.to.Fruit.Oday.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr.18.gam2)
hist(residuals(dtfr.18.gam2))
dtfr.18.out <-	summary(dtfr.18.gam2)

dtfr.18.gam3 <- glmer(Date.to.Fruit.Oday.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family = Gamma(link=log))
summary(dtfr.18.gam3)
hist(residuals(dtfr.18.gam3))
dtfr.18.out <-	summary(dtfr.18.gam3)

A15 <- AIC(dtfr.18.gam, dtfr.18.gam2, dtfr.18.gam3)

#####
#####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL AND DF##
##AND ALSO h2 TABLE ROW REFERENCE##
##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(dtfr.18.gam2))
vars
print(VarCorr(dtfr.18.gam2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(dtfr.18.gam2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(dtfr.18.gam2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(dtfr.18.gam2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(dtfr.18.gam2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(dtfr.18.gam2)[1] * model.matrix(dtfr.18.gam2)[, 1] +  fixef(dtfr.18.gam2)[2] * model.matrix(dtfr.18.gam2)[, 2] +
	fixef(dtfr.18.gam2)[3] * model.matrix(dtfr.18.gam2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(dtfr.18.gam2)[1] * model.matrix(dtfr.18.gam2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(dtfr.18.gam2)[2] * model.matrix(dtfr.18.gam2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(dtfr.18.gam2)[3] * model.matrix(dtfr.18.gam2)[, 3])

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
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtfr.18.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(dtfr.18.gam2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, custom.gamel = custom.functions)
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[43,1] <- "Date to Fruit"
h2[43,2] <- "2018"
h2[43,3] <- herit.gla$h2.obs
h2[43,4] <- ev.gla
h2[43,5] <- "GLA"

h2[44,1] <- "Date to Fruit"
h2[44,2] <- "2018"
h2[44,3] <- herit.mba$h2.obs
h2[44,4] <- ev.mba 
h2[44,5] <-	"MBA"

h2[45,1] <- "Date to Fruit"
h2[45,2] <- "2018"
h2[45,3] <- herit.pra$h2.obs
h2[45,4] <- ev.pra
h2[45,5] <- "PRA"
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

n.flr.18.mod <- glmer(Total.Flowers.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.18.mod))
nflr.mod.out <- summary(n.flr.18.mod)

n.flr.18.mod2 <- glmer(Total.Flowers.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.18.mod2))
nflr.mod.out2<- summary(n.flr.18.mod2)

n.flr.18.mod3 <- glmer(Total.Flowers.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.flr.18.mod3))
nflr.mod.out3<- summary(n.flr.18.mod3)

A16 <- AIC(n.flr.18.mod, n.flr.18.mod2, n.flr.18.mod3)

####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.flr.18.mod2))
vars
print(VarCorr(n.flr.18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.flr.18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(n.flr.18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(n.flr.18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(n.flr.18.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.flr.18.mod2)[1] * model.matrix(n.flr.18.mod2)[, 1] +  fixef(n.flr.18.mod2)[2] * model.matrix(n.flr.18.mod2)[, 2] +
	fixef(n.flr.18.mod2)[3] * model.matrix(n.flr.18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.flr.18.mod2)[1] * model.matrix(n.flr.18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.flr.18.mod2)[2] * model.matrix(n.flr.18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.flr.18.mod2)[3] * model.matrix(n.flr.18.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.18.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr18~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.flr.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[46,1] <- "Number of Flowers"
h2[46,2] <- "2018"
h2[46,3] <- herit.gla$h2.obs
h2[46,4] <- ev.gla
h2[46,5] <- "GLA"

h2[47,1] <- "Number of Flowers"
h2[47,2] <- "2018"
h2[47,3] <- herit.mba$h2.obs
h2[47,4] <- ev.mba 
h2[47,5] <-	"MBA"

h2[48,1] <- "Number of Flowers"
h2[48,2] <- "2018"
h2[48,3] <- herit.pra$h2.obs
h2[48,4] <- ev.pra
h2[48,5] <- "PRA"
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

n.fruit.18.mod <- glmer(No.Fruit.2018~Region + (1 | Population) + 
						 	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.18.mod))
frt.out <- summary(n.fruit.18.mod)

n.fruit.18.mod2 <- glmer(No.Fruit.2018~Region + (1 | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.18.mod2))
frt.out2 <- summary(n.fruit.18.mod2)

n.fruit.18.mod3 <- glmer(No.Fruit.2018~Region + (Region | Population) + 
						 	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						 family=neg.bin(theta = theta))
hist(residuals(n.fruit.18.mod3))
frt.out3 <- summary(n.fruit.18.mod3)

A17 <- AIC(n.fruit.18.mod, n.fruit.18.mod2, n.fruit.18.mod3)

####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(n.fruit.18.mod2))
vars
print(VarCorr(n.fruit.18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(n.fruit.18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(n.fruit.18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(n.fruit.18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(n.fruit.18.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(n.fruit.18.mod2)[1] * model.matrix(n.fruit.18.mod2)[, 1] +  fixef(n.fruit.18.mod2)[2] * model.matrix(n.fruit.18.mod2)[, 2] +
	fixef(n.fruit.18.mod2)[3] * model.matrix(n.fruit.18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(n.fruit.18.mod2)[1] * model.matrix(n.fruit.18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(n.fruit.18.mod2)[2] * model.matrix(n.fruit.18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(n.fruit.18.mod2)[3] * model.matrix(n.fruit.18.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.18.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr18~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(n.fruit.18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[49,1] <- "Number of Fruit"
h2[49,2] <- "2018"
h2[49,3] <- herit.gla$h2.obs
h2[49,4] <- ev.gla
h2[49,5] <- "GLA"

h2[50,1] <- "Number of Fruit"
h2[50,2] <- "2018"
h2[50,3] <- herit.mba$h2.obs
h2[50,4] <- ev.mba 
h2[50,5] <-	"MBA"

h2[51,1] <- "Number of Fruit"
h2[51,2] <- "2018"
h2[51,3] <- herit.pra$h2.obs
h2[51,4] <- ev.pra
h2[51,5] <- "PRA"
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
						  	(1 | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out <-	summary(seeds18.mod)
hist(residuals(seeds18.mod))

seeds18.mod2 <- glmer(sm.3~Region + (1 | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out2 <-	summary(seeds18.mod2)
seeds18.out2
hist(residuals(seeds18.mod2))
confint(seeds18.mod2)
seeds18.mod3 <- glmer(sm.3~Region + (Region | Population) + 
						  	(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
						  family = negative.binomial(theta = theta))
seeds18.out3 <-	summary(seeds18.mod3)
hist(residuals(seeds18.mod3))
A18 <- AIC(seeds18.mod, seeds18.mod2, seeds18.mod3)
####################EVERYTHING BELOW IS SAME FOR EACH TRAIT, JUST CHANGE MODEL, DF##
##AND distribution to match model, h2 TABLE ROW REFERENCE##

##pull coefficients: intercept and variance components for QGglmm##
vars <- as.data.frame(VarCorr(seeds18.mod2))
vars
print(VarCorr(seeds18.mod2), comp = "Variance")
vars
##Family.Unique variance in GL_alvar##
va.gla <-vars[2,4]
va.gla

##Family.Unique variance in MB_alvar##
va.mba <-vars[3,4]

##Family.Unique variance in Prairie##
va.pra <- vars[4,4]

##View latent-scale values region-specific mean##
##GL_alvar region mean (intercept)##
gla<-fixef(seeds18.mod2)['(Intercept)']*(nrow(dplyr::filter(flr18, flr18$Region =="GL_alvar"))/nrow(flr18))
gla
##Prairie region mean##
pra <-fixef(seeds18.mod2)['RegionPrairie']*(nrow(dplyr::filter(flr18, flr18$Region == "Prairie"))/nrow(flr18))
##MB_alvar region mean##
mba <- fixef(seeds18.mod2)['RegionMB_alvar']*(nrow(dplyr::filter(flr18, flr18$Region =="MB_alvar"))/nrow(flr18))

##look at model values to make sure mus makes sense##
fixef(seeds18.mod2)

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times family effect due to half-sibling design##
va <- 4*va.gla
va

##variance of Fixed effects (from design matrix, script adapted from Nakagawa Shielzeth 2013 S4)##
Fixed <- fixef(seeds18.mod2)[1] * model.matrix(seeds18.mod2)[, 1] +  fixef(seeds18.mod2)[2] * model.matrix(seeds18.mod2)[, 2] +
	fixef(seeds18.mod2)[3] * model.matrix(seeds18.mod2)[, 3]
##Calculation of the variance in fitted values
VarF <- var(Fixed)
##probably not used in model--should break out by region when running analysis##

##separate by region for analysis?##
##great lake alvar (intercept)##
varGLA <- var(fixef(seeds18.mod2)[1] * model.matrix(seeds18.mod2)[, 1])
##Manitoba alvar variance##
varMBA <- var(fixef(seeds18.mod2)[2] * model.matrix(seeds18.mod2)[, 2])
##prairie fixed effect##
varPRA <- var(fixef(seeds18.mod2)[3] * model.matrix(seeds18.mod2)[, 3])

##total variance in trait##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seeds18.mod2)
##Null model (No fixed effect) to get residual error variance##
#m0 <- glmer(No.Days.to.flr18~1 + (1 | Population) + 
#					(Region | Family.Unique) + (1 | Block.ID) + (1 | Sample.ID), data = flr18,
#				family=negbin.log(link=log), REML=F)
##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.gla+vars[7,4]+vars[8,4]+varGLA#+log(1 + 1/exp(as.numeric(fixef(m0))))
vars
vars[8,4]
##compare to model variance components
print(vars)
##view additive variance (half sib design so should be family unique x 4)##
va
##view vp (total variance)##
vp

##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams## mu = region level mean--gla or mba or pra
herit.gla <- QGparams(mu = gla, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.gla

##############MBA################
##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.mba
va

##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seeds18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.mba+vars[7,4]+vars[8,4]+varMBA#+log(1 + 1/exp(as.numeric(fixef(m0))))
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.mba <- QGparams(mu = mba, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.mba

##values are for variables from model, not yet converted to observation scale##
##additive variance NOTE: 4 times value due to half-sibling design##
va <- 4*va.pra
va
##total variance in trait Vp##
##old route--looking at iccs via sjstats##
#library(sjstats)
#icc(seeds18.mod2)

##total variance (including ve by log link function from null model)##
##		Va/family	Pop	Block		Vfixed			Ve	##
vp <- va.pra+vars[7,4]+vars[8,4]+varPRA#+log(1 + 1/exp(as.numeric(fixef(m0))))
print(vars)
##Latent-scale narrow-sense heritability##
lh2 <- va/vp
lh2
##Run QGparams to convert to observation scale (gives real heritability values)##
##put in QGparams##
herit.pra <- QGparams(mu = pra, var.a = va, var.p = vp, theta = theta, model = "negbin.log")
herit.pra

##Evolvability##
##Per Ned 6/26/19: == Va/mu^2##
ev.gla <- va.gla/(gla^2)
ev.mba <- va.mba/(mba^2)
ev.pra <- va.pra/(pra^2)

##Create a table to compile heritabilities## 
h2[52,1] <- "Seedmass (mg)"
h2[52,2] <- "2018"
h2[52,3] <- herit.gla$h2.obs
h2[52,4] <- ev.gla
h2[52,5] <- "GLA"

h2[53,1] <- "Seedmass (mg)"
h2[53,2] <- "2018"
h2[53,3] <- herit.mba$h2.obs
h2[53,4] <- ev.mba 
h2[53,5] <-	"MBA"

h2[54,1] <- "Seedmass (mg)"
h2[54,2] <- "2018"
h2[54,3] <- herit.pra$h2.obs
h2[54,4] <- ev.pra
h2[54,5] <- "PRA"
############################################

##EVALUATE THAT MODELS RAN BY VIEWING COMPARISON##
##IF OBJECT NOT FOUND, MODELS FAILED--ASSIGN TRAIT NA VALUE (FOR NOW)##
A1
A2
A3
A4
A5
#A6 #MODEL BREAK
#A7 #MODEL BREAK
A8
A9
A10 #MODEL 3 SEEMS BETTER?
A11
A12 #MODEL 3 SEEMS BETTER?
#A13 #MODEL BREAK
A14
A15
A16
A17
A18
h2[16,3:4]<-NA
h2[17,3:4]<-NA
h2[18,3:4]<-NA
h2[19,3:4]<-NA
h2[20,3:4]<-NA
h2[21,3:4]<-NA
h2[37,3:4]<-NA
h2[38,3:4]<-NA
h2[39,3:4]<-NA
#write.csv(h2, "Regional_heritabilities_and_evolvability_with_sampleID.csv", row.names = F)
#rm(list=ls())

##Ve table##
##Create a table to compile heritabilities## 
##is this how it should be done?##
#ve <- (va-(h2.obs*vp)/h2.obs)
##or is Sample.ID functional?##
col.classes2 = c("character", "character", "numeric")
col.names2 = c("Trait", "Year", "Individual variance")
Ve<-read.table(text = "",colClasses = col.classes2, col.names =col.names2)

vars1 <- as.data.frame(VarCorr(germ.mod2))
ve1 <-vars1[1,4]
Ve[1,1] <- "Germination"
Ve[1,2] <- "2015"
Ve[1,3] <- ve1

vars2 <- as.data.frame(VarCorr(TrueLeaf.mod2))
ve2 <-vars2[1,4]
Ve[2,1] <- " True Leaf"
Ve[2,2] <- "2015"
Ve[2,3] <- ve2

vars3 <- as.data.frame(VarCorr(dtff.16.gam2))
ve3 <-vars3[1,4]
Ve[3,1] <- "DTFF"
Ve[3,2] <- "2016"
Ve[3,3] <- ve3

vars4 <- as.data.frame(VarCorr(n.flr.16.mod2))
ve4 <-vars4[1,4]
Ve[4,1] <- "Number of Flowers"
Ve[4,2] <- "2016"
Ve[4,3] <- ve4

vars5 <- as.data.frame(VarCorr(n.fruit.16.mod))
ve5 <-vars5[1,4]
Ve[5,1] <- "Number of Fruit"
Ve[5,2] <- "2016"
Ve[5,3] <- ve5

Ve[6,1] <- "seedmass"
Ve[6,2] <- "2016"
Ve[6,3] <- 0

Ve[7,1] <- "DTFF"
Ve[7,2] <- "2017"
Ve[7,3] <- 0

vars8 <- as.data.frame(VarCorr(dtb.17.gam2))
ve8 <-vars8[1,4]
Ve[8,1] <- "Date to Bolt"
Ve[8,2] <- "2017"
Ve[8,3] <- ve8

vars9 <- as.data.frame(VarCorr(dtfr.17.mod2))
ve9 <-vars9[1,4]
Ve[9,1] <- "Date to Fruit"
Ve[9,2] <- "2017"
Ve[9,3] <- ve9

vars10 <- as.data.frame(VarCorr(n.flr.17.mod2))
ve10 <-vars10[1,4]
Ve[10,1] <- "Number of Flowers"
Ve[10,2] <- "2017"
Ve[10,3] <- ve10

vars11 <- as.data.frame(VarCorr(n.fruit.18.mod2))
ve11 <-vars11[1,4]
Ve[11,1] <- "Number of Fruit"
Ve[11,2] <- "2017"
Ve[11,3] <- ve11

vars12 <- as.data.frame(VarCorr(seed17.mod2))
ve12 <-vars12[1,4]
Ve[12,1] <- "Seedmass (mg)"
Ve[12,2] <- "2017"
Ve[12,3] <- ve12


Ve[13,1] <- "Date to first flower"
Ve[13,2] <- "2018"
Ve[13,3] <- 0

vars14 <- as.data.frame(VarCorr(dtb18.gam2))
ve14 <-vars14[1,4]
Ve[14,1] <- "Date to Bolt"
Ve[14,2] <- "2018"
Ve[14,3] <- 0

vars15 <- as.data.frame(VarCorr(dtfr.18.gam2))
ve15 <-vars15[1,4]
Ve[15,1] <- "Date to Fruit"
Ve[15,2] <- "2018"
Ve[15,3] <- ve15

vars16 <- as.data.frame(VarCorr(n.flr.18.mod2))
ve16 <-vars16[1,4]
Ve[16,1] <- "Number of Flowers"
Ve[16,2] <- "2018"
Ve[16,3] <- ve16

vars17 <- as.data.frame(VarCorr(n.fruit.18.mod2))
ve17 <-vars17[1,4]
Ve[17,1] <- "Number of Fruit"
Ve[17,2] <- "2018"
Ve[17,3] <- ve17

vars18 <- as.data.frame(VarCorr(seeds18.mod2))
ve18 <-vars18[1,4]
Ve[18,1] <- "Seedmass (mg)"
Ve[18,2] <- "2018"
Ve[18,3] <- ve18

write.csv(Ve, "Individual variance (Ve) of all traits.csv")
