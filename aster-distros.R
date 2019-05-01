setwd("C:/Users/Zebadiah/Onedrive/NDSU/r-mixed-models")
##load data##
# my old data: df2 <- read.csv('NV_CG_Experiment-mg.csv')
##cleanest data from mason's git##
df <- read.csv('NV_CG_Experiment2wdist2.csv')
##load libraries##
library(tidyverse);library(edgeR)
#new package##
#install.packages("fitdistrplus")
library(fitdistrplus)
#install.packages("blmeco")
library(blmeco)

##Days to Germination##
##############################
##lets check distributions (again)##
hist(df$No.Days.to.Germ)
##looks neg binomial##
class(df$No.Days.to.Germ)
range(df$No.Days.to.Germ)
##df w/out NA's for germination##
germ <- df[!is.na(df$No.Days.to.Germ),]
descdist(germ$No.Days.to.Germ, boot = 100)
f1g <- fitdist(germ$No.Days.to.Germ, "pois")
plot(f1g)
summary(f1g)
mean(germ$No.Days.to.Germ)
var(germ$No.Days.to.Germ)
##############################


##Number of Flowers 2016##
##############################
##Number of Flowers 2016##
flr.16 <- filter(df, Flower.Y.N.2016 >= 1)
nfl.hist<- ggplot(flr.16, aes(x=No.Flowers.2016)) + theme_bw() + 
	geom_histogram(aes(y=..density..),      
						binwidth=1,
						colour="black", fill="cyan") +
	geom_density(alpha=.2, fill="#FF6666") 
nfl.hist
hist(flr.16$No.Flowers.2016)
class(flr.16$No.Flowers.2016)
range(flr.16$No.Flowers.2016)
##df w/out NA's##
nfl16 <- flr.16[!is.na(flr.16$No.Flowers.2016),]
class(nfl16$No.Flowers.2016)
range(nfl16$No.Flowers.2016)
descdist(nfl16$No.Flowers.2016, boot = 100)
##NOTE--calculated from MASS, so comparison of AIC works##
f1g <- fitdistr(nfl16$No.Flowers.2016, "normal")
f2g <- fitdistr(nfl16$No.Flowers.2016, "poisson")
f3g <- fitdistr(nfl16$No.Flowers.2016, "negative binomial") ##estimate 1.2732
nfl16.nb<-f3g$estimate[1]

f4g <- fitdistr(nfl16$No.Flowers.2016, "gamma")
f5g <- fitdistr(nfl16$No.Flowers.2016, "exponential")
library(MASS)		
AIC(f1g,f2g,f3g,f4g,f5g)
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
##############################
##Conclusion: Negative binomial seems most appropriate##
summary(f3g)
##############################

##number fruit##
##############################
hist(df$No.Fruit.2016)
df1 <- filter(df, No.Fruit.2016 >= 0)
summary(df1$No.Fruit.2016)
nfl.hist<- ggplot(flr.16, aes(x=No.Fruit.2016)) + theme_bw() + 
	geom_histogram(aes(y=..density..),      
						binwidth=1,
						colour="black", fill="cyan") +
	geom_density(alpha=.2, fill="#FF6666") 
nfl.hist
hist(flr.16$No.Fruit.2016)
class(flr.16$No.Fruit.2016)
range(flr.16$No.Fruit.2016)
##df w/out NA's##
nfl16 <- flr.16[!is.na(flr.16$No.Fruit.2016),]
class(nfl16$No.Fruit.2016)
range(nfl16$No.Fruit.2016)
descdist(nfl16$No.Fruit.2016, boot = 100)
f1g <- fitdist(nfl16$No.Fruit.2016, "norm")
f2g <- fitdist(nfl16$No.Fruit.2016, "pois")
f3g <- fitdist(nfl16$No.Fruit.2016, "nbinom") ##parameter 0.1391973
nfr16.nb<-f3g$estimate[1]

f4g <- fitdist(nfl16$No.Fruit.2016, "lnorm")
##doesn't work##
f5g <- fitdist(nfl16$No.Fruit.2016, "exp")
?gofstat
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)
##Visual comparison of distribution##
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
summary(f3g)
##############################
##Conclusion: Looks fine w/ Pois or nbinom##
##############################


##seedmass 2016##
##############################
hist(df$sm)
##df w/out NA's##
seeds16 <- df[!is.na(df$sm),]
class(seeds16$sm)
range(seeds16$sm)
descdist(seeds16$sm, boot = 100)
f1g <- fitdist(seeds16$sm, "norm")
f2g <- fitdist(seeds16$sm, "pois")
f3g <- fitdist(seeds16$sm, "nbinom") ##estimate 68466523
sm16.nb <- f3g$estimate[1]

f4g <- fitdist(seeds16$sm, "lnorm")
##lnorm not working##
f5g <- fitdist(seeds16$sm, "exp")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: again either pois or nbinom--probably ZIP##
##############################


########2017 Season###########


##DTFF17##
##############################
#make date number?
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
hist(flr.17$DTFF.Ordinal.Day.2017)
flr17 <- flr.17[!is.na(flr.17$DTFF.Ordinal.Day.2017),]
f1g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "norm")
f2g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "pois")
f3g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "nbinom")
f4g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "lnorm")
f5g <- fitdist(flr17$DTFF.Ordinal.Day.2017, "exp")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
flr.17$DTFF.Ordinal.Day.2017 <-as.numeric(flr.17$DTFF.Ordinal.Day.2017)
##############################


##number flowers 2017##
##############################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
nfl.hist<- ggplot(flr.17, aes(x=Total.Flowers.2017)) + theme_bw() + 
	geom_histogram(aes(y=..density..),      
						binwidth=1,
						colour="black", fill="cyan") +
	geom_density(alpha=.2, fill="#FF6666") 
nfl.hist
hist(flr.17$Total.Flowers.2017)
class(flr.17$Total.Flowers.2017)
range(flr.17$Total.Flowers.2017)
##df w/out NA's for Total flowers##
nfl17 <- flr.17[!is.na(flr.17$Total.Flowers.2017),]
class(nfl17$Total.Flowers.2017)
range(nfl17$Total.Flowers.2017)
descdist(nfl17$Total.Flowers.2017, boot = 100)
f1g <- fitdist(nfl17$Total.Flowers.2017, "norm")
f2g <- fitdist(nfl17$Total.Flowers.2017, "pois")
f3g <- fitdist(nfl17$Total.Flowers.2017, "nbinom") ##estimate: 1.753497
nfl17.nb<-f3g$estimate[1]
f4g <- fitdist(nfl17$Total.Flowers.2017, "lnorm")
f5g <- fitdist(nfl17$Total.Flowers.2017, "exp")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
summary(f1g)
summary(f3g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: Looks most like nbinom, next closest is normal?##
##############################


##number fruit 2017##
##############################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
nfl.hist<- ggplot(flr.17, aes(x=No.Fruit.2017)) + theme_bw() + 
	geom_histogram(aes(y=..density..),      
						binwidth=1,
						colour="black", fill="cyan") +
	geom_density(alpha=.2, fill="#FF6666") 
nfl.hist
hist(flr.17$No.Fruit.2017)
class(flr.17$No.Fruit.2017)
range(flr.17$No.Fruit.2017)
##df w/out NA's for Total flowers##
nfl17 <- flr.17[!is.na(flr.17$No.Fruit.2017),]
class(nfl17$No.Fruit.2017)
range(nfl17$No.Fruit.2017)
descdist(nfl17$No.Fruit.2017, boot = 100)
f1g <- fitdist(nfl17$No.Fruit.2017, "norm")
f2g <- fitdist(nfl17$No.Fruit.2017, "pois")
f3g <- fitdist(nfl17$No.Fruit.2017, "nbinom") ##estimate 1.951956
nfr17.nb<-f3g$estimate[1]
f4g <- fitdist(nfl17$No.Fruit.2017, "lnorm")
##lnorm doesn't work##
f5g <- fitdist(nfl17$No.Fruit.2017, "exp")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
summary(f3g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: nbinom##
##############################


##Seedmass17##
##############################
flr.17 <- filter(df, Flower.Y.N.2017 >= 1)
flr.17 <- filter(df, sm.2 >=1)
#View(flr.17)
hist(df$sm.2)
summary(df$sm.2)
flr17 <- flr.17[!is.na(flr.17$sm.2),]
f1g <- fitdist(flr17$sm.2, "norm")
f2g <- fitdist(flr17$sm.2, "pois")
f3g <- fitdist(flr17$sm.2, "nbinom") ##estimate .9277568
sm17.nb<-f3g$estimate[1]
plot(f1g)
plot(f2g)
plot(f3g)
str(flr17)
summary(flr17$sm.2)					  
glimpse(flr17$sm.2)
mean(df$sm.2, na.rm =T)
summary(seed17.mod)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: nbinom##

########2018 Season###########

##DTFF18##
##############################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(df$DTFF.18.Oday)
hist(flr.18$DTFF.18.Oday)
flr18 <- flr.18[!is.na(flr.18$DTFF.18.Oday),]
f1g <- fitdist(flr18$DTFF.18.Oday, "norm")
f2g <- fitdist(flr18$DTFF.18.Oday, "pois")
f3g <- fitdist(flr18$DTFF.18.Oday, "nbinom")
f4g <- fitdist(flr18$DTFF.18.Oday, "lnorm")
f5g <- fitdist(flr18$DTFF.18.Oday, "exp")
plot(f1g)
plot(f2g)
plot(f3g)
plot(f4g)
plot(f5g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: Poisson##

##Total Flowers 2018##
##############################
str(df)
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
flr18 <- flr.18[!is.na(flr.18$Total.Flowers.2018),]
summary(flr.18$Total.Flowers.2018)
hist(flr.18$Total.Flowers.2018)
f1g <- fitdist(flr18$Total.Flowers.2018, "norm")
f2g <- fitdist(flr18$Total.Flowers.2018, "pois")
fg3 <- fitdist(log(flr18$Total.Flowers.2018), "norm")
f4g <- fitdist(flr18$Total.Flowers.2018, "nbinom") ##estimate 1.81233
nfl18.nb<-f4g$estimate[1]

plot(f1g)
plot(f2g)
plot(fg3)
plot(f4g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: nbinom or lnorm dist##

##No.Fruit.2018##
##############################
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
flr18 <- flr.18[!is.na(flr.18$No.Fruit.2018),]
summary(flr.18$No.Fruit.2018)
hist(flr.18$No.Fruit.2018)
f1g <- fitdist(flr18$No.Fruit.2018, "norm")
f2g <- fitdist(flr18$No.Fruit.2018, "pois")
fg3 <- fitdist(log(flr18$No.Fruit.2018), "norm")
f4g <- fitdist(flr18$No.Fruit.2018, "nbinom") ##estimate 2.511923
nfr18.nb<-f4g$estimate[1]
plot(f1g)
plot(f2g)
plot(fg3)
plot(f4g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Conclusion: nbinom dist##


##seedmass2018##
##############################
str(df)
summary(df$seedmass.2018.g.)
class(df$seedmass.2018.g.)
flr.18 <- filter(df, Flowering.Y.N.2018 >= 1)
hist(flr.18$seedmass.2018.g.)
flr18 <- flr.18[!is.na(flr.18$seedmass.2018.g.),]
f1g <- fitdist(flr18$seedmass.2018.g., "norm")
f2g <- fitdist(flr18$seedmass.2018.g., "pois")
fg3 <- fitdist(flr18$seedmass.2018.g., "lnorm")
f4g <- fitdist(flr18$seedmass.2018.g., "nbinom") ##estimate N/A
sm18.nb<-f4g$estimate[1]

summary(f3g)
plot(f1g)
plot(f2g)
plot(fg3)
plot(f4g)
fits <- list(f1g,f2g,f3g,f5g)
gofstat(fits)

##############################
##Says failed to converge for poisson and nbinom##

