##Aster-herit Visualizations!##
###############################
##Zebadiah.yoko@gmail.com##

df <- read.csv("Region-level_heritabilities_and_evolvability.csv")
library(tidyverse)

gla <- subset(df, Region == 'GLA')
mba <- subset(df, Region == 'MBA')
pra <- subset(df, Region == 'PRA')
y16 <- subset(df, Year != '2018' & Year !='2017')
y17 <- subset(df, Year == '2017')
y18 <- subset(df, Year == '2018')

##see if we can get panels a ggplot cheatsheet way##
gg <- ggplot(df, aes(x=Trait, y=Heritability, col = Region)) +
	geom_point() + facet_grid(.~Year)
gg				 

gg16 <- ggplot(y16, aes(x=Trait, y=Heritability, col = Region)) +
	geom_point(size=3) + theme_zpub()
gg16			 

gg17 <- ggplot(y17, aes(x=Trait, y=Heritability, col = Region)) +
	geom_point(size=3) + theme_zpub()
gg17	

gg18 <- ggplot(y18, aes(x=Trait, y=Heritability, col = Region)) +
	geom_point(size=3) + theme_zpub()
gg18	

ggsave('h2_2016.png', plot = gg16)
ggsave('h2_2017.png', plot = gg17)
ggsave('h2_2018.png', plot = gg18)

filter(df, Trait =="DTFF") %>% 
ggplot(aes(x=Year, y=Heritability, col = Region)) +geom_point()

##Fitness graphing data and code##
seeds <- read.csv('fitness_estimates.csv')
pd <- position_dodge(0.25) # move them .05 to the left and right

seeds$Year<-as.factor(seeds$Year)
sviz <-	ggplot(seeds, aes(x=Year, y=W, col = Region)) +
	geom_point(size=3, position = pd) + 
	geom_errorbar(aes(ymin=W-Stand.Err, ymax=W+Stand.Err), width=.1, position = pd ) +
	theme_zpub()+geom_hline(yintercept=0,linetype='dotted')
sviz

