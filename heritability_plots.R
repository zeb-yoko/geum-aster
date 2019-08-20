##script for visualizing trait heritabilities across years##
##zebadiah.yoko@gmail.com##
library(tidyverse)
##run theme_zpub function (separate script) to get graphs to work##

df <- read.csv("Region-level_heritabilities_evolvability_no_CI.csv")
df$Year <-as.factor(df$Year)
df$Region <- factor(df$Region, levels = c("PRA", "MBA", "GLA"))
##set colors	#Prairie   #GL ALvar  #MB Alvar##
col.esa <- c("#18563E", "#82BE42", "#FFC423")
col.grscl <- c("gray80", "black", "gray43")

##single-year traits##
##Germination##
filter(df, Trait =="Germination") %>% 
	ggplot(aes(x=Region, y=Heritability, col = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##Germination##
filter(df, Trait =="TrueLeaf") %>% 
	ggplot(aes(x=Region, y=Heritability, col = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##visualize multi-year traits##
##DTFF##
filter(df, Trait =="DTFF") %>% 
	ggplot(aes(x=Year, y=Heritability, col = Region, group = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##Date to Bolt##
filter(df, Trait =="Date to Bolt") %>% 
	ggplot(aes(x=Year, y=Heritability, col = Region, group = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##Number of Flowers##
filter(df, Trait =="Number of Flowers") %>% 
	ggplot(aes(x=Year, y=Heritability, col = Region, group = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##Number of Fruit##
filter(df, Trait =="Number of Fruit") %>% 
	ggplot(aes(x=Year, y=Heritability, col = Region, group = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)

##seedmass##
filter(df, Trait == "Seedmass (mg)") %>% 
	ggplot(aes(x=Year, y=Heritability, col = Region, group = Region)) +
	geom_point() +geom_line() + theme_zpub()+
	scale_color_manual("legend", values = col.esa)
