##Aster-herit Visualizations!##
###############################
##Zebadiah.yoko@gmail.com##

df <- read.csv("Region-level_heritabilities_and_evolvability.csv")

##Create a table to compile fitness estimates## 
col.classes = c("character", "numeric", "numeric")
col.names = c("Region", "Fitness", "SE")
fitness16 <- read.table(text = "",colClasses = col.classes, col.names = col.names)
fitness16[1,1] <- "GLA"
fitness16[1,2] <- 3.027
fitness16[1,3] <- 1.246
fitness16[2,1] <- "MBA"
fitness16[2,2] <- 0
fitness16[2,3] <- NA
fitness16[3,1] <- "PRA"
fitness16[3,2] <- 0
fitness16[3,3] <- NA


fitness16

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

gg16 <- ggplot(y16, aes(x=Trait, y=Heritability, col = Region)) +geom_point()
gg16			 

filter(df, Trait =="DTFF") %>% 
ggplot(aes(x=Year, y=Heritability, col = Region)) +geom_point()
