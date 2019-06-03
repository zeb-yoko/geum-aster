##Combine cleaned fitness/aster analysis data with phenology data##
##preserves NAs in data not included in Aster analysis##

##full dataset##
df <- read.csv("NV_CG_Experiment2wdist2.csv")
ncol(df) #93
library(tidyverse)
##remove columns in aster analysis## should remove 18 columns
df <- dplyr::select(df, -c(Germination.Y.N, Survival.Y.N, Survival.Y.N.2017, 
					  Flower.Y.N.2016, Flower.Y.N.2017, No.Flowers.2016, Total.Flowers.2017,
					  Fruit.Y.N.2016, Fruit.Y.N.2017, No.Fruit.2016, No.Fruit.2017,
					  sm, sm.2, Flowering.Y.N.2018, Total.Flowers.2018, Fruit.Y.N.2018, 
					  No.Fruit.2018, seedmass.2018.g.))
ncol(df) #75, good

##clean fitness data##
df2 <- read.csv("cleaned_data_for_aster.csv")

#check number of columns##
ncol(df2) # should be 18 selected + 6 shared in merge + new variables SURV2017, SURV2018, sm2017, sm2018, & sm.3)
str(df2) # good
##merge operation##
df.all <- merge(df, df2, by= c("Family.Unique", "Block.ID", "HabitatType",
				  "Region", "Population", "Dist.from.cg.km"), all =T)

##check number of columns##
ncol(df.all) #75+23 = 98
df.all
write.csv(df.all, "full_clean_geum_experiment.csv")
