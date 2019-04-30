# generate MLE of saturated model mean value parameter vector: mu
pout<- predict.aster(aout, se.fit=TRUE)


#exploring directions of recession
fred <- eigen(aout$fisher, symmetric = TRUE)
dor <- fred$vectors[ , fred$values == min(fred$values)]
names(dor) <- names(aout$coefficients)
dor <- zapsmall(dor / max(dor))
dor

modmat <- aout$modmat
dim(modmat)

modmat <- as.vector(modmat)
modmat <- matrix(modmat, ncol = length(dor))
dor.phi <- modmat %*% dor
dor.phi <- as.vector(dor.phi)


unique(dor.phi)

sum(dor.phi)

foo <- data.frame(Population = as.character(redata2016$Population),
                  id = redata2016$id,
                  varb = as.character(redata2016$varb),
                  resp = redata2016$resp, stringsAsFactors = FALSE)
foo <- foo[dor.phi == 1, ]

foo

unique(foo$Population) #population WA-BLK

unique(foo$varb) #seed mass

unique(foo$id) #these individuals

unique(foo$resp)# just zeros

#each individual in Population WA-BLK failed to set any seeds (all zeros)

# Various options exist: merge this pop with anohter, delete this pop from data
# and just "say" it has zero fitness, change the graphical model: in this case, make 
# fruit number terminal fitness node?

#looking at data, all failed to set fruit, too. But, some (2 inds) produced
#a single flower...maybe go with this?


#take a closer look at data
tapply(dat2$sm, dat2$Population, max) # hmmm 9 families set a max of zero seeds

tapply(dat2$No.Fruit.2016, dat2$Population, max)# same nine families failed to set >0 fruit, too

tapply(dat2$No.Flowers.2016, dat2$Population, max)# only three families did not produce >0 number of flowers 


#####################################################################
#  So, try No.Flowers.2016 as fitness node (i.e. truncate model)
#####################################################################

#set response variables -> these represent variables in graphical model
vars<- c("Germination.Y.N","Survival.Y.N", "Flower.Y.N.2016", "No.Flowers.2016")


#reshape data so that all response variables are located in a single vector in a new data
#set called "redata"
redata2016 <- reshape(dat2, varying = list(vars), direction = "long",timevar = "varb", times = as.factor(vars), v.names = "resp")


#Designation of fitness variable for 2016 data
fit <- grepl("No.Flowers.2016", as.character(redata2016$varb))
fit<- as.numeric(fit)

redata2016b$fit <- fit

#check
with(redata2016, sort(unique(as.character(varb)[fit == 0])))
with(redata2016, sort(unique(as.character(varb)[fit == 1])))


#add a variable "root" to redata files, where value is 1
redata2016b<- data.frame(redata2016, root=1)



#load aster package
library(aster)

#set graphical mode and dist. for fitness nodes
pred<- c(0,1,2,3)
fam<- c(1,1,1,2) #might want to play with these distributions, especially seedmass

#describe dist. of preds.
sapply(fam.default(), as.character)[fam]


aout<- aster(resp~varb + fit:(Population), pred, fam, varb, id, root, data=redata2016)

summary(aout, show.graph = TRUE) #same problem, this time pop "AB-RO" all failed to flower!




# so hard way, is to redisign the data, so that aster consideres the whole data set
# as a single individual

#But we have to re-specify the graphical model, which is difficult

#kill the graph model
outies <- dor.phi == 1
subdata <- redata2016[! outies, ]

#save real ID values
id <- subdata$id

#call it all one indiv
subdata$id <- 1

#remake the graph
idx <- seq(1, nrow(subdata))
varb <- as.character(subdata$varb)
pred <- rep(NA, length(idx))
fam <- rep(NA, length(idx))
pred[varb == "sm"] <- 0
fam[varb == "sm"] <- 1
head(idx[varb == "sm"])




sum(varb == "Germination.Y.N") == sum(varb == "Survival.Y.N")#true

pred[varb == "Germination.Y.N"] <- idx[varb == "Survival.Y.N"]
fam[varb == "Germination.Y.N."] <- 1

sum(varb == "Survival.Y.N") == sum(varb == "Flower.Y.N.2016")#true

pred[varb == "Survival.Y.N"] <- idx[varb == "Flower.Y.N.2016"]#true
fam[varb == "Survival.Y.N"] <- 1

sum(varb == "Flower.Y.N.2016") == sum(varb == "No.Flowers.2016")#true

pred[varb == "Flower.Y.N.2016"] <- idx[varb == "No.Flowers.2016"]
fam[varb == "Flower.Y.N.2016"] <- 1


sum(varb == "Flower.Y.N.2016") == sum(varb == "No.Flowers.2016")#true

pred[varb == "Flower.Y.N.2016"] <- idx[varb == "No.Flowers.2016"]
fam[varb == "Flower.Y.N.2016"] <- 1

sum(varb == "No.Flowers.2016") == sum(varb == "No.Fruit.2016")# true

pred[varb == "No.Flowers.2016"] <- idx[varb == "No.Fruit.2016"]
fam[varb == "No.Flowers.2016"] <- 2

sum(varb == "No.Fruit.2016") == sum(varb == "sm")# false

bar <- match(id[varb == "sm"], id[varb == "No.Fruit.2016"])
pred[varb == "sm"] <- idx[varb == "No.Fruit.2016"][bar]
fam[varb == "sm"] <- 2

#Now remake "varb"
subvarb <- paste(as.character(subdata$varb), id, sep = "")
subdata <- data.frame(subdata, subvarb = subvarb)

aout.sub <- aster(resp ~ varb + fit : Population,
                  +pred, fam, subvarb, id, root, data = subdata)
summary(aout.sub)