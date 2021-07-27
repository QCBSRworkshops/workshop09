## #Species community data frame (fish abundance): “DoubsSpe.csv”
## spe<- read.csv(file.choose(), row.names=1)
## spe<- spe[-8,] #Site number 8 contains no species and so row 8 (site 8) is removed. Be careful to
## #only run this command line once as you are overwriting "spe" each time.
## 
## #Environmental data frame: “DoubsEnv.csv”
## env<- read.csv(file.choose(), row.names=1)
## env<- env[-8,] #Remove corresponding abiotic data for site 8 (since removed from fish data).
## #Again, be careful to only run the last line once.

## names(spe) #see names of columns in spe
## dim(spe) #dimensions of spe; number of columns and rows
## str(spe) #displays internal structure of objects
## head(spe) #first few rows of the data frame
## summary(spe) #summary statistics for each column; min value, median value, max value, mean value etc.

## #Species distribution
## (ab <- table(unlist(spe))) #note that when you put an entire line of code in brackets like this, the output for that operation is displayed right away in the R console
## 
## barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

## sum(spe==0)

## sum(spe==0)/(nrow(spe)*ncol(spe))

## spe.pres <- colSums(spe>0) #compute the number of sites where each species is present.
## hist(spe.pres, main="Species occurrence", las=1, xlab="Frequency of occurrences", breaks=seq(0,30, by=5), col="grey")

## site.pres <- rowSums(spe>0) #number of species with a value greater than 0 in that site row
## hist(site.pres, main="Species richness", las=1, xlab="Frequency of sites", ylab="Number of species", breaks=seq(0,30, by=5), col="grey")

## names(env)
## dim(env)
## str(env)
## head(env)
## summary(env)
## pairs(env, main="Bivariate Plots of the Environmental Data" )

## env.z <- decostand(env, method="standardize")
## apply(env.z, 2, mean) # the data are now centered (means~0)
## apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)
