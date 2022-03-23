##Section: 01-preparing-for-the-workshop.R 

###Notice ###
#                                                                            #
#This is an automatically generated script based on the code chunks from the #
#book for this workshop.                                                     #
#                                                                            #
#It is minimally annotated to allow participants to provide their comments:  #
#a practice that we highly encourage.                                        #
#                                                                            #
#Note that the solutions to the challenges are also included in this script. #
#When solving the challenges by yourself, attempt to not scroll and peek at  #
#the solutions.                                                              #
#                                                                            #
#Happy coding!                                                               #


install.packages("vegan")
install.packages("gclus")
install.packages("ape")

library(vegan)
library(gclus)
library(ape)

source(file.choose()) # use coldiss.R which you have downloaded to your own directory


##Section: 02-introduction.R 




##Section: 03-data-exploration.R 

#Species community data frame (fish abundance): “DoubsSpe.csv”
spe<- read.csv(file.choose(), row.names=1)
spe<- spe[-8,] #Site number 8 contains no species and so row 8 (site 8) is removed. Be careful to
#only run this command line once as you are overwriting "spe" each time.

#Environmental data frame: “DoubsEnv.csv”
env<- read.csv(file.choose(), row.names=1)
env<- env[-8,] #Remove corresponding abiotic data for site 8 (since removed from fish data).
#Again, be careful to only run the last line once.

names(spe) #see names of columns in spe
dim(spe) #dimensions of spe; number of columns and rows
str(spe) #displays internal structure of objects
head(spe) #first few rows of the data frame
summary(spe) #summary statistics for each column; min value, median value, max value, mean value etc.

#Species distribution
(ab <- table(unlist(spe))) #note that when you put an entire line of code in brackets like this, the output for that operation is displayed right away in the R console

barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(5:0/5))

sum(spe==0)

sum(spe==0)/(nrow(spe)*ncol(spe))

spe.pres <- colSums(spe>0) #compute the number of sites where each species is present.
hist(spe.pres, main="Species occurrence", las=1, xlab="Frequency of occurrences", breaks=seq(0,30, by=5), col="grey")

site.pres <- rowSums(spe>0) #number of species with a value greater than 0 in that site row
hist(site.pres, main="Species richness", las=1, xlab="Frequency of sites", ylab="Number of species", breaks=seq(0,30, by=5), col="grey")

names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Bivariate Plots of the Environmental Data" )

env.z <- decostand(env, method="standardize")
apply(env.z, 2, mean) # the data are now centered (means~0)
apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)


##Section: 04-association-distances.R 

spe.db<-vegdist(spe, method="bray") # "bray" with presence-absence data is Sorensen dissimilarity
spe.dj<-vegdist(spe, method="jac") # Jaccard dissimilarity
spe.dg<-vegdist(spe, method="gower") # Gower dissimilarity
spe.db<-as.matrix(spe.db) #Put in matrix form (can visualize, write to .csv etc)

windows()
coldiss(spe.db, byrank=FALSE, diag=TRUE) # Heat map of Bray-Curtis dissimilarity
windows()
coldiss(spe.dj, byrank=FALSE, diag=TRUE) # Heat map of Jaccard dissimilarity
windows()
coldiss(spe.dg, byrank=FALSE, diag=TRUE) # Heat map of Gower dissimilarity

?dist # this function also compute dissimilarity matrix
env.de<-dist(env.z, method = "euclidean") # euclidean distance matrix of the standardized environmental variables
windows() #Creates a separate graphical window
coldiss(env.de, diag=TRUE)

(env.pearson<-cor(env)) # Pearson r linear correlation
round(env.pearson, 2) #Rounds the coefficients to 2 decimal points
(env.ken<-cor(env, method="kendall")) # Kendall tau rank correlation
round(env.ken, 2)

var.g1<-rnorm(30, 0, 1)
var.g2<-runif(30, 0, 5)
var.g3<-gl(3, 10)
var.g4<-gl(2, 5, 30)
(dat2<-data.frame(var.g1, var.g2, var.g3, var.g4))
str(dat2)
summary(dat2)

?daisy #This function can handle NAs in the data
(dat2.dg<-daisy(dat2, metric="gower"))
coldiss(dat2.dg)

spe.challenge<-spe[1:3,1:3] #”[1:3,” refers to rows 1 to 3 while “,1:3]” refers to the first 3 species columns (in #this case the three variables of interest)

(Abund.s1<-sum(spe.challenge[1,]))
(Abund.s2<-sum(spe.challenge[2,]))
(Abund.s3<-sum(spe.challenge[3,]))
#() around code will cause output to print right away in console

Spec.s1s2<-0
Spec.s1s3<-0
Spec.s2s3<-0
for (i in 1:3) {
  Spec.s1s2<-Spec.s1s2+abs(sum(spe.challenge[1,i]-spe.challenge[2,i]))
  Spec.s1s3<-Spec.s1s3+abs(sum(spe.challenge[1,i]-spe.challenge[3,i]))
  Spec.s2s3<-Spec.s2s3+abs(sum(spe.challenge[2,i]-spe.challenge[3,i])) }

(db.s1s2<-Spec.s1s2/(Abund.s1+Abund.s2)) #Site 1 compared to site 2
(db.s1s3<-Spec.s1s3/(Abund.s1+Abund.s3)) #Site 1 compared to site 3
(db.s2s3<-Spec.s2s3/(Abund.s2+Abund.s3)) #Site 2 compared to site 3

(spe.db.challenge<-vegdist(spe.challenge, method="bray"))

# Calculate the number of columns in your dataset
M<-ncol(spe.challenge)

# Calculate the species abundance differences between pairs of sites for each species
Spe1.s1s2<-abs(spe.challenge[1,1]-spe.challenge[2,1])
Spe2.s1s2<-abs(spe.challenge[1,2]-spe.challenge[2,2])
Spe3.s1s2<-abs(spe.challenge[1,3]-spe.challenge[2,3])
Spe1.s1s3<-abs(spe.challenge[1,1]-spe.challenge[3,1])
Spe2.s1s3<-abs(spe.challenge[1,2]-spe.challenge[3,2])
Spe3.s1s3<-abs(spe.challenge[1,3]-spe.challenge[3,3])
Spe1.s2s3<-abs(spe.challenge[2,1]-spe.challenge[3,1])
Spe2.s2s3<-abs(spe.challenge[2,2]-spe.challenge[3,2])
Spe3.s2s3<-abs(spe.challenge[2,3]-spe.challenge[3,3])

# Calculate the range of each species abundance between sites
Range.spe1<-max(spe.challenge[,1]) - min (spe.challenge[,1])
Range.spe2<-max(spe.challenge[,2]) - min (spe.challenge[,2])
Range.spe3<-max(spe.challenge[,3]) - min (spe.challenge[,3])

# Calculate the Gower dissimilarity
(dg.s1s2<-(1/M)*((Spe2.s1s2/Range.spe2)+(Spe3.s1s2/Range.spe3)))
(dg.s1s3<-(1/M)*((Spe2.s1s3/Range.spe2)+(Spe3.s1s3/Range.spe3)))
(dg.s2s3<-(1/M)*((Spe2.s2s3/Range.spe2)+(Spe3.s2s3/Range.spe3)))

# Compare your results
(spe.db.challenge<-vegdist(spe.challenge, method="gower"))

spe.pa<-decostand(spe, method="pa")

#Hellinger transformation
spe.hel<-decostand(spe, method="hellinger") #can also use method=”hell”

#Chi-square transformation
spe.chi<-decostand(spe, method="chi.square")

# Hellinger transformation
# First calculate the total species abundances by site
(site.totals=apply(spe, 1, sum))

# Scale species abundances by dividing them by site totals
(scale.spe<-spe/site.totals)

# Calculate the square root of scaled species abundances
(sqrt.scale.spe<-sqrt(scale.spe))

# Compare the results
sqrt.scale.spe
spe.hel
sqrt.scale.spe-spe.hel # or: sqrt.scale.spe/spe.hel

# Chi-square transformation
# First calculate the total species abundances by site
(site.totals<-apply(spe, 1, sum))

# Then calculate the square root of total species abundances
(sqrt.spe.totals<-sqrt(apply(spe, 2, sum)))

# Scale species abundances by dividing them by the site totals and the species totals
scale.spe2<-spe
for (i in 1:nrow(spe)) {
  for (j in 1:ncol(spe)) {
   (scale.spe2[i,j]=scale.spe2[i,j]/(site.totals[i]*sqrt.spe.totals[j]))   }}

#Adjust the scale abundance species by multiplying by the square root of the species matrix total
(adjust.scale.spe2<-scale.spe2*sqrt(sum(rowSums(spe))))

# Compare the results
adjust.scale.spe2
spe.chi
adjust.scale.spe2-spe.chi # or: adjust.scale.spe2/spe.chi

spe.dhel<-vegdist(spe.hel,method="euclidean") #generates the distance matrix from Hellinger transformed data

#See difference between the two matrices
head(spe.hel)# Hellinger-transformed species data
head(spe.dhel)# Hellinger distances among sites

#Faire le groupement à liens simples
#Perform single linkage clustering
spe.dhel.single<-hclust(spe.dhel, method="single")
plot(spe.dhel.single)

#Perform complete linkage clustering
spe.dhel.complete<-hclust(spe.dhel, method="complete")
plot(spe.dhel.complete)

#Perform Ward minimum variance clustering
spe.dhel.ward<-hclust(spe.dhel, method="ward.D2")
plot(spe.dhel.ward)

#Re-plot the dendrogram by using the square roots of the fusion levels
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward)
plot(spe.dhel.ward, hang=-1) # hang=-1 aligns all objets on the same line


##Section: 05-unconstrained-ordination.R 




##Section: 06-principal-component-analysis.R 

#Run the PCA using the rda() function (NB: rda() is used for both PCA and RDA)
spe.h.pca <- rda(spe.hel)

#Extract the results
summary(spe.h.pca) #overall results

summary(spe.h.pca, display=NULL) # only eigenvalues and their contribution to the variance
eigen(cov(spe.hel)) # also compute the eigenvalues

spe.scores <- scores(spe.h.pca, display="species", choices=c(1,2)) # species scores on the first two PCA axes
site.scores <- scores(spe.h.pca, display="sites", choices=c(1,2)) # sites scores on the first two PCA axes
#Note: if you don’t specify the number of principal components to extract (e.g. choices=c(1,2) or choices=c(1:2) then all of the scores will be extracted for all of the principal components.

# Identify the significant axis using the Kaiser-Guttman criterion
ev <- spe.h.pca$CA$eig
ev[ev>mean(ev)]
n <- length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

#Run the PCA
env.pca <- rda(env.z) # or rda(env, scale=TRUE)

#Extract the results
summary(env.pca)
summary(env.pca, scaling=2)

# Identify the significant axis using the Kaiser-Guttman criterion
ev <- env.pca$CA$eig
ev[ev>mean(ev)]
n <- length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

plot(spe.h.pca)

plot(spe.h.pca, type="n") #produces a blank biplot with nothing displayed but the axes
points(spe.h.pca, dis="sp", col="blue") #points are added for the species (columns) (dis=)
#use text() instead of points() if you want the labels
points(spe.h.pca, dis="sites", col="red") #points are added for the sites (rows)

#Biplot of the PCA on transformed species data (scaling 1)
windows()
plot(spe.h.pca)
windows()
biplot(spe.h.pca)
windows()
plot(spe.h.pca, scaling=1, type="none", # scaling 1 = distance biplot :
                                        # distances among objects in the biplot approximate their Euclidean distances
                                        # but angles among descriptor vectors DO NOT reflect their correlation
     xlab = c("PC1 (%)", round((spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig))*100,2)), #this comes from the summary
     ylab = c("PC2 (%)", round((spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig))*100,2)))
points(scores(spe.h.pca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(spe.h.pca, display="species", choices=c(1), scaling=1),
     scores(spe.h.pca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spe.h.pca, display="species", scaling=1)),
     col="red", cex=0.8)

#Biplot of the PCA on the environmental variables (scaling 2)
windows()
plot(env.pca)
windows()
plot(env.pca, scaling=2, type="none", # scaling 2 = correlation biplot :
                                      # distances among abjects in the biplot DO NOT approximate their Euclidean distances
                                      # but angles among descriptor vectors reflect their correlation
     xlab = c("PC1 (%)", round((env.pca$CA$eig[1]/sum(env.pca$CA$eig))*100,2)),
     ylab = c("PC2 (%)", round((env.pca$CA$eig[2]/sum(env.pca$CA$eig))*100,2)),
     xlim = c(-1,1), ylim=c(-1,1))
points(scores(env.pca, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="darkgreen", cex=1.2)
text(scores(env.pca, display="species", choices=c(1), scaling=2),
     scores(env.pca, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(env.pca, display="species", scaling=2)),
     col="red", cex=0.8)

Sites_scores_Env_Axis1<- scores(env.pca, display="sites", choices=c(1), scaling=2)
spe$ANG
plot( Sites_scores_Env_Axis1, spe$TRU)
summary(lm(spe$TRU~Sites_scores_Env_Axis1))
abline(lm(spe$TRU~Sites_scores_Env_Axis1))

mite.spe <- mite #mite data is from the vegan package

#Hellinger transformation of mite data and PCA
mite.spe.hel <- decostand(mite.spe, method="hellinger")
mite.spe.h.pca <- rda(mite.spe.hel)

#What are the significant axes?
ev <- mite.spe.h.pca$CA$eig
ev[ev>mean(ev)]
n <- length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

#Output summary/results
summary(mite.spe.h.pca, display=NULL)
windows()

#Plot the biplot
plot(mite.spe.h.pca, scaling=1, type="none",
     xlab=c("PC1 (%)", round((mite.spe.h.pca$CA$eig[1]/sum(mite.spe.h.pca$CA$eig))*100,2)),
     ylab=c("PC2 (%)", round((mite.spe.h.pca$CA$eig[2]/sum(mite.spe.h.pca$CA$eig))*100,2)))
points(scores(mite.spe.h.pca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.h.pca, display="species", choices=c(1), scaling=1),
     scores(mite.spe.h.pca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.h.pca, display="species", scaling=1)),
     col="red", cex=0.8)


##Section: 07-correspondence-analysis.R 

#Run the CA using the cca() function (NB: cca() is used for both CA and CCA)
spe.ca <- cca(spe)

# Identify the significant axes
ev<-spe.ca$CA$eig
ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

summary(spe.h.pca) #overall results
summary(spe.h.pca, diplay=NULL)# only axis eigenvalues and contribution

par(mfrow=c(1,2))
##scaling 1
plot(spe.ca, scaling=1, type="none", main='CA - biplot scaling 1', xlab=c("CA1 (%)", round((spe.ca$CA$eig[1]/sum(spe.ca$CA$eig))*100,2)),
ylab=c("CA2 (%)", round((spe.ca$CA$eig[2]/sum(spe.ca$CA$eig))*100,2)))

points(scores(spe.ca, display="sites", choices=c(1,2), scaling=1), pch=21, col="black", bg="steelblue", cex=1.2)

text(scores(spe.ca, display="species", choices=c(1), scaling=1),
     scores(spe.ca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spe.ca, display="species", scaling=1)),col="red", cex=0.8)

##scaling 2
plot(spe.ca, scaling=1, type="none", main='CA - biplot scaling 2', xlab=c("CA1 (%)", round((spe.ca$CA$eig[1]/sum(spe.ca$CA$eig))*100,2)),
     ylab=c("CA2 (%)", round((spe.ca$CA$eig[2]/sum(spe.ca$CA$eig))*100,2)), ylim=c(-2,3))

points(scores(spe.ca, display="sites", choices=c(1,2), scaling=2), pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(spe.ca, display="species", choices=c(1), scaling=2),
     scores(spe.ca, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(spe.ca, display="species", scaling=2)),col="red", cex=0.8)

# CA on mite species
mite.spe.ca<-cca(mite.spe)

#What are the significant axes ?
ev<-mite.spe.ca$CA$eig
ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

#Output summary/results
summary(mite.spe.ca, display=NULL)

#Plot the biplot
windows()
plot(mite.spe.ca, scaling=1, type="none",
     xlab=c("PC1 (%)", round((mite.spe.ca$CA$eig[1]/sum(mite.spe.ca$CA$eig))*100,2)),
     ylab=c("PC2 (%)", round((mite.spe.ca$CA$eig[2]/sum(mite.spe.ca$CA$eig))*100,2)))
points(scores(mite.spe.ca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.ca, display="species", choices=c(1), scaling=1),
     scores(mite.spe.ca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.ca, display="species", scaling=1)),
     col="red", cex=0.8)


##Section: 08-principal-coordinate-analysis.R 

#Using cmdscale()
?cmdscale
cmdscale(dist(spe.hel), k=(nrow(spe)-1), eig=TRUE)

#Using pcoa()
?pcoa
spe.h.pcoa <- pcoa(dist(spe.hel))

# Extract the results
spe.h.pcoa

#Construct the biplot
biplot.pcoa(spe .h.pcoa, spe.hel, dir.axis2=-1)

spe.bray.pcoa <- pcoa(spe.db) #where spe.db is the species dissimilarity matrix using Bray-Curtis.
spe.bray.pcoa
biplot.pcoa(spe.bray.pcoa, spe.hel, dir.axis2=-1)
#Note that the distance measure chosen strongly influences the results.

mite.spe.h.pcoa <- pcoa(dist(mite.spe.hel))
mite.spe.h.pcoa
windows()
biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel, dir.axis2=-1)


##Section: 09-non-metric-multidimensional-scaling.R 

# Run the NMDS
spe.nmds<-metaMDS(spe, distance='bray', k=2)

#Extract the results
spe.nmds

#Assess the goodness of fit and draw a Shepard plot
spe.nmds$stress
stressplot(spe.nmds, main='Shepard plot')

# Construct the biplot
windows()
plot(spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(spe.nmds$stress, 3)),
     xlab=c("NMDS1"),
     ylab=c("NMDS2"))
points(scores(spe.nmds, display="sites", choices=c(1,2)),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(spe.nmds, display="species", choices=c(1)),
     scores(spe.nmds, display="species", choices=c(2)),
     labels=rownames(scores(spe.nmds, display="species")),
     col="red", cex=0.8)

mite.spe.nmds<-metaMDS(mite.spe, distance='bray', k=2)
#Extract the results
mite.spe.nmds

#Assess the goodness of fit
mite.spe.nmds$stress
stressplot(mite.spe.nmds, main='Shepard plot')

#Construct the biplot
windows()
plot(mite.spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(mite.spe.nmds$stress, 3)),
     xlab=c("NMDS1"),
     ylab=c("NMDS2"))
points(scores(mite.spe.nmds, display="sites", choices=c(1,2)),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.nmds, display="species", choices=c(1)),
     scores(mite.spe.nmds, display="species", choices=c(2)),
     labels=rownames(scores(mite.spe.nmds, display="species")),
     col="red", cex=0.8)


##Section: 10-final-considerations.R 




##Section: 11-references.R 




