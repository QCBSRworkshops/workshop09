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


library (ade4)
data (doubs)

spe <- doubs$fish
env <- doubs$env

library (codep)
data (Doubs)

spe <- Doubs.fish
env <- Doubs.env

list.of.packages <- c("ape", "ade4", "codep", "gclus", "vegan", "GGally", "PlaneGeometry", "remotes", "matlib")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE) 
  print(paste0("The following package was installed:", new.packages)) 
} else if(length(new.packages) == 0) {
  print("All packages were already installed previously")
}

# Load all required libraries at once
lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE)

# source(file.choose()) # use coldiss.R which you have downloaded to your own directory


##Section: 02-introduction.R 




##Section: 03-data-exploration.R 

spe <- read.csv("data/doubsspe.csv", 
                row.names = 1) 

env <- read.csv("data/doubsenv.csv", 
                row.names = 1)

library (ade4)
data (doubs)

spe <- doubs$fish
env <- doubs$env

library (codep)
data (Doubs)

spe <- Doubs.fish
env <- Doubs.env

head(spe)[, 1:8]

str(spe)

# Try some of these!

names(spe)   # names of objects
dim(spe)     # dimensions

str(spe)     # structure of objects
summary(spe) # summary statistics

head(spe)    # first 6 rows

str(env)

summary(env) # summary statistics


##Section: 04-association-distances.R 

dist(spe)

class(dist(spe))

str(dist(spe))

as.matrix(dist(spe))

dim(as.matrix(dist(spe)))

as.matrix(dist(spe))[1:3, 1:3]

#setting up the plot
xlim <- c(0,4)
ylim <- c(-1,5)
par(mar=c(1,1,1,1)+.1)
plot(xlim, 
     ylim, type="n", 
     xlab="X1", 
     ylab="X2", 
     asp=1)
grid()
# define some vectors
a=c(4, 0)
b=c(0, 4)
# plot the vectors
vectors(a, labels="D(y21, y11)", pos.lab=3, frac.lab=.5, col="grey")
vectors(a + b, labels="D1(x1,x2)", pos.lab=4, frac.lab=.5, col="red")
# vector a+b starting from a is equal to b.
vectors(a + b, labels="D(y12, y22)", pos.lab=4, frac.lab=.5, origin=a, col="grey")

points(x = 4, y = 0, type = "p")
text(x=4.1, y=-0.2, labels="")

points(x = 0, y = 0, type = "p")
text(x=-0.1, y=-0.2, labels="x1")

points(x = 4, y = 4, type = "p")
text(x=4.1, y=4.2, labels="x2")

spe.D.Euclid <- dist(x = spe,  
                     method = "euclidean")

is.euclid(spe.D.Euclid)

Y.hmm <- data.frame(
  y1 = c(0, 0, 1),
  y2 = c(4, 1, 0),
  y3 = c(8, 1, 0))

Y.hmm.DistEu <- dist(x = Y.hmm,  
                     method = "euclidean")

as.matrix(Y.hmm.DistEu)

spe.D.Ch <- vegdist(spe,  
                    method = "chord")

as.matrix(spe.D.Ch)[1:3, 1:3]

Y.hmm.DistCh <- vegdist(Y.hmm,  
                    method = "chord")

as.matrix(Y.hmm.DistCh)

as.matrix(Y.hmm.DistEu)

Y.hmm

spe.D.Jac <- vegdist(spe, 
                     method = "jaccard",
                     binary = TRUE)

spe.D.Sor <- vegdist(spe, 
                     method = "bray",
                     binary = TRUE)

spe.db.pa <- vegdist(spe, 
                      method = "bray",
                      binary = FALSE)
spe.db <- as.matrix(spe.db.pa)

# coldiss() function
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012
#

"coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
  require(gclus)
  
  if (max(D)>1) D <- D/max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1-D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow=c(1,2), pty="s")
  
  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank= TRUE		equal-sized classes
# byrank= FALSE		equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Example:
# coldiss(spe.dj, nc=9, byrank=F, diag=T)


coldiss(spe.D.Jac)

spe[1:6, 1:6]

spe.pa <- decostand(spe, method = "pa")
spe.pa[1:6, 1:6]

spe.total <- decostand(spe, 
                       method = "total")
spe.total[1:5, 1:6]

spe.total <- decostand(spe, 
                       method = "hellinger")
spe.total[1:5, 1:6]

?decostand
env.z <- decostand(env, method = "standardize")

apply(env.z, 2, mean)
apply(env.z, 2, sd)

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

# Demonstration of a cluster dendrogram
spe.hel<-decostand(spe, method="hellinger")
spe.dhel <- vegdist(spe.hel,method="euclidean")
spe.dhel.ward <- hclust(spe.dhel, method="ward.D2")
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward, hang=-1) # hang=-1 aligns all objets on the same line


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




