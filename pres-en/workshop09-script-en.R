list.of.packages <- c("ape", "gclus", "vegan", "GGally", "PlaneGeometry", "remotes")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE) 
  print(paste0("The following package was installed:", new.packages)) 
} else if(length(new.packages) == 0) {
    print("All packages were already installed previously")
}

# Load all required libraries at once
lapply(list.of.packages, require, character.only = TRUE, quietly = TRUE)

install.packages(c("ape",
                   "gclus",
                   "vegan",
                   "GGally",
                   "PlaneGeometry",
                   "remotes"))

spe <- read.csv("data/doubsspe.csv", row.names = 1) # this line will vary depending on where you saved the dataset
spe <-  spe[-8,] # remove site with no data
# for larger datasets:
row_sub = apply(spe, 1, function(row) any(row !=0 ))
# find rows where there are not any values that equal zero
spe <- spe[row_sub,]
# only keep rows that have no zeros  

env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # remove site with no data

names(spe) # Names of objects
dim(spe) # dimensions
str(spe) # structure of objects
summary(spe) # summary statistics
head(spe) # first 6 rows

head(spe) # first 6 rows

par(mar = c(4,4,.5,.5), cex = 1.5)
ab <- table(unlist(spe))
barplot(ab, las = 1, col = grey(5:0/5),
        xlab = "Abundance class", ylab = "Frequency")

sum(spe == 0)

sum(spe == 0)/(nrow(spe)*ncol(spe))

par(mar = c(4,4,1,.5), cex = 1.5)
site.pre <- rowSums(spe > 0)
barplot(site.pre, main = "Species richness",
        xlab = "Sites", ylab = "Number of species",
        col = "grey ", las = 1)

library(vegan)
spec.pa <- decostand(spe, method = "pa")

spec.hel <- decostand(spe, method = "hellinger")
spec.chi <- decostand(spe, method = "chi.square")

spe.pa <- decostand(spe,method = "log")

names(env) # Names of objects
dim(env) # dimensions
str(env) # structure of objects
summary(env) # summary statistics
head(env) # first 6 rows

head(env) # first 6 rows



library(GGally)
ggpairs(env, main="Bivariate Plots of the Environmental Data")

?decostand
env.z <- decostand(env, method = "standardize")

apply(env.z, 2, mean)
apply(env.z, 2, sd)

?vegdist

spe.db.pa <- vegdist(spe, method = "bray")
spe.db <- as.matrix(spe.db.pa)

#         1         2         3         4         5         6         7         9
# 1  0.0000000 0.6000000 0.6842105 0.7500000 0.8918919 0.7500000 0.6842105 1.0000000
# 2  0.6000000 0.0000000 0.1428571 0.3333333 0.6956522 0.3939394 0.1428571 0.6923077
# 3  0.6842105 0.1428571 0.0000000 0.1891892 0.6800000 0.2972973 0.1250000 0.7333333
#...


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

# the code for the coldiss() function is in the workshop script.
coldiss(spe.db.pa)

# Demonstration of a cluster dendrogram
spe.hel<-decostand(spe, method="hellinger")
spe.dhel<-vegdist(spe.hel,method="euclidean")
spe.dhel.ward<-hclust(spe.dhel, method="ward.D2")
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward, hang=-1) # hang=-1 aligns all objets on the same line


x <- matrix(abs(rnorm(1:40)),nrow=5,byrow=TRUE)
x.hel <- decostand(x,method="hellinger")
x.dhel <- vegdist(x.hel, method="euclidean")
x.dhel.single <- hclust(x.dhel, method="single")
x.dhel.complete <- hclust(x.dhel, method="complete")
plot <- par(mfrow=c(1,2),mar=c(5,2,4,2))
plot(as.dendrogram(x.dhel.single), horiz=TRUE, main="Single linkage clustering")
plot(as.dendrogram(x.dhel.complete), horiz=TRUE, main = "Complete linkage clustering")

par(mar=c(.5,3.8,2,.5), cex = 1.5)
spe.dhe1 <- vegdist(spec.hel, method = "euclidean")
spe.dhe1.single <- hclust(spe.dhe1, method = "single")
plot(spe.dhe1.single)

par(mfrow=c(1,2), mar=c(.5,2.5,1.5,2.5), cex=1)
spe.dhe1 <- vegdist(spec.hel, method = "euclidean")
spe.dhe1.complete <- hclust(spe.dhe1, method = "complete")
plot(spe.dhe1.single, main="Single linkage clustering", hang =-1)
plot(spe.dhe1.complete, main="Complete linkage clustering", hang=-1)

par(mar=c(.5,2.5,1.5,.5), cex = 1)
spe.dhel.ward <- hclust(spe.dhe1, method = "ward.D2")
spe.dhel.ward$height <- sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward, hang = -1) # hang = -1 aligns objects at the same level

par(mar=c(.5,3.8,2,.5), cex = 1.5)
plot(spe.dhel.ward, hang = -1) # hang = -1 aligns objects at the same level

x <- rnorm(5000, mean = 0, sd = 1)
y <- rnorm(5000, mean = 0, sd = 1)
z <- rnorm(5000, mean = 0, sd = 1)

xyz <- data.frame(x, y, z)

GGally::ggpairs(xyz) 

cov(xyz)

round(cov(xyz), digits = 5)

library(PlaneGeometry)

P <- c(0, 0) 
w <- c(1, 0) 

ratio <- 1 
angle <- 30

shear <- Shear$new(P, 
                   w, 
                   ratio, 
                   angle)

wt <- ratio * c(-w[2], w[1])

Q <- P + w
R <- Q + wt
S <- P + wt
A <- shear$transform(P)
B <- shear$transform(Q)
C <- shear$transform(R)
D <- shear$transform(S)


plot(0, 0, type = "n", asp = 1, xlim = c(0,1), ylim = c(0,2))

lines(rbind(P, Q, R, S, P), 
      lwd = 2) # unit square

lines(rbind(A, B, C, D, A), 
      lwd = 2, 
      col = "blue") # image by the shear

arrows(x0 = A[1], 
       y0 = A[2], 
       x1 = B[1],
       y1 = B[2],
       col = "red",
       lwd = 2)

arrows(x0 = A[1], 
       y0 = A[2], 
       x1 = D[1],
       y1 = D[2],
       col = "purple",
       lwd = 2)

data(varechem)

str(varechem)

data(varechem)

# Step 1 
Y <- varechem[, 1:2]

head(Y)

data(varechem)

# Step 1 
Y <- varechem[, 1:2]

# Step 2
Y_std <- as.matrix(scale(Y))

head(Y_std)

round(apply(Y_std, 2, mean))
round(apply(Y_std, 2, sd))


data(varechem)

# Step 1 
Y <- varechem[, 1:2]

# Step 2
Y_std <- as.matrix(scale(Y))

# Step 3
(Y_R <- cov(Y_std))

data(varechem)

# Step 1 
Y <- varechem[, 1:2]

# Step 2
Y_std <- as.matrix(scale(Y))

# Step 3
Y_R <- cov(Y_std)

# Step 4 
(Eigenvalues <- eigen(Y_R)$values)

(Eigenvectors <- eigen(Y_R)$vectors)

Eigenvectors. <- as.data.frame(Eigenvectors)
row.names(Eigenvectors.) <- c("P", "N")
colnames(Eigenvectors.) <- c("PC1", "PC2")
Eigenvectors.[, 1] <- Eigenvectors.[, 1]*-1
Y_std <- as.data.frame(Y_std)

plot(N ~ P, 
     col = as.factor(rownames(Y_std)), 
     main="Distances to PC1", 
     pch = 19, 
     xlim=c(-2.2, 2.2), 
     ylim = c(-2.2,2.2), 
     data = as.data.frame(Y_std))

abline(v=0 , h=0, 
       col = "dark gray")

#Overlap pertinent evectors

abline(0, 
       Eigenvectors[2, 1]/Eigenvectors[1, 1],
       col='purple')
# abline(0, 
#        Eigenvectors[1, 2]/Eigenvectors[2, 2],
#        col='orange')

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[1]*Eigenvectors[1, 1],
       y1 = Eigenvalues[1]*Eigenvectors[2, 1],
       col = "purple",
       lwd = 2)

# arrows(x0 = 0, 
#        y0 = 0, 
#        x1 = Eigenvalues[2]*Eigenvectors[1,2], 
#        y1 = Eigenvalues[2]*Eigenvectors[2, 2],
#        col = "orange", 
#        lwd = 2)

# Plot the lines from first evector to points

line1 <- c(0, 
           Eigenvectors[2, 1]/Eigenvectors[1, 1])

perp.segment.coord <- function(x0, y0, line1){
  #finds endpoint for a perpendicular segment from the point (x0,y0) to the line1
  a <- line1[1]  #intercept
  b <- line1[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}

ss <- perp.segment.coord(Y_std$P, 
                         Y_std$N, 
                         line1)
# do.call(segments, ss)
# which is the same as:

segments(x0 = ss$x0, 
         x1 = ss$x1, 
         y0 = ss$y0, 
         y1 = ss$y1, 
         col = 'purple')

points(N ~ P, 
       col = as.factor(rownames(Y_std)), 
       pch = 19,
       data = Y_std)
with(Y_std,
     text(N ~ P, 
          labels = as.factor(rownames(Y_std)),
                     pos = 1, 
          cex=1.4))

Eigenvectors. <- as.data.frame(Eigenvectors)
row.names(Eigenvectors.) <- c("P", "N")
colnames(Eigenvectors.) <- c("PC1", "PC2")
Eigenvectors.[, 1] <- Eigenvectors.[, 1]*-1
Y_std <- as.data.frame(Y_std)

plot(N ~ P, 
     col = as.factor(rownames(Y_std)), 
     main="Distances to PC2", 
     pch = 19, 
     xlim=c(-2.2, 2.2), 
     ylim = c(-2.2,2.2), 
     data = as.data.frame(Y_std))

abline(v=0 , h=0, 
       col = "dark gray")

#Overlap pertinent evectors

abline(0, 
       Eigenvectors[2, 1]/Eigenvectors[1, 1],
       col='purple')
abline(0, 
       Eigenvectors[1, 2]/Eigenvectors[2, 2],
       col='orange')

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[1]*Eigenvectors[1, 1],
       y1 = Eigenvalues[1]*Eigenvectors[2, 1],
       col = "purple",
       lwd = 2)

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[2]*Eigenvectors[1,2], 
       y1 = Eigenvalues[2]*Eigenvectors[2, 2],
       col = "orange", 
       lwd = 2)


line2 <- c(0, 
           Eigenvectors[1, 2]/Eigenvectors[1, 1])

perp.segment.coord <- function(x0, y0, line2){
  a <- line2[1]  #intercept
  b <- line2[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}

ss <- perp.segment.coord(Y_std$P, 
                         Y_std$N, 
                         line2)

segments(x0 = ss$x0, 
         x1 = ss$x1, 
         y0 = ss$y0, 
         y1 = ss$y1, 
         col = 'orange')

points(N ~ P, 
       col = as.factor(rownames(Y_std)), 
       pch = 19,
       data = Y_std)

with(Y_std,
     text(N ~ P, 
          labels = as.factor(rownames(Y_std)),
                     pos = 1, 
          cex=1.4)
     )

sum(diag(cov(Y_std)))
sum(eigen(cov(Y_std))$values)

# Step 1 
Y <- varechem[, 1:2]

# Step 2
Y_std <- as.matrix(scale(Y))

# Step 3
Y_R <- cov(Y_std)

# Step 4 
Eigenvalues <- eigen(Y_R)$values
Eigenvectors <- eigen(Y_R)$vectors

# Step 5
F_PrComps <- Y_std %*% Eigenvectors
head(F_PrComps)

Eigenvectors. <- as.data.frame(Eigenvectors)
row.names(Eigenvectors.) <- c("P", "N")
colnames(Eigenvectors.) <- c("PC1", "PC2")
Eigenvectors.[, 1] <- Eigenvectors.[, 1]*-1
Y_std <- as.data.frame(Y_std)

op <- par(mfrow = c(2, 1),     # 2x2 layout
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0)    # axis label at 2 rows distance, tick labels at 1 row
    )       

plot(N ~ P, 
     col = as.factor(rownames(Y_std)),
     pch = 19, 
     xlim=c(-2.2, 2.2), 
     ylim = c(-2.2,2.2), 
     data = as.data.frame(Y_std))

abline(v=0 , h=0, 
       col = "dark gray")

#Overlap pertinent evectors

abline(0, 
       Eigenvectors[2, 1]/Eigenvectors[1, 1],
       col='purple')
# abline(0, 
#        Eigenvectors[1, 2]/Eigenvectors[2, 2],
#        col='orange')

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[1]*Eigenvectors[1, 1],
       y1 = Eigenvalues[1]*Eigenvectors[2, 1],
       col = "purple",
       lwd = 2)

# arrows(x0 = 0, 
#        y0 = 0, 
#        x1 = Eigenvalues[2]*Eigenvectors[1,2], 
#        y1 = Eigenvalues[2]*Eigenvectors[2, 2],
#        col = "orange", 
#        lwd = 2)

# Plot the lines from first evector to points

line1 <- c(0, 
           Eigenvectors[2, 1]/Eigenvectors[1, 1])

perp.segment.coord <- function(x0, y0, line1){
  a <- line1[1]  #intercept
  b <- line1[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}

ss <- perp.segment.coord(Y_std$P, 
                         Y_std$N, 
                         line1)
# do.call(segments, ss)
# which is the same as:

segments(x0 = ss$x0, 
         x1 = ss$x1, 
         y0 = ss$y0, 
         y1 = ss$y1, 
         col = 'purple')

points(N ~ P, 
       col = as.factor(rownames(Y_std)), 
       pch = 19,
       data = Y_std)
with(Y_std,
     text(N ~ P, 
          labels = as.factor(rownames(Y_std)),
                     pos = 1, 
          cex=1.4))


plot(N ~ P, 
     col = as.factor(rownames(Y_std)),
     pch = 19, 
     xlim=c(-2.2, 2.2), 
     ylim = c(-2.2,2.2), 
     data = as.data.frame(Y_std))

abline(v=0 , h=0, 
       col = "dark gray")

#Overlap pertinent evectors

abline(0, 
       Eigenvectors[2, 1]/Eigenvectors[1, 1],
       col='purple')
abline(0, 
       Eigenvectors[1, 2]/Eigenvectors[2, 2],
       col='orange')

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[1]*Eigenvectors[1, 1],
       y1 = Eigenvalues[1]*Eigenvectors[2, 1],
       col = "purple",
       lwd = 2)

arrows(x0 = 0, 
       y0 = 0, 
       x1 = Eigenvalues[2]*Eigenvectors[1,2], 
       y1 = Eigenvalues[2]*Eigenvectors[2, 2],
       col = "orange", 
       lwd = 2)


line2 <- c(0, 
           Eigenvectors[1, 2]/Eigenvectors[1, 1])

perp.segment.coord <- function(x0, y0, line2){
  a <- line2[1]  #intercept
  b <- line2[2]  #slope
  x1 <- (x0 + b * y0 - a * b)/(1 + b^2)
  y1 <- a + b * x1
  list(x0 = x0, y0 = y0, 
       x1 = x1, y1 = y1)
}

ss <- perp.segment.coord(Y_std$P, 
                         Y_std$N, 
                         line2)

segments(x0 = ss$x0, 
         x1 = ss$x1, 
         y0 = ss$y0, 
         y1 = ss$y1, 
         col = 'orange')

points(N ~ P, 
       col = as.factor(rownames(Y_std)), 
       pch = 19,
       data = Y_std)

with(Y_std,
     text(N ~ P, 
          labels = as.factor(rownames(Y_std)),
                     pos = 1, 
          cex=1.4)
     )

title(xlab = "N",
      ylab = "P",
      outer = TRUE, line = 3)

par(op)

score <- as.data.frame(F_PrComps)

colnames(score) <- c("PC1", "PC2")

op <- par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
          mar = c(1, 1, 0, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0)    # axis label at 2 rows distance, tick labels at 1 row
)


plot(PC2 ~ PC1, 
     col = as.factor(rownames(score)), 
     pch = 19, 
     xlim = c(-2.2, 2.2), ylim = c(-2.2,2.2), xlab='PC1', ylab='PC2',data = score)

abline(h = 0, col = 'purple')
abline(v = 0, col='orange')


perp.segment.horiz <- function(x0, y0){
  x1 <- x0
  y1 <- 0
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}

ss1 <- perp.segment.horiz(score[,1], score[,2])

segments(x0 = ss1$x0, x1 = ss1$x1, y0 = ss1$y0, y1 = ss1$y1, col='purple')


points(PC2 ~ PC1, col=as.factor(rownames(score)), pch = 19, xlab='V1', ylab='V2',data=score)
with(score,text(PC2 ~ PC1, labels=as.factor(rownames(score)), pos = 3, cex=1.4))


plot(PC2 ~ PC1, col=as.factor(rownames(score)), 
     pch = 19, xlim=c(-2.2, 2.2), ylim = c(-2.2,2.2),
     xlab='PC1', ylab='PC2',data=score)

abline(h = 0, col = 'purple')
abline(v = 0, col ='orange')


perp.segment.vert <- function(x0, y0){
  x1 <- 0
  y1 <- y0
  
  list(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}

ss1a <- perp.segment.vert(score[,1], score[,2])
segments(x0 = ss1a$x0, x1 = ss1a$x1, y0 = ss1a$y0, y1 = ss1a$y1, col='orange')


points(PC2 ~ PC1, col=as.factor(rownames(score)), pch = 19, xlab='V1', ylab='V2',data=score)

with(score,text(PC2 ~ PC1, labels=as.factor(rownames(score)), pos = 3, cex=1.4))

title(xlab = "PC1",
      ylab = "PC2",
      outer = TRUE, 
      line = 3)

par(op)

data(varechem)

Y <- varechem[, 1:2] 
Y_std <- as.matrix(scale(Y))
Y_R <- cov(Y_std)

Eigenvalues <- eigen(Y_R)$values
Eigenvectors <- eigen(Y_R)$vectors

F_PrComps <- Y_std %*% Eigenvectors

head(F_PrComps)

PCA_prcomp <- prcomp(Y, 
                     center = TRUE, 
                     scale = TRUE)

# or PCA_prcomp <- prcomp(Y_std)

head(PCA_prcomp$x)

data(varechem)

Y <- varechem[, 1:2] 
Y_std <- as.matrix(scale(Y))
Y_R <- cov(Y_std)

Eigenvalues <- eigen(Y_R)$values
Eigenvectors <- eigen(Y_R)$vectors

F_PrComps <- Y_std %*% Eigenvectors

head(F_PrComps)

PCA_princomp <- princomp(Y_std)

head(PCA_princomp$scores)

data(varechem)

Y <- varechem[, 1:2] 
Y_std <- as.matrix(scale(Y))
Y_R <- cov(Y_std)

Eigenvalues <- eigen(Y_R)$values
Eigenvectors <- eigen(Y_R)$vectors

F_PrComps <- Y_std %*% Eigenvectors

head(F_PrComps)

PCA_vegan_rda <- rda(Y_std)

scores(PCA_vegan_rda, 
       display = "sites", 
       scaling = 1,
       choices = seq_len(PCA_vegan_rda$CA$rank),
       const = sqrt(PCA_vegan_rda$tot.chi * (nrow(PCA_vegan_rda$CA$u) - 1)))[1:5, ]

spe.h.pca <- rda(spec.hel)

# summary(spe.h.pca)

paste(capture.output(summary(spe.h.pca))[5:8])

paste(capture.output(summary(spe.h.pca))[c(12:16, 21:24)])

sum(spe.h.pca$CA$eig)

paste(capture.output(summary(spe.h.pca))[c(26:29, 31:32, 34:40, 63:64, 66:72)])

paste(capture.output(summary(spe.h.pca))[c(32, 34:40)])

paste(capture.output(summary(spe.h.pca))[c(64, 66:72)])

scores(spe.h.pca,
       display = "species" or "sites")

ev <- spe.h.pca$CA$eig
# ev[ev > mean(ev)]

par(mar=c(4,4,2.5,.5), cex = 1.5)
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue",
       lwd = 2, col = "red3" , bty = "n")

head(bstick(spe.h.pca))

screeplot(spe.h.pca, 
          bstick = TRUE, type = "lines")

env.pca <- rda(env.z)
# summary(env.pca, scaling  = 2)

ev <- env.pca$CA$eig

ev[ev>mean(ev)]

par(mar=c(4,4,2.5,.5), cex = 1.5)
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue",
       lwd = 2, col = "red3" , bty = "n")

par(mar=c(4,4, 0.1,0.1), cex = 1.5)
plot(spe.h.pca)

par(mar = c(4,4,0.05,0.05), cex = 1.2)
biplot(spe.h.pca)

par(mar = c(4,4,0.05,0.05), cex = 1.2)
biplot(spe.h.pca, scaling = 1)

par(mar = c(4,4,0.05,0.05), cex = 1.2)
biplot(spe.h.pca, scaling = 2)

data(mite)

mite.spe.hel <- decostand(mite, 
                          method = "hellinger")

mite.spe.h.pca <- rda(mite.spe.hel)

ev <- mite.spe.h.pca$CA$eig
ev[ev>mean(ev)]
n <- length(ev)
barplot(ev, main = "Eigenvalues",
        col = "grey", las = 2)
abline(h = mean(ev),
       col = "red3", lwd = 2)
legend("topright",
       "Average eigenvalue",
       lwd = 2,
       col = "red3", bty = "n")

par(mar=c(4,4,2,1), cex = 1.2)
ev <- mite.spe.h.pca$CA$eig
n <- length(ev)
barplot(ev, main = "Eigenvalues", 
        col = "grey", las = 2)
abline(h = mean(ev),
       col = "red3", lwd = 2)
legend("topright", 
       "Average eigenvalue", 
       lwd = 2, col = "red3", 
       bty = "n")

par(mar = c(4,4,0.05,0.05), cex = 1.5)
biplot(mite.spe.h.pca, 
       col = c("red3", "grey15"))

library(ape)
spe.h.pcoa <- pcoa(dist(spec.hel))
summary(spe.h.pcoa)

head(spe.h.pcoa$values)

head(spe.h.pcoa$vectors)[, 1:5]

biplot.pcoa(spe.h.pcoa, spec.hel)

spe.bray.pcoa <- pcoa(spe.db.pa)

spe.bray.pcoa$values$Eigenvalues

spe.bray.pcoa <- pcoa(spe.db.pa, 
                      correction = "cailliez")

spe.bray.pcoa$values$Corr_eig

par(mar=c(3,3,.5,1), cex = 1.2)
biplot.pcoa(spe.bray.pcoa)

mite.spe <- mite
mite.spe.hel <- decostand(mite.spe, method = "hellinger")

mite.spe.h.pcoa <- pcoa(dist(mite.spe.hel))

biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel)

spe.nmds <- metaMDS(spe, distance = 'bray', k = 2)

spe.nmds <- metaMDS(spe, distance = 'bray', k = 2)

spe.nmds$stress
stressplot(spe.nmds, main = "Shepard plot")

plot(spe.nmds, type = "none",
     main = paste("NMDS/Bray - Stress =",
                  round(spe.nmds$stress, 3)),
     xlab = c("NMDS1"), ylab = "NMDS2")

points(scores(spe.nmds, display = "sites",
              choiches = c(1,2),
              pch = 21,
              col = "black",
              g = "steelblue",
              cex = 1.2))
text(scores(spe.nmds, display = "species", choices = c(1)),
            scores(spe.nmds, display = "species", choices = c(2)),
            labels = rownames(scores(spe.nmds, display = "species")),
            col = "red", cex = 0.8)

plot(spe.nmds, type = "none",
     main = paste("NMDS/Bray - Stress =",
                  round(spe.nmds$stress, 3)),
     xlab = c("NMDS1"), ylab = "NMDS2")

points(scores(spe.nmds, display = "sites",
              choiches = c(1,2),
              pch = 21,
              col = "black",
              g = "steelblue",
              cex = 1.2))
text(scores(spe.nmds, display = "species", choices = c(1)),
            scores(spe.nmds, display = "species", choices = c(2)),
            labels = rownames(scores(spe.nmds, display = "species")),
            col = "red", cex = 0.8)

plot(spe.h.pca, scaling  = 1,
     type = "none",
     xlab = c("PC1 (%)",
              round(spe.h.pca$CA$eig[1]/sum(spe.h.pca$CAeig)*100,2)),
     ylab = c("PC2 (%)", round(spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig)*100,2)))
points(scores(spe.h.pca, display = "sites", choices = c(1,2), scaling = 1),
       pch=21, col = "black", bg = "steelblue" , cex  = 1.2)
text(scores(spe.h.pca, display = "species", choices = 1, scaling = 1),
     scores(spe.h.pca, display = "species", choices = 2, scaling = 1),
     labels = rownames(scores(spe.h.pca, display = "species", scaling = 1)),
     col = "red", cex = 0.8)
spe.cs <- scores(spe.h.pca, choices = 1:2, scaling = 1 , display = "sp")
arrows(0, 0, spe.cs[,1], spe.cs[,2], length = 0)

par(mar = c(4,4,0.05,0.05), cex = 1.2)
plot(spe.h.pca, scaling  = 1, type = "none",
     xlab = c("PC1 (%)", round(spe.h.pca$CA$eig[1]/sum(spe.h.pca$CAeig)*100,2)),
     ylab = c("PC2 (%)", round(spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig)*100,2)))
points(scores(spe.h.pca, display = "sites", choices = c(1,2), scaling = 1),
       pch=21, col = "black", bg = "steelblue" , cex  = 1.2)
text(scores(spe.h.pca, display="species", choices=c(1), scaling = 1),
     scores(spe.h.pca, display = "species", choices = c(2), scaling = 1),
     labels=rownames(scores(spe.h.pca, display = "species", scaling = 1)),
     col = "red", cex = 0.8)
spe.cs <- scores(spe.h.pca, choices = 1:2, scaling = 1 , display = "sp")
arrows(0,0,spe.cs[,1], spe.cs[,2], length = 0)

spe.ca <- cca(spe[-8,])
# only take columns which rowsums are > than 0.


summary(spe.ca)

mite.spe <- mite

mite.spe.ca <- cca(mite.spe)

ev <- mite.spe.ca$CA$eig
ev[ev > mean(ev)]
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue", lwd = 2, col = red3, bty = "n")

par(mar=c(4,4,2,1), cex = 1.2)
ev <- mite.spe.ca$CA$eig
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue", lwd = 2, col = "red3", bty = "n")
