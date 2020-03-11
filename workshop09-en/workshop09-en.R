## ----setup2, echo = FALSE-----------------------------------------------------
## Setup for your presentation
knitr::opts_chunk$set(
  eval = TRUE,
  cache = TRUE,
  comment = "#",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width=5, fig.height=5, fig.retina=3,
  fig.align = 'center'
)
options(repos=structure(c(CRAN="http://cran.r-project.org")))


## ----install_pkgs, echo = FALSE, results = "asis"-----------------------------
cat(
  qcbsRworkshops::first_slides(9, c('ape', 'gclus', 'vegan'))
)


## -----------------------------------------------------------------------------
spe <- read.csv("data/doubsspe.csv", row.names = 1)
spe <-  spe[-8,] # remove site with no data


## -----------------------------------------------------------------------------
env <- read.csv("data/doubsenv.csv", row.names = 1)
env <- env[-8,] # remove site with no data


## ---- eval = F----------------------------------------------------------------
## names(spe) # Names of objects
## dim(spe) # dimensions
## str(spe) # structure of objects
## summary(spe) # summary statistics
## head(spe) # first 6 rows


## ---- echo = F----------------------------------------------------------------
head(spe) # first 6 rows


## ---- fig.width=8, echo = -1--------------------------------------------------
par(mar = c(4,4,.5,.5), cex = 1.5)
ab <- table(unlist(spe))
barplot(ab, las = 1, col = grey(5:0/5),
        xlab = "Abundance class", ylab = "Frequency")


## -----------------------------------------------------------------------------
sum(spe == 0)


## -----------------------------------------------------------------------------
sum(spe == 0)/(nrow(spe)*ncol(spe))


## ---- fig.width=10, fig.height=5, echo=-1-------------------------------------
par(mar = c(4,4,1,.5), cex = 1.5)
site.pre <- rowSums(spe > 0)
barplot(site.pre, main = "Species richness",
        xlab = "Sites", ylab = "Number of species",
        col = "grey ", las = 1)


## -----------------------------------------------------------------------------
library(vegan)
spec.pa <- decostand(spe, method = "pa")


## -----------------------------------------------------------------------------
spec.hel <- decostand(spe, method = "hellinger")
spec.chi <- decostand(spe, method = "chi.square")


## -----------------------------------------------------------------------------
spe.pa <- decostand(spe,method = "log")


## ---- eval = F----------------------------------------------------------------
## names(env) # Names of objects
## dim(env) # dimensions
## str(env) # structure of objects
## summary(env) # summary statistics
## head(env) # first 6 rows


## ---- echo = T----------------------------------------------------------------
head(env) # first 6 rows


## ---- eval = F----------------------------------------------------------------
## pairs(env, main = "Bivariate Plots of the Environmental Data")


## ---- eval = -1---------------------------------------------------------------
?decostand
env.z <- decostand(env, method = "standardize")


## -----------------------------------------------------------------------------
apply(env.z, 2, mean)
apply(env.z, 2, sd)


## ---- eval = F----------------------------------------------------------------
## ?vegdist


## -----------------------------------------------------------------------------
spe.db.pa <- vegdist(spe, method = "bray")


## ---- fig.width=7, echo = -1--------------------------------------------------
par(mar=c(.5,3.8,2,.5), cex = 1.5)
spe.dhe1 <- vegdist(spec.hel, method = "euclidean")
spe.dhe1.single <- hclust(spe.dhe1, method = "single")
plot(spe.dhe1.single)


## ---- fig.width=9, echo = -1--------------------------------------------------
par(mar=c(.5,3.8,2,.5), cex = 1.5)
spe.dhel.ward <- hclust(spe.dhe1, method = "ward.D2")
spe.dhel.ward$height <- sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward, hang = -1) # hang = -1 aligns objects at the same level


## ---- fig.width=9, echo = F---------------------------------------------------
par(mar=c(.5,3.8,2,.5), cex = 1.5)
plot(spe.dhel.ward, hang = -1) # hang = -1 aligns objects at the same level


## -----------------------------------------------------------------------------
spe.h.pca <- rda(spec.hel)

summary(spe.h.pca)


## -----------------------------------------------------------------------------
summary(spe.h.pca, display = NULL)


## -----------------------------------------------------------------------------
eigen(cov(spec.hel))


## -----------------------------------------------------------------------------
spe.scores <- scores(spe.h.pca,
                     display = "species",
                     choices = c(1,2))


## -----------------------------------------------------------------------------
site.scores <- scores(spe.h.pca,
                      display = "sites",
                      choices = c(1,2))


## -----------------------------------------------------------------------------
ev <- spe.h.pca$CA$eig


## -----------------------------------------------------------------------------
ev[ev>mean(ev)]


## ---- echo = -1, fig.width=10, fig.height = 5.5-------------------------------
par(mar=c(4,4,2.5,.5), cex = 1.5)
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue",
       lwd = 2, col = "red3" , bty = "n")


## -----------------------------------------------------------------------------
env.pca <- rda(env.z)
summary(env.pca, scaling  = 2) # default


## -----------------------------------------------------------------------------
ev <- env.pca$CA$eig


## -----------------------------------------------------------------------------
ev[ev>mean(ev)]


## ---- echo = -1, fig.width=8, fig.height = 5----------------------------------
par(mar=c(4,4,2.5,.5), cex = 1.5)
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue",
       lwd = 2, col = "red3" , bty = "n")


## ---- echo = -1---------------------------------------------------------------
par(mar=c(4,4,.1,.1), cex = 1.5)
plot(spe.h.pca)


## ---- echo = -1, fig.height=4, fig.width=4.5----------------------------------
par(mar = c(4,4,0.05,0.05), cex = 1.2)
biplot(spe.h.pca)


## ---- eval = F----------------------------------------------------------------
## plot(spe.h.pca, scaling  = 1, type = "none",
##      xlab = c("PC1 (%)", round(spe.h.pca$CA$eig[1]/sum(spe.h.pca$CAeig)*100,2)),
##      ylab = c("PC2 (%)", round(spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig)*100,2)))
## points(scores(spe.h.pca, display = "sites", choices = c(1,2), scaling = 1),
##        pch=21, col = "black", bg = "steelblue" , cex  = 1.2)
## text(scores(spe.h.pca, display = "species", choices = 1, scaling = 1),
##      scores(spe.h.pca, display = "species", choices = 2, scaling = 1),
##      labels = rownames(scores(spe.h.pca, display = "species", scaling = 1)),
##      col = "red", cex = 0.8)
## spe.cs <- scores(spe.h.pca, choices = 1:2, scaling = 1 , display = "sp")
## arrows(0, 0, spe.cs[,1], spe.cs[,2], length = 0)


## ---- echo = F, fig.width=7,fig.height=7--------------------------------------
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


## ---- eval = F----------------------------------------------------------------
## install.packages("devtools")
## require("devtools")
## install_github("ggvegan", "gavinsimpson")
## require("ggvegan")
## autoplot()


## ---- eval = F----------------------------------------------------------------
## require(rgl)
## require(vegan3d)
## ordirgl(spe.h.pca)


## -----------------------------------------------------------------------------
data(mite)


## -----------------------------------------------------------------------------
mite.spe.hel <- decostand(mite, method = "hellinger")

mite.spe.h.pca <- rda(mite.spe.hel)


## ---- eval = F----------------------------------------------------------------
## ev <- mite.spe.h.pca$CA$eig
## ev[ev>mean(ev)]
## n <- length(ev)
## barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
## abline(h = mean(ev), col = "red3", lwd = 2)
## legend("topright", "Average eigenvalue", lwd = 2, col = "red3", bty = "n")


## ---- echo = F, fig.width=10, fig.height=7------------------------------------
par(mar=c(4,4,2,1), cex = 1.2)
ev <- mite.spe.h.pca$CA$eig
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue", lwd = 2, col = "red3", bty = "n")


## ---- echo = -1, fig.height=6.5, fig.width=7----------------------------------
par(mar = c(4,4,0.05,0.05), cex = 1.5)
biplot(mite.spe.h.pca, col = c("red3", "grey15"))


## -----------------------------------------------------------------------------
spe.ca <- cca(spe[-8,])
# only take columns which rowsums are > than 0.



## -----------------------------------------------------------------------------
summary(spe.ca)


## -----------------------------------------------------------------------------
mite.spe <- mite


## -----------------------------------------------------------------------------
mite.spe.ca <- cca(mite.spe)


## ---- eval = F----------------------------------------------------------------
## ev <- mite.spe.ca$CA$eig
## ev[ev > mean(ev)]
## n <- length(ev)
## barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
## abline(h = mean(ev), col = "red3", lwd = 2)
## legend("topright", "Average eigenvalue", lwd = 2, col = red3, bty = "n")


## ---- echo = F, fig.width=10, fig.height=7------------------------------------
par(mar=c(4,4,2,1), cex = 1.2)
ev <- mite.spe.ca$CA$eig
n <- length(ev)
barplot(ev, main = "Eigenvalues", col = "grey", las = 2)
abline(h = mean(ev), col = "red3", lwd = 2)
legend("topright", "Average eigenvalue", lwd = 2, col = "red3", bty = "n")


## ---- include = FALSE---------------------------------------------------------
library(ape)


## ---- eval = FALSE------------------------------------------------------------
## ?cmdscale
## library(ape)
## ?pcoa


## -----------------------------------------------------------------------------
spe.h.pcoa <- pcoa(dist(spec.hel))
summary(spe.h.pcoa)


## ---- fig.height=5, fig.width=8-----------------------------------------------
biplot.pcoa(spe.h.pcoa, spec.hel)


## -----------------------------------------------------------------------------
spe.bray.pcoa <- pcoa(spe.db.pa)
# spe.bray.pcoa


## -----------------------------------------------------------------------------
spe.bray.pcoa <- pcoa(spe.db.pa, correction = "cailliez")
# spe.bray.pcoa


## ---- fig.width=6, fig.height=5.5, echo = -1----------------------------------
par(mar=c(3,3,.5,1), cex = 1.2)
biplot.pcoa(spe.bray.pcoa)


## -----------------------------------------------------------------------------
mite.spe.hel <- decostand(mite.spe, method = "hellinger")


## -----------------------------------------------------------------------------
mite.spe.h.pcoa <- pcoa(dist(mite.spe.hel))


## ---- fig.width=6, fig.height=6-----------------------------------------------
biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel)


## -----------------------------------------------------------------------------
spe.nmds <- metaMDS(spe, distance = 'bray', k = 2)


## -----------------------------------------------------------------------------
spe.nmds$stress
stressplot(spe.nmds, main = "Shepard plot")


## ---- eval = F----------------------------------------------------------------
## spe.nmds <- metaMDS(spe, distance = 'bray', k = 2)
## spe.nmds$stress
## stressplot(spe.nmds, main = "Shepard plot")


## ---- eval = F----------------------------------------------------------------
## plot(spe.nmds, type = "none",
##      main = paste("NMDS/Bray - Stress =",
##                   round(spe.nmds$stress, 3)),
##      xlab = c("NMDS1"), ylab = "NMDS2")
## 
## points(scores(spe.nmds, display = "sites",
##               choiches = c(1,2),
##               pch = 21,
##               col = "black",
##               g = "steelblue",
##               cex = 1.2))
## text(scores(spe.nmds, display = "species", choices = c(1)),
##             scores(spe.nmds, display = "species", choices = c(2)),
##             labels = rownames(scores(spe.nmds, display = "species")),
##             col = "red", cex = 0.8)

