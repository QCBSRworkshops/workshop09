##Section: 01-preparation-pour-l-atelier.R 

###Avis ###
#                                                                            #
#Ceci est un script généré automatiquement basé sur les morceaux de code du  #
#livre pour cet atelier.                                                     #
#                                                                            #
#Il est minimalement annoté pour permettre aux participants de fournir leurs #
#commentaires : une pratique que nous encourageons vivement.                 #
#                                                                            #
#Notez que les solutions aux défis sont également incluses dans ce script.   #
#Lorsque vous résolvez les défis par vous-méme, essayez de ne pas parcourir  #
#le code et de regarder les solutions.                                       #
#                                                                            #
#Bon codage !                                                               #


library (ade4)
data (doubs)

spe <- doubs$fish
env <- doubs$env

library (codep)
data (Doubs)

spe <- Doubs.fish
env <- Doubs.env

list.of.packages <- c("ape", "ade4", "codep", 
                      "gclus", "vegan", "GGally", 
                      "PlaneGeometry", "remotes", 
                      "matlib",
                      "MASS")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0) {
  install.packages(new.packages, dependencies = TRUE) 
  print(paste0("The following package was installed:", new.packages)) 
} else if(length(new.packages) == 0) {
  print("All packages were already installed previously")
}

# Chargement de toutes les bibliothèques nécessaires en une seule fois
invisible(lapply(list.of.packages, library, character.only = TRUE, quietly = TRUE))

# source(file.choose()) #  utilisez coldiss.R que vous avez téléchargé dans votre propre répertoire


##Section: 02-introduction-fr.R 




##Section: 03-exploration-des-donnees.R 

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

# Essayez-en quelques-uns !

names(spe) # noms des objets
dim(spe) # dimensions

str(spe) # structure des objets
summary(spe) # statistiques sommaires

head(spe) # 6 premières lignes

str(env)

summary(env) # statistiques sommaires


##Section: 04-association-distances.R 

dist(spe)

class(dist(spe))

as.matrix(dist(spe))

str(dist(spe))

dim(as.matrix(dist(spe)))

as.matrix(dist(spe))[1:6, 1:6]

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

# Créer une matrice
x <- matrix(rnorm(100*3), ncol = 3)

# Calculer la matrice de covariance et son inverse généralisé
cov_mat <- cov(x)
cov_inv <- MASS::ginv(cov_mat)

# Calculer la distance de Mahalanobis en utilisant l'inverse généralisé
mah_dist <- mahalanobis(x, 
                        colMeans(x), 
                        cov_inv)

# Imprimer la distance de Mahalanobis
mah_dist

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

# Utilisation :
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)

# Si D n'est pas une matrice de dissimilarité (max(D) > 1), alors D est divisée par max(D)
# nc nombre de couleurs (classes)
# byrank = TRUE classes de taille égale
# byrank = FALSE intervalles de longueur égale
# diag = TRUE imprime les étiquettes des objets également sur la diagonale

# Exemple :
# coldiss(spe.dj, nc=9, byrank=F, diag=T)

coldiss(spe.D.Jac)

# obtenir l'ordre des lignes et des colonnes
order_spe.D.Jac <- hclust(spe.D.Jac, method = "complete")$order

# réorganiser la matrice pour produire une figure ordonnée par similarités
order_spe.D.Jac_matrix <- as.matrix(spe.D.Jac)[order_spe.D.Jac, order_spe.D.Jac]

# converts to data frame
molten_spe.D.Jac <- reshape2::melt(
  as.matrix(order_spe.D.Jac_matrix)
  )

# créer un objet ggplot
ggplot(data = molten_spe.D.Jac, 
       aes(x = Var1, y = Var2, 
           fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal()

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

?dist

# matrice de distance euclidienne des variables environnementales standardisées
env.de <- dist(env.z, method = "euclidean")

windows() # Créer une fenêtre graphique séparée
coldiss(env.de, diag=TRUE)

(env.pearson<-cor(env)) # Calcul du r de Pearson entre les variables
round(env.pearson, 2) # Arrondit les coefficients à 2 points décimaux
(env.ken <- cor(env, method="kendall")) # Corrélation de rang du tau de Kendall
round(env.ken, 2)

var.g1 <- rnorm(30, 0, 1)
var.g2 <- runif(30, 0, 5)
var.g3 <- gl(3, 10)
var.g4 <- gl(2, 5, 30)

(dat2 <- data.frame(var.g1, var.g2, var.g3, var.g4))

str(dat2)
summary(dat2)


`?daisy #Cette fonction peut gérer les NA dans les données
(dat2.dg <- daisy(dat2, metric="gower"))
coldiss(dat2.dg)`

spe.challenge <- spe[1:3,1:3] #”[1:3,” refers to rows 1 to 3 while “,1:3]” refers to the first 3 species columns (in #this case the three variables of interest)

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

(db.s1s2<-Spec.s1s2/(Abund.s1+Abund.s2)) #Site 1 comparé au site 2
(db.s1s3<-Spec.s1s3/(Abund.s1+Abund.s3)) #Site 1 comparé au site 3
(db.s2s3<-Spec.s2s3/(Abund.s2+Abund.s3)) #Site 2 comparé au site 3

(spe.db.challenge<-vegdist(spe.challenge, method="bray"))

# Calculez le nombre de colonnes dans votre jeu de données
M<-ncol(spe.challenge)

# Calculez les différences d'abondance entre les paires de sites pour chaque espèce
Spe1.s1s2<-abs(spe.challenge[1,1]-spe.challenge[2,1])
Spe2.s1s2<-abs(spe.challenge[1,2]-spe.challenge[2,2])
Spe3.s1s2<-abs(spe.challenge[1,3]-spe.challenge[2,3])
Spe1.s1s3<-abs(spe.challenge[1,1]-spe.challenge[3,1])
Spe2.s1s3<-abs(spe.challenge[1,2]-spe.challenge[3,2])
Spe3.s1s3<-abs(spe.challenge[1,3]-spe.challenge[3,3])
Spe1.s2s3<-abs(spe.challenge[2,1]-spe.challenge[3,1])
Spe2.s2s3<-abs(spe.challenge[2,2]-spe.challenge[3,2])
Spe3.s2s3<-abs(spe.challenge[2,3]-spe.challenge[3,3])

# Calculer l'étendue de l'abondance de chaque espèce entre les sites
Range.spe1<-max(spe.challenge[,1]) - min (spe.challenge[,1])
Range.spe2<-max(spe.challenge[,2]) - min (spe.challenge[,2])
Range.spe3<-max(spe.challenge[,3]) - min (spe.challenge[,3])

# Calculer la dissimilarité de Gower
(dg.s1s2<-(1/M)*((Spe2.s1s2/Range.spe2)+(Spe3.s1s2/Range.spe3)))
(dg.s1s3<-(1/M)*((Spe2.s1s3/Range.spe2)+(Spe3.s1s3/Range.spe3)))
(dg.s2s3<-(1/M)*((Spe2.s2s3/Range.spe2)+(Spe3.s2s3/Range.spe3)))

# Comparez vos résultats
(spe.db.challenge<-vegdist(spe.challenge, method="gower"))

# Demonstration of a cluster dendrogram

spe.hel<-decostand(spe, method="hellinger")
spe.dhel <- vegdist(spe.hel,method="euclidean")
spe.dhel.ward <- hclust(spe.dhel, method="ward.D2")
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)

plot(spe.dhel.ward, hang=-1) # hang=-1 aligns all objets on the same line

# générer des échantillons de données
set.seed(123)
x <- matrix(rnorm(20), ncol = 2)

# effectuer un clustering agglomératif à lien unique
hc <- hclust(dist(x), method = "single")

# tracer le dendrogramme
plot(hc, 
     main = "Dendrogramme du regroupement agglomératif à lien unique",
     hang = -1)

# effectuer un clustering agglomératif à liens complets
hc <- hclust(dist(x), method = "complete")

# tracer le dendrogramme
plot(hc, 
     main = "Dendrogramme de l'agglomération de liens complète",
     hang = -1)

# effectuer un clustering par la méthode des groupes de paires non pondérés avec moyenne arithmétique
hc <- hclust(dist(x), method = "average")

# tracer le dendrogramme
plot(hc, 
     main = "Dendrogramme de la méthode du groupe de paires non pondéré avec la moyenne arithmétique \nAgglomerative Clustering",
     hang = -1)

# exécuter la méthode des groupes de paires pondérés avec la moyenne arithmétique du clustering
hc <- hclust(dist(x), method = "mcquitty")

# trace le dendrogramme
plot(hc, 
     main = "Dendrogramme de la méthode des groupes de paires pondérées avec le clustering agglomératif à moyenne arithmétique",
     hang = -1)

# effectuer le regroupement à variance minimale de Ward
hc <- hclust(dist(x), method = "ward.D")

# tracer le dendrogramme
plot(hc, 
     main = "Dendrogramme de \nWard's minimum variance Agglomerative Clustering",
     hang = -1)

# effectuer le regroupement à variance minimale de Ward
hc <- hclust(dist(x), method = "ward.D2")

# tracer le dendrogramme
plot(hc, 
     main = "Dendrogramme de \nWard's D2 minimum variance Agglomerative Clustering",
     hang = -1)

# Génère la matrice de distance à partir des données transformées de Hellinger
spe.dhel <- vegdist(spe.hel, method="euclidean") 

# Voir la différence entre les deux matrices
head(spe.hel)# Données d'espèces transformées par Hellinger
head(spe.dhel)# distances de Hellinger entre les sites

spe.dhel.single <- hclust(spe.dhel, method="single")

plot(spe.dhel.single)

spe.dhel.complete <- hclust(spe.dhel, method="complete")
plot(spe.dhel.complete)

# Effectuer le clustering de variance minimale de Ward
spe.dhel.ward <- hclust(spe.dhel, method="ward.D2")
plot(spe.dhel.ward)

# Reprendre le dendrogramme en utilisant les racines carrées des niveaux de fusion
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward)
plot(spe.dhel.ward, hang=-1) # hang=-1 aligne tous les objets sur la même ligne


##Section: 05-ordination-non-contrainte.R 




##Section: 06-analyse-des-composantes-principales.R 

library(knitr)
opts_chunk$set(fig.align = 'center')

# Chargement du paquet datasets
library(datasets)

# Charger le jeu de données varechem
data(varechem)

# Sélectionner les données
(data <- varechem[, 1:2])


data_std <- scale(data)

cov_matrix <- cov(data_std)

eigen_decomp <- eigen(cov_matrix)
Eigenvalues <- eigen_decomp$values
Eigenvectors <- eigen_decomp$vectors

F_PrComps <- data_std %*% Eigenvectors

head(F_PrComps)

Eigenvectors. <- as.data.frame(Eigenvectors)
row.names(Eigenvectors.) <- c("P", "N")
colnames(Eigenvectors.) <- c("PC1", "PC2")
Eigenvectors.[, 1] <- Eigenvectors.[, 1]*-1
Y_std <- as.data.frame(data_std)

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

PCA_princomp <- princomp(Y_std)

head(PCA_princomp$scores)

PCA_vegan_rda <- rda(Y_std)

scores(PCA_vegan_rda, 
       display = "sites", 
       scaling = 1,
       choices = seq_len(PCA_vegan_rda$CA$rank),
       const = sqrt(PCA_vegan_rda$tot.chi * (nrow(PCA_vegan_rda$CA$u) - 1)))[1:5, ]

spe.h.pca <- rda(spe.hel)

# summary(spe.h.pca)

sum(spe.h.pca$CA$eig)

paste(capture.output(summary(spe.h.pca))[c(26:29, 31:32, 34:40, 63:64, 66:72)])

paste(capture.output(summary(spe.h.pca))[c(32, 34:40)])

paste(capture.output(summary(spe.h.pca))[c(64, 66:72)])

scores(spe.h.pca,
       display = "species" or "sites")

ev <- spe.h.pca$CA$eig
# ev[ev > mean(ev)]

par(mar=c(1,4,2.5,.5), cex = 1.5)
n <- length(ev)
barplot(ev, main = "", col = "grey", las = 2)
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
biplot(spe.h.pca, scaling = 2)

par(mar = c(4,4,0.05,0.05), cex = 1.2)
biplot(spe.h.pca, scaling = 1)

data(mite)

mite.spe.hel <- decostand(mite, 
                          method = "hellinger")

mite.spe.h.pca <- rda(mite.spe.hel)

ev <- mite.spe.h.pca$CA$eig
ev[ev>mean(ev)]
n <- length(ev)
barplot(ev, main = "Valeurs propres",
        col = "grey", las = 2)
abline(h = mean(ev),
       col = "red3", lwd = 2)
legend("topright",
       "Valeur propre moyenne",
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
       "Valeur propre moyenne", 
       lwd = 2, col = "red3", 
       bty = "n")

par(mar = c(4,4,0.05,0.05), cex = 1.5)
biplot(mite.spe.h.pca, 
       col = c("red3", "grey15"))


##Section: 07-analyse-de-correspondence.R 

#Effectuer une CA à l'aide de la fonction cca() (NB: cca() est utilisée à la fois pour les CA et CCA)
spe.ca <- cca(spe[ -8, ])
 
# Identifier les axes significatifs
ev<-spe.ca$CA$eig

ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

# summary(spe.h.pca) 
# summary(spe.h.pca, diplay=NULL)

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

# CA 
mite.spe <- mite
mite.spe.ca<-cca(mite.spe)

# Quels sont les axes importants? 
ev<-mite.spe.ca$CA$eig
ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

# Résultats  
# summary(mite.spe.ca, display=NULL)

# Biplot
plot(mite.spe.ca, scaling=1, type="none",
     xlab=c("PC1 (%)", round((mite.spe.ca$CA$eig[1]/sum(mite.spe.ca$CA$eig))*100,2)),
     ylab=c("PC2 (%)", round((mite.spe.ca$CA$eig[2]/sum(mite.spe.ca$CA$eig))*100,2)))
points(scores(mite.spe.ca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.ca, display="species", choices=c(1), scaling=1),
     scores(mite.spe.ca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.ca, display="species", scaling=1)),
     col="red", cex=0.8)


##Section: 08-analyse-des-coordonnees-principales.R 

# En utilisant la fonction cmdscale()
# ?cmdscale
# cmdscale(dist(spe.hel), k=(nrow(spe)-1), eig=TRUE)

# En utilisant la fonction pcoa()
# ?pcoa
spe.h.pcoa<-pcoa(dist(spe.hel))

# Extraction des résultats
spe.h.pcoa 

# Représentation graphique
biplot.pcoa(spe.h.pcoa, spe.hel, dir.axis2 = -1)

spe.bray.pcoa<-pcoa(spe.db) # il s'agit de la matrice de distances de Bray-Curtis qu'on a générée plus tôt 
spe.bray.pcoa

biplot.pcoa(spe.bray.pcoa, spe.hel, dir.axis2 = -1)
# Le choix d'une mesure de distance est très important car ça influence les résultats! 

mite.spe.h.pcoa<-pcoa(dist(mite.spe.hel))
mite.spe.h.pcoa

biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel, dir.axis2=-1)


##Section: 09-echelle-multidimensionnelle-non-metrique.R 

# NMDS
spe.nmds<-metaMDS(spe[ -8, ], distance='bray', k=2)
 
#Extraction des résultats
spe.nmds

#Évaluation de la qualité de l'ajustement et construction du diagramme de Shepard
spe.nmds$stress
stressplot(spe.nmds, main='Shepard plot')

# Construction du biplot

plot(spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(spe.nmds$stress, 3)),
     xlab=c("NMDS1"),
     ylab=c("NMDS2"))
points(scores(spe.nmds, display="sites", choices=c(1,2)),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(spe.nmds, display="species", choices=c(1)),
     scores(spe.nmds, display="species", choices=c(2)),
     labels=rownames(scores(spe.nmds, display="species")),
     col="red", cex=0.8)

#NMDS
mite.spe.nmds<-metaMDS(mite.spe, distance='bray', k=2)

#Extraction des résultats
mite.spe.nmds

#Évaluation de la qualité de l'ajustement
mite.spe.nmds$stress
stressplot(mite.spe.nmds, main='Shepard plot')

#Construction du biplot

plot(mite.spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(mite.spe.nmds$stress, 3)),
     xlab=c("NMDS1"),
     ylab=c("NMDS2"))
points(scores(mite.spe.nmds, display="sites", choices=c(1,2)),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.nmds, display="species", choices=c(1)),
     scores(mite.spe.nmds, display="species", choices=c(2)),
     labels=rownames(scores(mite.spe.nmds, display="species")),
     col="red", cex=0.8)


##Section: 10-considerations-finales.R 




##Section: 11-references-fr.R 




