##Section: 01-preparation-pour-l-atelier.R 

install.packages("vegan")
install.packages("gclus")
install.packages("ape")

library(vegan)
library(gclus)
library(ape)

source(file.choose()) #coldiss.R


##Section: 02-introduction-fr.R 




##Section: 03-exploration-des-donnees.R 

#Matrice d'abondances d'espèces: “DoubsSpe.csv”
spe<- read.csv(file.choose(), row.names=1)
spe<- spe[-8,] #Pas d'espèces dans le site 8, supprimer site 8.
#Exécutez cette ligne une seule fois.

#L'environnement: “DoubsEnv.csv”
env<- read.csv(file.choose(), row.names=1)
env<- env[-8,] #Supprimer site 8 puisqu'on l'a retiré de la matrice d'abondances.
#Exécutez cette ligne une seule fois.

names(spe) # Les noms des colonnes
dim(spe) # Le nombre de lignes et de colonnes.
str(spe) # La structure interne de la matrice.
head(spe) # Les premières lignes.
summary(spe) # Les statistiques descriptives.

(ab<-table(unlist(spe))) #Les parenthèses signifient que la sortie s'affiche immédiatement

barplot(ab, las=1, xlab=”Abundance class”, ylab=”Frequency”, col=grey(5:0/5))

sum(spe==0)

sum(spe==0)/(nrow(spe)*ncol(spe))

spe.pres<- colSums(spe>0) # Somme des sites où chaque espèce est présente.
hist(spe.pres, main=”Cooccurrence des espèces”, las=1, xlab=”Fréquence”, breaks=seq(0,30, by=5), col=”grey”)

site.pres<- rowSums(spe>0) # Nombre d'espèces présentes dans chaque site
hist(site.pres, main=”Richesse en espèces”, las=1, xlab=”Fréquence des sites”, ylab=”Nombre d'espèces”, breaks=seq(0,30, by=5), col=”grey”)

names(env)
dim(env)
str(env)
head(env)
summary(env)
pairs(env, main="Données environnementales" )

env.z<-decostand(env, method="standardize")
apply(env.z, 2, mean) # Les données sont maintenant centrées (moyennes~0)...
apply(env.z, 2, sd)   # et réduites (écart-type=1)


##Section: 04-association-distances.R 

spe.db<-vegdist(spe, method="bray") # distance de Bray (avec des données de présence-absence, correspond à Sorensen)
spe.dj<-vegdist(spe, method="jac") # distance de Jaccard
spe.dg<-vegdist(spe, method="gower") # distance de Gower
spe.db<-as.matrix(spe.db) # réarranger en format matrice (pour visualisation, ou pour exporter en .csv)

windows()
coldiss(spe.db, byrank=FALSE, diag=TRUE) # Carte des points chauds Bray-Curtis
windows()
coldiss(spe.dj, byrank=FALSE, diag=TRUE) # Carte des points chauds Jaccard
windows()
coldiss(spe.dg, byrank=FALSE, diag=TRUE) # Carte des points chauds Gower

env.de<-dist(env.z, method = "euclidean") # matrice de distances euclidiennes des données env. standardisées
windows() # crée une nouvelle fenêtre graphique
coldiss(env.de, diag=TRUE)

(env.pearson<-cor(env)) # coefficient r de corrélation de Pearson
round(env.pearson, 2)  # arrondit les coefficients à deux décimales
(env.ken<-cor(env, method="kendall")) # coefficient tau de corrélation de rang de Kendall
round(env.ken, 2)

var.g1<-rnorm(30, 0, 1)
var.g2<-runif(30, 0, 5)
var.g3<-gl(3, 10)
var.g4<-gl(2, 5, 30)
(dat2<-data.frame(var.g1, var.g2, var.g3, var.g4))
str(dat2)
summary(dat2)

?daisy #Cette fonction peut gérer la présence de NA dans les données
(dat2.dg<-daisy(dat2, metric="gower"))
coldiss(dat2.dg)

spe.challenge<-spe[1:3,1:3] # les 3 premières lignes et 3 premières espèces (colonnes)

(Abund.s1<-sum(spe.challenge[1,]))
(Abund.s2<-sum(spe.challenge[2,]))
(Abund.s3<-sum(spe.challenge[3,]))

Spec.s1s2<-0
Spec.s1s3<-0
Spec.s2s3<-0
for (i in 1:3) {
  Spec.s1s2<-Spec.s1s2+abs(sum(spe.challenge[1,i]-spe.challenge[2,i]))
  Spec.s1s3<-Spec.s1s3+abs(sum(spe.challenge[1,i]-spe.challenge[3,i]))
  Spec.s2s3<-Spec.s2s3+abs(sum(spe.challenge[2,i]-spe.challenge[3,i])) }

(db.s1s2<-Spec.s1s2/(Abund.s1+Abund.s2)) #1 comparé à 2
(db.s1s3<-Spec.s1s3/(Abund.s1+Abund.s3)) #1 comparé à 3
(db.s2s3<-Spec.s2s3/(Abund.s2+Abund.s3)) #2 comparé à 3

(spe.db.challenge<-vegdist(spe.challenge, method="bray"))

# Calculer le nombre de colonnes
M<-ncol(spe.challenge)

# Calculer les différences d'abondance de chaque espèce entre paires de sites
Spe1.s1s2<-abs(spe.challenge[1,1]-spe.challenge[2,1])
Spe2.s1s2<-abs(spe.challenge[1,2]-spe.challenge[2,2])
Spe3.s1s2<-abs(spe.challenge[1,3]-spe.challenge[2,3])
Spe1.s1s3<-abs(spe.challenge[1,1]-spe.challenge[3,1])
Spe2.s1s3<-abs(spe.challenge[1,2]-spe.challenge[3,2])
Spe3.s1s3<-abs(spe.challenge[1,3]-spe.challenge[3,3])
Spe1.s2s3<-abs(spe.challenge[2,1]-spe.challenge[3,1])
Spe2.s2s3<-abs(spe.challenge[2,2]-spe.challenge[3,2])
Spe3.s2s3<-abs(spe.challenge[2,3]-spe.challenge[3,3])

# Calculer l'étendue d'abondance de chaque espèces parmi les sites
Range.spe1<-max(spe.challenge[,1]) - min (spe.challenge[,1])
Range.spe2<-max(spe.challenge[,2]) - min (spe.challenge[,2])
Range.spe3<-max(spe.challenge[,3]) - min (spe.challenge[,3])

# Calculer la distance de Gower
(dg.s1s2<-(1/M)*((Spe2.s1s2/Range.spe2)+(Spe3.s1s2/Range.spe3)))
(dg.s1s3<-(1/M)*((Spe2.s1s3/Range.spe2)+(Spe3.s1s3/Range.spe3)))
(dg.s2s3<-(1/M)*((Spe2.s2s3/Range.spe2)+(Spe3.s2s3/Range.spe3)))

# Vérifier vos résultats
(spe.db.challenge<-vegdist(spe.challenge, method="gower"))

spe.pa<-decostand(spe, method="pa")

#La transformation Hellinger
spe.hel<-decostand(spe, method="hellinger") # vous pouvez aussi simplement écrire "hel"

#Transformation de chi-carré
spe.chi<-decostand(spe, method="chi.square")

# Hellinger
# Calculer l'abondance des espèces par site
(site.totals=apply(spe, 1, sum))

# Réduire les abondances d'espèces en les divisant par les totaux par sites
(scale.spe<-spe/site.totals)

# Calculer la racine carrée des abondances d'espèces réduites
(sqrt.scale.spe<-sqrt(scale.spe))

# Comparer les résultats
sqrt.scale.spe
spe.hel
sqrt.scale.spe-spe.hel # ou: sqrt.scale.spe/spe.hel

# Chi-carré
# Premièrement calculer le total des abondances d'espèces par site
(site.totals<-apply(spe, 1, sum))

# Ensuite calculer la racine carrée du total des abondances d'espèces
(sqrt.spe.totals<-sqrt(apply(spe, 2, sum)))

# Réduire les abondances d'espèces en les divisant par les totaux par sites et les totaux par espèces
scale.spe2<-spe
for (i in 1:nrow(spe)) {
  for (j in 1:ncol(spe)) {
   (scale.spe2[i,j]=scale.spe2[i,j]/(site.totals[i]*sqrt.spe.totals[j]))   }}

#Ajuster les abondances en les multipliant par la racine carrée du total de la matrice des espèces
(adjust.scale.spe2<-scale.spe2*sqrt(sum(rowSums(spe))))

#Vérifier les résultats
adjust.scale.spe2
spe.chi
adjust.scale.spe2-spe.chi # or: adjust.scale.spe2/spe.chi

spe.dhel<-vegdist(spe.hel,method="euclidean") #crée une matrice de distances Hellinger à partir des données d’abondance transformées

#Pour voir la différence entre les deux types d’objets
head(spe.hel)# données d’abondances transformées Hellingerhead(spe.dhel)# matrice de distances de Hellinger entre les sites

#Faire le groupement à liens simples
spe.dhel.single<-hclust(spe.dhel, method="single")
plot(spe.dhel.single)

#Faire le groupement à liens complet
spe.dhel.complete<-hclust(spe.dhel, method="complete")
plot(spe.dhel.complete)

#Faire le groupement de Ward
spe.dhel.ward<-hclust(spe.dhel, method="ward.D2")
plot(spe.dhel.ward)

#Refaire le dendrogramme en utilisant la racine carrée des distances
spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
plot(spe.dhel.ward)
plot(spe.dhel.ward, hang=-1) # hang=-1 permet d’afficher les objets sur la même ligne


##Section: 05-ordination-non-contrainte.R 




##Section: 06-analyse-des-composantes-principales.R 

#Exécuter la PCA avec la fonction rda()- cette fonction calcule à la fois des PCA et des RDA
spe.h.pca<-rda(spe.hel)

#Extraire les résultats
summary(spe.h.pca)

summary(spe.h.pca, display=NULL) # seulement les valeurs propres
eigen(cov(spe.hel)) # vous pouvez aussi trouver les valeurs propres par cette ligne de code

spe.scores<-scores(spe.h.pca, display="species", choices=c(1,2)) # scores des espèces selon les premier et deuxième axes
site.scores<-scores(spe.h.pca, display="sites", choices=c(1,2)) # scores des sites selon les premier et deuxième axes
#Remarque: si vous ne spécifiez pas le nombre de composantes principales à l'aide de choices = c (1,2)
#(ou choices = c (1: 2)), les scores selon toutes les composantes principales seront extraits.

# Identification des axes significatifs de la PCA à l'aide du critère de Kaiser-Guttman
ev<-spe.h.pca$CA$eig
ev[ev>mean(ev)]
n<-length(ev)
bsm<-data.frame(j=seq(1:n), p=0)
bsm$p[1]=1/n
for (i in 2:n) {
  bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
bsm$p=100*bsm$p/n
bsm
barplot(ev, main="valeurs propres", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "moyenne des valeurs propres", lwd=1, col=2, bty="n")

#Exécuter la PCA
env.pca<-rda(env.z) # ou rda(env, scale=TRUE)

#Extraction des résultats
summary(env.pca)
summary(env.pca, scaling=2)

ev<-env.pca$CA$eig
ev[ev>mean(ev)]
n<-length(ev)
bsm<-data.frame(j=seq(1:n), p=0)
bsm$p[1]=1/n
for (i in 2:n) {
  bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
bsm$p=100*bsm$p/n
bsm
barplot(ev, main="valeurs propres", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "moyenne des valeurs propres", lwd=1, col=2, bty="n")

plot(spe.h.pca)

plot(spe.h.pca, type=”n”) # Produit une figure vierge
points(spe.h.pca, dis=”sp”, col=”blue”) # ajoute les points correspondant aux espèces
#utilizer text() plutôt que points() si vous préférez que les codes des espèces s'affichent (nom des colonnes)
points(spe.h.pca, dis=”sites”, col=”red”) # ajoute les points correspondant aux sites

#Scaling 1
windows()
plot(spe.h.pca)
windows()
biplot(spe.h.pca)
windows()
# scaling 1 = distance biplot :
# distances entre les objets est une approximation de leur distance euclidienne
# les angles entre les descripteurs ne réflètent PAS leur corrélation
plot(spe.h.pca, scaling=1, type="none",
     xlab<-c("PC1 (%)", round((spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig))*100,2)),
     ylab<-c("PC2 (%)", round((spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig))*100,2)))
points(scores(spe.h.pca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(spe.h.pca, display="species", choices=c(1), scaling=1),
     scores(spe.h.pca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(spe.h.pca, display="species", scaling=1)),
     col="red", cex=0.8)

#Scaling 2
windows()
plot(env.pca)
windows()
# scaling 2 = graphique de corrélations :
# les distances entre les objets ne sont PAS des approximations de leur distance euclidienne
# les angles entres les descripteurs reflètent leur corrélation
plot(env.pca, scaling=2, type="none",
     xlab<-c("PC1 (%)", round((env.pca$CA$eig[1]/sum(env.pca$CA$eig))*100,2)),
     ylab<-c("PC2 (%)", round((env.pca$CA$eig[2]/sum(env.pca$CA$eig))*100,2)),
     xlim<-c(-1,1), ylim=c(-1,1))
points(scores(env.pca, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="darkgreen", cex=1.2)
text(scores(env.pca, display="species", choices=c(1), scaling=2),
     scores(env.pca, display="species", choices=c(2), scaling=2),
     labels<-rownames(scores(env.pca, display="species", scaling=2)),
     col="red", cex=0.8)

Sites_scores_Env_Axis1<- scores(env.pca, display="sites", choices=c(1), scaling=2)
spe$ANG
plot( Sites_scores_Env_Axis1, spe$TRU)
summary(lm(spe$TRU~Sites_scores_Env_Axis1))
abline(lm(spe$TRU~Sites_scores_Env_Axis1))

mite.spe<-data(mite) # données disponibles dans vegan

# Transformation de Hellinger
mite.spe.hel<-decostand(mite.spe, method="hellinger")
mite.spe.h.pca<-rda(mite.spe.hel)

# Quels sont les axes significatifs?
ev<-mite.spe.h.pca$CA$eig
ev[ev>mean(ev)]
n<-length(ev)
bsm<-data.frame(j=seq(1:n), p=0)
bsm$p[1]=1/n
for (i in 2:n) {
  bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
bsm$p=100*bsm$p/n
bsm
barplot(ev, main="Valeurs propres", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Moyenne des valeurs propres", lwd=1, col=2, bty="n")

# Résultats
summary(mite.spe.h.pca, display=NULL)
windows()

# Représentation graphique de la PCA
plot(mite.spe.h.pca, scaling=1, type="none",
     xlab=c("PC1 (%)", round((mite.spe.h.pca$CA$eig[1]/sum(mite.spe.h.pca$CA$eig))*100,2)),
     ylab=c("PC2 (%)", round((mite.spe.h.pca$CA$eig[2]/sum(mite.spe.h.pca$CA$eig))*100,2)))
points(scores(mite.spe.h.pca, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
text(scores(mite.spe.h.pca, display="species", choices=c(1), scaling=1),
     scores(mite.spe.h.pca, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(mite.spe.h.pca, display="species", scaling=1)),
     col="red", cex=0.8)


##Section: 07-analyse-de-correspondence.R 

#Effectuer une CA à l'aide de la fonction cca() (NB: cca() est utilisée à la fois pour les CA et CCA)
spe.ca <- cca(spe)

# Identifier les axes significatifs
ev<-spe.ca$CA$eig
ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

summary(spe.h.pca)
summary(spe.h.pca, diplay=NULL)

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
mite.spe.ca<-cca(mite.spe)

# Quels sont les axes importants?
ev<-mite.spe.ca$CA$eig
ev[ev>mean(ev)]
n=length(ev)
barplot(ev, main="Eigenvalues", col="grey", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

# Résultats
summary(mite.spe.ca, display=NULL)

# Biplot
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


##Section: 08-analyse-des-coordonnees-principales.R 

# En utilisant la fonction cmdscale()
?cmdscale
cmdscale(dist(spe.hel), k=(nrow(spe)-1), eig=TRUE)

# En utilisant la fonction pcoa()
?pcoa
spe.h.pcoa<-pcoa(dist(spe.hel))

# Extraction des résultats
spe.h.pcoa

# Représentation graphique
biplot.pcoa(spe .h.pcoa, spe.hel, dir.axis2=-1)

spe.bray.pcoa<-pcoa(spe.db) # il s'agit de la matrice de distances de Bray-Curtis qu'on a générée plus tôt
spe.bray.pcoa
biplot.pcoa(spe.bray.pcoa, spe.hel, dir.axis2=-1)
# Le choix d'une mesure de distance est très important car ça influence les résultats!

mite.spe.h.pcoa<-pcoa(dist(mite.spe.hel))
mite.spe.h.pcoa
windows()
biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel, dir.axis2=-1)


##Section: 09-echelle-multidimensionnelle-non-metrique.R 

# NMDS
spe.nmds<-metaMDS(spe, distance='bray', k=2)

#Extraction des résultats
spe.nmds

#Évaluation de la qualité de l'ajustement et construction du diagramme de Shepard
spe.nmds$stress
stressplot(spe.nmds, main='Shepard plot')

# Construction du biplot
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

#NMDS
mite.spe.nmds<-metaMDS(mite.spe, distance='bray', k=2)

#Extraction des résultats
mite.spe.nmds

#Évaluation de la qualité de l'ajustement
mite.spe.nmds$stress
stressplot(mite.spe.nmds, main='Shepard plot')

#Construction du biplot
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


##Section: 10-considerations-finales.R 




##Section: 11-references-fr.R 




