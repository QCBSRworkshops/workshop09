## spe.db<-vegdist(spe, method="bray") # distance de Bray (avec des données de présence-absence, correspond à Sorensen)
## spe.dj<-vegdist(spe, method="jac") # distance de Jaccard
## spe.dg<-vegdist(spe, method="gower") # distance de Gower
## spe.db<-as.matrix(spe.db) # réarranger en format matrice (pour visualisation, ou pour exporter en .csv)

## windows()
## coldiss(spe.db, byrank=FALSE, diag=TRUE) # Carte des points chauds Bray-Curtis
## windows()
## coldiss(spe.dj, byrank=FALSE, diag=TRUE) # Carte des points chauds Jaccard
## windows()
## coldiss(spe.dg, byrank=FALSE, diag=TRUE) # Carte des points chauds Gower

## env.de<-dist(env.z, method = "euclidean") # matrice de distances euclidiennes des données env. standardisées
## windows() # crée une nouvelle fenêtre graphique
## coldiss(env.de, diag=TRUE)

## (env.pearson<-cor(env)) # coefficient r de corrélation de Pearson
## round(env.pearson, 2)  # arrondit les coefficients à deux décimales
## (env.ken<-cor(env, method="kendall")) # coefficient tau de corrélation de rang de Kendall
## round(env.ken, 2)

## var.g1<-rnorm(30, 0, 1)
## var.g2<-runif(30, 0, 5)
## var.g3<-gl(3, 10)
## var.g4<-gl(2, 5, 30)
## (dat2<-data.frame(var.g1, var.g2, var.g3, var.g4))
## str(dat2)
## summary(dat2)

## ?daisy #Cette fonction peut gérer la présence de NA dans les données
## (dat2.dg<-daisy(dat2, metric="gower"))
## coldiss(dat2.dg)

## spe.challenge<-spe[1:3,1:3] # les 3 premières lignes et 3 premières espèces (colonnes)

## (Abund.s1<-sum(spe.challenge[1,]))
## (Abund.s2<-sum(spe.challenge[2,]))
## (Abund.s3<-sum(spe.challenge[3,]))

## Spec.s1s2<-0
## Spec.s1s3<-0
## Spec.s2s3<-0
## for (i in 1:3) {
##   Spec.s1s2<-Spec.s1s2+abs(sum(spe.challenge[1,i]-spe.challenge[2,i]))
##   Spec.s1s3<-Spec.s1s3+abs(sum(spe.challenge[1,i]-spe.challenge[3,i]))
##   Spec.s2s3<-Spec.s2s3+abs(sum(spe.challenge[2,i]-spe.challenge[3,i])) }

## (db.s1s2<-Spec.s1s2/(Abund.s1+Abund.s2)) #1 comparé à 2
## (db.s1s3<-Spec.s1s3/(Abund.s1+Abund.s3)) #1 comparé à 3
## (db.s2s3<-Spec.s2s3/(Abund.s2+Abund.s3)) #2 comparé à 3

## (spe.db.challenge<-vegdist(spe.challenge, method="bray"))

## # Calculer le nombre de colonnes
## M<-ncol(spe.challenge)
## 
## # Calculer les différences d'abondance de chaque espèce entre paires de sites
## Spe1.s1s2<-abs(spe.challenge[1,1]-spe.challenge[2,1])
## Spe2.s1s2<-abs(spe.challenge[1,2]-spe.challenge[2,2])
## Spe3.s1s2<-abs(spe.challenge[1,3]-spe.challenge[2,3])
## Spe1.s1s3<-abs(spe.challenge[1,1]-spe.challenge[3,1])
## Spe2.s1s3<-abs(spe.challenge[1,2]-spe.challenge[3,2])
## Spe3.s1s3<-abs(spe.challenge[1,3]-spe.challenge[3,3])
## Spe1.s2s3<-abs(spe.challenge[2,1]-spe.challenge[3,1])
## Spe2.s2s3<-abs(spe.challenge[2,2]-spe.challenge[3,2])
## Spe3.s2s3<-abs(spe.challenge[2,3]-spe.challenge[3,3])
## 
## # Calculer l'étendue d'abondance de chaque espèces parmi les sites
## Range.spe1<-max(spe.challenge[,1]) - min (spe.challenge[,1])
## Range.spe2<-max(spe.challenge[,2]) - min (spe.challenge[,2])
## Range.spe3<-max(spe.challenge[,3]) - min (spe.challenge[,3])
## 
## # Calculer la distance de Gower
## (dg.s1s2<-(1/M)*((Spe2.s1s2/Range.spe2)+(Spe3.s1s2/Range.spe3)))
## (dg.s1s3<-(1/M)*((Spe2.s1s3/Range.spe2)+(Spe3.s1s3/Range.spe3)))
## (dg.s2s3<-(1/M)*((Spe2.s2s3/Range.spe2)+(Spe3.s2s3/Range.spe3)))
## 
## # Vérifier vos résultats
## (spe.db.challenge<-vegdist(spe.challenge, method="gower"))

## spe.pa<-decostand(spe, method="pa")

## #La transformation Hellinger
## spe.hel<-decostand(spe, method="hellinger") # vous pouvez aussi simplement écrire "hel"
## 
## #Transformation de chi-carré
## spe.chi<-decostand(spe, method="chi.square")

## # Hellinger
## # Calculer l'abondance des espèces par site
## (site.totals=apply(spe, 1, sum))
## 
## # Réduire les abondances d'espèces en les divisant par les totaux par sites
## (scale.spe<-spe/site.totals)
## 
## # Calculer la racine carrée des abondances d'espèces réduites
## (sqrt.scale.spe<-sqrt(scale.spe))
## 
## # Comparer les résultats
## sqrt.scale.spe
## spe.hel
## sqrt.scale.spe-spe.hel # ou: sqrt.scale.spe/spe.hel
## 
## # Chi-carré
## # Premièrement calculer le total des abondances d'espèces par site
## (site.totals<-apply(spe, 1, sum))
## 
## # Ensuite calculer la racine carrée du total des abondances d'espèces
## (sqrt.spe.totals<-sqrt(apply(spe, 2, sum)))
## 
## # Réduire les abondances d'espèces en les divisant par les totaux par sites et les totaux par espèces
## scale.spe2<-spe
## for (i in 1:nrow(spe)) {
##   for (j in 1:ncol(spe)) {
##    (scale.spe2[i,j]=scale.spe2[i,j]/(site.totals[i]*sqrt.spe.totals[j]))   }}
## 
## #Ajuster les abondances en les multipliant par la racine carrée du total de la matrice des espèces
## (adjust.scale.spe2<-scale.spe2*sqrt(sum(rowSums(spe))))
## 
## #Vérifier les résultats
## adjust.scale.spe2
## spe.chi
## adjust.scale.spe2-spe.chi # or: adjust.scale.spe2/spe.chi

## spe.dhel<-vegdist(spe.hel,method="euclidean") #crée une matrice de distances Hellinger à partir des données d’abondance transformées
## 
## #Pour voir la différence entre les deux types d’objets
## head(spe.hel)# données d’abondances transformées Hellingerhead(spe.dhel)# matrice de distances de Hellinger entre les sites

## #Faire le groupement à liens simples
## spe.dhel.single<-hclust(spe.dhel, method="single")
## plot(spe.dhel.single)
## 
## #Faire le groupement à liens complet
## spe.dhel.complete<-hclust(spe.dhel, method="complete")
## plot(spe.dhel.complete)

## #Faire le groupement de Ward
## spe.dhel.ward<-hclust(spe.dhel, method="ward.D2")
## plot(spe.dhel.ward)
## 
## #Refaire le dendrogramme en utilisant la racine carrée des distances
## spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
## plot(spe.dhel.ward)
## plot(spe.dhel.ward, hang=-1) # hang=-1 permet d’afficher les objets sur la même ligne
