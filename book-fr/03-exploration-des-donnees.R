## #Matrice d'abondances d'espèces: “DoubsSpe.csv”
## spe<- read.csv(file.choose(), row.names=1)
## spe<- spe[-8,] #Pas d'espèces dans le site 8, supprimer site 8.
## #Exécutez cette ligne une seule fois.
## 
## #L'environnement: “DoubsEnv.csv”
## env<- read.csv(file.choose(), row.names=1)
## env<- env[-8,] #Supprimer site 8 puisqu'on l'a retiré de la matrice d'abondances.
## #Exécutez cette ligne une seule fois.

## names(spe) # Les noms des colonnes
## dim(spe) # Le nombre de lignes et de colonnes.
## str(spe) # La structure interne de la matrice.
## head(spe) # Les premières lignes.
## summary(spe) # Les statistiques descriptives.

## (ab<-table(unlist(spe))) #Les parenthèses signifient que la sortie s'affiche immédiatement
## 
## barplot(ab, las=1, xlab=”Abundance class”, ylab=”Frequency”, col=grey(5:0/5))

## sum(spe==0)

## sum(spe==0)/(nrow(spe)*ncol(spe))

## spe.pres<- colSums(spe>0) # Somme des sites où chaque espèce est présente.
## hist(spe.pres, main=”Cooccurrence des espèces”, las=1, xlab=”Fréquence”, breaks=seq(0,30, by=5), col=”grey”)

## site.pres<- rowSums(spe>0) # Nombre d'espèces présentes dans chaque site
## hist(site.pres, main=”Richesse en espèces”, las=1, xlab=”Fréquence des sites”, ylab=”Nombre d'espèces”, breaks=seq(0,30, by=5), col=”grey”)

## names(env)
## dim(env)
## str(env)
## head(env)
## summary(env)
## pairs(env, main="Données environnementales" )

## env.z<-decostand(env, method="standardize")
## apply(env.z, 2, mean) # Les données sont maintenant centrées (moyennes~0)...
## apply(env.z, 2, sd)   # et réduites (écart-type=1)
