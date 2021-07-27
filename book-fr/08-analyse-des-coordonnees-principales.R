## # En utilisant la fonction cmdscale()
## ?cmdscale
## cmdscale(dist(spe.hel), k=(nrow(spe)-1), eig=TRUE)
## 
## # En utilisant la fonction pcoa()
## ?pcoa
## spe.h.pcoa<-pcoa(dist(spe.hel))
## 
## # Extraction des résultats
## spe.h.pcoa
## 
## # Représentation graphique
## biplot.pcoa(spe .h.pcoa, spe.hel, dir.axis2=-1)

## spe.bray.pcoa<-pcoa(spe.db) # il s'agit de la matrice de distances de Bray-Curtis qu'on a générée plus tôt
## spe.bray.pcoa
## biplot.pcoa(spe.bray.pcoa, spe.hel, dir.axis2=-1)
## # Le choix d'une mesure de distance est très important car ça influence les résultats!

## mite.spe.h.pcoa<-pcoa(dist(mite.spe.hel))
## mite.spe.h.pcoa
## windows()
## biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel, dir.axis2=-1)
