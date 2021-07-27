## # NMDS
## spe.nmds<-metaMDS(spe, distance='bray', k=2)
## 
## ### Extraction des résultats
## spe.nmds
## 
## ### Évaluation de la qualité de l'ajustement et construction du diagramme de Shepard
## spe.nmds$stress
## stressplot(spe.nmds, main='Shepard plot')
## 
## # Construction du biplot
## windows()
## plot(spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(spe.nmds$stress, 3)),
##      xlab=c("NMDS1"),
##      ylab=c("NMDS2"))
## points(scores(spe.nmds, display="sites", choices=c(1,2)),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(spe.nmds, display="species", choices=c(1)),
##      scores(spe.nmds, display="species", choices=c(2)),
##      labels=rownames(scores(spe.nmds, display="species")),
##      col="red", cex=0.8)

## ### NMDS
## mite.spe.nmds<-metaMDS(mite.spe, distance='bray', k=2)
## 
## ### Extraction des résultats
## mite.spe.nmds
## 
## ### Évaluation de la qualité de l'ajustement
## mite.spe.nmds$stress
## stressplot(mite.spe.nmds, main='Shepard plot')
## 
## ### Construction du biplot
## windows()
## plot(mite.spe.nmds, type="none", main=paste('NMDS/Bray - Stress=', round(mite.spe.nmds$stress, 3)),
##      xlab=c("NMDS1"),
##      ylab=c("NMDS2"))
## points(scores(mite.spe.nmds, display="sites", choices=c(1,2)),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(mite.spe.nmds, display="species", choices=c(1)),
##      scores(mite.spe.nmds, display="species", choices=c(2)),
##      labels=rownames(scores(mite.spe.nmds, display="species")),
##      col="red", cex=0.8)
