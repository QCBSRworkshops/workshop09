## #Exécuter la PCA avec la fonction rda()- cette fonction calcule à la fois des PCA et des RDA
## spe.h.pca<-rda(spe.hel)
## 
## #Extraire les résultats
## summary(spe.h.pca)

## summary(spe.h.pca, display=NULL) # seulement les valeurs propres
## eigen(cov(spe.hel)) # vous pouvez aussi trouver les valeurs propres par cette ligne de code

## spe.scores<-scores(spe.h.pca, display="species", choices=c(1,2)) # scores des espèces selon les premier et deuxième axes
## site.scores<-scores(spe.h.pca, display="sites", choices=c(1,2)) # scores des sites selon les premier et deuxième axes
## #Remarque: si vous ne spécifiez pas le nombre de composantes principales à l'aide de choices = c (1,2)
## #(ou choices = c (1: 2)), les scores selon toutes les composantes principales seront extraits.

## # Identification des axes significatifs de la PCA à l'aide du critère de Kaiser-Guttman
## ev<-spe.h.pca$CA$eig
## ev[ev>mean(ev)]
## n<-length(ev)
## bsm<-data.frame(j=seq(1:n), p=0)
## bsm$p[1]=1/n
## for (i in 2:n) {
##   bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
## bsm$p=100*bsm$p/n
## bsm
## barplot(ev, main="valeurs propres", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "moyenne des valeurs propres", lwd=1, col=2, bty="n")

## #Exécuter la PCA
## env.pca<-rda(env.z) # ou rda(env, scale=TRUE)
## 
## #Extraction des résultats
## summary(env.pca)
## summary(env.pca, scaling=2)

## ev<-env.pca$CA$eig
## ev[ev>mean(ev)]
## n<-length(ev)
## bsm<-data.frame(j=seq(1:n), p=0)
## bsm$p[1]=1/n
## for (i in 2:n) {
##   bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
## bsm$p=100*bsm$p/n
## bsm
## barplot(ev, main="valeurs propres", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "moyenne des valeurs propres", lwd=1, col=2, bty="n")

## plot(spe.h.pca)

## plot(spe.h.pca, type=”n”) # Produit une figure vierge
## points(spe.h.pca, dis=”sp”, col=”blue”) # ajoute les points correspondant aux espèces
## #utilizer text() plutôt que points() si vous préférez que les codes des espèces s'affichent (nom des colonnes)
## points(spe.h.pca, dis=”sites”, col=”red”) # ajoute les points correspondant aux sites

## #Scaling 1
## windows()
## plot(spe.h.pca)
## windows()
## biplot(spe.h.pca)
## windows()
## # scaling 1 = distance biplot :
## # distances entre les objets est une approximation de leur distance euclidienne
## # les angles entre les descripteurs ne réflètent PAS leur corrélation
## plot(spe.h.pca, scaling=1, type="none",
##      xlab<-c("PC1 (%)", round((spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig))*100,2)),
##      ylab<-c("PC2 (%)", round((spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig))*100,2)))
## points(scores(spe.h.pca, display="sites", choices=c(1,2), scaling=1),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(spe.h.pca, display="species", choices=c(1), scaling=1),
##      scores(spe.h.pca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(spe.h.pca, display="species", scaling=1)),
##      col="red", cex=0.8)

## #Scaling 2
## windows()
## plot(env.pca)
## windows()
## # scaling 2 = graphique de corrélations :
## # les distances entre les objets ne sont PAS des approximations de leur distance euclidienne
## # les angles entres les descripteurs reflètent leur corrélation
## plot(env.pca, scaling=2, type="none",
##      xlab<-c("PC1 (%)", round((env.pca$CA$eig[1]/sum(env.pca$CA$eig))*100,2)),
##      ylab<-c("PC2 (%)", round((env.pca$CA$eig[2]/sum(env.pca$CA$eig))*100,2)),
##      xlim<-c(-1,1), ylim=c(-1,1))
## points(scores(env.pca, display="sites", choices=c(1,2), scaling=2),
##        pch=21, col="black", bg="darkgreen", cex=1.2)
## text(scores(env.pca, display="species", choices=c(1), scaling=2),
##      scores(env.pca, display="species", choices=c(2), scaling=2),
##      labels<-rownames(scores(env.pca, display="species", scaling=2)),
##      col="red", cex=0.8)

## Sites_scores_Env_Axis1<- scores(env.pca, display="sites", choices=c(1), scaling=2)
## spe$ANG
## plot( Sites_scores_Env_Axis1, spe$TRU)
## summary(lm(spe$TRU~Sites_scores_Env_Axis1))
## abline(lm(spe$TRU~Sites_scores_Env_Axis1))

## mite.spe<-data(mite) # données disponibles dans vegan

## # Transformation de Hellinger
## mite.spe.hel<-decostand(mite.spe, method="hellinger")
## mite.spe.h.pca<-rda(mite.spe.hel)
## 
## # Quels sont les axes significatifs?
## ev<-mite.spe.h.pca$CA$eig
## ev[ev>mean(ev)]
## n<-length(ev)
## bsm<-data.frame(j=seq(1:n), p=0)
## bsm$p[1]=1/n
## for (i in 2:n) {
##   bsm$p[i]=bsm$p[i-1]+(1/(n=1-i))}
## bsm$p=100*bsm$p/n
## bsm
## barplot(ev, main="Valeurs propres", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Moyenne des valeurs propres", lwd=1, col=2, bty="n")
## 
## # Résultats
## summary(mite.spe.h.pca, display=NULL)
## windows()
## 
## # Représentation graphique de la PCA
## plot(mite.spe.h.pca, scaling=1, type="none",
##      xlab=c("PC1 (%)", round((mite.spe.h.pca$CA$eig[1]/sum(mite.spe.h.pca$CA$eig))*100,2)),
##      ylab=c("PC2 (%)", round((mite.spe.h.pca$CA$eig[2]/sum(mite.spe.h.pca$CA$eig))*100,2)))
## points(scores(mite.spe.h.pca, display="sites", choices=c(1,2), scaling=1),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(mite.spe.h.pca, display="species", choices=c(1), scaling=1),
##      scores(mite.spe.h.pca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(mite.spe.h.pca, display="species", scaling=1)),
##      col="red", cex=0.8)
