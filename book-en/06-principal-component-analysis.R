## #Run the PCA using the rda() function (NB: rda() is used for both PCA and RDA)
## spe.h.pca <- rda(spe.hel)
## 
## #Extract the results
## summary(spe.h.pca) #overall results

## summary(spe.h.pca, display=NULL) # only eigenvalues and their contribution to the variance
## eigen(cov(spe.hel)) # also compute the eigenvalues

## spe.scores <- scores(spe.h.pca, display="species", choices=c(1,2)) # species scores on the first two PCA axes
## site.scores <- scores(spe.h.pca, display="sites", choices=c(1,2)) # sites scores on the first two PCA axes
## #Note: if you don’t specify the number of principal components to extract (e.g. choices=c(1,2) or choices=c(1:2) then all of the scores will be extracted for all of the principal components.

## # Identify the significant axis using the Kaiser-Guttman criterion
## ev <- spe.h.pca$CA$eig
## ev[ev>mean(ev)]
## n <- length(ev)
## barplot(ev, main="Eigenvalues", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

## #Run the PCA
## env.pca <- rda(env.z) # or rda(env, scale=TRUE)
## 
## #Extract the results
## summary(env.pca)
## summary(env.pca, scaling=2)

## # Identify the significant axis using the Kaiser-Guttman criterion
## ev <- env.pca$CA$eig
## ev[ev>mean(ev)]
## n <- length(ev)
## barplot(ev, main="Eigenvalues", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

## plot(spe.h.pca)

## plot(spe.h.pca, type="n") #produces a blank biplot with nothing displayed but the axes
## points(spe.h.pca, dis="sp", col="blue") #points are added for the species (columns) (dis=)
## #use text() instead of points() if you want the labels
## points(spe.h.pca, dis="sites", col="red") #points are added for the sites (rows)

## #Biplot of the PCA on transformed species data (scaling 1)
## windows()
## plot(spe.h.pca)
## windows()
## biplot(spe.h.pca)
## windows()
## plot(spe.h.pca, scaling=1, type="none", # scaling 1 = distance biplot :
##                                         # distances among objects in the biplot approximate their Euclidean distances
##                                         # but angles among descriptor vectors DO NOT reflect their correlation
##      xlab = c("PC1 (%)", round((spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig))*100,2)), #this comes from the summary
##      ylab = c("PC2 (%)", round((spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig))*100,2)))
## points(scores(spe.h.pca, display="sites", choices=c(1,2), scaling=1),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(spe.h.pca, display="species", choices=c(1), scaling=1),
##      scores(spe.h.pca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(spe.h.pca, display="species", scaling=1)),
##      col="red", cex=0.8)

## #Biplot of the PCA on the environmental variables (scaling 2)
## windows()
## plot(env.pca)
## windows()
## plot(env.pca, scaling=2, type="none", # scaling 2 = correlation biplot :
##                                       # distances among abjects in the biplot DO NOT approximate their Euclidean distances
##                                       # but angles among descriptor vectors reflect their correlation
##      xlab = c("PC1 (%)", round((env.pca$CA$eig[1]/sum(env.pca$CA$eig))*100,2)),
##      ylab = c("PC2 (%)", round((env.pca$CA$eig[2]/sum(env.pca$CA$eig))*100,2)),
##      xlim = c(-1,1), ylim=c(-1,1))
## points(scores(env.pca, display="sites", choices=c(1,2), scaling=2),
##        pch=21, col="black", bg="darkgreen", cex=1.2)
## text(scores(env.pca, display="species", choices=c(1), scaling=2),
##      scores(env.pca, display="species", choices=c(2), scaling=2),
##      labels=rownames(scores(env.pca, display="species", scaling=2)),
##      col="red", cex=0.8)

## Sites_scores_Env_Axis1<- scores(env.pca, display="sites", choices=c(1), scaling=2)
## spe$ANG
## plot( Sites_scores_Env_Axis1, spe$TRU)
## summary(lm(spe$TRU~Sites_scores_Env_Axis1))
## abline(lm(spe$TRU~Sites_scores_Env_Axis1))

## mite.spe <- mite #mite data is from the vegan package

## #Hellinger transformation of mite data and PCA
## mite.spe.hel <- decostand(mite.spe, method="hellinger")
## mite.spe.h.pca <- rda(mite.spe.hel)
## 
## #What are the significant axes?
## ev <- mite.spe.h.pca$CA$eig
## ev[ev>mean(ev)]
## n <- length(ev)
## barplot(ev, main="Eigenvalues", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
## 
## #Output summary/results
## summary(mite.spe.h.pca, display=NULL)
## windows()
## 
## #Plot the biplot
## plot(mite.spe.h.pca, scaling=1, type="none",
##      xlab=c("PC1 (%)", round((mite.spe.h.pca$CA$eig[1]/sum(mite.spe.h.pca$CA$eig))*100,2)),
##      ylab=c("PC2 (%)", round((mite.spe.h.pca$CA$eig[2]/sum(mite.spe.h.pca$CA$eig))*100,2)))
## points(scores(mite.spe.h.pca, display="sites", choices=c(1,2), scaling=1),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(mite.spe.h.pca, display="species", choices=c(1), scaling=1),
##      scores(mite.spe.h.pca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(mite.spe.h.pca, display="species", scaling=1)),
##      col="red", cex=0.8)
