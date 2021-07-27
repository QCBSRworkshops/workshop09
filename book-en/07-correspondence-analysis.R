## #Run the CA using the cca() function (NB: cca() is used for both CA and CCA)
## spe.ca <- cca(spe)
## 
## # Identify the significant axes
## ev<-spe.ca$CA$eig
## ev[ev>mean(ev)]
## n=length(ev)
## barplot(ev, main="Eigenvalues", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")

## summary(spe.h.pca) #overall results
## summary(spe.h.pca, diplay=NULL)# only axis eigenvalues and contribution

## par(mfrow=c(1,2))
## #### scaling 1
## plot(spe.ca, scaling=1, type="none", main='CA - biplot scaling 1', xlab=c("CA1 (%)", round((spe.ca$CA$eig[1]/sum(spe.ca$CA$eig))*100,2)),
## ylab=c("CA2 (%)", round((spe.ca$CA$eig[2]/sum(spe.ca$CA$eig))*100,2)))
## 
## points(scores(spe.ca, display="sites", choices=c(1,2), scaling=1), pch=21, col="black", bg="steelblue", cex=1.2)
## 
## text(scores(spe.ca, display="species", choices=c(1), scaling=1),
##      scores(spe.ca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(spe.ca, display="species", scaling=1)),col="red", cex=0.8)
## 
## #### scaling 2
## plot(spe.ca, scaling=1, type="none", main='CA - biplot scaling 2', xlab=c("CA1 (%)", round((spe.ca$CA$eig[1]/sum(spe.ca$CA$eig))*100,2)),
##      ylab=c("CA2 (%)", round((spe.ca$CA$eig[2]/sum(spe.ca$CA$eig))*100,2)), ylim=c(-2,3))
## 
## points(scores(spe.ca, display="sites", choices=c(1,2), scaling=2), pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(spe.ca, display="species", choices=c(1), scaling=2),
##      scores(spe.ca, display="species", choices=c(2), scaling=2),
##      labels=rownames(scores(spe.ca, display="species", scaling=2)),col="red", cex=0.8)

## # CA on mite species
## mite.spe.ca<-cca(mite.spe)
## 
## #What are the significant axes ?
## ev<-mite.spe.ca$CA$eig
## ev[ev>mean(ev)]
## n=length(ev)
## barplot(ev, main="Eigenvalues", col="grey", las=2)
## abline(h=mean(ev), col="red")
## legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
## 
## #Output summary/results
## summary(mite.spe.ca, display=NULL)
## 
## #Plot the biplot
## windows()
## plot(mite.spe.ca, scaling=1, type="none",
##      xlab=c("PC1 (%)", round((mite.spe.ca$CA$eig[1]/sum(mite.spe.ca$CA$eig))*100,2)),
##      ylab=c("PC2 (%)", round((mite.spe.ca$CA$eig[2]/sum(mite.spe.ca$CA$eig))*100,2)))
## points(scores(mite.spe.ca, display="sites", choices=c(1,2), scaling=1),
##        pch=21, col="black", bg="steelblue", cex=1.2)
## text(scores(mite.spe.ca, display="species", choices=c(1), scaling=1),
##      scores(mite.spe.ca, display="species", choices=c(2), scaling=1),
##      labels=rownames(scores(mite.spe.ca, display="species", scaling=1)),
##      col="red", cex=0.8)
