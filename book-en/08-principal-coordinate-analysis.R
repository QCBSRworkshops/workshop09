## #Using cmdscale()
## ?cmdscale
## cmdscale(dist(spe.hel), k=(nrow(spe)-1), eig=TRUE)
## 
## #Using pcoa()
## ?pcoa
## spe.h.pcoa <- pcoa(dist(spe.hel))
## 
## # Extract the results
## spe.h.pcoa
## 
## #Construct the biplot
## biplot.pcoa(spe .h.pcoa, spe.hel, dir.axis2=-1)

## spe.bray.pcoa <- pcoa(spe.db) #where spe.db is the species dissimilarity matrix using Bray-Curtis.
## spe.bray.pcoa
## biplot.pcoa(spe.bray.pcoa, spe.hel, dir.axis2=-1)
## #Note that the distance measure chosen strongly influences the results.

## mite.spe.h.pcoa <- pcoa(dist(mite.spe.hel))
## mite.spe.h.pcoa
## windows()
## biplot.pcoa(mite.spe.h.pcoa, mite.spe.hel, dir.axis2=-1)
