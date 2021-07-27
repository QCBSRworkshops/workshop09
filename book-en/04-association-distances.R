## spe.db<-vegdist(spe, method="bray") # "bray" with presence-absence data is Sorensen dissimilarity
## spe.dj<-vegdist(spe, method="jac") # Jaccard dissimilarity
## spe.dg<-vegdist(spe, method="gower") # Gower dissimilarity
## spe.db<-as.matrix(spe.db) #Put in matrix form (can visualize, write to .csv etc)

## windows()
## coldiss(spe.db, byrank=FALSE, diag=TRUE) # Heat map of Bray-Curtis dissimilarity
## windows()
## coldiss(spe.dj, byrank=FALSE, diag=TRUE) # Heat map of Jaccard dissimilarity
## windows()
## coldiss(spe.dg, byrank=FALSE, diag=TRUE) # Heat map of Gower dissimilarity

## ?dist # this function also compute dissimilarity matrix
## env.de<-dist(env.z, method = "euclidean") # euclidean distance matrix of the standardized environmental variables
## windows() #Creates a separate graphical window
## coldiss(env.de, diag=TRUE)

## (env.pearson<-cor(env)) # Pearson r linear correlation
## round(env.pearson, 2) #Rounds the coefficients to 2 decimal points
## (env.ken<-cor(env, method="kendall")) # Kendall tau rank correlation
## round(env.ken, 2)

## var.g1<-rnorm(30, 0, 1)
## var.g2<-runif(30, 0, 5)
## var.g3<-gl(3, 10)
## var.g4<-gl(2, 5, 30)
## (dat2<-data.frame(var.g1, var.g2, var.g3, var.g4))
## str(dat2)
## summary(dat2)

## ?daisy #This function can handle NAs in the data
## (dat2.dg<-daisy(dat2, metric="gower"))
## coldiss(dat2.dg)

## spe.challenge<-spe[1:3,1:3] #”[1:3,” refers to rows 1 to 3 while “,1:3]” refers to the first 3 species columns (in #this case the three variables of interest)

## (Abund.s1<-sum(spe.challenge[1,]))
## (Abund.s2<-sum(spe.challenge[2,]))
## (Abund.s3<-sum(spe.challenge[3,]))
## #() around code will cause output to print right away in console

## Spec.s1s2<-0
## Spec.s1s3<-0
## Spec.s2s3<-0
## for (i in 1:3) {
##   Spec.s1s2<-Spec.s1s2+abs(sum(spe.challenge[1,i]-spe.challenge[2,i]))
##   Spec.s1s3<-Spec.s1s3+abs(sum(spe.challenge[1,i]-spe.challenge[3,i]))
##   Spec.s2s3<-Spec.s2s3+abs(sum(spe.challenge[2,i]-spe.challenge[3,i])) }

## (db.s1s2<-Spec.s1s2/(Abund.s1+Abund.s2)) #Site 1 compared to site 2
## (db.s1s3<-Spec.s1s3/(Abund.s1+Abund.s3)) #Site 1 compared to site 3
## (db.s2s3<-Spec.s2s3/(Abund.s2+Abund.s3)) #Site 2 compared to site 3

## (spe.db.challenge<-vegdist(spe.challenge, method="bray"))

## # Calculate the number of columns in your dataset
## M<-ncol(spe.challenge)
## 
## # Calculate the species abundance differences between pairs of sites for each species
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
## # Calculate the range of each species abundance between sites
## Range.spe1<-max(spe.challenge[,1]) - min (spe.challenge[,1])
## Range.spe2<-max(spe.challenge[,2]) - min (spe.challenge[,2])
## Range.spe3<-max(spe.challenge[,3]) - min (spe.challenge[,3])
## 
## # Calculate the Gower dissimilarity
## (dg.s1s2<-(1/M)*((Spe2.s1s2/Range.spe2)+(Spe3.s1s2/Range.spe3)))
## (dg.s1s3<-(1/M)*((Spe2.s1s3/Range.spe2)+(Spe3.s1s3/Range.spe3)))
## (dg.s2s3<-(1/M)*((Spe2.s2s3/Range.spe2)+(Spe3.s2s3/Range.spe3)))
## 
## # Compare your results
## (spe.db.challenge<-vegdist(spe.challenge, method="gower"))

## spe.pa<-decostand(spe, method="pa")

## #Hellinger transformation
## spe.hel<-decostand(spe, method="hellinger") #can also use method=”hell”
## 
## #Chi-square transformation
## spe.chi<-decostand(spe, method="chi.square")

## # Hellinger transformation
## # First calculate the total species abundances by site
## (site.totals=apply(spe, 1, sum))
## 
## # Scale species abundances by dividing them by site totals
## (scale.spe<-spe/site.totals)
## 
## # Calculate the square root of scaled species abundances
## (sqrt.scale.spe<-sqrt(scale.spe))
## 
## # Compare the results
## sqrt.scale.spe
## spe.hel
## sqrt.scale.spe-spe.hel # or: sqrt.scale.spe/spe.hel
## 
## # Chi-square transformation
## # First calculate the total species abundances by site
## (site.totals<-apply(spe, 1, sum))
## 
## # Then calculate the square root of total species abundances
## (sqrt.spe.totals<-sqrt(apply(spe, 2, sum)))
## 
## # Scale species abundances by dividing them by the site totals and the species totals
## scale.spe2<-spe
## for (i in 1:nrow(spe)) {
##   for (j in 1:ncol(spe)) {
##    (scale.spe2[i,j]=scale.spe2[i,j]/(site.totals[i]*sqrt.spe.totals[j]))   }}
## 
## #Adjust the scale abundance species by multiplying by the square root of the species matrix total
## (adjust.scale.spe2<-scale.spe2*sqrt(sum(rowSums(spe))))
## 
## # Compare the results
## adjust.scale.spe2
## spe.chi
## adjust.scale.spe2-spe.chi # or: adjust.scale.spe2/spe.chi

## spe.dhel<-vegdist(spe.hel,method="euclidean") #generates the distance matrix from Hellinger transformed data
## 
## #See difference between the two matrices
## head(spe.hel)# Hellinger-transformed species data
## head(spe.dhel)# Hellinger distances among sites

## #Faire le groupement à liens simples
## #Perform single linkage clustering
## spe.dhel.single<-hclust(spe.dhel, method="single")
## plot(spe.dhel.single)
## 
## #Perform complete linkage clustering
## spe.dhel.complete<-hclust(spe.dhel, method="complete")
## plot(spe.dhel.complete)

## #Perform Ward minimum variance clustering
## spe.dhel.ward<-hclust(spe.dhel, method="ward.D2")
## plot(spe.dhel.ward)
## 
## #Re-plot the dendrogram by using the square roots of the fusion levels
## spe.dhel.ward$height<-sqrt(spe.dhel.ward$height)
## plot(spe.dhel.ward)
## plot(spe.dhel.ward, hang=-1) # hang=-1 aligns all objets on the same line
