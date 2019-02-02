library(scatterplot3d)
library(sqldf)    #for compare
library(vegan)    #for dist
library(MASS)     #for NMDS
library(graphics) #for graphics
library(rgl)      #for 3D plots
library(cluster)  #for hclust
library(pvclust)  #for p values in clusters
library(A2R)      #for dendrograms!!!
#library(foreign)  #for DBF

#required by A2R
library(trimcluster)
library(prabclus)
library(MASS)
library(cluster)
library(mclust)
library(flexmix)
library(modeltools)
library(stats4)
library(multcomp)
library(mvtnorm)

#############
# FUNCTIONS #
#############

#distance function
d <- function(comm, method) 
  betadiver(x=comm, method="sim")

#############
# Data code #
#############
setwd("C:/Data/Moss Biogeo")
OGU <- read.table("input/OGUused.txt", sep="\t", header=T)

###############################
# Filter species subset here! #
###############################
Dist <- read.table("input/AllSpecDist3.txt", sep="\t", header=T, na.strings="nix")
#Only species that are legitimate OR
Dist <- Dist[Dist$STATUS == "Legitimate",]
#Only Species of rank 4
#Dist <- Dist[Dist$RANK == 4,]

#kick out some OGUs that cause problems
Dist <- Dist[Dist$AREACODE != "AR-AN",]
Dist <- Dist[Dist$AREACODE != "AR-NO",]
Dist <- Dist[Dist$AREACODE != "BD",]

#Only OGUs of 100 or more species
SpecNum <- aggregate(Dist$AREACODE, by=list(Dist$AREACODE), FUN=length)
colnames(SpecNum) <- c("AREACODE", "SPECNUM")
Dist <- merge(Dist, SpecNum, by="AREACODE")
Dist <- Dist[Dist$SPECNUM >= 100,]
Dist <- as.data.frame(as.matrix(Dist[c(1,2)]))

ktabs <- table(Dist$AREACODE,Dist$NAMEID)
#ktabg <- table(Dist$AreaCode,Dist$Genus)
#ktabf <- table(Dist$AreaCode,Dist$Family)

#countspec <- aggregate(Dist,by=list(Dist$AREACODE),FUN=length)[,1:2]
#colnames(countspec) <- c("AreaCode", "SpeciesNumber")
#Regions <- merge(countspec, OGU, by.x="AreaCode", by.y="AreaCode", all.x=TRUE, all.y=FALSE)

#test species
#test <- scatterplot3d(SpecR4$points[,1], SpecR4$points[,2], SpecR4$points[,3], color="red")
#test$points3d(SpecR4$species[,1], SpecR4$species[,2], SpecR4$species[,3])
#test$points3d(SpecR4$points[,1], SpecR4$points[,2], SpecR4$points[,3], col="red", pch=16, cex=1.5)

##################
# NMS Ordination #
##################

SpecR4 <- metaMDS(ktabs, distance="mixed", distfun=d, zerodist="add", k=3, trymax=100, autotransform=T, wascores=T, expand=T, noshare=0.1, pc=T)
#GenR4 <- metaMDS(ktabg, distance="mixed", distfun=d, zerodist="add", k=3, trymax=100, autotransform=T, wascores=T, expand=T, noshare=0.1, pc=T)
#FamR4 <- metaMDS(ktabf, distance="mixed", distfun=d, zerodist="add", k=3, trymax=100, autotransform=T, wascores=T, expand=T, noshare=0.1, pc=T)

scatterplot3d(SpecR4$points[,1],SpecR4$points[,2],SpecR4$points[,3], xlab="NMS1 (red)", ylab="NMS2 (green)", zlab="NMS3 (blue)", box=FALSE, pch=16, type="h", cex.symbols=1, lty.hplot=3)

###############
# Adjust Axes #
###############

#configure ordination plot orientation and scale axes here!!!
    x <- SpecR4$points[,3]
    y <- SpecR4$points[,2]
    z <- SpecR4$points[,1]
    x <- x-min(x)
    x <- x*255/max(x)
    y <- y-min(y)
    y <- y*255/max(y)
    z <- z-min(z)
    z <- z*255/max(z)
    x <- (x*-1)+255
    y <- (y*-1)+255
    z <- (z*-1)+255
    NMDS <- as.data.frame(cbind(x,y,z))
    rownames(NMDS) <- names(SpecR4$points[,3])
    
    rm(x,y,z)

plot3d(NMDS$y,NMDS$x,NMDS$z, col=rgb(NMDS$x,NMDS$y,NMDS$z, maxColorValue=255), pch=16, cex.symbols=16, size=10)

NMDS[,1] <- NMDS[,1]*-1+255
NMDS[,2] <- NMDS[,2]*-1+255


####################
# Cluster Analysis #
####################

similaritytable <- betadiver(ktabs, method="sim")
fullclusters <- hclust(similaritytable, method="average")
h3c <- cutree(fullclusters, k=2)
h5c <- cutree(fullclusters, k=4)
h7c <- cutree(fullclusters, k=7)
h9c <- cutree(fullclusters, k=9)
h11c <- cutree(fullclusters, k=11)
clusters <- cbind(h3c,h5c,h7c,h9c,h11c)

dend <- cut(as.dendrogram(fullclusters), h=0.67)$upper

tree <- hclust(similaritytable, method="average", members=h11c)
plot(tree)

plot(dend, horiz=TRUE)
ggdendrogram(dend, rotate=TRUE)

#RGB color for cluster
clustercolor <- cbind(NMDS,clusters,rownames(NMDS))
meancolorvalues <- aggregate(clustercolor[,1:3], FUN="mean", by=list(h11c))
colnames(meancolorvalues) <- c("h11c", "R", "G", "B")
MDS <- merge(clustercolor, meancolorvalues, by="h11c")
colnames(MDS) <- c("h11c","x","y","z", "h3c", "h5c", "h7c", "h9c", "AreaCode", "R","G","B")

##############
# Dendrogram #
##############
memb <- cutree(fullclusters, k=11)
cent <- NULL
for(h in 1:11){
  cent <- rbind(cent, colMeans(ktabs[memb == h, , drop = FALSE]))
  }
cutclusters <- hclust(betadiver(cent, method="sim"), method="average", members=table(memb))

plot(cutclusters)
#A2Rplot(cutclusters)

#Species richness of clusters
memb2 <- as.data.frame(memb)
memb2 <- cbind(rownames(memb2), memb2)
colnames(memb2) <- c("AREACODE", "CLUSTER")

specclust <- merge(Dist, memb2, by="AREACODE")
specclust <- aggregate(specclust, by=list(specclust$CLUSTER, specclust$NAMEID), FUN=length)
specclust2 <- specclust[,c(1,2)]

clustrich <- aggregate(specclust2, by=list(specclust[,1]), FUN=length)



#p values for clusters
#pvalu <- pvclust(similaritytable, method.hclust="average")


####

#******
# UPGMA mit flashClust
#
distmat.upgma.flash <-hclust(distmat.bray,method="average") #UPGMA
distmat.upgma.flash.coph <- cophenetic(distmat.upgma.flash)
cor(distmat.bray,distmat.upgma.flash) # cophenetic r = 0.7374355



#******
# Plotte vereinfachtes Dendrogramm
#

library(maptree)
library(cluster)

#clipped.tree <- clip.clust(hclust(betadiver(ktabs, method="sim"), method="average"), ktabs, k=11)
draw.clust(cutclusters, ktabs, pch="")


Tree$order.lab<-as.matrix(spp[,2]) # this reconstructs the object as needed for 'clip.clust' commmand
Tree$data<-spp[,2:3] # 2:3 this reconstructs the object as needed for 'clip.clust' commmand
group<-group.clust(Tree,k=n) # creates a reasonable order of cluster groups based on the pruned tree (left to right)
draw.clust(clip.clust(Tree,k=n),size=6,pch=19,col=my.pal) # prunes tree at specified height and plots it


#write output file
write.table(MDS, "output/MDS.csv", sep="\t")


################
# Backup point #
################

setwd("C:/Data/Moss Biogeo")
#MDS <- read.table("output/MDS.csv", sep="\t")
#OGU <- read.table("input/OGUused.txt", sep="\t", header=T)

#Plot the clusters as map
#backupmapdata <- mapdata
mapdata <- merge(OGU, MDS, by="AreaCode")

#plot first cluster group
plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=mapdata$h3c)
legend(x="bottomleft", legend=unique(mapdata$h3c), fill=unique(mapdata$h3c))
#plot second cluster group
plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=mapdata$h5c)
legend(x="bottomleft", legend=unique(mapdata$h5c), fill=unique(mapdata$h5c))
#plot third cluster group
plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=mapdata$h7c)
legend(x="bottomleft", legend=unique(mapdata$h7c), fill=unique(mapdata$h7c))
#plot fourth cluster group
plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=mapdata$h11c)
legend(x="bottomleft", legend=unique(mapdata$h11c), fill=unique(mapdata$h9c))

plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=rgb(mapdata$x,mapdata$y,mapdata$z, maxColorValue=255))
plot(mapdata$Longitude, mapdata$Latitude, pch=20, cex=sqrt(mapdata$Area/100000000000), col=rgb(mapdata$R,mapdata$G,mapdata$B, maxColorValue=255))
#identify(mapdata$Longitude, mapdata$Latitude)

plot(MDS[,2],MDS[,3],pch=20, col=MDS$h11c)
#identify(MDS[,2],MDS[,3])

MDS$x <- MDS$x*-1+255
MDS$y <- MDS$y*-1+255
MDS$z <- MDS$z*-1+255

#3D plots for NMDS and Clusters
plot3d(MDS$y,MDS$x,MDS$z, col=rgb(MDS$x,MDS$y,MDS$z, maxColorValue=255), pch=16, cex.symbols=16, size=10)
plot3d(MDS[,3],MDS[,4],MDS[,2], col=rgb(MDS$x,MDS$y,MDS$z, maxColorValue=255), pch=16, cex.symbols=16, size=10)

plot3d(MDS$x,MDS$y,MDS$z, col=rgb(MDS$R,MDS$G,MDS$B, maxColorValue=255), pch=16, cex.symbols=16, size=10)
plot3d(MDS$x,MDS$y,MDS$z, col=MDS$h11c, pch=16, cex.symbols=16, size=10)


##########################
# Plot NMDS and Clusters #
##########################

tiff(file="output/NMDS plots.tiff", compression="none", units="cm", width=24, height=19.7, res=600)
  par(mfcol=c(3,3), oma=c(2,2,2,2), xpd=TRUE)
  s3d <- scatterplot3d(MDS[,2],MDS[,3]*-1+255,MDS[,4], y.ticklabs=c(300,250,200,150,100,50,0), color=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255),
    xlab="NMDS axis 1", ylab="", zlab="NMDS axis 3", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
    cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=16, type="h", cex.symbols=1, lty.hplot=3, mar=c(1,1,1,2))
    s3d$points3d(MDS[,2],MDS[,3]*-1+255,MDS[,4], pch=21, cex=1, bg=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255))
    text(9.5, 1, "NMDS axis 2", cex=0.9, srt=40)
    text(s3d$xyz.convert(-20,0,365), labels="(a)", cex=2)

  s3d <- scatterplot3d(MDS[,3],MDS[,4],MDS[,2], color=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255),
    xlab="NMDS axis 2", ylab="", zlab="NMDS axis 1", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
    cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=16, type="h", cex.symbols=1, lty.hplot=3, mar=c(1,1,1,2))
    s3d$points3d(MDS[,3],MDS[,4],MDS[,2], pch=21, cex=1, bg=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255))
    s3d$points3d(subset(MDS,h11c==3)[,3],subset(MDS,h11c==3)[,4],subset(MDS,h11c==3)[,2], pch=21, cex=1, bg=rgb(subset(MDS,h11c==3)[,2],subset(MDS,h11c==3)[,3],subset(MDS,h11c==3)[,4], maxColorValue=255))
    text(9.5, 1, "NMDS axis 3", cex=0.9, srt=40)
    text(s3d$xyz.convert(-20,0,365), labels="(b)", cex=2)

  s3d <- scatterplot3d(MDS[,4],MDS[,2],MDS[,3]*-1+255, z.ticklabs=c(300,250,200,150,100,50,0), color=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255),
    xlab="NMDS axis 3", ylab="", zlab="NMDS axis 2", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
    cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=16, type="h", cex.symbols=1, lty.hplot=3, mar=c(1,1,1,2))
    s3d$points3d(MDS[,4],MDS[,2],MDS[,3]*-1+255, pch=21, cex=1, bg=rgb(MDS[,2],MDS[,3],MDS[,4], maxColorValue=255))
    text(9.5, 1, "NMDS axis 1", cex=0.9, srt=40)
    text(s3d$xyz.convert(-20,0,365), labels="(c)", cex=2)

#Clusterplots
  #Plot1
     s3d <- scatterplot3d(MDS[,2],MDS[,3]*-1+255,MDS[,4], y.ticklabs=c(300,250,200,150,100,50,0), color="grey", bg="grey",
         xlab="NMDS axis 1", ylab="", zlab="NMDS axis 3", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
         cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3, mar=c(1,3,1,0))
         text(9.5, 1, "NMDS axis 2", cex=0.9, srt=40)
         text(s3d$xyz.convert(-20,0,365), labels="(d)", cex=2)
MDS[,3] <- MDS[,3]*-1+255
    #cluster1
     points(s3d$xyz.convert(subset(MDS, h9c == 1)[,c(2,3,4)]), pch=21, col="black", bg="blue", cex=1)
     polygon(s3d$xyz.convert(subset(MDS, h9c==1)[chull(s3d$xyz.convert
        (subset(MDS, h9c==1)[,2], subset(MDS, h9c==1)[,3], subset(MDS, h9c==1)[,4])),c(2,3,4)]),
        col="blue", density=40, border="black")
    #cluster4&5
     points(s3d$xyz.convert(subset(MDS, h9c==4|h9c==5)[,c(2,3,4)]), pch=21, col="black", bg="darkturquoise", cex=1)
     polygon(s3d$xyz.convert(subset(MDS, h9c==4|h9c==5)[chull(s3d$xyz.convert
         (subset(MDS, h9c==4|h9c==5)[,2], subset(MDS, h9c==4|h9c==5)[,3], subset(MDS, h9c==4|h9c==5)[,4])),c(2,3,4)]),
         col="darkturquoise", density=40, border="black")
MDS[,3] <- MDS[,3]*-1+255

  #Plot2
    s3d <- scatterplot3d(MDS[,3],MDS[,4],MDS[,2], color="grey", bg="grey",
        xlab="NMDS axis 2", ylab="", zlab="NMDS axis 1", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
        cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3, mar=c(1,3,1,0))
        text(9.5, 1, "NMDS axis 3", cex=0.9, srt=40)
        text(s3d$xyz.convert(-20,0,365), labels="(e)", cex=2)
    #cluster3
      points(s3d$xyz.convert(subset(MDS, h11c==3)[,c(3,4,2)]), pch=21, col="black", bg="orange", cex=1)
      polygon(s3d$xyz.convert(subset(MDS, h11c==3)[chull(s3d$xyz.convert
      (subset(MDS, h11c==3)[,3], subset(MDS, h11c==3)[,4], subset(MDS, h11c==3)[,2])),c(3,4,2)]),
      col="orange", density=40, border="black")
    #cluster8
      points(s3d$xyz.convert(subset(MDS, h11c==8)[,c(3,4,2)]), pch=21, col="black", bg="maroon4", cex=1)
      polygon(s3d$xyz.convert(subset(MDS, h11c==8)[chull(s3d$xyz.convert
      (subset(MDS, h11c==8)[,3], subset(MDS, h11c==8)[,4], subset(MDS, h11c==8)[,2])),c(3,4,2)]),
      col="maroon4", density=40, border="black")
    #cluster(7)2
      points(s3d$xyz.convert(subset(MDS, h7c==2)[,c(3,4,2)]), pch=21, col="black", bg="forestgreen", cex=1)
      polygon(s3d$xyz.convert(subset(MDS, h7c==2)[chull(s3d$xyz.convert
      (subset(MDS, h7c==2)[,3], subset(MDS, h7c==2)[,4], subset(MDS, h7c==2)[,2])),c(3,4,2)]),
      col="forestgreen", density=40, border="black")

  #Plot3

    s3d <- scatterplot3d(MDS[,4],MDS[,2],MDS[,3]*-1+255, color="grey", bg="grey", z.ticklabs=c(300,250,200,150,100,50,0),
        xlab="NMDS axis 3", ylab="", zlab="NMDS axis 2", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
        cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3, mar=c(1,3,1,0))
        text(9.5, 1, "NMDS axis 1", cex=0.9, srt=40)
        text(s3d$xyz.convert(-20,0,365), labels="(f)", cex=2)
MDS[,3] <- MDS[,3]*-1+255
    #cluster4
        points(s3d$xyz.convert(subset(MDS, h11c==4)[,c(4,2,3)]), pch=21, col="black", bg="lightgreen", cex=1)
        polygon(s3d$xyz.convert(subset(MDS, h11c==4)[chull(s3d$xyz.convert
        (subset(MDS, h11c==4)[,4], subset(MDS, h11c==4)[,2], subset(MDS, h11c==4)[,3])),c(4,2,3)]),
        col="lightgreen", density=40, border="black")
    #cluster5
        points(s3d$xyz.convert(subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,c(4,2,3)]), pch=21, col="black", bg="turquoise", cex=1)
        polygon(s3d$xyz.convert(subset(MDS, h11c==5 & AreaCode!="AU-QLD")[chull(s3d$xyz.convert
        (subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,4], subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,2], subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,3])),c(4,2,3)]),
        col="turquoise", density=40, border="black")
    #cluster6
        points(s3d$xyz.convert(subset(MDS, h11c==6)[,c(4,2,3)]), pch=21, col="black", bg="mediumorchid2", cex=1)
        polygon(s3d$xyz.convert(subset(MDS, h11c==6)[chull(s3d$xyz.convert
        (subset(MDS, h11c==6)[,4], subset(MDS, h11c==6)[,2], subset(MDS, h11c==6)[,3])),c(4,2,3)]),
        col="mediumorchid2", density=40, border="black")
    #Queensland
        points(s3d$xyz.convert(subset(MDS, AreaCode=="AU-QLD")[,c(4,2,3)]), pch=21, col="black", bg="turquoise", cex=1)
        text(s3d$xyz.convert(subset(MDS, AreaCode=="AU-QLD")[,c(4,2,3)]), labels="Queensland", col="black", cex=0.7, pos=1)
    #Hawaii
        points(s3d$xyz.convert(subset(MDS, AreaCode=="I-HW")[,c(4,2,3)]), pch=21, col="black", bg="yellow", cex=1)
        text(s3d$xyz.convert(subset(MDS, AreaCode=="I-HW")[,c(4,2,3)]), labels="Hawaii", col="black", cex=0.7, pos=3)
    #New Caldedonia
        points(s3d$xyz.convert(subset(MDS, h11c==11)[,c(4,2,3)]), pch=21, col="black", bg="darkgreen", cex=1)
        text(s3d$xyz.convert(subset(MDS, h11c==11)[,c(4,2,3)]), labels="New Caledonia", col="black", cex=0.7, pos=3)

MDS[,3] <- MDS[,3]*-1+255

dev.off()



######################
# Plots for talk IAB #
######################

tiff(file="output/talkplots/NMDS1.tiff", compression="none", units="cm", width=10, height=10, res=600)
s3d <- scatterplot3d(MDS[,2],MDS[,3],MDS[,4], color=rgb(MDS[,2],MDS[,3]*-1+255,MDS[,4], maxColorValue=255),
xlab="", ylab="", zlab="", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=16, type="h", cex.symbols=0.75, lty.hplot=3)
s3d$points3d(MDS[,2],MDS[,3],MDS[,4], pch=21, cex=0.75, bg=rgb(MDS[,2],MDS[,3]*-1+255,MDS[,4], maxColorValue=255))
dev.off()

tiff(file="output/talkplots/NMDS2.tiff", compression="none", units="cm", width=10, height=10, res=600)
s3d <- scatterplot3d(MDS[,2],MDS[,3],MDS[,4], color="grey", bg="grey",
xlab="", ylab="", zlab="", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3)
#MDS[,3] <- MDS[,3]*-1+255
#cluster1
points(s3d$xyz.convert(subset(MDS, h9c == 1)[,c(2,3,4)]), pch=21, col="black", bg="blue", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h9c==1)[chull(s3d$xyz.convert
(subset(MDS, h9c==1)[,2], subset(MDS, h9c==1)[,3], subset(MDS, h9c==1)[,4])),c(2,3,4)]),
col="blue", density=40, border="black")
#cluster4&5
points(s3d$xyz.convert(subset(MDS, h9c==4|h9c==5)[,c(2,3,4)]), pch=21, col="black", bg="darkturquoise", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h9c==4|h9c==5)[chull(s3d$xyz.convert
(subset(MDS, h9c==4|h9c==5)[,2], subset(MDS, h9c==4|h9c==5)[,3], subset(MDS, h9c==4|h9c==5)[,4])),c(2,3,4)]),
col="darkturquoise", density=40, border="black")
dev.off()

tiff(file="output/talkplots/NMDS3.tiff", compression="none", units="cm", width=10, height=10, res=600)
s3d <- scatterplot3d(MDS[,2],MDS[,3],MDS[,4], color="grey", bg="grey",
xlab="", ylab="", zlab="", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3)
#cluster3
points(s3d$xyz.convert(subset(MDS, h11c==3)[,c(2,3,4)]), pch=21, col="black", bg="orange", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h11c==3)[chull(s3d$xyz.convert
(subset(MDS, h11c==3)[,2], subset(MDS, h11c==3)[,3], subset(MDS, h11c==3)[,4])),c(2,3,4)]),
col="orange", density=40, border="black")
#cluster8
points(s3d$xyz.convert(subset(MDS, h11c==8)[,c(2,3,4)]), pch=21, col="black", bg="maroon4", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h11c==8)[chull(s3d$xyz.convert
(subset(MDS, h11c==8)[,2], subset(MDS, h11c==8)[,3], subset(MDS, h11c==8)[,4])),c(2,3,4)]),
col="maroon4", density=40, border="black")
#cluster(7)2
points(s3d$xyz.convert(subset(MDS, h7c==2)[,c(2,3,4)]), pch=21, col="black", bg="forestgreen", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h7c==2)[chull(s3d$xyz.convert
(subset(MDS, h7c==2)[,2], subset(MDS, h7c==2)[,3], subset(MDS, h7c==2)[,4])),c(2,3,4)]),
col="forestgreen", density=40, border="black")
dev.off()

#Plot3
tiff(file="output/talkplots/NMDS4.tiff", compression="none", units="cm", width=10, height=10, res=600)
s3d <- scatterplot3d(MDS[,2],MDS[,3],MDS[,4], color="grey", bg="grey",
xlab="", ylab="", zlab="", xlim=c(0,255), ylim=c(0,255), zlim=c(0,255),
cex.lab=0.6, cex.axis=0.6, box=FALSE, pch=21, type="p", cex.symbols=1, lty.hplot=3)
#cluster4
points(s3d$xyz.convert(subset(MDS, h11c==4)[,c(2,3,4)]), pch=21, col="black", bg="lightgreen", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h11c==4)[chull(s3d$xyz.convert
(subset(MDS, h11c==4)[,2], subset(MDS, h11c==4)[,3], subset(MDS, h11c==4)[,4])),c(2,3,4)]),
col="lightgreen", density=40, border="black")
#cluster5
points(s3d$xyz.convert(subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,c(2,3,4)]), pch=21, col="black", bg="turquoise", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h11c==5 & AreaCode!="AU-QLD")[chull(s3d$xyz.convert
(subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,2], subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,3], subset(MDS, h11c==5 & AreaCode!="AU-QLD")[,4])),c(2,3,4)]),
col="turquoise", density=40, border="black")
#cluster6
points(s3d$xyz.convert(subset(MDS, h11c==6)[,c(2,3,4)]), pch=21, col="black", bg="mediumorchid2", cex=1)
polygon(s3d$xyz.convert(subset(MDS, h11c==6)[chull(s3d$xyz.convert
(subset(MDS, h11c==6)[,2], subset(MDS, h11c==6)[,3], subset(MDS, h11c==6)[,4])),c(2,3,4)]),
col="mediumorchid2", density=40, border="black")
#Queensland
points(s3d$xyz.convert(subset(MDS, AreaCode=="AU-QLD")[,c(2,3,4)]), pch=21, col="black", bg="turquoise", cex=1)
text(s3d$xyz.convert(subset(MDS, AreaCode=="AU-QLD")[,c(2,3,4)]), labels="Queensland", col="black", cex=0.7, pos=1)
#Hawaii
points(s3d$xyz.convert(subset(MDS, AreaCode=="I-HW")[,c(2,3,4)]), pch=21, col="black", bg="yellow", cex=1)
text(s3d$xyz.convert(subset(MDS, AreaCode=="I-HW")[,c(2,3,4)]), labels="Hawaii", col="black", cex=0.7, pos=3)
#New Caldedonia
points(s3d$xyz.convert(subset(MDS, h11c==11)[,c(2,3,4)]), pch=21, col="black", bg="darkgreen", cex=1)
text(s3d$xyz.convert(subset(MDS, h11c==11)[,c(2,3,4)]), labels="New Caledonia", col="black", cex=0.7, pos=3)
dev.off()
