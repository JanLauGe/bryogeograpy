setwd("D:/Arbeit/Data")
AllSpec <- read.table("D:/Arbeit/Data/AllSpec.txt", sep="\t", header=T, dec=",")

Status <- cbind(AllSpec[,1:2])
Status <- unique(Status)
SumStat <- summary(Status[,2])

  x <- barplot(SumStat)
  y <- SumStat
  barplot(SumStat, col=c("red", "orange", "darkgreen", "grey"),
      border=T, xlab="Taxonomischer Status", ylab="Anzahl Arten",
      names.arg=c("illegitim", "ung?ltig", "legitim", "Keine Angabe"),
      ylim=c(0,10000))
  text(x,y+300, labels=SumStat)


Rank <- cbind(AllSpec[,1],AllSpec[,3])
Rank <- as.data.frame(unique(Rank))
Rank[,2] <- as.factor(Rank[,2])
SumRank <- summary(Rank[,2])

  x <- barplot(SumRank)
  y <- SumRank

par(bty="n", col.axis="white")
  barplot(SumRank, col=c("red", "orange", "yellow", "green", "blue", "grey"),
      border=T, xlab="Taxonomischer Status", ylab="Anzahl Arten",
      names.arg=c("Rang 0", "Rang 1", "Rang 2", "Rang 3", "Rang 4", "NA"),
      ylim=c(0,6000), xlim=c(0,8))
  text(x,y+300, labels=SumRank)
par()
