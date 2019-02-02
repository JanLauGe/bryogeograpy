##################################
# R Script global moss diversity #
##################################
# Jan Laurens Geffert            #
# 2.2.2011                       #
##################################

require(lmtest)   #for coeftest
require(sandwich) #for coeftest
require(foreign)  #for dbf
require(plotrix)  #for corner.label 
require(rgl)      #for colors
require(Hmisc)
require(BiodiversityR)
setwd("D:/Arbeit/Data")

#read revised table of SpecDiv and predictor variables for OGU subset
Data <- read.table("D:/Arbeit/Data/Modell/GUs.txt", sep="\t", header=T, dec=",", na.strings="", row.names="Name")
detach(Data)
attach(Data)
str(Data)

##########################
# Exploratory Analysis   #
##########################

############
# Barplots #
############

#Biom-Barplot as pdf
pdf("bp.pdf")
  par(mar = c(15,6,4,2) + 0.1)
    plot(Div~Biom, xlab="", ylab="Dokumentierte Artenzahl", axes=F,
    col=c(rgb(0,255,0, max=255), rgb(178,229,25, max=255), rgb(173,255,47, max=255), rgb(100,255,138, max=255), rgb(0,200,15, max=255), rgb(60,179,113, max=255), rgb(240,234,150, max=255), rgb(255,215,0, max=255), rgb(210,180,140, max=255), rgb(75,225,170, max=255), rgb(255,0,0, max=255), rgb(240,128,128, max=255)))
    axis(1, at=1:12, labels=F)
    axis(2)
    text(x=1:12,y=-100, srt = 45, adj = 1, labels=c("Tropische Feuchtwälder", "Tropische Trockenwälder", "Tropische Nadelwälder", "Gemäßigte Mischwälder", "Gemäßigte Nadelwälder", "Boreale Wälder", "Tropisches Grasland", "Gemäßigtes Grasland", "Montanes Grasland", "Arktische Tundra", "Mediterrane Sklerophyllwälder", "Wüsten und Trockengebiete"), xpd = TRUE)
    box(lty=1, col="black")
    
  par(mar = c(15,6,4,2) + 0.1)
    plot(StdDiv~Biom, xlab="", ylab="Standardisierte Artenzahl", axes=F,
    col=c(rgb(0,255,0, max=255), rgb(178,229,25, max=255), rgb(173,255,47, max=255), rgb(100,255,138, max=255), rgb(0,200,15, max=255), rgb(60,179,113, max=255), rgb(240,234,150, max=255), rgb(255,215,0, max=255), rgb(210,180,140, max=255), rgb(75,225,170, max=255), rgb(255,0,0, max=255), rgb(240,128,128, max=255)))
    axis(1, at=1:12, labels=F)
    axis(2)
    text(x=1:12,y=-100, srt = 45, adj = 1, labels=c("Tropische Feuchtwälder", "Tropische Trockenwälder", "Tropische Nadelwälder", "Gemäßigte Mischwälder", "Gemäßigte Nadelwälder", "Boreale Wälder", "Tropisches Grasland", "Gemäßigtes Grasland", "Montanes Grasland", "Arktische Tundra", "Mediterrane Sklerophyllwälder", "Wüsten und Trockengebiete"), xpd = TRUE)
    box(lty=1, col="black")
dev.off()   

#Continent and floristic kingdom barplots as pdf
  pdf("bp2.pdf")
  par(font=2, ps=17, mar=c(4,4,1,1), mgp=c(2.5,1,0))
  layout(matrix(c(1:6),3,2, byrow=T))
    plot(Div~Kont, xlab="", ylab="Dokumentierte Artenzahl", col="grey")
    plot(Div~King, xlab="", ylab="", col="grey")
    plot(StdDiv~Kont, xlab="Kontinente", ylab="Standardisierte Artenzahl", col="grey")
    plot(StdDiv~King, xlab="Florenreich", ylab="", col="grey")
dev.off()

#Area Classes
  pArea <- Area/1000
  brks <- c(0,100,200,seq(300, 3000, 300),3500)
  ArrClass <- cut(pArea, breaks=brks)
  par(mar = c(10,6,4,2) + 0.1)
  plot(ArrClass, xlab="Gebietsgröße in 1000km²", ylab="Anzahl Gebietseinheiten", axes=F)
  mtext(brks, side=1, at=) 

#Correlation table
cortab <- as.matrix(Data[9:49])
  cr <- rcorr(cortab, type="pearson")
  cr2 <- cbind(cr$r, cr$n, cr$P)
  write.dbf(cr2, file="D:/Arbeit/Data/correl.dbf")

#Latitude Plots
pdf("lat.pdf")
  par(font=2, ps=17)
  plot(Div~absLat, xlim=c(0,70), ylim=c(0,900), xlab="Absoluter Breitengrad", 
       ylab="Dokumentierte Artenzahl", cex=log2(Area)/10)
       lines(lowess(Div~absLat, f=0.5), col="red")

  plot(StdDiv~absLat,xlim=c(0,70), ylim=c(0,900), xlab="Absoluter Breitengrad", 
       ylab="Standardisierte Artenzahl", cex=log2(Area)/10)
       lines(lowess(StdDiv~absLat, f=0.5), col="red")
dev.off()
    
#####################
# Correlation Plots #    
#####################

pdf("cor1.pdf", paper="a4")
  layout(matrix(c(1:15),5,3, byrow=T))
    par(mar=c(4,3,1,1), mgp=c(2.5,1,0))
    for(i in 16:30)
        {correl <- cor.test(Div,Data[,i])
        Signif <- symnum(correl$p.value, corr = FALSE, na = FALSE,
              cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols=c("***", "**", "*", ".", " "))
        plot(Div~Data[,i], ylab="", xlab=colnames(Data[i]))
        lines(lowess(Div~Data[,i], f=0.5), col="red")
        corner.label(label=Signif, x=-1, y=1, col.lab=2)}                 
dev.off()
 
pdf("cor2.pdf", paper="a4")
  layout(matrix(c(1:15),5,3, byrow=T))
    par(mar=c(4,3,1,1), mgp=c(2.5,1,0))
    for(i in c(14, 15, 31:43))
        {correl <- cor.test(Div,Data[,i])
        Signif <- symnum(correl$p.value, corr = FALSE, na = FALSE,
              cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols=c("***", "**", "*", ".", " "))
        plot(Div~Data[,i], ylab="", xlab=colnames(Data[i]))
        lines(lowess(Div~Data[,i], f=0.5), col="red")
        corner.label(label=Signif, x=-1, y=1, col.lab=2)}                 
dev.off() 
 
pdf("cor3.pdf", paper="a4")
  layout(matrix(c(1:15),5,3, byrow=T))
    par(mar=c(4,3,1,1), mgp=c(2.5,1,0))
    for(i in 44:49)
        {correl <- cor.test(Div,Data[,i])
        Signif <- symnum(correl$p.value, corr = FALSE, na = FALSE,
              cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols=c("***", "**", "*", ".", " "))
        plot(Div~Data[,i], ylab="", xlab=colnames(Data[i]))
        lines(lowess(Div~Data[,i], f=0.5), col="red")
        corner.label(label=Signif, x=-1, y=1, col.lab=2)}                 
dev.off()
 

################  
##Species-Area #
################
#
#    plot(log(Div)~log(Area))
#      abline(lm(log(Div)~log(Area)))
#
#  #Outlier in Arten-Fläche-Plots
#    SpecArea <- cbind(Code, Div, Area, log(Div), log(Area))
#    colnames(SpecArea) <- c("SCode", "SDiv", "SArea", "logDiv", "logArea")
#    SpecArea <- as.data.frame(SpecArea)
#    attach(SpecArea)
#  
#    plot(logDiv ~ logArea, xlab="Log Area [km²]", ylab="Log Species number")
#    abline(lm(logDiv ~ logArea))
#  #compute conf intervals (default 95%, vary with "level = xxx")
#    confin <- predict(lm(logDiv ~ logArea), interval="confidence", level=0.95, type="response")
#  #compute standard error
#    stander <- predict(lm(logDiv ~ logArea), interval="confidence", level=0.95, type="response", se.fit=T)
#
#  interv <-cbind(SpecArea, stander$fit ,stander$se.fit)
#    interv<-interv[order(interv[,5]),] # bring in order
#    lines(interv[,5], interv[,6], lwd=1, lty=2) # plot upper 95% conf int
#    lines(interv[,5], interv[,6], lwd=1, lty=2) # plot lower 95% conf int
#    lines(interv[,5], (interv[,6]+ 4*interv[,9]), lwd=1, lty=4) #plot fit + 4x standard error
#    lines(interv[,5], (interv[,6]- 4*interv[,9]), lwd=1, lty=4) #plot fit - 4x standard error
#
## Formel und R² für die Regression mit allen Datenpunkten in Grafik schreiben lassen
#  text(2.1, 2.99 , "all records, dotted line", pos=4)
#    Beschr <- paste('Log S = ',round((coef(ModEU)[1]),3),' + ',round((coef(ModEU)[2]),3),'Log A')
#  text(2.1, 2.95 , Beschr, pos=4, cex = 0.8)
#  R2 <- paste ('R² = ',round(summary(ModEU)$r.squared, 2),'***')
#  text(2.1, 2.9, R2, pos = 4, cex = 0.8)
#
## Ausreißer, die mehr als 4 standard errors nach unten abweichen, in Grafik beschriften lassen
#  MatrEUKomb<- cbind(MatrEU, conf[,2:5])
#  MatrEUAusr<-subset(MatrEUKomb[(MatrEUKomb$LogS < (MatrEUKomb$fit - (4 *MatrEUKomb[,38]))),] )
#    text(MatrEUAusr$LogA+0.01, MatrEUAusr$LogS, MatrEUAusr$Unit, pos=4, cex = 0.8)
#
## Regression noch einmal ohne Ausreißer rechnen und plotten
#  MatrEUrest<-subset(MatrEUKomb[(MatrEUKomb$LogS >= (MatrEUKomb$fit - (4 *MatrEUKomb[,38]))),] )
#  summary(ModEUrest <- lm(MatrEUrest$LogS ~ MatrEUrest$LogA)); AIC(ModEUrest) # LogS, LogA
#   abline(ModEUrest, lty=1)
#  text(2.1, 2.8 , "without outlier, solid line", pos=4)
#    Beschr <- paste('Log S = ',round((coef(ModEUrest)[1]),3),' + ',round((coef(ModEUrest)[2]),3),'Log A')
#  text(2.1, 2.75 , Beschr, pos=4, cex = 0.8)
#  R2 <- paste ('R² = ',round(summary(ModEUrest)$r.squared, 2),'***')
#  text(2.1, 2.7, R2, pos = 4, cex = 0.8)
#dev.off()
#
#
#
##glmulti to test for good combination of predictor variables
#testmodel <- glmulti(y="Div", xr=c("Area", "Alt", "Wet", "WetMax", "Cloud", "Bio15", "Bio18",
#"GGBiom", "EcoNum", "Frahm", "PET", "AET", "Treec"), data=Data,
#maxsize= 5, level=1, method="h", fitfunction=glm, name="mossglm", intercept=T, crit="aic", report=T, plotty=T)
#

#######################################
# Modelling species richness for OGUs #
# based on environmental predictors   #
#######################################

  model <- glm(Div~log(Area)+Alt+WetMax, data=Data, family=quasipoisson())

  model1 <- glm(Div~log(Area)+Alt+WetMax+Bio15, data=Data, family=quasipoisson())
  model2 <- glm(Div~log(Area)+Alt+WetMax+Bio18, data=Data, family=quasipoisson())
  model3 <- glm(Div~log(Area)+Alt+WetMax+Bio15+Snow, data=Data, family=quasipoisson())
  model4 <- glm(Div~log(Area)+Alt+WetMax+Bio18+Snow, data=Data, family=quasipoisson())
  model5 <- glm(Div~log(Area)+Alt+WetMax+Snow, data=Data, family=quasipoisson()) 
  model6 <- glm(Div~Frahm+log(Area)+Alt+WetMax+PET, data=Data, family=quasipoisson())
  model7 <- glm(Div~Frahm+log(Area)+Alt+WetMax+Dry, data=Data, family=quasipoisson())
  model8 <- glm(Div~Frahm+log(Area)+Alt+WetMax+Snow, data=Data, family=quasipoisson())
  model9 <- glm(Div~Frahm+log(Area)+Alt+WetMax+PET+Dry+Snow, data=Data, family=quasipoisson())
 
#Test for residual correlation
  residual <- residuals.glm(model1, type="response")  
    rescor <- cor(residual,Data[,7:49])
    rescor <- as.data.frame(rescor)
    sort(rescor)
  
#Write residuals to DBF
    pred1 <- predict.glm(model1, newdata=Data, type="response")
    pred2 <- predict.glm(model2, newdata=Data, type="response")
    
    rn <- cbind(Data, pred1, pred2, residual1, residual2)
    plot(rn$Lon, rn$Lat, cex=abs(residual1)/100)
    
    write.dbf(rn, file="D:/Arbeit/Data/residuals3.dbf")  
  
#Evaluate model fit
    summary.glm(model1)
    
    anova.glm(model1)
    coeftest(model1, vcov = sandwich)
    deviancepercentage(model1, test="Chisq", data=Data, digits=5)

    plot(model <- glm(Div~Frahm+log(Area)+Alt+WetMax+Dry, data=Data, family=quasipoisson()))
    sort(residuals.glm(model, "response"))
        
    plot(residual~Div)
    abline(lm(residual~Div))
    identify(residual~Div, labels=Code)


    
########################
# Prediction for Grid: #
########################

  grid100 <- read.table("D:/Arbeit/Data/Modell/100Grid/GRID.txt", sep="\t", header=T, dec=",", na.strings="")
  SpecRich1 <- predict.glm(model1, newdata=grid100, type="response")
  SpecRich2 <- predict.glm(model2, newdata=grid100, type="response")
  SpecRich3 <- predict.glm(model3, newdata=grid100, type="response")
  SpecRich4 <- predict.glm(model4, newdata=grid100, type="response")
  SpecRich5 <- predict.glm(model5, newdata=grid100, type="response")
  SpecRich6 <- predict.glm(model6, newdata=grid100, type="response")
  SpecRich7 <- predict.glm(model7, newdata=grid100, type="response")
  SpecRich8 <- predict.glm(model8, newdata=grid100, type="response")
  SpecRich9 <- predict.glm(model9, newdata=grid100, type="response")
    
#Write results   
  SpecRich <- cbind(grid100$Id, round(SpecRich1), round(SpecRich2), round(SpecRich3), round(SpecRich4))
  , round(SpecRich5), round(SpecRich6), round(SpecRich7), round(SpecRich8), round(SpecRich9))
  colnames(SpecRich) <- c("Id", "AAW", "FAAW", "AAWP", "AAWD", "AAWS", "FAAWP","FAAWD", "FAAWS", "FAAWPDS")
 
  Results <- cbind(grid100$Id, round(SpecRich1), round(SpecRich2))
  colnames(Results) <- c("Id", "AAWDS", "AAWDPS")
  
  write.dbf(SpecRich, file="D:/Arbeit/Data/predictions_30.dbf")
