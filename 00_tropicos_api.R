##################################################
# Script to check for accepted names on Tropicos #
# using API URL and XML extraction               #
# ---------------------------------------------- #
# Jan Laurens Geffert, 8.4.2011                  #
##################################################

require(XML)    #for XML structure
require(RCurl)  #for getURL

#Define function to generate URLs and extract XML content
getacnames <- function(species)
    {TropURL1 <- paste("http://services.tropicos.org/Name", sep="")
    TropURL2 <- paste("AcceptedNames?apikey=28fe69d0-3bb5-4782-9454-a5ce15f426c9&format=xml", sep="")

    results <- data.frame(SynID=character(0), SynName=character(0), SynFullName=character(0), SynFam=character(0), AccID=character(0), AccName=character(0), AccFullName=character(0), AccFam=character(0), RefID=character(0))

    for (i in 1:length(species))
        {doc <- xmlTreeParse(getURL(URLencode(paste(TropURL1,species[i],TropURL2, sep="/"))))
        r <- xmlRoot(doc)
        
        for (j in 1:xmlSize(r))   
            {if (xmlName(r[[j]]) == "Synonym")
                {if (xmlName(r[[j]][[1]]) == "SynonymName")
                    {hit <- data.frame(xmlValue(r[[1]][[1]][[1]][[1]]), xmlValue(r[[1]][[1]][[2]][[1]]), xmlValue(r[[1]][[1]][[3]][[1]]), xmlValue(r[[1]][[1]][[4]][[1]]),
                    xmlValue(r[[1]][[2]][[1]][[1]]), xmlValue(r[[1]][[2]][[2]][[1]]), xmlValue(r[[1]][[2]][[3]][[1]]), xmlValue(r[[1]][[2]][[4]][[1]]),
                    xmlValue(r[[1]][[3]][[1]][[1]]))
                    colnames(hit) <- c("SynID", "SynName", "SynFullName", "SynFam", "AccID", "AccName", "AccFullName", "AccFam", "RefID")
            results <- rbind(results, hit)}}}}
        return(results)}

#Actual extraction
#based on "speclist" vector with NameIds
#split into 1k chunks for security and management

setwd("C:/Users/Laurens/Desktop/Tropicos Check")
species <- scan("Speclist.txt", sep=";", what="list")    
    
syns1 <- getacnames(species[1:1000])
syns2 <- getacnames(species[1001:2000])
syns3 <- getacnames(species[2001:3000])
syns4 <- getacnames(species[3001:4000])
syns5 <- getacnames(species[4001:5000])
syns6 <- getacnames(species[5001:6000])
syns7 <- getacnames(species[6001:7000])
syns8 <- getacnames(species[7001:8000])
syns9 <- getacnames(species[8001:9000])
syns10 <- getacnames(species[9001:10000])
syns11 <- getacnames(species[10001:11000])
syns12 <- getacnames(species[11001:12000])
syns13 <- getacnames(species[12001:13000])
syns14 <- getacnames(species[13001:14000])

syns <- rbind(syns1,syns2,syns3,syns4,syns5,syns6,syns7,syns8,syns9,syns10,syns11,syns12,syns13,syns14)

write.table(syns, "results.txt", sep="\t", quote=F)