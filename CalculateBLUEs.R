library(dplyr)
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(sommer)

df1<-read.csv("../input/PhenotypesRaw.csv")

#Remove any lines without a PI identifier
df2<-subset(df1,PI!="NA")

#We want our genotype and phenotype files to match, in this script we filter the phenotype file for only lines
#in the genotype file, and generate a list of lines in the genotype file not present in the phenotype:

#Import the list of of genotypes:
LinesInHmp<-read.delim("../input/geno/hmpHeader.tsv",header = FALSE)
LinesInHmp2<-as.data.frame(t(LinesInHmp[12:369]))

#While we're thinking about lines and their names...
#If you peak at the data sets you see the nominclature in the genotype file contains underscores while the phenotype file does not (Typical)
#Lets just add those underscores in the phenotype file:

df2$PI<-sub("PI","PI_",df2$PI)

#Remove lines from the phenotype file that aren't in the genotype:
df3<-df2[which(df2$PI %in% LinesInHmp2$V1),]

#Now let's make and export a New list of all the lines in the phenotype file:
LinesInPheno<-unique(df3$PI)
write.table(LinesInPheno,"LinesInPheno.csv",col.names=FALSE,row.names=FALSE,quote=FALSE)

#extract list of traits for testing:
Traits<-as.data.frame(names(df3[6:10]))
names(Traits)<-"Traits"

#We'll use year as a random effect to calculate BLUES, we'll need to change it from a numeric to a character
df3$Year<-as.character(df3$Year)

#Calculate BLUEs for single trait:
fitB <-mmer(Phenols~PI, random=~Year, rcov=~units, data=df3)
#Let's look at the fit of our model:
summary(fitB)

#Use the predict.mmer function to calculate BLUEs:

BLUEs2<-predict.mmer(object = fitB,classify = "PI")
BLUEs3<-BLUEs2$pvals[,2:3]

#This will be the combined file for all of our BLUEs:
BLUEsdf<-as.data.frame(BLUEs3)
names(BLUEsdf)<-c("PI","Phenols")
BLUEsdf$PI<-BLUEs3$PI

#A function to calculate BLUEs (useful when screeing a ton of traits
getBLUE <- function(trait){
  f <- formula(paste0(trait, '~PI'))
  fit <- mmer(f , random=~Year, rcov=~units, data=df3)
  BLUEs2<-predict.mmer(object = fit,classify = "PI")
  BLUEs3<-BLUEs2$pvals[,3]
  names(BLUEs3)<-paste0(trait)
  return(BLUEs3)
}

#A for loop to apply our function to our list of traits

for (i in Traits$Traits){
  BLUEsdf<-cbind(BLUEsdf,getBLUE(i))
}

names(BLUEsdf)<-c("PI","Phenols",Traits$Traits)

#Export your BLUEs:
write.table(BLUEsdf,"Phenotypes.csv",row.names=FALSE,quote=FALSE,sep = ",")

#df3[13:15] are 3 additional traits that won't work if try to apply the function to them. This is because for some values of "PI" there is missing data (No data for any Reps)
#Extra credit if you fix these and add them to your final phenotype file.

