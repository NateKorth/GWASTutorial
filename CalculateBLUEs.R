library(sommer)

df1<-read.csv("../input/nir_SC_Compiled_Rhodes2014_Korth2020.csv")

#Remove any lines without a PI identifier
df2<-subset(df1,PI!="NA")

#We want our genotype and phenotype files to match, in this script we filter the phenotype file for only lines
#in the genotype file, and generate a list of lines in the genotype file not present in the phenotype:

#Import the list of of genotypes:
LinesInHmp<-read.delim("../input/geno/hmpHeader.tsv",header = FALSE)
LinesInHmp2<-as.data.frame(t(LinesInHmp[12:369]))
names(LinesInHmp2)<-"PI"

#While we're thinking about lines and their names...
#If you peak at the data sets you see the nominclature in the genotype file contains underscores while the phenotype file does not (Typical)
#Lets just add those underscores in the phenotype file:

df2$PI<-sub("PI","PI_",df2$PI)

#Remove lines from the phenotype file that aren't in the genotype:
df3<-df2[which(df2$PI %in% LinesInHmp2$PI),]

#Now let's make and export a New list of all the lines in the phenotype file:
#You may have generated the table in an earlier step
LinesInPheno<-unique(df3$PI)
#write.table(LinesInPheno,"LinesInPheno.csv",col.names=FALSE,row.names=FALSE,quote=FALSE)

#extract list of traits for testing:
Traits<-as.data.frame(names(df3[6:9]))
names(Traits)<-"Traits"

#We'll use year as a random effect to calculate BLUES, we'll need to change it from a numeric to a character
df3$Year<-as.character(df3$Year)

#Calculate BLUEs for single trait:
#TODO
fitB <-mmer(Tannins~PI, random=~(rand1+rand2+rand3), rcov=~units, data=df3)
#Let's look at the fit of our model:
summary(fitB)

#Use the predict.mmer function to calculate BLUEs:

BLUEs2<-predict.mmer(object = fitB,classify = "PI")
BLUEs3<-BLUEs2$pvals[,2:3]

#This will be the combined file for all of our BLUEs:
BLUEsdf<-as.data.frame(BLUEs3)
names(BLUEsdf)<-c("PI","Tannins")
BLUEsdf$PI<-BLUEs3$PI

#A function to calculate BLUEs (useful when screeing a ton of traits)
#TODO add the proper random effects to the model
getBLUE <- function(trait){
  f <- formula(paste0(trait, '~PI'))
  fit <- mmer(f , random=~rand1+rand2+rand3, rcov=~units, data=df3)
  BLUEs2<-predict.mmer(object = fit,classify = "PI")
  BLUEs3<-BLUEs2$pvals[,3]
  names(BLUEs3)<-paste0(trait)
  return(BLUEs3)
}

#A for loop to apply our function to our list of traits

for (i in Traits$Traits){
  BLUEsdf<-cbind(BLUEsdf,getBLUE(i))
}

names(BLUEsdf)<-c("PI","Tannins",Traits$Traits)

#Sort the BLUEs by lines in the genotype file:
LinesInPheno<-as.data.frame(LinesInPheno)
names(LinesInPheno)<-"PI"
LinesInHmp3<-LinesInHmp2[which(LinesInHmp2$PI %in% LinesInPheno$PI),]
BLUEsdf<-BLUEsdf[match(LinesInHmp3,BLUEsdf$PI),]

#Export your BLUEs:
write.table(BLUEsdf,"Phenotypes.csv",row.names=FALSE,quote=FALSE,sep = ",")

#df3[10:15] are additional traits that won't work if try to apply the function to them. This is because for some values of "PI" there is missing data (No data for any Reps) and or some traits don't have value for a specific environment or year...
#Extra credit if you fix these, add them to your final phenotype file and conduct GWAS on them.
