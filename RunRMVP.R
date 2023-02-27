library(rMVP)
phenotype1<-read.csv("../input/Phenotypes.csv",head=TRUE)
genotype<-attach.big.matrix("../input/geno/SAP.geno.desc")
map<-read.table("../input/geno/SAP.geno.map", head = TRUE)
kinship<-attach.big.matrix("../input/geno/SAP.kin.desc")
covariates_PC<-bigmemory::as.matrix(attach.big.matrix("../input/geno/SAP.pc.desc"))

#Run a single GWAS
phenotype2<-phenotype1[1:2]
imMVP <- MVP(phe=phenotype2,geno=genotype,map=map,K=kinship,CV.FarmCPU=covariates_PC,priority="speed",maxLoop=10,ncpus=64,method.bin="FaST-LMM",method=c("FarmCPU","MLM"))


#A loop to run all phenotypes

#for(i in 2:ncol(phenotype1)){
#   imMVP <- MVP(phe=phenotype1[,c(1,i)],geno=genotype,map=map,K=kinship,
#    CV.FarmCPU=covariates_PC,priority="speed",maxLoop=10,ncpus=64,
#    method.bin="FaST-LMM",method=c("FarmCPU","MLM")
#  )
#  gc()
#}
