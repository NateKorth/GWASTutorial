# GWAS Tutorial
This tutorial is specifically designed for students of Complex Biosystems 852 at the university of Nebraska, it relies on a basic knowledge of Unix and R.
In this tutorial we will be analyzing biochemical features of Sorghum biocolor in a diversity panel and scanning for associated SNPs

## Step1: Setup
First things first make a new project folder in your work directory 
```
cd $WORK
mkdir GWASTutorial
mkdir GWASTutorial/Scripts
mkdir GWASTutorial/Output
mkdir GWASTutorial/input
mkdir GWASTutorial/input/Genotype
```

For BLUE calculation we'll need a bit of computational resourses (best not to do this on the work node) 
Politely ask HCC for some resourses
```
srun --nodes=1 --ntasks-per-node=4 --mem-per-cpu=2024 --qos=short --pty $SHELL
```


Download the genotype data from the following location and transfer to Genotype folder:
```
https://figshare.com/ndownloader/files/22670489

#An easy way to do this is (once you've downladed the file) use scp from your command prompt on your computer:
scp .\Downloads\SAP_imputed.hmp.zip USERNAME@crane.unl.edu:/work/GROUP/USERNAME/GWASTutorial/input/Genotype
```

This is an imputed SNP file from: Miao, C., Xu, Y., Liu, S., Schnable, P.S. and Schnable, J.C., 2020. Increased power and accuracy of causal locus identification in time series genome-wide association in sorghum. Plant physiology, 183(4), pp.1898-1909.

We're going to use the pre-imputed file for the sake of time, for more on imputation see Beagle: http://faculty.washington.edu/browning/beagle/beagle.html

Once the genotype is in your Genotype folder, let's generate a list of all the lines we have genotype informations for:
```
head -1 SAP_imputed.hmp > hmpHeader.tsv
```
Let's convert our genotype file to vcf, so in later steps we can use vcftools to work with our data
```
ml tassel/5.2
run_pipeline.pl -Xms512m -Xmx10G -fork1 -h SAP_imputed.hmp -export -exportType VCF

```
next grab the phenotype file:
make sure it ends up in your input folder
```
curl -o PhenotypesRaw.csv https://raw.githubusercontent.com/NateKorth/GWASTutorial/main/nir_SC_Compiled_Rhodes2014_Korth2020.csv
```
We'll need the lines in our Genotype file and Phenotype file to match. There's one hundred ways to do this, we'll do it in R (You can use the command line to access R as shown below, or if you prefer Rstudio feel free to use that).
First we'll remove lines in the phenotype file that are not in the genotype file 
```
cd $WORK/GWASTutorial/Scripts
ml R/4.0
R
df1<-read.csv("../input/PhenotypesRaw.csv")

#Remove any lines without a PI identifier
df2<-subset(df1,PI!="NA")

#Import the list of of genotypes:
LinesInHmp<-read.delim("../input/Genotype/hmpHeader.tsv",header = FALSE)
LinesInHmp2<-as.data.frame(t(LinesInHmp[12:369]))
names(LinesInHmp2)<-"PI"

#While we're thinking about lines and their names...
#If you peak at the data sets you see the nominclature in the genotype file contains underscores while the phenotype file does not (Typical)
#Lets just add those underscores in the phenotype file:

df2$PI<-sub("PI","PI_",df2$PI)

#Remove lines from the phenotype file that aren't in the genotype:
df3<-df2[which(df2$PI %in% LinesInHmp2$PI),]
```

Use df3 for downstream BLUE calculation or output a .csv file that contains the phenotypes filtered.

Now output a list of lines in the genotype file that are not in the phenotype file
```
LinesInPheno<-unique(df3$PI)
write.table(LinesInPheno,"LinesInPheno.csv",col.names=FALSE,row.names=FALSE,quote=FALSE)
```
##Step2: Calculate BLUEs
### Calculate BLUEs using Sommer package
Still within R on HCC

Install the R packages we need to calculate BLUEs and run GWAS we'll design and test a linear model for for GWAS in R using the sommer package 
```
>install.packages("sommer")
>install.packages("rMVP")
```
Next we'll design and test a linear model looking at different factors we may need to include
```
library(sommer)
#We'll use some of the "numerical" data as a random effect including year, to calculate BLUES, we'll need to change it from a numeric to a character
df3$Year<-as.character(df3$Year)
df3$Rep<-as.character(df3$Rep)

#Let's start with Tannins
#Fit a the most complex Linear model you can imagine:
##Format: model <- mmer (Phenotype~FixedEffect, random=~RandomEffects, data=YourDataFrame)

fit <-mmer(Tannins~PI, random=~Year+Env+Rep+Year:Env+Year:Rep+Env:Rep+PI:Year+PI:Env, rcov=~units, data=df3)
#Let's look at the fit of our model:
summary(fit)

#Look at the AIC and BIC values for goodness of fit, and the varComp for how much variance is explained by each Random effect
#Rep has a pretty small percentage of variance, let's remove it and see how it effects our Fit (when you remove Rep you should also remove interaction terms)

fit2 <-mmer(Tannins~PI, random=~Year+Env+Year:Env+PI:Year+PI:Env, rcov=~units, data=df3)
summary(fit2)

#okay the AIC and BIC are slightly larger, that's okay a small increase percentage wise.
#Notice something weird? All the variance-covariance components for Year are the same as Year:Env
#Why is that? There's two environments and 4 years, but it's unbalanced data for the NE Env only 1 year is represented, the other 3 years are all present in Texas
#This is causing the interaction effect to be equivilent to the year effect, let's remove it from our model
#For a very similar reason the the gxe interactions are also confounded, let's remove the term that's lower and see what happens:

fit3 <-mmer(Tannins~PI, random=~Year+Env+PI:Year, rcov=~units, data=df3)
summary(fit3)

#Okay this is looking better (not perfect but things are rarely perfect when working with real biological data)
```

Edit the R script: CalculateBLUEs.R to reflect the model you've come up with (using command nano)
Batch the script as a job using RunR_1.sh

For more details on he sommer package see: https://cran.r-project.org/web/packages/sommer/vignettes/v3.sommer.qg.pdf

## Step3: Genotype filtering

We need to make sure the Lines in the phenotype file and genotype file match, using the "LinesInPheno.csv" we just generated in R, now filter the genotype file:
```
#We're starting to get to bigger jobs, it might be wise to start batching these via .slurm files
ml vcftools
vcftools --vcf SAP_imputed.vcf --out SAP_imputed_Filter1 --keep ../../Scripts/LinesInPheno.csv --recode
```
Next we're going to filter the genotype file minor allele frequency:
```
ml vcftools
vcftools --vcf SAP_imputed_Filter1.recode.vcf --out SAP_imputed_Filter2 --maf 0.05 --recode
```
filter snps with high heterozygosity using bcftools
```
ml bcftools
bcftools filter SAP_imputed_Filter2.recode.vcf --exclude 'F_PASS(GT=="het") > 0.1' -o SAP_imputed_Filter3.vcf
```
## Step.4 GWAS
Almost! We'll be using rMVP, first we'll need to format the genotype data specifically for the program by running the Rscript: PrepGenoForMVP.R
This is getting big enough that we should batch this job in a slurm file (RunR_2.sh)

This script will prepare your genotype file, a map file, a kinship matrix, and the first 3 principle components describing each sorghum genotype. Check your output you should have binary files and .desc (description) files that contain information that describes the location and contents of each outputted file 
```
sbatch RunR_2.sh
```
For more details on rMVP see: 
https://github.com/xiaolei-lab/rMVP

Now we're ready to start a GWAS the basic script to use is RunRMVP.R, which will run mlm and FarmCPU and automatically make figures for you. 
#TODO make a slurm file called RunR_3.sh to batch on hcc, allocate 40 GB on 1 node and this will run pretty fast. 
Have the script submit the R file and move all the output files to the output folder 

Use your favorite file transfer system to move all the output to your computer to view the files.

## Step.5 RMIP
We can use all the same input data from the last step, just make a slurm file to run RunRMIP.R and you'll have output
This script doesn't contain any fancy visualization so we'll have to do that manually, here's some basic code to do it in R:

```
ml R/4.0
R
install.packages("tidyverse)

GWAS1<-read.csv(YOUROUTPUTFILEHERE)

data_cum <- GWAS1 %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHR, bp_add)

GWAS1 <- GWAS1 %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

axis_set <- GWAS1 %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

GWAS2<-cbind(GWAS1,GWASM)

ggplot(GWAS5, aes(x=bp_cum,y=FREQ))+ geom_point(),size=2)+
  scale_x_continuous(label = axis_set2$CHR, breaks = axis_set4$center) + scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_size_continuous(range = c(0.5,3)) + ylab("RMIP Frequency") + xlab(NULL)+unique(length(axis_set3$CHR))))+
  theme_minimal() + theme(panel.border = element_blank(),panel.grid.major.x = element_blank(),plot.margin=unit(c(1,0,0,0),"cm"), 
  panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle = 0, size = 16, vjust = 0.5))+)+ labs(fill='Trait Type')+
  geom_hline(yintercept=10, linetype="dashed",color="black")
```






