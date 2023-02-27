# GWAS Tutorial
This tutorial is specifically designed for students of Complex Biosystems 852 at the university of Nebraska, it relies on a basic knowledge of Unix and R.
In this tutorial we will be analyzing biochemical features of Sorghum biocolor in a diversity panel and scanning for associated SNPs

## Step1: Setup
First things first make a new project folder in your work directory 
```
cd $WORK
mkdir GWASTutorial
mkdir GWASTutorial/input
mkdir GWASTutorial/input/Genotype
mkdir GWASTutorial/Scripts
mkdir GWASTutorial/Output

```
Next we'll install the R packages we'll need for this pipeline in R on HCC:
```
ml R/4.0
R
>install.packages("sommer")
>install.packages("rMVP")
```

Download the genotype data from the following location and transfer to Genotype folder:
```
https://figshare.com/ndownloader/files/22670489
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
curl -o PhenotypesRaw.csv https://raw.githubusercontent.com/NateKorth/GWASTutorial/main/nir_SC_Compiled_Rhodes2014.csv
```

### Calculate BLUEs using Sommer package

Use the R script: CalculateBLUEs.R to calculate BLUEs for the phenotypes
Within this script you'll import the list of genotypes to filter out any lines we don't have genetic information for.
And you'll generate a list of sorghum lines in the phenotype to filter the genotype.

You can do this is your local R, on HCC in R, or batch it as a job using RunR_1.sh (But I encourage you to first look through and see what the script is doing)

For more details on he sommer package see: https://cran.r-project.org/web/packages/sommer/vignettes/v3.sommer.qg.pdf

### Genotype file filtering

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
bcftools filter SAP_imputed_Filter2.recode.vcf --exclude 'F_PASS(GT=="het") < 0.1' -o SAP_imputed_Filter3.vcf
```
### GWAS
Almost! We'll be using rMVP, first we'll need to format the genotype data specifically for the program by running the Rscript: PrepGenoForMVP.R
This is getting big enough that we should batch this job in a slurm file (RunR_2.sh)

This script will prepare your genotype file, a map file, a kinship matrix, and the first 3 principle components describing each sorghum genotype. Check your output you should have binary files and .desc (description) files that contain information that describes the location and contents of each outputted file 
```
sbatch RunR_2.sh
```
For more details on rMVP see: 
https://github.com/xiaolei-lab/rMVP

Now we're ready to start a GWAS the basic script to use is RunRMVP.R, which will run mlm and FarmCPU and automatically make figures for you. 
#TODO make a slurm file called RunR_3.sh to batch it to hcc 






