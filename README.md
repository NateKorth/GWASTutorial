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

```

Download the genotype data from the following location and transfer to Genotype folder:
```
https://figshare.com/ndownloader/files/22670489
```

This is an imputed SNP file from: Miao, C., Xu, Y., Liu, S., Schnable, P.S. and Schnable, J.C., 2020. Increased power and accuracy of causal locus identification in time series genome-wide association in sorghum. Plant physiology, 183(4), pp.1898-1909.

Once the genotype is in your Genotype folder, let's generate a list of all the lines we have genotype informations for:
```
head -1 SAP_imputed.hmp > hmpHeader.tsv
```

next grab the phenotype file:
make sure it ends up in your input folder
```
curl -o PhenotypesRaw.csv https://raw.githubusercontent.com/NateKorth/GWASTutorial/main/nir_SC_Compiled_Rhodes2014.csv
```
### Calculate BLUPs using Sommer package
First download the sommer package in R v4.0
```
ml R/4.0
R
>install.packages("sommer")
```
Use the R script: CalculateBLUEs.R to calculate BLUEs for the phenotypes
Within this script you'll import the list of genotypes to filter out any lines we don't have genetic information for.
And you'll generate a list of sorghum lines in the phenotype to filter the genotype.

We need to make sure the Lines in the phenotype file and genotype file match 





