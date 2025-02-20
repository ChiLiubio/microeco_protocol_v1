
# allow more waiting time to download each package
options(timeout = 1000)


# install microeco package 
install.packages("microeco")

# install file2meco package for the file reading
install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())

# invoke some visualization functions
install.packages("ggpubr")

# arrange and align several plots on one page
install.packages("aplot")

# multiple facets in barplot
install.packages("ggh4x")

# UniFrac distance
install.packages("GUniFrac")

# for some statistics
install.packages("agricolae")

# generate significant letters in the post hoc test
install.packages("FSA")
install.packages("rcompanion")

# evaluate R2 in mixed-effects model
install.packages("performance")

# visualize correlation network of metabolomic data in R
install.packages("ggraph")

# machine learning
## feature selection
install.packages("Boruta")
## split samples
install.packages("rsample")
# random forest model
install.packages("randomForest")
# Support Vector Machine model
install.packages("kernlab")
# eXtreme Gradient Boosting model
install.packages("xgboost")
# determine significance of features in random forest model
install.packages("rfPermute")
# ROC curve
install.packages("multiROC")
# for the analysis of multiple models
install.packages("caretEnsemble")


# differential test methods
install.packages("MicrobiomeStat", repos = BiocManager::repositories())
BiocManager::install("metagenomeSeq")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("Maaslin2")

# for data conversion
BiocManager::install("phyloseq")
BiocManager::install("microbiome")

# some tree-like plot 
BiocManager::install("ggtree")

# PLS-DA
BiocManager::install("ropls")

# use SparseDOSSA2 package to simulate community data
# It is not available for Bioconductor version '3.19' when run: BiocManager::install("SparseDOSSA2")
# Install the package from source in GitHub https://github.com/biobakery/SparseDOSSA2
# first install some dependent packages
install.packages("ks")
install.packages("huge")
install.packages("future.apply")
install.packages("truncnorm")
install.packages("Rmpfr")
install.packages("extdata/SparseDOSSA2_0.99.2.tar.gz", repos = NULL, type = "source")

# for rarefaction curve
install.packages("extdata/mecodev_0.2.0.tar.gz", repos = NULL, type = "source")


######################################################################
######################################################################
# Test whether each package has been correctly installed
library(microeco)
if(packageVersion("microeco") < '1.13.0'){stop("Minimum version of microeco package should be 1.13.0! Current version is ", packageVersion("microeco"), " ! Please reinstall it!")}
library(file2meco)
if(packageVersion("file2meco") < '0.9.0'){stop("Minimum version of file2meco package should be 0.9.0! Current version is ", packageVersion("file2meco"), " ! Please reinstall it!")}
library(ANCOMBC)
if(packageVersion("ANCOMBC") < '2.8.1'){stop("Minimum version of ANCOMBC package should be 2.8.1! Current version is ", packageVersion("ANCOMBC"), " ! Please reinstall it from Bioconductor!")}
library(ggpubr)
library(aplot)
library(ggh4x)
library(cowplot)
library(GUniFrac)
library(agricolae)
library(FSA)
library(rcompanion)
library(performance)
library(ggraph)
library(Boruta)
library(rsample)
library(randomForest)
library(rfPermute)
library(multiROC)
library(kernlab)
library(xgboost)
library(caretEnsemble)
library(MicrobiomeStat)
library(metagenomeSeq)
library(ALDEx2)
library(edgeR)
library(DESeq2)
library(Maaslin2)
library(phyloseq)
library(microbiome)
library(ggtree)
library(ropls)
library(SparseDOSSA2)
library(mecodev)










