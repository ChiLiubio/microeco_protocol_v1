

# install microeco package 
install.packages("microeco")

# install file2meco package for the file reading
install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())

# generate some combined graph
install.packages("aplot")

# multiple facets in barplot
install.packages("ggh4x")

# save plot to local file
install.packages("cowplot")

# statistics
install.packages("agricolae")

# generate significant letters
install.packages("FSA")
install.packages("rcompanion")

# evaluate R2 in mixed model
install.packages("performance")


# machine learning
install.packages("Boruta")
install.packages("rsample")
install.packages("randomForest")
install.packages("rfPermute")
install.packages("multiROC")
install.packages("kernlab")
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

# some tree-like plot 
BiocManager::install("ggtree")


# generate simulated community data
# not available for Bioconductor version '3.19'
BiocManager::install("SparseDOSSA2")




######################################################################
######################################################################
# Trouble shooting for package downloading: The internet speed is too slow, or repository connection fails due to the firewall
# please find a mirror name in https://cran.r-project.org/mirrors.html
mirror = "https://mirrors.tuna.tsinghua.edu.cn"
options("repos" = c(CRAN = paste0(mirror, "/CRAN/")))
options(BioC_mirror = paste0(mirror, "/bioconductor/"))
######################################################################
######################################################################

