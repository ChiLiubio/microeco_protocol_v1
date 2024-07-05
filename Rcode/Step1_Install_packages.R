

# install CRAN packages
install.packages("microeco")


# machine learning
install.packages("caretEnsemble")


# install Bioconductor packages
install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
install.packages("MicrobiomeStat", repos = BiocManager::repositories())
BiocManager::install("metagenomeSeq")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")



BiocManager::install("SparseDOSSA2")








# Trouble shooting for package downloading: slow or repository connection failure
# find a mirror in https://cran.r-project.org/mirrors.html
mirror = "https://mirrors.tuna.tsinghua.edu.cn"
options("repos" = c(CRAN = paste0(mirror, "/CRAN/")))
options(BioC_mirror = paste0(mirror, "/bioconductor/"))



