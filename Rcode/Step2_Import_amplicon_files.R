## 
## QIIME2-qza files to microtable object
## 


###################################################
# load packages
library(microeco)
library(magrittr)
library(file2meco)
library(readxl)
###################################################
# create output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage2_amplicon_microtable"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
###################################################

# 16S QIIME2-qza directory path
tmp <- "Input/1.Amplicon/QIIME2_qza"

# use a temporary name (tmp_) for convenience
# If phylogenetic tree or representative fasta is not available for user, please set: phylo_tree = NULL and rep_fasta = NULL
tmp_microtable <- qiime2meco(file.path(tmp, "feature_table_filter.qza"),
	taxonomy_table = file.path(tmp, "taxonomy.qza"),
	phylo_tree = file.path(tmp, "rooted_tree.qza"),
	rep_fasta = file.path(tmp, "sequences_filter.qza")
	)

# read sample information into R
tmp_sample <- as.data.frame(read_excel("Input/1.Amplicon/Sample_info.xlsx", col_names = TRUE), stringsAsFactors = FALSE)
# row names are necessary and should be in accordance with those in feature abundance table
rownames(tmp_sample) <- tmp_sample[, 1]
# assign the table into sample_table of tmp_microtable
tmp_microtable$sample_table <- tmp_sample

# tidy all the data inside the object to ensure that samples and features correspond across different data
tmp_microtable$tidy_dataset()

# use clone function to copy the object completely
amplicon_16S_microtable_raw <- clone(tmp_microtable)

# use save function to save the microtable object to output directory
save(amplicon_16S_microtable_raw, file = file.path(output_dir, "amplicon_16S_microtable_raw.RData"), compress = TRUE)







