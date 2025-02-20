## Amplicon sequencing data analysis
## 
## Import QIIME2 .qza format files to generate microtable object
## 


###################################################
# load packages
library(microeco) # data analysis
library(file2meco) # for data importing
library(magrittr) # for pipe operator
library(readxl) # import data from Excel spreadsheet
###################################################
# create output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage2_amplicon_microtable"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
###################################################

# 16S QIIME2-qza directory path
tmp <- "Input/1.Amplicon/QIIME2_qza"

# In this example, the QIIME2 analysis data and sample information are processed separately, meaning that the QIIME2 analysis results are read first, and then the sample information table located in the Excel file is read and integrated.
#	If users need to read them together, they can provide the path of the sample information table to the sample_table parameter of the qiime2meco function, which supports four types of data input:
#		1) q2-type tab seperated file of QIIME2; 2) comma seperated file with the suffix csv or tab seperated file with suffix tsv or txt; 3) Excel type file with the suffix xlsx or xls; 4) data.frame object in R
# If the sample names in the sequence analysis do not match the names in the sample information table, a table can be provided to the match_table parameter for name replacement. For specific instructions, please refer to the help documentation of the qiime2meco function.
# use a temporary name (tmp_) for convenience
tmp_microtable <- qiime2meco(
	feature_table = file.path(tmp, "feature_table_filter.qza"),
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
# The benefits of saving data to a local file in .RData format: 1) backing up data to a local folder; 2) facilitating the quick reloading of data into R in the future analysis.


##################################################################
##################################################################
# Other examples of data reading with QIIME2 qza files
# all files input, including Sample_info.xlsx
tmp_microtable <- qiime2meco(
	feature_table = file.path(tmp, "feature_table_filter.qza"),
	sample_table = "Input/1.Amplicon/Sample_info.xlsx",
	taxonomy_table = file.path(tmp, "taxonomy.qza"),
	phylo_tree = file.path(tmp, "rooted_tree.qza"),
	rep_fasta = file.path(tmp, "sequences_filter.qza")
	)
tmp_microtable$tidy_dataset()

# If phylogenetic tree or representative fasta is not available, please delete those parameters or set: phylo_tree = NULL, rep_fasta = NULL
tmp_microtable2 <- qiime2meco(
	feature_table = file.path(tmp, "feature_table_filter.qza"),
	sample_table = "Input/1.Amplicon/Sample_info.xlsx",
	taxonomy_table = file.path(tmp, "taxonomy.qza")
	)
tmp_microtable2$tidy_dataset()

# Provided sample_table is a data.frame object when using the qiime2meco function
tmp_sample <- as.data.frame(read_excel("Input/1.Amplicon/Sample_info.xlsx", col_names = TRUE), stringsAsFactors = FALSE)
rownames(tmp_sample) <- tmp_sample[, 1]

tmp_microtable3 <- qiime2meco(
	feature_table = file.path(tmp, "feature_table_filter.qza"),
	sample_table = tmp_sample,
	taxonomy_table = file.path(tmp, "taxonomy.qza")
	)
tmp_microtable3$tidy_dataset()


##################################################################
##################################################################
# Other examples of data reading with QIIME2 qza-converted files (This is a small dataset, intended solely for showing the steps)
# use microtable$new function to construct a microtable object

# feature abundance
feature_table <- read.table("Input/1.Amplicon/QIIME2_converted/feature_table.tsv", sep = "\t", skip = 1, comment.char = "$", header = TRUE, row.names = 1)

# taxonomy table
taxonomy_table_raw <- read.table("Input/1.Amplicon/QIIME2_converted/taxonomy.tsv", sep = "\t", skip = 1, comment.char = "$", header = FALSE)
# split the merged taxonomy to generate a taxonomic table
# please adjust split parameter when the connector symbol is not ";"
taxonomy_table <- file2meco:::split_assignments(taxonomy_table_raw[, 2], ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), split = ";") %>% as.data.frame
rownames(taxonomy_table) <- taxonomy_table_raw[, 1]

# phylogenetic tree
tree <- ape::read.tree("Input/1.Amplicon/QIIME2_converted/tree.nwk")

# sequences
sequences <- Biostrings::readDNAStringSet("Input/1.Amplicon/QIIME2_converted/dna-sequences.fasta")

# construct microtable object
test <- microtable$new(sample_table = tmp_sample, otu_table = feature_table, tax_table = taxonomy_table, phylo_tree = tree, rep_fasta = sequences)

# tidy all the data inside the object to ensure that samples and features correspond across different data
test$tidy_dataset()

# tidy taxonomy information to unify the naming of taxonomy
test %<>% tidy_taxonomy

##################################################################





