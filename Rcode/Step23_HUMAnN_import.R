## 
## Import HUMAnN analysis results
## 


######################################################
# load packages
library(microeco)
library(file2meco)
library(magrittr)
library(readxl)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/2.Metagenome/Stage8_HUMAnN"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}

# use sample information table in amplicon_16S_microtable object (stage2)
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################

# File 'match_table.xlsx' is used to replace the sample names in abundance file to make them same with those in provided sample metadata.
# If the samples names in provided sample metadata is totally same with those in feature abundance table, this file is not needed.
match_table_path <- "./Input/2.Metagenome/match_table.xlsx"


# MetaCyc
file_path <- "./Input/2.Metagenome/HUMAnN3/Metacyc_pathabundance_joint_table.tsv"

# generate microtable object
tmp_microtable <- humann2meco(file_path, db = "MetaCyc", match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table)
# remove the pathways classified into "unclassified" class
tmp_microtable$tax_table %<>% subset(Superclass1 != "unclassified")
# prune all the data
tmp_microtable$tidy_dataset()

MetaCyc_microtable <- clone(tmp_microtable)
save(MetaCyc_microtable, file = file.path(output_dir, "MetaCyc_microtable.RData"))



# KEGG abundance file
file_path <- "./Input/2.Metagenome/HUMAnN3/KEGG_pathabundance_joint_table.tsv"

# generate microtable object
tmp_microtable <- humann2meco(file_path, db = "KEGG", match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table)
# remove the pathways classified into "unclassified" class
tmp_microtable$tax_table %<>% subset(Level.1 != "unclassified")
# prune all the data
tmp_microtable$tidy_dataset()

KEGG_microtable <- clone(tmp_microtable)
save(KEGG_microtable, file = file.path(output_dir, "KEGG_microtable.RData"))














