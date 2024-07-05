


######################################################
# load packages
library(microeco)
library(file2meco)
library(magrittr)
library(readxl)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/2.Metagenome/StageⅨ_Bracken"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
######################################################

file_path <- "./Input/2.Metagenome/Bracken/bracken_abundance_table.txt"
match_table_path <- "./Input/2.Metagenome/match_table.xlsx"

# use sample information table of amplicon sequencing data
input_path <- "Output/1.Amplicon/StageⅡ_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in StageⅡ !")
}
load(input_path)

# generate microtable object
# rel = FALSE: original abundances in taxa_abund list
tmp_microtable <- mpa2meco(file_path, match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table, rel = FALSE, use_level = "s__")
tmp_microtable$tidy_dataset()

Bracken_microtable_relFALSE <- clone(tmp_microtable)
save(Bracken_microtable_relFALSE, file = file.path(output_dir, "Bracken_microtable_relFALSE.RData"))


# rel = TRUE: relative abundances in taxa_abund list
# sel_same = 1: select the first Kingdom name when multiple names are found (e.g., k__Eukaryota|k__Fungi). 1 is the default selection.
tmp_microtable <- mpa2meco(file_path, match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table, rel = TRUE, use_level = "s__", sel_same = 1)
tmp_microtable$tidy_dataset()

Bracken_microtable_relTRUE_K1 <- clone(tmp_microtable)
save(Bracken_microtable_relTRUE_K1, file = file.path(output_dir, "Bracken_microtable_relTRUE_K1.RData"))

# sel_same = 2: select the second Kingdom name when multiple names are found (e.g., k__Eukaryota|k__Fungi)
tmp_microtable <- mpa2meco(file_path, match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table, rel = TRUE, use_level = "s__", sel_same = 2)
tmp_microtable$tidy_dataset()

Bracken_microtable_relTRUE_K2 <- clone(tmp_microtable)
save(Bracken_microtable_relTRUE_K2, file = file.path(output_dir, "Bracken_microtable_relTRUE_K2.RData"))

# sel_same = 3: select both Kingdom names when multiple names are found (e.g., k__Eukaryota|k__Fungi)
tmp_microtable <- mpa2meco(file_path, match_table = match_table_path, sample_table = amplicon_16S_microtable$sample_table, rel = TRUE, use_level = "s__", sel_same = 3)
tmp_microtable$tidy_dataset()

Bracken_microtable_relTRUE_K3 <- clone(tmp_microtable)
save(Bracken_microtable_relTRUE_K3, file = file.path(output_dir, "Bracken_microtable_relTRUE_K3.RData"))






