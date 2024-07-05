

######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/3.Metabolome/StageⅩ_Metabolome"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
######################################################

file_path <- "./Input/3.Metabolome/Metabolome.xlsx"

tmp_feature <- as.data.frame(read_excel(file_path, col_names = TRUE), stringsAsFactors = FALSE)
# delete those unknown and unidentified
tmp_feature %<>% .[!grepl("unknown", .$Peak), ] %>% .[!grepl("Analyte", .$Peak), ]
# assign rownames
rownames(tmp_feature) <- tmp_feature$Peak
# delete extra columns
tmp_feature <- tmp_feature[, -c(1:2)]
# assign 0 to NA, i.e. those not detected
tmp_feature[is.na(tmp_feature)] <- 0

# construct a taxonomic table
tmp_tax <- data.frame(class = rownames(tmp_feature))
rownames(tmp_tax) <- tmp_tax[, 1]



# use sample information table of amplicon sequencing data
input_path <- "Output/1.Amplicon/StageⅡ_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in StageⅡ !")
}
load(input_path)


# use sample_table in 16S data
tmp_mt <- microtable$new(otu_table = tmp_feature, tax_table = tmp_tax, sample_table = amplicon_16S_microtable$sample_table)
tmp_mt$tidy_dataset()
# calculate abundance
tmp_mt$cal_abund(rel = FALSE)

# save the object to RData
metab_microtable <- clone(tmp_mt)
save(metab_microtable, file = file.path(output_dir, "metab_microtable.RData"))





