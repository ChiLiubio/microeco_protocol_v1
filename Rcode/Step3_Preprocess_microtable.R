

###################################################
# load packages
library(microeco)
library(magrittr)
###################################################
# check output directory
output_dir <- "./Output/1.Amplicon/StageⅡ_amplicon_microtable"
if(! dir.exists(output_dir)){
	stop("Please first run Step2_Import_amplicon_files.R!")
}
###################################################

# load raw microtable object generated in Step2
input_path <- file.path(output_dir, "amplicon_16S_microtable_raw.RData")
load(input_path)

# clone an object for the convenience of operation
tmp_microtable <- clone(amplicon_16S_microtable_raw)

# filter the ASV not classified in the Kingdom Archaea or Bacteria
tmp_microtable$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

# remove the lines containing mitochondria or chloroplast regardless of word case in the tax_table
tmp_microtable$filter_pollution(taxa = c("mitochondria", "chloroplast"))

# prune all the data inside the object
tmp_microtable$tidy_dataset()

# we use ASV as the prefix to generate new feature names instead of raw character strings
tmp_microtable$rename_taxa(newname_prefix = "ASV_")

# have a look on the range of sample sequencing depths
tmp_microtable$sample_sums() %>% range

# assign factor levels to make the elements ordered in the following statistics and visualization
tmp_microtable$sample_table$Compartment %<>% factor(., levels = c("Bulk soil", "Rhizosphere", "Endophyte"))
tmp_microtable$sample_table$Fertilization %<>% factor(., levels = c("CK", "NPK", "NPKS"))
tmp_microtable$sample_table$Cropping %<>% factor(., levels = c("CC", "RC"))
tmp_microtable$sample_table$Group %<>% factor(., levels = c("CC-CK", "RC-CK", "CC-NPK", "RC-NPK", "CC-NPKS", "RC-NPKS"))

# use clone function to copy the object completely
amplicon_16S_microtable <- clone(tmp_microtable)

# use save function to save the microtable object to output directory
save(amplicon_16S_microtable, file = file.path(output_dir, "amplicon_16S_microtable.RData"), compress = TRUE)







