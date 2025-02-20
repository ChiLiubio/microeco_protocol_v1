## 
## Preprocess the raw microtable object
## 


###################################################
# load packages
library(microeco)
library(magrittr)
###################################################
# check output directory
output_dir <- "./Output/1.Amplicon/Stage2_amplicon_microtable"
# load raw microtable object RData generated in last step
input_path <- file.path(output_dir, "amplicon_16S_microtable_raw.RData")
if(!file.exists(input_path)){
	stop("Please first run the script in last step !")
}
load(input_path)
###################################################
# clone an object for the convenience of operation
tmp_microtable <- clone(amplicon_16S_microtable_raw)

# filter the ASV not classified in the Kingdom Archaea or Bacteria
# If only bacteria are needed, please use: subset(Kingdom == "k__Bacteria")
tmp_microtable$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

# Optional; remove the lines containing mitochondria or chloroplast regardless of word case in the tax_table
tmp_microtable$filter_pollution(taxa = c("mitochondria", "chloroplast"))

# tidy all the data inside the object
tmp_microtable$tidy_dataset()

# Optional
# use ASV as the prefix to generate new feature names instead of raw character strings
tmp_microtable$rename_taxa(newname_prefix = "ASV_")

# have a look on the range of sequencing depths across samples
tmp_microtable$sample_sums() %>% range

# Optional step
# If the sequencing data was performed in batches, 
# we recommend using the ConQuR package (https://www.nature.com/articles/s41467-022-33071-9) to correct the tmp_microtable$otu_table to remove the batch effects.
# The data here did not undergo correction because there was no batch effect.
# Pseudo code like (ConQuR repository: https://github.com/wdl2459/ConQuR):
# library(ConQuR)
# taxa <- tmp_microtable$otu_table %>% t %>% as.data.frame
# batchid <- tmp_microtable$sample_table[, 'batchid']
# taxa_corrected <- ConQuR(tax_tab = taxa, batchid = batchid, batch_ref="0") # see https://wdl2459.github.io/ConQuR/ConQuR.Vignette.html for more details
# tmp_microtable$otu_table <- taxa_corrected %>% t %>% as.data.frame


# Optional
# assign factor levels to make the elements ordered in the following statistics and visualization
tmp_microtable$sample_table$Compartment %<>% factor(., levels = c("Bulk soil", "Rhizosphere", "Endophyte"))
tmp_microtable$sample_table$Fertilization %<>% factor(., levels = c("CK", "NPK", "NPKS"))
tmp_microtable$sample_table$Cropping %<>% factor(., levels = c("CC", "RC"))
tmp_microtable$sample_table$Group %<>% factor(., levels = c("CC-CK", "RC-CK", "CC-NPK", "RC-NPK", "CC-NPKS", "RC-NPKS"))

# use clone function to copy the object completely
amplicon_16S_microtable <- clone(tmp_microtable)

# save the microtable object to the output folder
save(amplicon_16S_microtable, file = file.path(output_dir, "amplicon_16S_microtable.RData"), compress = TRUE)







