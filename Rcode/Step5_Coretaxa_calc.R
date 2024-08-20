## 
## Perform core taxa analysis for each compartment
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage3_coremicrobiome"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load rarefied data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
# first check whether data path exists
if(! file.exists(input_path)){
	stop("Please first run all the scripts in Stage2 !")
}
load(input_path)
######################################################


# use temporary name tmp_microtable for convenience
tmp_microtable <- clone(amplicon_16S_microtable_rarefy)


# define core taxa: occurrence frequency 50%; mean relative abundance 0.05%
freq <- 0.5
abund <- 0.0005

# select core taxa for each compartment
S <- clone(tmp_microtable)
S$sample_table %<>% .[.$Compartment != "Bulk soil", ]
S$tidy_dataset()
S$filter_taxa(rel_abund = abund, freq = freq)

RS <- clone(tmp_microtable)
RS$sample_table %<>% .[.$Compartment != "Rhizosphere", ]
RS$tidy_dataset()
RS$filter_taxa(rel_abund = abund, freq = freq)

R <- clone(tmp_microtable)
R$sample_table %<>% .[.$Compartment != "Endophyte", ]
R$tidy_dataset()
R$filter_taxa(rel_abund = abund, freq = freq)

# merge the core ASV names into a data.frame object and save it into a file
res <- rbind(data.frame(compartment = "S", ASV = S$taxa_names()), data.frame(compartment = "RS", ASV = RS$taxa_names()), data.frame(compartment = "R", ASV = R$taxa_names()))
# write the data.frame object containing core ASV names in each compartment to the output directory
write.csv(res, file.path(output_dir, "Coretaxa_calc_compartments.csv"))

# also save each microtable object to directory
save(S, file = file.path(output_dir, "Coretaxa_Bulk.RData"), compress = TRUE)
save(RS, file = file.path(output_dir, "Coretaxa_Rhizosphere.RData"), compress = TRUE)
save(R, file = file.path(output_dir, "Coretaxa_Endophyte.RData"), compress = TRUE)







