##
## 
## Calculate alpha diversity and save the results to a comma separated file
## 
## 


###########################
# load packages
library(microeco)
library(magrittr)
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage4_AlphaDiversity"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
###########################


tmp_microtable <- clone(amplicon_16S_microtable_rarefy)

# calculate alpha diversity
tmp_microtable$cal_alphadiv()
write.csv(tmp_microtable$alpha_diversity, file.path(output_dir, "AlphaDiv_metrics.csv"))

save(tmp_microtable, file = file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData"), compress = TRUE)






