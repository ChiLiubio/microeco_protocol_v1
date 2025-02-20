##
## 
## Calculate alpha diversity and save the data
## 
## 


###########################
# load packages
library(microeco)
library(magrittr)
library(mecodev)
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage4_AlphaDiversity"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
# first check whether the data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
# load amplicon_16S_microtable for rarefaction curve
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
load(input_path)
###########################

# rarefaction curve based on the Shannon-Weiver diversity
# use trans_rarefy class in mecodev package
tmp <- trans_rarefy$new(amplicon_16S_microtable, alphadiv = c("Shannon"), depth = c(0, 10, 50, 500, 1000, 2000, 4000, 6000, 10000, 15000, 20000, 25000, 30000))
# show_samplename: add the sample labels; For the colors of different groups, please use color parameter, such as: color = "Group"
g1 <- tmp$plot_rarefy(show_samplename = TRUE, color_values = rep("grey50", 100), show_legend = FALSE)

cowplot::save_plot(file.path(output_dir, "AlphaDiv_Rarefactioncurve_Shannon.png"), g1, base_aspect_ratio = 1.4, dpi = 300, base_height = 5)




tmp_microtable <- clone(amplicon_16S_microtable_rarefy)

# calculate alpha diversity table based on the the rarefied abundance table
# Available options of measures argument include "Observed", "Coverage", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher" and "Pielou"
# If the user want to calculate Faith's phylogenetic diversity, please add parameter: PD = TRUE
tmp_microtable$cal_alphadiv()
# Export the diversity table of calculation results to a local file on the computer
write.csv(tmp_microtable$alpha_diversity, file.path(output_dir, "AlphaDiv_metrics.csv"))
# Description of the result file: 
# 	The first column represents the sample names
#	The "Shannon" column denotes the Shannon–Weaver index
#	Simpson and InvSimpson are two variants of Simpson's index. Simpson: 1−D; InvSimpson: 1/D, D=∑pi^2, pi is the proportional abundance of species i.
#	"Coverage" represents good's coverage.
# 	"Pielou" denotes the Pielou evenness index.
# 	"Chao1" and "ACE" reprenset Chao1 and ACE indexes, respectively.


# save the microtable object with the alpha diversity
save(tmp_microtable, file = file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData"), compress = TRUE)





