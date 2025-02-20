## 
## Correlation between alpha diversity and environmental variables
## 


###########################
# load packages
library(microeco)
library(magrittr)
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage4_AlphaDiversity"

# load data
input_path <- file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData")
if(! file.exists(input_path)){
	stop("Please first run the script in step7 !")
}
load(input_path)
###########################



# select samples from bulk soil
tmp_microtable_bulksoil <- clone(tmp_microtable)
tmp_microtable_bulksoil$sample_table %<>% .[.$Compartment == "Bulk soil", ]
tmp_microtable_bulksoil$tidy_dataset()

# create a trans_env object for the following correlation analysis
# env_cols parameter to select the columns in sample_table of microtable object
tmp <- trans_env$new(dataset = tmp_microtable_bulksoil, env_cols = 9:21)

# calculate the correlation; default method: Pearson correlation. Please use cor_method to select others, e.g., Spearman correlation (cor_method = "spearman")
tmp$cal_cor(add_abund_table = tmp_microtable_bulksoil$alpha_diversity)
# The 'AdjPvalue' column represents the p value after adjustment. 'Significance': *: P < 0.05; **: P < 0.01; ***: P < 0.001.
write.csv(tmp$res_cor, file.path(output_dir, "AlphaDiv_Env_bulksoil_correlation.csv"))

# heatmap on all diversity indexes
g1 <- tmp$plot_cor()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Env_bulksoil_correlation.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)








