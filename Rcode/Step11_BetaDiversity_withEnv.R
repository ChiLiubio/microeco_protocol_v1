## 
## Beta diversity analysis involving environmental factors
## 


###########################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage5_BetaDiversity"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
###########################

# select rarefied data for rhizosphere soil
tmp_microtable_rarefy_rhizo <- clone(amplicon_16S_microtable_rarefy)
tmp_microtable_rarefy_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rarefy_rhizo$tidy_dataset()

measure <- "bray"
tmp_microtable_rarefy_rhizo$cal_betadiv(method = measure)

#########################################################################################
## Beta diversity methods with environmental variables
#########################################################################################

# dbRDA: distance-based RDA
method <- "dbRDA"
# standardize = TRUE: standardize variables
t1 <- trans_env$new(dataset = tmp_microtable_rarefy_rhizo, env_cols = 7:19, standardize = TRUE)
t1$cal_ordination(method = method, use_measure = "bray")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
# save the R2 of model to the directory
capture.output(t1$res_ordination_R2, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_rawoutput.txt"))
# save fitted R2 and significance of each factor to the directory
capture.output(t1$res_ordination_envfit, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_rawoutput.txt"), append = TRUE)
# save the significances of factors in the model to the directory
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_termssig.csv"))
# save the significances of axes in the model to the directory
write.csv(t1$res_ordination_axis, file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_axissig.csv"))

# transform raw results for visualization
t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
# save the transformed data of scores and arrows to the directory
write.csv(t1$res_ordination_trans$df_sites, file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows, file.path(output_dir, "BetaDiv_rarefy_rhizo_dbRDA_trans_arrow.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_dbRDA.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)



# RDA at Genus level
method <- "RDA"
t1 <- trans_env$new(dataset = tmp_microtable_rarefy_rhizo, env_cols = 7:19, standardize = TRUE)
t1$cal_ordination(method = method, taxa_level = "Genus")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
# save statistical results to local file
capture.output(t1$res_ordination_R2, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_rawoutput.txt"))
capture.output(t1$res_ordination_envfit, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_rawoutput.txt"), append = TRUE)
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_termssig.csv"))
write.csv(t1$res_ordination_axis, file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_axissig.csv"))

# transform raw results for visualization
t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1, min_perc_tax = 0.2, max_perc_tax = 1)
write.csv(t1$res_ordination_trans$df_sites, file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows, file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_trans_arrowenv.csv"))
write.csv(t1$res_ordination_trans$df_arrows_spe, file.path(output_dir, "BetaDiv_rarefy_rhizo_RDA_Genus_trans_arrowspe.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_RDA_Genus.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)




# CCA at Genus level; same steps with RDA at Genus level
method <- "CCA"
t1 <- trans_env$new(dataset = tmp_microtable_rarefy_rhizo, env_cols = 7:19, standardize = TRUE)
t1$cal_ordination(method = method, taxa_level = "Genus")
# get the significance of the terms
t1$cal_ordination_anova()
# fit factors onto the ordination to get R2 for each factor
t1$cal_ordination_envfit()
# save statistical results to local file
capture.output(t1$res_ordination_R2, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_rawoutput.txt"))
capture.output(t1$res_ordination_envfit, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_rawoutput.txt"), append = TRUE)
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_termssig.csv"))
write.csv(t1$res_ordination_axis, file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_axissig.csv"))

t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1, min_perc_tax = 0.2, max_perc_tax = 1)
write.csv(t1$res_ordination_trans$df_sites, file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows, file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_trans_arrowenv.csv"))
write.csv(t1$res_ordination_trans$df_arrows_spe, file.path(output_dir, "BetaDiv_rarefy_rhizo_CCA_Genus_trans_arrowspe.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group")
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_CCA_Genus.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)




##################################################
# Mantel test

t1 <- trans_env$new(dataset = tmp_microtable_rarefy_rhizo, env_cols = 7:19, standardize = TRUE)
t1$cal_mantel(use_measure = "bray")
# save the results to the directory; 'p.adjusted' column: p value adjustment. The meaning of 'Significance' is "*: P < 0.05, **: P < 0.01; ***: P < 0.001".
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_rhizo_bray.csv"))
# partial mantel test
t1$cal_mantel(partial_mantel = TRUE, use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_rhizo_bray_partial.csv"))
# mantel test for different groups
t1$cal_mantel(by_group = "Cropping", use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_rhizo_bray_byCropping.csv"))
# partial mantel test for different groups
t1$cal_mantel(by_group = "Cropping", partial_mantel = TRUE, use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_rhizo_bray_byCropping_partial.csv"))



#########################################
# correlation
t1 <- trans_env$new(dataset = tmp_microtable_rarefy_rhizo, env_cols = 7:19, standardize = TRUE)

g1 <- t1$plot_scatterfit(
    x = tmp_microtable_rarefy_rhizo$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)], y = "pH",
	type = "cor", group = "Cropping", group_order = c("CC", "RC"),
	point_size = 3, point_alpha = 0.6, line_se = TRUE, line_size = 1.5, shape_values = c(16, 17),
    y_axis_title = "Euclidean distance of pH", x_axis_title = "Bray-Curtis distance", size = 5
)

cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_scatterfit_pH_Bray.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)






