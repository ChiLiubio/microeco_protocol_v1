## 
## Correlation analysis between differential taxa and environmental factors
## 


######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)

# select data for rhizosphere soil
tmp_microtable_rhizo <- clone(amplicon_16S_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()
# get the relative abundance data: taxa_abund list stored in the object
tmp_microtable_rhizo$cal_abund(rel = TRUE)
######################################################



######################################################

formula <- "Cropping+Fertilization"
taxlevel <- "Genus"

# when no random effect is found, "glmm_beta" represents generalized linear model with a family function of beta distribution, and is totally same with "betareg" (beta regression)
tmp_transdiff <- trans_diff$new(dataset = tmp_microtable_rhizo, method = "glmm_beta", formula = formula, taxa_level = taxlevel, filter_thres = 0.002)

# select the feature with significance in rotational cropping and two fertilization treatments
select_taxa <- tmp_transdiff$res_diff %>% .[! is.na(.$Estimate), ] %>% .[.$Factors != "(Intercept)", ] %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% unique
# further select the feature with standard taxonomic names; \\d is a regular expression usage, meaning any number from 0 to 9
select_taxa %<>% .[!grepl("\\d", .)]

# use trans_env class to perform correlation between taxa and environmental variables
tmp_transenv_rhizo <- trans_env$new(dataset = tmp_microtable_rhizo, env_cols = 7:19, standardize = TRUE)
# use_data = "other" means provide customized taxa names via other_taxa parameter
tmp_transenv_rhizo$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv_rhizo$res_cor, file.path(output_dir, "Diff_abund_env_cor_betareg.csv"))

# visualization
g1 <- tmp_transenv_rhizo$plot_cor(cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_cluster_betareg.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)

# When the cluster plot is added, the output object is aplot class, not ggplot class. 
# If you want to change the plot like ggplot2 usage, you should access the object like managing a list object.
g1[[1]]
# optimize the x axis text name
g1[[1]] <- g1[[1]] + ggplot2::scale_x_discrete(labels = c(NH4 = expression(NH[4]^'+'-N), NO3 = expression(NO[3]^'-'-N)))
g1
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_cluster_betareg_newxaxistext.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)


# also visualize the differential test result
# select factors and taxa that will be visualized in the figure
tmp <- tmp_transdiff$res_diff %>% .[.$Taxa %in% select_taxa, ] %>% .[!.$Factors %in% c("(Intercept)"), ]
# assign back to object
tmp_transdiff$res_diff <- tmp
write.csv(tmp_transdiff$res_diff, file.path(output_dir, "Diff_abund_diff_betareg.csv"))

# when the formula is found, the plot_diff_bar function can automatically employ the heatmap instead of bar plot.
g1 <- tmp_transdiff$plot_diff_bar(heatmap_cell = "Estimate", heatmap_lab_fill = "Betareg\nCoef", cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_diff_betareg.png"), g1, base_aspect_ratio = 1, dpi = 300, base_height = 8)















