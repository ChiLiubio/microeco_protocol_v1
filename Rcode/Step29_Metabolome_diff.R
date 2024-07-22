## 
## Differential test of metabolome data across groups
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)
######################################################
output_dir <- "Output/3.Metabolome/Stage10_Metabolome"
# load data
input_path <- file.path(output_dir, "metab_microtable.RData")
if(! file.exists(input_path)){
	stop("Please first run the script in the step27!")
}
load(input_path)

# load the features selected in the last step
input_path <- file.path(output_dir, "Classification_rf_select_features.RData")
if(! file.exists(input_path)){
	stop("Please first run the script in the step28!")
}
load(input_path)

######################################################

tmp_microtable <- clone(metab_microtable)
# filter features
tmp_microtable$otu_table %<>% .[rownames(.) %in% tmp_sel_features, ]
tmp_microtable$tidy_dataset()
# regenerate the taxa_abund list for the differential test
tmp_microtable$cal_abund(rel = FALSE)



# Linear regression using log-transformed data
t1 <- trans_diff$new(dataset = tmp_microtable, method = "lm", formula = "Compartment + Cropping + Fertilization", taxa_level = "class", transformation = "log")

write.csv(t1$res_diff, file.path(output_dir, "Metabolome_diff_lm_log.csv"))

tmp_table <- t1$res_diff
tmp_table$Significance %<>% gsub("ns", "", .)
tmp_table %<>% .[.$Factors != "(Intercept)", ]
tmp_table$Factors %<>% gsub("CompartmentRhizosphere", "Compartment: Rhizosphere", .)
tmp_table$Factors %<>% gsub("FertilizationNPKS", "Fertilization: NPKS", .)
tmp_table$Factors %<>% gsub("FertilizationNPK", "Fertilization: NPK", .)
tmp_table$Factors %<>% gsub("CroppingRC", "Cropping: RC", .)


t1$res_diff <- tmp_table
g1 <- t1$plot_diff_bar(filter_feature = c("ns", "", "*"), heatmap_cell = "Estimate", heatmap_lab_fill = "Estimate", cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Metabolome_diff_lm_log_filter_nonextremesig.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 8)

g1 <- t1$plot_diff_bar(heatmap_cell = "Estimate", heatmap_lab_fill = "Estimate", cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Metabolome_diff_lm_log.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 8)


###########################################################
# One-way ANOVA for rhizosphere soil using log-transformed data

tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = "anova", group = "Group", taxa_level = "class", transformation = "log")

write.csv(t1$res_diff, file.path(output_dir, "Metabolome_diff_ANOVA_log_Group.csv"))




