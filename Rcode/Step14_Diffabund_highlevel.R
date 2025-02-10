## 
## Differential abundance test at high taxonomic level
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


tmp_microtable <- clone(amplicon_16S_microtable)
tmp_microtable$cal_abund(rel = TRUE)

# select data for rhizosphere soil
tmp_microtable_rhizo <- clone(amplicon_16S_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()
# get the relative abundance data: taxa_abund list stored in the object
tmp_microtable_rhizo$cal_abund(rel = TRUE)
######################################################



##################################################################
##  Genus level     relative abundance      single factor    ####
##################################################################

group <- "Fertilization"
taxlevel <- "Genus"

# LEfSe <doi:10.1186/gb-2011-12-6-r60>
method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, filter_thres = 0.001)
# In the result, LDA means the score of linear discriminant analysis
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_singlefactor_lefse_LDAbar.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_singlefactor_lefse_abundsig.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


# Wilcoxon test with arcsine transformation
method <- "wilcox"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_singlefactor_wilcox_abundsig.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


# one-way ANOVA with arcsine transformation
method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
# In the result, 'Letter' represents the post hoc test for one-way anova. To change the default 'duncan.test' method, please add anova_post_test parameter, which will be passed to trans_alpha class.
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_singlefactor_ANOVA_abundsig.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


# Dunn's Kruskal-Wallis Multiple Comparisons with arcsine transformation
method <- "KW_dunn"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_singlefactor_KWdunn_abundsig.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


## Beta regression
## We did not run beta regression because the glmm_beta method that follows is exactly the same as beta regression when there are no random effects. Moreover, the glmm_beta method is more robust
# method <- "betareg"
# tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = group, taxa_level = taxlevel, filter_thres = 0.001)
## "Estimate" and "Std.Error" represent the fitted coefficient and its standard error, respectively. Zvalue is the statistic.
# write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))

# generalized linear mixed-effects model with the beta distribution family
# when no random effect is provided, "glmm_beta" represents generalized linear model with a family function of beta distribution and is totally same with "betareg" (beta regression)
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = group, taxa_level = taxlevel, filter_thres = 0.001)
# R2 represents the variance that the factors account for in the model. "Estimate" and "Std.Error" represent the fitted coefficient and its standard error, respectively.
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))






##################################################################
##  Genus level     relative abundance      multiple factor   ####
##################################################################

formula <- "Compartment+Cropping+Fertilization"
taxlevel <- "Genus"

method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))

# linear regression with arcsine transformation
method <- "lm"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))


# generalized linear mixed-effects model with the beta distribution family
# when no random effect is provided, "glmm_beta" represents generalized linear model with a family function of beta distribution and is totally same with "betareg" (beta regression)
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))

# Figure 4c

tmp$res_diff %<>% .[.$Factors != "(Intercept)", ]
tmp$res_diff %<>% .[.$Factors != "Model", ]
tmp$res_diff %<>% .[.$Factors != "Compartment", ]
tmp$res_diff %<>% .[.$Factors != "Cropping", ]
tmp$res_diff %<>% .[.$Factors != "Fertilization", ]

# filter those without extreme significance in treatments
tmp_sel <- c()
for(x in unique(tmp$res_diff$Taxa)){
	tmp_table <- tmp$res_diff[tmp$res_diff$Taxa == x, ]
	tmp_table <- tmp_table[tmp_table$Factors %in% c("CroppingRC", "FertilizationNPK", "FertilizationNPKS"), ]
	if(any(grepl("**", tmp_table$Significance, fixed = TRUE))){
		tmp_sel <- c(tmp_sel, x)
	}else{
		next
	}
}

tmp$res_diff %<>% .[tmp$res_diff$Taxa %in% tmp_sel, ]

# adjust the group names
tmp$res_diff$Factors %<>% gsub("CompartmentRhizosphere", "Compartment: Rhizosphere", .)
tmp$res_diff$Factors %<>% gsub("CompartmentEndophyte", "Compartment: Endophyte", .)
tmp$res_diff$Factors %<>% gsub("CroppingRC", "Cropping: RC", .)
tmp$res_diff$Factors %<>% gsub("FertilizationNPK", "Fertilization: NPK", .)
tmp$res_diff$Factors %<>% gsub("FertilizationNPKS", "Fertilization: NPKS", .)

# remove the genera without standard names
tmp_sel <- tmp$res_diff$Taxa %>% grepl("\\d+", .)
tmp$res_diff %<>% .[! tmp_sel, ]
# delete the prefix
tmp$res_diff$Taxa %<>% gsub(".*g__", "", .)

tmp$res_diff$Taxa[tmp$res_diff$Taxa == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
tmp$res_diff$Taxa[tmp$res_diff$Taxa == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"

g2 <- tmp$plot_diff_bar(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef", xtext_angle = 30, xtext_size = 16, filter_feature = "",
	text_x_order = c("Compartment: Rhizosphere", "Compartment: Endophyte", "Cropping: RC", "Fertilization: NPK", "Fertilization: NPKS")) + 
	theme(legend.position = "left")

cowplot::save_plot(file.path(output_dir, "Diff_abund_test_Genus_multifactor_glmmbeta.png"), g2, base_aspect_ratio = 1.3, dpi = 300, base_height = 9)





##############################################################################
##  Multiple taxonomic level     relative abundance      single factor    ####
##############################################################################

group <- "Fertilization"
taxlevel <- "all"

method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, alpha = 0.01, lefse_subgroup = NULL, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_multitax_singlefactor_lefse_LDAbar.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 5)
cowplot::save_plot(file.path(output_dir, "Diff_abund_test_multitax_singlefactor_lefse_cladogram.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 9)


method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))

method <- "KW_dunn"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))



##############################################################################
##  Multiple taxonomic level     relative abundance      multiple factor  ####
##############################################################################

formula <- "Cropping+Fertilization"
taxlevel <- "all"

# ANOVA with arcsine transformation
method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))

# linear regression with arcsine transformation
method <- "lm"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))

# generalized linear mixed-effects model with the beta distribution family
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.005)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))


##################################################################################################
##  Random effects   Multiple or single high taxonomic level     relative abundance      ##############
##################################################################################################

formula <- "Cropping + Fertilization + (1|Compartment)"
taxlevel <- "Genus"

# linear mixed-effects model with arcsine transformation
method <- "lme"
tmp <- trans_diff$new(dataset = amplicon_16S_microtable, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.005)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_mixedeff_Genus_", method, ".csv")))

# generalized linear mixed-effects model with the beta distribution family
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = amplicon_16S_microtable, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.005)
# In R2, Conditional R2 means the R2 explained by both the fixed and random effects. Marginal R2 represents the R2 explained by the fixed effects. 
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_mixedeff_Genus_", method, ".csv")))





