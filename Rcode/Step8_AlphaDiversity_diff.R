## 
## Differential test for alpha diversity
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
	stop("Please first run the script of last step !")
}
load(input_path)
###########################



# Compare two cropping treatments for each compartment
measure <- "Shannon"

# If the by_group parameter is not needed, please delete it
tmp <- trans_alpha$new(dataset = tmp_microtable, group = "Cropping", by_group = "Compartment")
tmp$cal_diff(method = "wilcox", measure = measure)
# In the result, the P.adj column is the p value after adjustment with FDR method for method = "wilcox". To change the adjustment method, please use p_adjust_method parameter.
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_bycompart_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_bycompart_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_alpha(measure = measure, plot_type = "errorbar", add_line = TRUE, line_type = 2)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_bycompart_wilcox_errorbar_line.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


# Compare three fertilization treatments for each compartment
measure <- "Shannon"

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "Fertilization", by_group = "Compartment")

tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure, y_increase = 0.05)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

tmp$cal_diff(method = "anova", measure = measure)
# The letters in the 'Letter' column are the result of post hoc test for one-way anova. To change the default 'duncan.test' method, please use anova_post_test parameter.
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_anova.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_anova_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

tmp$cal_diff(method = "KW_dunn", measure = measure, KW_dunn_letter = TRUE)
# The file format is same with one-way ANOVA. The 'MonoLetter' column is for aligning letters to facilitate the viewing of results.
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_KW_dunn.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_KW_dunn_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)




# two-way anova: Cropping + Fertilization + Cropping:Fertilization
tmp <- trans_alpha$new(dataset = tmp_microtable, by_group = "Compartment")

# for single measure
measure <- "Shannon"
tmp$cal_diff(method = "anova", measure = measure, formula = "Cropping*Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova.csv"))

## for all measures
tmp <- trans_alpha$new(dataset = tmp_microtable, by_group = "Compartment")
tmp$cal_diff(method = "anova", formula = "Cropping*Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova_allmeasures.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova_allmeasures_heatmap.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 7)



# paired t-test or paired Wilcoxon test
measure <- "Shannon"
# parameter by_ID is used to provide the identities of paired data
tmp <- trans_alpha$new(dataset = tmp_microtable, group = "Cropping", by_group = "Compartment", by_ID = "Plant_ID")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_bycompart_wilcox_paireddata.csv"))
tmp$cal_diff(method = "t.test", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_bycompart_ttest_paireddata.csv"))



# Linear model
tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "lm", formula = "Compartment+Cropping+Fertilization")

tmp$res_diff %<>% .[.$Factors != "(Intercept)", ]
tmp$res_diff$Factors %<>% gsub("CompartmentRhizosphere", "Compartment: Rhizosphere", .)
tmp$res_diff$Factors %<>% gsub("CompartmentEndophyte", "Compartment: Endophyte", .)
tmp$res_diff$Factors %<>% gsub("CroppingRC", "Cropping: RC", .)
tmp$res_diff$Factors %<>% gsub("FertilizationNPK", "Fertilization: NPK", .)
tmp$res_diff$Factors %<>% gsub("FertilizationNPKS", "Fertilization: NPKS", .)

# Figure 2e
g1 <- tmp$plot_alpha(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef", xtext_angle = 50, 
	text_x_order = c("Compartment: Rhizosphere", "Compartment: Endophyte", "Cropping: RC", "Fertilization: NPK", "Fertilization: NPKS")) + 
	theme(legend.position = "left")
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_Compartment_lm_heatmap.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)


# use errorbar plot to show one index
tmp$res_diff %<>% .[.$Measure == "Shannon", ]

# Figure 2f
g1 <- tmp$plot_alpha(coefplot_sig_pos = 3)
g1 <- g1 + theme(axis.text.y = element_text(size = 13, color = "black"), axis.title.x = element_text(size = 13, color = "black"))

cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_Compartment_lm_Shannon_errorbar.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)


# linear mixed-effects model
# for the formula usage, please see https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html
tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "lme", formula = "Cropping+Fertilization+(1|Compartment)")
# "Estimate" and "Std.Error" columns represent the fitted coefficient and its standard error, respectively.
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lme.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lme.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)



##########################################################################
##########################################################################
# select one compartment: Endophyte

tmp_microtable_endophyte <- clone(tmp_microtable)
tmp_microtable_endophyte$sample_table %<>% .[.$Compartment == "Endophyte", ]
tmp_microtable_endophyte$tidy_dataset()

# linear regression
tmp <- trans_alpha$new(dataset = tmp_microtable_endophyte)
tmp$cal_diff(method = "lm", formula = "Cropping+Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lm_Endophyte.csv"))
g1 <- tmp$plot_alpha(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef")
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lm_Endophyte.png"), g1, base_aspect_ratio = 1, dpi = 300, base_height = 7)

# one-way anova
measure <- "Shannon"
tmp <- trans_alpha$new(dataset = tmp_microtable_endophyte, group = "Group")
tmp$cal_diff(method = "anova")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Group_Endophyte_anova.csv"))
g1 <- tmp$plot_alpha(measure = measure, add_sig_text_size = 5)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Group_Endophyte_anova.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 5)

# Kruskal-Wallis test with Dunn's multiple comparisons
tmp <- trans_alpha$new(dataset = tmp_microtable_endophyte, group = "Fertilization")
tmp$cal_diff(method = "KW_dunn", KW_dunn_letter = TRUE)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_KW_dunn_Endophyte.csv"))
# Figure 2c
g1 <- tmp$plot_alpha(plot_type = "ggdotplot", fill = "Fertilization", alpha = 0.3, size = 2, y_increase = 0.4, add = "mean_se", xtext_angle = 0, add_sig_text_size = 5)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_KW_dunn_Endophyte.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# Wilcoxon test for fertilization treatments in each cropping treatment group
tmp <- trans_alpha$new(dataset = tmp_microtable_endophyte, group = "Fertilization", by_group = "Cropping")
tmp$cal_diff(method = "wilcox")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_byCropping_wilcox_Endophyte.csv"))
# Figure 2d
g1 <- tmp$plot_alpha(plot_type = "ggviolin", fill = "Fertilization", alpha = 0.3, y_start = 0.3, y_increase = 0.15, add = "mean_se", xtext_angle = 0)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_byCropping_wilcox_Endophyte.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)



