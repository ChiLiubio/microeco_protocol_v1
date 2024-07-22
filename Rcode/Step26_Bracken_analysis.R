## 
## Statistical analysis and visualization on imported abundances
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
######################################################
output_dir <- "Output/2.Metagenome/Stage9_Bracken"
# load data
input_path <- file.path(output_dir, "Bracken_microtable_relFALSE.RData")
if(!file.exists(input_path)){
	stop("Please first run the script in last step !")
}
load(input_path)
input_path <- file.path(output_dir, "Bracken_microtable_relTRUE_K1.RData")
load(input_path)
input_path <- file.path(output_dir, "Bracken_microtable_relTRUE_K2.RData")
load(input_path)
input_path <- file.path(output_dir, "Bracken_microtable_relTRUE_K3.RData")
load(input_path)
######################################################

# original abundance visualization
tmp_microtable <- clone(Bracken_microtable_relFALSE)

# use_percentage = FALSE: show the original abundance
tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom", use_percentage = FALSE)
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE, bar_full = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum", use_percentage = FALSE)
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE, bar_full = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus", use_percentage = FALSE)
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE, bar_full = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_rawabund.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)

# show the mean of original abundance for each group
tmp_microtable_group <- tmp_microtable$merge_samples("Group")

tmp <- microtable$new(tmp_microtable$taxa_abund$Genus, sample_table = tmp_microtable$sample_table)
tmp_2 <- tmp$merge_samples("Group")
tmp_microtable_group$taxa_abund$Genus <- tmp_2$otu_table

tmp <- trans_abund$new(dataset = tmp_microtable_group, taxrank = "Genus", use_percentage = FALSE)
g1 <- tmp$plot_bar(others_color = "grey70", bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + ylab("Raw abundance")
cowplot::save_plot(file.path(output_dir, "Bracken_barplot_rawabund_groupmean.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


# relative abundance visualization
tmp_microtable <- clone(Bracken_microtable_relTRUE_K1)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom")
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum")
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus")
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)




# recalculate relative abundance based on the species abundance (otu_table) for Bracken_microtable_relTRUE_K1
tmp_microtable$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom")
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum")
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus")
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K1.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)


# recalculate relative abundance based on the species abundance (otu_table) for Bracken_microtable_relTRUE_K2
tmp_microtable <- clone(Bracken_microtable_relTRUE_K2)

tmp_microtable$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom")
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum")
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus")
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K2.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)


# recalculate relative abundance based on the species abundance (otu_table) for Bracken_microtable_relTRUE_K3
tmp_microtable <- clone(Bracken_microtable_relTRUE_K3)

tmp_microtable$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom")
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum")
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus")
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K3.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)



####################################################################

# PCA at species level; all samples
tmp_microtable <- clone(Bracken_microtable_relTRUE_K1)

t1 <- trans_beta$new(dataset = tmp_microtable)

t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
rownames(t1$res_ordination$loading) %<>% gsub(".*__", "", .)
# write the results to local file
write.csv(t1$res_ordination$scores, file.path(output_dir, "Bracken_PCA_Species_Score.csv"))
write.csv(t1$res_ordination$loading, file.path(output_dir, "Bracken_PCA_Species_Loading.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Compartment", loading_arrow = TRUE, loading_text_italic = FALSE)
cowplot::save_plot(file.path(output_dir, "Bracken_PCA_Species.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


# Rhizosphere
tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

t1 <- trans_beta$new(dataset = tmp_microtable_rhizo)

t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
rownames(t1$res_ordination$loading) %<>% gsub(".*__", "", .)
write.csv(t1$res_ordination$scores, file.path(output_dir, "Bracken_PCA_Species_Rhizosphere_Score.csv"))
write.csv(t1$res_ordination$loading, file.path(output_dir, "Bracken_PCA_Species_Rhizosphere_Loading.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", loading_arrow = TRUE, loading_text_italic = FALSE)
g1 <- g1 + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, color = "grey50", linetype = 2) + geom_hline(yintercept = 0, color = "grey50", linetype = 2)
cowplot::save_plot(file.path(output_dir, "Bracken_PCA_Species_Rhizosphere.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)



####################################################################
# differential test
tmp_microtable <- clone(Bracken_microtable_relTRUE_K1)
tmp_microtable$cal_abund(rel = TRUE)

tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

# perform lefse
tmp <- trans_diff$new(tmp_microtable_rhizo, method = "lefse", group = "Cropping", p_adjust_method = "none")
write.csv(tmp$res_diff, file.path(output_dir, "Bracken_lefse_Rhizosphere_Cropping.csv"))
g1 <- tmp$plot_diff_bar(threshold = 2, width = 0.7)
cowplot::save_plot(file.path(output_dir, "Bracken_lefse_Rhizosphere_Cropping.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)


# multiple factors test; ANOVA; non-relative abundance with log trans
tmp <- trans_diff$new(dataset = tmp_microtable, method = "anova", formula = "Cropping*Fertilization*Compartment", taxa_level = "Phylum", transformation = "AST")
# By default, the x and y axis text order come from the clustering.
# To adjust the text order, please use parameters text_x_order or text_y_order
g1 <- tmp$plot_diff_bar(text_x_order = tmp$res_diff$Factors[1:7])
cowplot::save_plot(file.path(output_dir, "Bracken_anova_multiway_AST.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 8)


# multiple factors test; beta regression
# select genera with high relative abundance
t1 <- trans_diff$new(dataset = tmp_microtable, method = "glmm_beta", formula = "Cropping + Fertilization + Compartment", taxa_level = "Genus", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, "Bracken_betareg_Genus.csv"))

t1$res_diff %<>% .[.$Factors != "(Intercept)", ]
t1$res_diff %<>% .[!is.na(.$Estimate), ]
# further filter features
tmp_table <- t1$res_diff
tmp_feature_1 <- tmp_table %>% .[.$Factors == "CroppingRC" & grepl("**", .$Significance, fixed = TRUE), ] %>% .$Taxa
tmp_feature_2 <- tmp_table %>% .[.$Factors == "FertilizationNPK" & grepl("**", .$Significance, fixed = TRUE), ] %>% .$Taxa
t1$res_diff %<>% .[.$Taxa %in% c(tmp_feature_1, tmp_feature_2), ]

g1 <- t1$plot_diff_bar(filter_feature = "", heatmap_cell = "Estimate", heatmap_lab_fill = "Estimate")
cowplot::save_plot(file.path(output_dir, "Bracken_betareg_Genus_extremsig.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 10)












