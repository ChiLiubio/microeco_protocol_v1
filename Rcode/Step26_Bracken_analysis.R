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

# show the mean of original abundance (in the taxa_abund list of microtable object) for each group
tmp_microtable_group <- tmp_microtable$merge_samples("Group")
# extract the abundance and generate a new microtable object instead of runing cal_abund function
tmp <- microtable$new(tmp_microtable$taxa_abund$Species, sample_table = tmp_microtable$sample_table)
tmp_2 <- tmp$merge_samples("Group")
tmp_microtable_group$taxa_abund$Species <- tmp_2$otu_table
# thus the data for visualization comes from the original taxa_abund list of microtable object, not the calculated abundance from otu_table
tmp <- trans_abund$new(dataset = tmp_microtable_group, taxrank = "Species", use_percentage = FALSE)
# Figure 6e
g1 <- tmp$plot_bar(others_color = "grey70", bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + ylab("Raw abundance")
cowplot::save_plot(file.path(output_dir, "Bracken_barplot_rawabund_groupmean_Species.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


######################################################
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

# show the mean of relative abundance (in the taxa_abund list of microtable object) for each group
tmp_microtable_group <- tmp_microtable$merge_samples("Group")
# extract the abundance and generate a new microtable object instead of runing cal_abund function
tmp <- microtable$new(tmp_microtable$taxa_abund$Species, sample_table = tmp_microtable$sample_table)
tmp_2 <- tmp$merge_samples("Group")
tmp_microtable_group$taxa_abund$Species <- tmp_2$otu_table
# Figure 6f
tmp <- trans_abund$new(dataset = tmp_microtable_group, taxrank = "Species")
g1 <- tmp$plot_bar(others_color = "grey70", bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + theme(legend.text = element_text(size = 11))
cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_groupmean_Species.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


######################################################
# recalculate relative abundance based on the species abundance (otu_table) for Bracken_microtable_relTRUE_K1
tmp_microtable <- clone(Bracken_microtable_relTRUE_K1)
tmp_microtable$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Kingdom")
g1 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Phylum")
g2 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

tmp <- trans_abund$new(dataset = tmp_microtable, taxrank = "Genus")
g3 <- tmp$plot_bar(others_color = "grey70", facet = c("Compartment", "Fertilization", "Cropping"), xtext_keep = FALSE)

p1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 1, nrow = 3, heights = c(1, 1, 1))

cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K1.png"), p1, base_aspect_ratio = 1.1, dpi = 300, base_height = 11)

# merge samples for each group and recalculate the relative abundance
tmp_microtable_group <- tmp_microtable$merge_samples("Group")
tmp_microtable_group$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable_group, taxrank = "Kingdom")
tmp$data_abund %<>% .[! .$Taxonomy %in% c("Bacteria", "Archaea"), ]
# Figure 6g
g1 <- tmp$plot_bar(others_color = "grey70", bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + theme(legend.text = element_text(size = 11))
cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K1_groupmean_Kingdom.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


######################################################
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



######################################################
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


tmp_microtable_group <- tmp_microtable$merge_samples("Group")
tmp_microtable_group$cal_abund()

tmp <- trans_abund$new(dataset = tmp_microtable_group, taxrank = "Kingdom")
tmp$data_abund %<>% .[! .$Taxonomy %in% c("Bacteria", "Archaea"), ]
# Figure 6h
g1 <- tmp$plot_bar(others_color = "grey70", bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + theme(legend.text = element_text(size = 11))
cowplot::save_plot(file.path(output_dir, "Bracken_barplot_relabund_recal_K3_groupmean_Kingdom.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


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




########################################################################################################################
########################################################################################################################
########################################################################################################################
# Optional
# merge the 16S QIIME2 result and Bracken result into one for the comparison
library(tibble)
library(dplyr)
library(microeco)
library(magrittr)

# load all the required data
# load 16S data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
load(input_path)
# load Bracken data
input_path <- "Output/2.Metagenome/Stage9_Bracken/Bracken_microtable_relFALSE.RData"
load(input_path)

# calculate the relative abundance at Genus level, i.e. the total sum scaling
# first filter out the samples that are not used in the metagenomic sequencing
amplicon_16S_microtable_rarefy_select <- clone(amplicon_16S_microtable_rarefy)
amplicon_16S_microtable_rarefy_select$sample_table %<>% .[rownames(.) %in% Bracken_microtable_relFALSE$sample_names(), ]
amplicon_16S_microtable_rarefy_select$tidy_dataset()
amplicon_16S_microtable_rarefy_select$cal_abund(select_cols = "Genus")
# filter out the taxa that donot belong to prokaryotes to make two data comparable
Bracken_microtable_relFALSE_prok <- clone(Bracken_microtable_relFALSE)
Bracken_microtable_relFALSE_prok$tax_table %<>% .[.$Kingdom %in% c("k__Archaea", "k__Bacteria"), ]
Bracken_microtable_relFALSE_prok$tidy_dataset()
Bracken_microtable_relFALSE_prok$cal_abund(select_cols = "Genus")

# new otu_table
tmp_16S <- amplicon_16S_microtable_rarefy_select$taxa_abund$Genus
colnames(tmp_16S) %<>% paste0("amplicon_", .)
tmp_16S %<>% rownames_to_column
tmp_Bracken <- Bracken_microtable_relFALSE_prok$taxa_abund$Genus	
colnames(tmp_Bracken) %<>% paste0("metagenome_", .)
tmp_Bracken %<>% rownames_to_column
# merge two new tables
tmp_otutable <- full_join(tmp_16S, tmp_Bracken, by = c("rowname" = "rowname"))
rownames(tmp_otutable) <- tmp_otutable[, 1]
tmp_otutable <- tmp_otutable[, -1]
# assign 0 to NA
tmp_otutable[is.na(tmp_otutable)] <- 0

# new tax_table
# select Kingdom and Genus as the example
tmp_16S <- amplicon_16S_microtable_rarefy_select$tax_table[, c(1, 6)] %>% .[.$Genus != "g__", ] %>% unique
tmp_Bracken <- Bracken_microtable_relFALSE_prok$tax_table[, c(1, 6)] %>% .[.$Genus != "g__", ] %>% unique
tmp_taxtable <- rbind.data.frame(tmp_16S, tmp_Bracken) %>% unique
# must assign row names to the table
rownames(tmp_taxtable) <- tmp_taxtable$Genus


# create microtable object
merged_16S_Kraken <- microtable$new(otu_table = tmp_otutable, tax_table = tmp_taxtable)
# complement sample_table
tmp_sampletable <- merged_16S_Kraken$sample_table
colnames(tmp_sampletable)[2] <- "Method"
# generate the Method column for the further comparison
tmp_sampletable$Method %<>% gsub("_.*", "", .)
# add original metadata into the new sample table
tmp_sampletable$raw_SampleID <- tmp_sampletable$SampleID
tmp_sampletable$raw_SampleID %<>% gsub(".*_", "", .)
tmp_sampletable <- left_join(tmp_sampletable, amplicon_16S_microtable_rarefy_select$sample_table, by = c("raw_SampleID" = "SampleID"))
# add row names
rownames(tmp_sampletable) <- tmp_sampletable[, 1]
# assign the table back to the microtable object
merged_16S_Kraken$sample_table <- tmp_sampletable

# generate abundance; rel = FALSE means the original abundance in the otu_table
merged_16S_Kraken$tidy_dataset()
merged_16S_Kraken$cal_abund(rel = FALSE)

# trans_abund
t1 <- trans_abund$new(dataset = merged_16S_Kraken, taxrank = "Genus", ntaxa = 20)
t1$plot_bar(others_color = "grey70", facet = "Method", xtext_size = 5, xtext_angle = 30)



















