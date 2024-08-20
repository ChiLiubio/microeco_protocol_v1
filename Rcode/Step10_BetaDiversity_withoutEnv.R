## 
## Beta diversity analysis that do not refer to environmental factors
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
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
load(input_path)
###########################


# select data for rhizosphere soil
tmp_microtable_rhizo <- clone(amplicon_16S_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

# select rarefied data for rhizosphere soil
tmp_microtable_rarefy_rhizo <- clone(amplicon_16S_microtable_rarefy)
tmp_microtable_rarefy_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rarefy_rhizo$tidy_dataset()


# calculate Bray-Curtis and UniFrac dissimilarity for rarefied rhizosphere data
measure <- "bray"
tmp_microtable_rarefy_rhizo$cal_betadiv(method = measure, unifrac = TRUE)
# save the Bray-Curtis dissimilarity matrix
write.csv(tmp_microtable_rarefy_rhizo$beta_diversity[[measure]], file.path(output_dir, "BetaDiv_rarefy_rhizo_bray.csv"))
# save the weighted UniFrac dissimilarity matrix
write.csv(tmp_microtable_rarefy_rhizo$beta_diversity$wei_unifrac, file.path(output_dir, "BetaDiv_rarefy_rhizo_weightedUniFrac.csv"))
# save the unweighted UniFrac dissimilarity matrix
write.csv(tmp_microtable_rarefy_rhizo$beta_diversity$unwei_unifrac, file.path(output_dir, "BetaDiv_rarefy_rhizo_unweightedUniFrac.csv"))

# calculate Aitchison dissimilarity for non-rarefied rhizosphere data
measure <- "aitchison"
tmp_microtable_rhizo$cal_betadiv(method = measure)
write.csv(tmp_microtable_rhizo$beta_diversity[[measure]], file.path(output_dir, "BetaDiv_rarefy_rhizo_aitchison.csv"))


###########################
# PCoA based on Bray-Curtis distance with rarefied data
measure <- "bray"
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Group", measure = measure)
t1$cal_ordination(method = "PCoA", ncomp = 2)
# 'PCo1' column means the score of first coordinate axis; 'PCo2' column represents the score of second coordinate axis; Other columns are the metadata
write.csv(t1$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_rarefy_PCoA_Bray_Score.csv"))
# visualization
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_PCoA_Bray.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# PCoA based on Aitchison distance with raw data
measure <- "aitchison"
t2 <- trans_beta$new(dataset = tmp_microtable_rhizo, group = "Group", measure = measure)
t2$cal_ordination(method = "PCoA", ncomp = 2)
write.csv(t2$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_nonrarefy_PCoA_Aitchison_Score.csv"))

g2 <- t2$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_nonrarefy_PCoA_Aitchison.png"), g2, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


# NMDS based on Bray-Curtis distance with rarefied data
measure <- "bray"
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Group", measure = measure)
t1$cal_ordination(method = "NMDS")
# 'MDS1' column means the score of first coordinate axis; 'MDS2' column represents the score of second coordinate axis; Other columns are the metadata
write.csv(t1$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_rarefy_NMDS_Bray_Score.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"), NMDS_stress_pos = c(1.1, 1.4), NMDS_stress_text_prefix = "Stress: ")
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_NMDS_Bray.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# NMDS based on Aitchison distance with raw data
measure <- "aitchison"
t2 <- trans_beta$new(dataset = tmp_microtable_rhizo, group = "Group", measure = measure)
t2$cal_ordination(method = "NMDS")
write.csv(t1$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_nonrarefy_NMDS_Aitchison_Score.csv"))

g2 <- t2$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"), NMDS_stress_pos = c(1.1, 1.3), NMDS_stress_text_prefix = "Stress: ")
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_nonrarefy_NMDS_Aitchison.png"), g2, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)



# PCA at genus level
# first generate a microtable object with genera as features
tmp_microtable_rarefy_rhizo_genus <- tmp_microtable_rarefy_rhizo$merge_taxa(taxa = "Genus")
# remove unidentified genera
tmp_microtable_rarefy_rhizo_genus$tax_table %<>% .[.$Genus != "g__", ]
# also delete the genera with number in the names
tmp_microtable_rarefy_rhizo_genus$tax_table %<>% .[!grepl("\\d+", .$Genus), ]
tmp_microtable_rarefy_rhizo_genus$tidy_dataset()
rownames(tmp_microtable_rarefy_rhizo_genus$tax_table) <- rownames(tmp_microtable_rarefy_rhizo_genus$otu_table) <- gsub("g__", "", tmp_microtable_rarefy_rhizo_genus$tax_table$Genus)


t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo_genus)
t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
write.csv(t1$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_rarefy_PCA_Genus_Score.csv"))
# Columns PC1-PC3 represent the loadings in each axis. 'dist' is sum of squares for loadings of PC1 and PC2 and used to order the features.
write.csv(t1$res_ordination$loading, file.path(output_dir, "BetaDiv_Rhizo_rarefy_PCA_Genus_Loading.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", loading_arrow = TRUE)
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_PCA_Genus.png"), g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)


# DCA
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo_genus)
t1$cal_ordination(method = "DCA", scale_species = TRUE)
write.csv(t1$res_ordination$scores, file.path(output_dir, "BetaDiv_Rhizo_rarefy_DCA_Genus_Score.csv"))
write.csv(t1$res_ordination$loading, file.path(output_dir, "BetaDiv_Rhizo_rarefy_DCA_Genus_Loading.csv"))

g1 <- t1$plot_ordination(plot_color = "Group", loading_arrow = TRUE)
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_DCA_Genus.png"), g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)


# group distance transformation and boxplot

# group distances within Cropping treatments
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Cropping", measure = "bray")
t1$cal_group_distance(within_group = TRUE, by_group = "Fertilization")
# manipulate res_group_distance to remove the distances between same group
t1$res_group_distance %<>% .[!.$Fertilization %in% c("CK vs CK", "NPK vs NPK", "NPKS vs NPKS"), ]
# save the converted distance values (long-format table) to local directory
write.csv(t1$res_group_distance, file.path(output_dir, "BetaDiv_Rhizo_rarefy_bray_Cropping_within.csv"))
t1$cal_group_distance_diff(method = "wilcox")
# save the differential test results to directory
write.csv(t1$res_group_distance_diff, file.path(output_dir, "BetaDiv_Rhizo_rarefy_bray_Cropping_within_diff.csv"))
g1 <- t1$plot_group_distance()
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_bray_Cropping_within_boxplot.png"), g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# group distances within Fertilization treatments
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Fertilization", measure = "bray")
t1$cal_group_distance(within_group = TRUE, by_group = "Cropping")
t1$cal_group_distance_diff(method = "anova")
g1 <- t1$plot_group_distance()
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_bray_Fertilization_within_boxplot.png"), g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# group distances between different Fertilization treatments
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Fertilization", measure = "bray")
t1$cal_group_distance(within_group = FALSE, by_group = "Cropping")
t1$cal_group_distance_diff(method = "anova")
g1 <- t1$plot_group_distance()
cowplot::save_plot(file.path(output_dir, "BetaDiv_Rhizo_rarefy_bray_Fertilization_between_boxplot.png"), g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)


# one-way perMANOVA for all groups
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_manova(manova_all = TRUE, group = "Group")
write.csv(t1$res_manova, file.path(output_dir, "BetaDiv_rarefy_rhizo_perMANOVA_oneway_all.csv"))

# one-way perMANOVA for paired groups
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_manova(manova_all = FALSE, group = "Group")
write.csv(t1$res_manova, file.path(output_dir, "BetaDiv_rarefy_rhizo_perMANOVA_oneway_paired.csv"))

# one-way perMANOVA for paired groups with constraints
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_manova(manova_all = FALSE, group = "Fertilization", by_group = "Cropping")
write.csv(t1$res_manova, file.path(output_dir, "BetaDiv_rarefy_rhizo_perMANOVA_oneway_bygroup.csv"))

# two-way perMANOVA for Cropping and Fertilization
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_manova(manova_set = "Cropping*Fertilization")
write.csv(t1$res_manova, file.path(output_dir, "BetaDiv_rarefy_rhizo_perMANOVA_twoway.csv"))


# ANOSIM overall test for all groups
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_anosim(paired = FALSE, group = "Group")
write.csv(t1$res_anosim, file.path(output_dir, "BetaDiv_rarefy_rhizo_ANOSIM_all.csv"))

# paired test
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_anosim(paired = TRUE, group = "Group")
write.csv(t1$res_anosim, file.path(output_dir, "BetaDiv_rarefy_rhizo_ANOSIM_paired.csv"))

# paired test with group constraints
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, measure = "bray")
t1$cal_anosim(paired = TRUE, group = "Fertilization", by_group = "Cropping")
write.csv(t1$res_anosim, file.path(output_dir, "BetaDiv_rarefy_rhizo_ANOSIM_bygroup.csv"))


# permdisp method: test the dispersion of variances
t1 <- trans_beta$new(dataset = tmp_microtable_rarefy_rhizo, group = "Group", measure = "bray")
t1$cal_betadisper()
# The result res_betadisper is not data.frame class, we directly save the printed output to txt file
capture.output(t1$res_betadisper, file = file.path(output_dir, "BetaDiv_rarefy_rhizo_permdisp.txt"))












