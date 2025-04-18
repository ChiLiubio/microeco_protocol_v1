## 
## Statistical analysis and visualization on metabolic pathways
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
######################################################
# check the data
output_dir <- "Output/2.Metagenome/Stage8_HUMAnN"
# load data
input_path <- file.path(output_dir, "MetaCyc_microtable.RData")
if(!file.exists(input_path)){
	stop("Please first run the script in the last step to import files!")
}
load(input_path)

# load(file.path(output_dir, "KEGG_microtable.RData"))
######################################################

tmp_microtable <- clone(MetaCyc_microtable)

# select rhizosphere data
tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

# merge samples into one accroding to the groups
tmp_microtable_rhizo_group <- tmp_microtable_rhizo$merge_samples("Group")


####################################
# Abundance visualization

# select pathway
# rel = FALSE: the abundance of enrichment taxa is not converted to relative abundance
tmp_microtable$cal_abund(select_cols = 1:3, rel = FALSE)
# donot use percentage because the abundance is RPK, not relative abundance
tmp <- trans_abund$new(tmp_microtable, taxrank = "Superclass1", ntaxa = 10, use_percentage = FALSE)
# bar_full = FALSE: show original abundance instead of normalized 0-1
g1 <- tmp$plot_bar(facet = c("Compartment", "Fertilization", "Cropping"), bar_full = FALSE, xtext_size = 4) + ylab("Abundance (RPK)")
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)

# select both pathway and taxa
tmp_microtable$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
# delete_taxonomy_lineage = FALSE: show the original names in front of target level
tmp <- trans_abund$new(tmp_microtable, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = FALSE)
g1 <- tmp$plot_bar(facet = c("Compartment", "Fertilization", "Cropping"), xtext_size = 4)
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1_Phylum.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)


# visualization for the merged samples of rhizosphere data
# rel = FALSE: the abundance of enrichment taxa is not converted to relative abundance
tmp_microtable_rhizo_group$cal_abund(select_cols = 1:3, rel = FALSE)
# donot use percentage as the abundance is RPK, not relative abundance
tmp <- trans_abund$new(tmp_microtable_rhizo_group, taxrank = "Superclass1", ntaxa = 10, use_percentage = FALSE)
# Figure 6a
g1 <- tmp$plot_bar(bar_full = FALSE, xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + ylab("Abundance (RPK)") + theme(legend.text = element_text(size = 10.5))
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1_rhizo_merge.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)


# select both pathway and taxa
tmp_microtable_rhizo_group$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
# delete_taxonomy_lineage = FALSE: remain the original names in front of target level
tmp <- trans_abund$new(tmp_microtable_rhizo_group, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = FALSE)
# Figure 6b
g1 <- tmp$plot_bar(xtext_size = 15, xtext_angle = 30, barwidth = 0.618) + theme(legend.text = element_text(size = 10.5))
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1_rhizo_merge_Phylum.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)





######################################################
# Differential test of features based on the rhizosphere data

# pathway
# select all the three pathway levels; calculate relative abundance
tmp_microtable_rhizo$cal_abund(select_cols = 1:3, rel = TRUE)
tmp <- trans_diff$new(tmp_microtable_rhizo, method = "lefse", group = "Fertilization")
# The format of result is same with that in amplicon sequencing data part
write.csv(tmp$res_diff, file.path(output_dir, "MetaCyc_lefse_pathway_Fertilization.csv"))
# Figure 6c
g1 <- tmp$plot_diff_bar(use_number = 1:20, width = 0.618) + theme(legend.position = "top", legend.justification = c("right", "center"), legend.location = "plot")
cowplot::save_plot(file.path(output_dir, "MetaCyc_lefse_Pathway_Fertilization_rhizo.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# taxa
# select all the taxonomic levels
tmp_microtable_rhizo$cal_abund(select_cols = 4:10, rel = TRUE)
# p_adjust_method = "none" shut down the p value adjustment
tmp <- trans_diff$new(tmp_microtable_rhizo, method = "lefse", group = "Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "MetaCyc_lefse_Taxa_Fertilization.csv"))
# Figure 6d
g1 <- tmp$plot_diff_bar(threshold = 3, width = 0.618) + theme(legend.position = "top", legend.justification = c("right", "center"), legend.location = "plot")
cowplot::save_plot(file.path(output_dir, "MetaCyc_lefse_Taxa_Fertilization_rhizo.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)


# multiple factors test on pathway; ANOVA; non-relative abundance with log transformation
tmp_microtable_rhizo$cal_abund(select_cols = 1:3, rel = FALSE)
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = "anova", formula = "Cropping*Fertilization", taxa_level = "Superclass1", transformation = "log")
# By default, the x and y axis text order come from the clustering.
# To adjust the text order, please use parameters text_x_order or text_y_order
g1 <- tmp$plot_diff_bar(text_x_order = tmp$res_diff$Factors[1:7])
cowplot::save_plot(file.path(output_dir, "MetaCyc_anova_multiway_Superclass1_log_rhizo.png"), g1, base_aspect_ratio = 1.4, dpi = 300, base_height = 5)





###################################
# PCA at pathway level
# prepare data
tmp_microtable_pathway <- tmp_microtable$merge_taxa(taxa = "pathway")
rownames(tmp_microtable_pathway$tax_table) <- rownames(tmp_microtable_pathway$otu_table) <- tmp_microtable_pathway$tax_table$pathway
# perform PCA
t1 <- trans_beta$new(dataset = tmp_microtable_pathway)
t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
write.csv(t1$res_ordination$scores, file.path(output_dir, "MetaCyc_PCA_pathway_Score.csv"))
# Columns PC1-PC3 represent the loadings in each axis. 'dist' is sum of squares for loadings of PC1 and PC2 and used to order the features.
write.csv(t1$res_ordination$loading, file.path(output_dir, "MetaCyc_PCA_pathway_Loading.csv"))
# loading_arrow = TRUE: show the arrows of loadings
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Compartment", loading_arrow = TRUE, loading_text_italic = FALSE)
cowplot::save_plot(file.path(output_dir, "MetaCyc_PCA_pathway.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)



###################################
# functional diversity
tmp_microtable_pathway <- tmp_microtable$merge_taxa(taxa = "pathway")
tmp_microtable_pathway$cal_alphadiv()

write.csv(tmp_microtable_pathway$alpha_diversity, file.path(output_dir, "MetaCyc_functional_alphadiversity.csv"))






