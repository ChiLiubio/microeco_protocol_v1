
######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
######################################################
# create an output directory if it does not exist
output_dir <- "Output/2.Metagenome/StageⅧ_HUMAnN"
if(! dir.exists(output_dir)){
	stop("Please first run the last step to import files!")
}
# load data
input_path <- file.path(output_dir, "MetaCyc_microtable.RData")
load(input_path)
input_path <- file.path(output_dir, "KEGG_microtable.RData")
load(input_path)
######################################################

tmp_microtable <- clone(MetaCyc_microtable)



##################
# Abundance visualization

# rel = FALSE means the abundance of enrichment taxa is not converted to relative abundance
tmp_microtable$cal_abund(select_cols = 1:3, rel = FALSE)
# donot use percentage as the abundance is RPK, not relative abundance
tmp <- trans_abund$new(tmp_microtable, taxrank = "Superclass1", ntaxa = 10, use_percentage = FALSE)
# bar_full = FALSE show original abundance instead of normalized 0-1
g1 <- tmp$plot_bar(facet = c("Compartment", "Fertilization", "Cropping"), bar_full = FALSE, xtext_size = 4) + ylab("Abundance (RPK)")
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)

# select both pathway and taxa
tmp_microtable$cal_abund(select_cols = c("Superclass1", "Phylum", "Genus"), rel = TRUE)
# delete_taxonomy_lineage = FALSE: remain the original names in front of target level
tmp <- trans_abund$new(tmp_microtable, taxrank = "Phylum", ntaxa = 10, delete_taxonomy_lineage = FALSE)
g1 <- tmp$plot_bar(facet = c("Compartment", "Fertilization", "Cropping"), xtext_size = 4)
cowplot::save_plot(file.path(output_dir, "MetaCyc_barplot_Superclass1_Phylum.png"), g1, base_aspect_ratio = 1.8, dpi = 300, base_height = 6)




###################################
# Differential test of features

# functional biomarker
# select three pathway levels; calculate relative abundance
tmp_microtable$cal_abund(select_cols = 1:3, rel = TRUE)
tmp <- trans_diff$new(tmp_microtable, method = "lefse", group = "Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "MetaCyc_lefse_pathway_Fertilization.csv"))
g1 <- tmp$plot_diff_bar(use_number = 1:20)
cowplot::save_plot(file.path(output_dir, "MetaCyc_lefse_Pathway_Fertilization.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# taxonomic biomarker
# select all the taxonomic levels
tmp_microtable$cal_abund(select_cols = 4:10, rel = TRUE)
# p_adjust_method = "none" shut down the p value adjustment
tmp <- trans_diff$new(tmp_microtable, method = "lefse", group = "Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "MetaCyc_lefse_Taxa_Fertilization.csv"))
g1 <- tmp$plot_diff_bar(threshold = 1)
cowplot::save_plot(file.path(output_dir, "MetaCyc_lefse_Taxa_Fertilization.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)


# multiple factors test on pathway; ANOVA; non-relative abundance with log trans
tmp_microtable$cal_abund(select_cols = 1:3, rel = FALSE)
tmp <- trans_diff$new(dataset = tmp_microtable, method = "anova", formula = "Cropping*Fertilization*Compartment", taxa_level = "Superclass1", transformation = "log")
# By default, the x and y axis text order come from the clustering.
# To adjust the text order, please use parameters text_x_order or text_y_order
g1 <- tmp$plot_diff_bar(text_x_order = tmp$res_diff$Factors[1:7])
cowplot::save_plot(file.path(output_dir, "MetaCyc_anova_multiway_Superclass1_log.png"), g1, base_aspect_ratio = 1.4, dpi = 300, base_height = 5)






###################################
# PCA at pathway level
# prepare data
tmp_microtable_pathway <- tmp_microtable$merge_taxa(taxa = "pathway")
rownames(tmp_microtable_pathway$tax_table) <- rownames(tmp_microtable_pathway$otu_table) <- tmp_microtable_pathway$tax_table$pathway
# perform PCA
t1 <- trans_beta$new(dataset = tmp_microtable_pathway)
t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Compartment", loading_arrow = TRUE, loading_text_italic = FALSE)
cowplot::save_plot(file.path(output_dir, "MetaCyc_PCA_pathway.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)












