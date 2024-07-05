

######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/3.Metabolome/StageⅩ_Metabolome"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load metabolome microtable object
input_path <- file.path(output_dir, "metab_microtable.RData")
load(input_path)
######################################################
set.seed(123)

tmp_metab_microtable <- clone(metab_microtable)
tmp_metab_microtable$filter_taxa(freq = 0.2)

######################################################

# load amplicon 16S microtable object
input_path <- "./Output/1.Amplicon/StageⅡ_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in StageⅡ !")
}
load(input_path)

# select data for rhizosphere soil
tmp_amplicon_microtable <- clone(amplicon_16S_microtable)

formula <- "Cropping+Fertilization+Compartment"
taxlevel <- "Genus"

# perform differential test to select significant genera
# use beta regression method (glmm_beta method without random effects) as an example to show the correlation between differential taxa and environmental variables
tmp_betareg <- trans_diff$new(dataset = tmp_amplicon_microtable, method = "glmm_beta", formula = formula, taxa_level = taxlevel, filter_thres = 0.0005)

# only select those taxa significant in RC compared to CC
select_taxa <- tmp_betareg$res_diff %>% .[.$Factors == "CroppingRC", ] %>% 
	.[grepl("*", .$Significance, fixed = TRUE), ] %>% 
	.$Taxa %>% unique %>% .[!grepl("\\d", .)]


# load metabolites selected in the previous step
load(file.path(output_dir, "Classification_rf_select_features.RData"))

tmp2 <- tmp_metab_microtable$otu_table %>% t %>% as.data.frame
tmp2 %<>% .[, colnames(.) %in% tmp_sel_features]


# perform correlation analysis for different compartments
tmp_amplicon_microtable_rhizo <- clone(tmp_amplicon_microtable)
tmp_amplicon_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_amplicon_microtable_rhizo$tidy_dataset()

tmp_transenv <- trans_env$new(dataset = tmp_amplicon_microtable_rhizo, add_data = tmp2, standardize = FALSE)
tmp_transenv$cal_cor(cor_method = "spearman", use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor, file.path(output_dir, "Metabolome_Genera_cor_spearman_Rhizosphere.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both", cluster_height_rows = 0.3, cluster_height_cols = 0.15)
cowplot::save_plot(file.path(output_dir, "Metabolome_Genera_cor_spearman_Rhizosphere.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)



tmp_amplicon_microtable_bulk <- clone(tmp_amplicon_microtable)
tmp_amplicon_microtable_bulk$sample_table %<>% .[.$Compartment == "Bulk soil", ]
tmp_amplicon_microtable_bulk$tidy_dataset()

tmp_transenv <- trans_env$new(dataset = tmp_amplicon_microtable_bulk, add_data = tmp2, standardize = FALSE)
tmp_transenv$cal_cor(cor_method = "spearman", use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor, file.path(output_dir, "Metabolome_Genera_cor_spearman_Bulk.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both", cluster_height_rows = 0.3, cluster_height_cols = 0.15)
cowplot::save_plot(file.path(output_dir, "Metabolome_Genera_cor_spearman_Bulk.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)



