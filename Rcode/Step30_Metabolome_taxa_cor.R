## 
## Correlation analysis between metabolites and genera
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)
######################################################
output_dir <- "Output/3.Metabolome/Stage10_Metabolome"
# load metabolome microtable object
input_path <- file.path(output_dir, "metab_microtable.RData")
if(!file.exists(input_path)){
	stop("Please first run the script in step27 !")
}
load(input_path)

# load metabolites selected in the previous step
input_path <- file.path(output_dir, "Classification_rf_select_features.RData")
if(! file.exists(input_path)){
	stop("Please first run the scripts in step28 !")
}
load(input_path)

# load amplicon 16S microtable object
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)

######################################################
set.seed(123)

tmp_metab_microtable <- clone(metab_microtable)
tmp_metab_microtable$filter_taxa(freq = 0.2)


# first perform differential test to select significant genera
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


tmp_sel_metab_table <- tmp_metab_microtable$otu_table %>% t %>% as.data.frame
tmp_sel_metab_table %<>% .[, colnames(.) %in% tmp_sel_features]


# perform correlation analysis for different compartments
tmp_amplicon_microtable_rhizo <- clone(tmp_amplicon_microtable)
tmp_amplicon_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_amplicon_microtable_rhizo$tidy_dataset()

tmp_transenv <- trans_env$new(dataset = tmp_amplicon_microtable_rhizo, add_data = tmp_sel_metab_table, standardize = FALSE)
tmp_transenv$cal_cor(cor_method = "spearman", use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor, file.path(output_dir, "Metabolome_Genera_cor_spearman_Rhizosphere.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both", cluster_height_rows = 0.3, cluster_height_cols = 0.15)
cowplot::save_plot(file.path(output_dir, "Metabolome_Genera_cor_spearman_Rhizosphere.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)



tmp_amplicon_microtable_bulk <- clone(tmp_amplicon_microtable)
tmp_amplicon_microtable_bulk$sample_table %<>% .[.$Compartment == "Bulk soil", ]
tmp_amplicon_microtable_bulk$tidy_dataset()

tmp_transenv <- trans_env$new(dataset = tmp_amplicon_microtable_bulk, add_data = tmp_sel_metab_table, standardize = FALSE)
tmp_transenv$cal_cor(cor_method = "spearman", use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor, file.path(output_dir, "Metabolome_Genera_cor_spearman_Bulk.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both", cluster_height_rows = 0.3, cluster_height_cols = 0.15)
cowplot::save_plot(file.path(output_dir, "Metabolome_Genera_cor_spearman_Bulk.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)


#########################################
# correlation network

# rhizosphere
tmp_amplicon_microtable_rhizo_genus <- tmp_amplicon_microtable_rhizo$merge_taxa("Genus")
tmp_amplicon_microtable_rhizo_genus$filter_taxa(rel_abund = 0.0005)
tmp_amplicon_microtable_rhizo_genus$tax_table %<>% .[.$Genus != "g__", ] %>% .[!grepl("\\d", .$Genus), ]
tmp_amplicon_microtable_rhizo_genus$tidy_dataset()
rownames(tmp_amplicon_microtable_rhizo_genus$otu_table) <- rownames(tmp_amplicon_microtable_rhizo_genus$tax_table) <- gsub("g__", "", tmp_amplicon_microtable_rhizo_genus$tax_table$Genus)

tmp_rhizo <- trans_network$new(dataset = tmp_amplicon_microtable_rhizo_genus, cor_method = "spearman", add_data = tmp_sel_metab_table)

# create network
tmp_rhizo$cal_network(COR_p_thres = 0.01, COR_cut = 0.8)
# calculate modules
tmp_rhizo$cal_module()
# calculate network attributes
tmp_rhizo$cal_network_attr()
# node properties
tmp_rhizo$get_node_table()
write.csv(tmp_rhizo$res_node_table, file.path(output_dir, "network_rhizo_metabolites_genus_spearman_nodeproperties.csv"))
# save network to gexf format file which can be opened by Gephi software directly
tmp_rhizo$save_network(file.path(output_dir, "network_rhizo_metabolites_genus_spearman.gexf"))
g1 <- tmp_rhizo$plot_network(method = "ggraph")
cowplot::save_plot(file.path(output_dir, "network_rhizo_metabolites_genus_spearman_ggraph.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 8)



# bulk soil
tmp_amplicon_microtable_bulk_genus <- tmp_amplicon_microtable_bulk$merge_taxa("Genus")
tmp_amplicon_microtable_bulk_genus$filter_taxa(rel_abund = 0.0005)
tmp_amplicon_microtable_bulk_genus$tax_table %<>% .[.$Genus != "g__", ] %>% .[!grepl("\\d", .$Genus), ]
tmp_amplicon_microtable_bulk_genus$tidy_dataset()
rownames(tmp_amplicon_microtable_bulk_genus$otu_table) <- rownames(tmp_amplicon_microtable_bulk_genus$tax_table) <- gsub("g__", "", tmp_amplicon_microtable_bulk_genus$tax_table$Genus)

tmp_bulk <- trans_network$new(dataset = tmp_amplicon_microtable_bulk_genus, cor_method = "spearman", add_data = tmp_sel_metab_table)

# create network
tmp_bulk$cal_network(COR_p_thres = 0.01, COR_cut = 0.8)
# calculate modules
tmp_bulk$cal_module()
# calculate network attributes
tmp_bulk$cal_network_attr()
# node properties
tmp_bulk$get_node_table()
write.csv(tmp_bulk$res_node_table, file.path(output_dir, "network_bulk_metabolites_genus_spearman_nodeproperties.csv"))
# save network to gexf file
tmp_bulk$save_network(file.path(output_dir, "network_bulk_metabolites_genus_spearman.gexf"))
g1 <- tmp_bulk$plot_network(method = "ggraph")
cowplot::save_plot(file.path(output_dir, "network_bulk_metabolites_genus_spearman_ggraph.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 8)







