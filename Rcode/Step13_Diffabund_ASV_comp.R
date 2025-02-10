

######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	stop("Please first run the script in the last step!")
}

######################################################
# load tmp_microtable_rhizo object saved in the last step
load(file.path(output_dir, "tmp_microtable_rhizo.RData"))

# load differential test objects saved in the last step
load(file.path(output_dir, "ASV_Fertilization_metagenomeSeq.RData"))
load(file.path(output_dir, "ASV_Fertilization_ALDEx2_t.RData"))
load(file.path(output_dir, "ASV_Fertilization_ALDEx2_kw.RData"))
load(file.path(output_dir, "ASV_Fertilization_DESeq2.RData"))
load(file.path(output_dir, "ASV_Fertilization_edgeR.RData"))
load(file.path(output_dir, "ASV_Fertilization_ancombc2.RData"))
load(file.path(output_dir, "ASV_Fertilization_linda.RData"))
load(file.path(output_dir, "ASV_Fertilization_maaslin2.RData"))
load(file.path(output_dir, "ASV_Fertilization_wilcox_GMPR.RData"))
load(file.path(output_dir, "ASV_Fertilization_wilcox_Wrench.RData"))
######################################################

# Then analyze intersection of significant feature among all the methods using above results
# use list to store all the results
diff_taxa_list <- list()
# select the significant ASVs for each method
diff_taxa_list[["metagenomeSeq"]] <- tmp_metagenomeSeq$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["ALDEx2_t"]] <- tmp_ALDEx2_t$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["ALDEx2_kw"]] <- tmp_ALDEx2_kw$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["DESeq2"]] <- tmp_DESeq2$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["edgeR"]] <- tmp_edgeR$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["ancombc2"]] <- tmp_ancombc2$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["linda"]] <- tmp_linda$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["maaslin2"]] <- tmp_maaslin2$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["GMPR_wilcox"]] <- tmp_GMPR_wilcox$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_list[["Wrench_wilcox"]] <- tmp_Wrench_wilcox$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)

# generate a data.frame object for different methods and features: 1 represents positive; 0 denotes negative
# extract all feature names
all_taxa_names <- rownames(tmp_microtable_rhizo$otu_table)
# generate a new list with feature vectors (0 and 1) for each method
tmp <- lapply(diff_taxa_list, function(x){
	vec <- match(all_taxa_names, x)
	names(vec) <- all_taxa_names
	vec[is.na(vec)] <- 0
	vec[vec > 0] <- 1
	vec
})
# convert the list to data.frame object
tmp1 <- do.call(rbind, tmp) %>% t %>% as.data.frame

# create microtable object
tmp_mtobj <- microtable$new(otu_table = tmp1, tax_table = tmp_microtable_rhizo$tax_table)

# use trans_venn class to analyze the features intersections
tmp_transvennobj <- trans_venn$new(tmp_mtobj, ratio = "numratio", name_joint = "-")
# only show the elements with a relative large number
tmp_transvennobj$data_summary %<>% .[.$Counts > 4, ]

# Figure 4a
g1 <- tmp_transvennobj$plot_bar(sort_samples = FALSE)
cowplot::save_plot(file.path(output_dir, "diff_methods_singlefactor_threegroups_venn_bar.png"), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 7)



# further analyze the taxonomic composition in the intersection sets
tmp_transvennobj_mt <- tmp_transvennobj$trans_comm(use_frequency = TRUE)
tmp_transvennobj_mt$otu_table %<>% .[, colnames(.) %in% rownames(tmp_transvennobj$data_summary)]
tmp_transvennobj_mt$tidy_dataset()
# calculate the relative abundance, i.e. ratio here
tmp_transvennobj_mt$cal_abund()
# save the proportion data at Genus level to the directory
write.csv(tmp_transvennobj_mt$taxa_abund$Genus, file.path(output_dir, "diff_methods_singlefactor_threegroups_venn_comp_Genus.csv"))

# visualize the proportion data at Genus level
library(tibble)
tmp_x <- tmp_transvennobj$data_summary %>% rownames_to_column %>% .[order(.$Counts, decreasing = TRUE), ] %>% .$rowname

tmp_transvennobj_mt_abund <- trans_abund$new(dataset = tmp_transvennobj_mt, taxrank = "Genus", ntaxa = 10, use_percentage = TRUE)
# Figure 4a
g2 <- tmp_transvennobj_mt_abund$plot_bar(bar_full = FALSE, legend_text_italic = T, xtext_angle = 50, order_x = tmp_x) + ylab("Ratio (%)") + 
	theme(legend.position = "left", plot.margin = unit(c(0, 0, 0, 4), "cm"))
cowplot::save_plot(file.path(output_dir, "diff_methods_singlefactor_threegroups_venn_comp_Genus.png"), g2, base_aspect_ratio = 1.5, dpi = 300, base_height = 8)











