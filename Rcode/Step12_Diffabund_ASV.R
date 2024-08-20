## 
## Differential abundance test at ASV level
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
######################################################
# preprocess data
# select data for rhizosphere soil
tmp_microtable_rhizo <- clone(amplicon_16S_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()
# get the relative abundance data: taxa_abund list stored in the object
tmp_microtable_rhizo$cal_abund(rel = TRUE)

# delete the ASV with very low abundance
tmp_microtable_rhizo$filter_taxa(rel_abund = 0.0001)

# As an example, we delete the genera without clear taxonomic information to reduce running time at ASV level
tmp_microtable_rhizo$tax_table %<>% .[.$Genus != "g__", ]
tmp_microtable_rhizo$tax_table %<>% .[!grepl("\\d+", .$Genus), ]
tmp_microtable_rhizo$tidy_dataset()

save(tmp_microtable_rhizo, file = file.path(output_dir, "tmp_microtable_rhizo.RData"), compress = TRUE)


#############################################
##  ASV level   single factor  2 groups    ##
#############################################

group <- "Cropping"
taxlevel <- "ASV"

# metastat
method <- "metastat"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# metagenomeSeq
method <- "metagenomeSeq"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# ALDEx2_t
method <- "ALDEx2_t"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# In result, 'wi.ep': expected p-value of the Wilcoxon Rank Sum test; 'wi.eBH': corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# ALDEx2_kw
method <- "ALDEx2_kw"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# In result, 'kw.ep': expected p-value of the Kruskal-Wallis test for each feature; 'kw.eBH': corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# DESeq2
method <- "DESeq2"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# lfcSE gives the standard error of the log2FoldChange; stat is the Wald statistic
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# edgeR
# based on exactTest function of edgeR package
method <- "edgeR"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# For results, logFC: log2-fold-change; logCPM: average log2-counts-per-million
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# ancombc2
method <- "ancombc2"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# In the table; lfc: log fold changes; se: standard errors (SEs); W: test statistics
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# linda
method <- "linda"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
# In the table, log2FoldChange: bias-corrected coefficients; lfcSE: standard errors of the coefficients; stat: log2FoldChange / lfcSE
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# maaslin2
tmp_dir <- "tmp_maaslin2"
dir.create(tmp_dir)
method <- "maaslin2"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, standardize = FALSE, fixed_effects = group, 
	tmp_input_maaslin2 = file.path(tmp_dir, "maaslin2_tmp_input"), tmp_output_maaslin2 = file.path(tmp_dir, "maaslin2_tmp_output"))
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, ".csv")))

# GMPR + wilcox
tmp <- trans_norm$new(tmp_microtable_rhizo)
norm_obj <- tmp$norm(method = "GMPR")
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)

method <- "wilcox"
tmp <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, "_GMPR.csv")))

# Wrench + wilcox
tmp <- trans_norm$new(tmp_microtable_rhizo)
norm_obj <- tmp$norm(method = "Wrench", condition = group)
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)

method <- "wilcox"
tmp <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_twogroups_", method, "_Wrench.csv")))





#############################################
##  ASV level   single factor  3 groups  ##
#############################################

group <- "Fertilization"
taxlevel <- "ASV"

# metagenomeSeq
method <- "metagenomeSeq"
tmp_metagenomeSeq <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_metagenomeSeq$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, ".csv")))
save(tmp_metagenomeSeq, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)


# ALDEx2_t
method <- "ALDEx2_t"
tmp_ALDEx2_t <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ALDEx2_t$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, ".csv")))
save(tmp_ALDEx2_t, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# ALDEx2_kw
method <- "ALDEx2_kw"
tmp_ALDEx2_kw <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
save(tmp_ALDEx2_kw, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# DESeq2
method <- "DESeq2"
tmp_DESeq2 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
save(tmp_DESeq2, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# edgeR
method <- "edgeR"
tmp_edgeR <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
save(tmp_edgeR, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# ancombc2
method <- "ancombc2"
tmp_ancombc2 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ancombc2$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, ".csv")))
save(tmp_ancombc2, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# linda
method <- "linda"
tmp_linda <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_linda$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, ".csv")))
save(tmp_linda, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# maaslin2
tmp_dir <- "tmp_maaslin2"
dir.create(tmp_dir)
method <- "maaslin2"
tmp_maaslin2 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, standardize = FALSE, fixed_effects = group, 
	tmp_input_maaslin2 = file.path(tmp_dir, "maaslin2_tmp_input"), tmp_output_maaslin2 = file.path(tmp_dir, "maaslin2_tmp_output"))
write.csv(tmp_maaslin2$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, ".csv")))
save(tmp_maaslin2, file = file.path(output_dir, paste0("ASV_Fertilization_", method, ".RData")), compress = TRUE)

# GMPR + wilcox
tmp <- trans_norm$new(tmp_microtable_rhizo)
norm_obj <- tmp$norm(method = "GMPR")
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)

method <- "wilcox"
tmp_GMPR_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_GMPR_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, "_GMPR.csv")))
save(tmp_GMPR_wilcox, file = file.path(output_dir, paste0("ASV_Fertilization_", method, "_GMPR.RData")), compress = TRUE)

# Wrench + wilcox
tmp <- trans_norm$new(tmp_microtable_rhizo)
norm_obj <- tmp$norm(method = "Wrench", condition = group)
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)

method <- "wilcox"
tmp_Wrench_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_Wrench_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_singlefactor_threegroups_", method, "_Wrench.csv")))
save(tmp_Wrench_wilcox, file = file.path(output_dir, paste0("ASV_Fertilization_", method, "_Wrench.RData")), compress = TRUE)



#############################################
##  ASV level   multiple factors  ##
#############################################

formula <- "Cropping+Fertilization"
taxlevel <- "ASV"

# DESeq2
method <- "DESeq2"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_multifactor_", method, ".csv")))

# ancombc2
method <- "ancombc2"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = NULL, fix_formula = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_multifactor_", method, ".csv")))

# linda
method <- "linda"
tmp <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_test_ASV_multifactor_", method, ".csv")))





