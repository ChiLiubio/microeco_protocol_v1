

##################################################################
##  Genus level     relative abundance      single factor    ####
##################################################################

group <- "Fertilization"
taxlevel <- "Genus"

# LEfSe
method <- "lefse"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))

# Wilcoxon test with arcsine transformation
method <- "wilcox"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))

# one-way ANOVA with arcsine transformation
method <- "anova"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))

# Dunn's Kruskal-Wallis Multiple Comparisons with arcsine transformation
method <- "KW_dunn"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))

# Beta regression
method <- "betareg"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = group, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_singlefactor_", method, ".csv")))




##################################################################
##  Genus level     relative abundance      multiple factor   ####
##################################################################

formula <- "Cropping+Fertilization"
taxlevel <- "Genus"

method <- "anova"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))

method <- "lm"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))

method <- "betareg"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))

method <- "glmm_beta"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_Genus_multifactor_", method, ".csv")))




##############################################################################
##  Multiple taxonomic level     relative abundance      single factor    ####
##############################################################################

group <- "Fertilization"
taxlevel <- "all"

method <- "lefse"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, alpha = 0.01, lefse_subgroup = NULL, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))

method <- "anova"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))

method <- "KW_dunn"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, group = group, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_singlefactor_", method, ".csv")))



##############################################################################
##  Multiple taxonomic level     relative abundance      multiple factor  ####
##############################################################################

formula <- "Cropping+Fertilization"
taxlevel <- "all"

method <- "anova"
t1 <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))

method <- "lm"
t1_lm <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))

method <- "betareg"
t1_betareg <- trans_diff$new(dataset = tmp_microtable_rhizo, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.001)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_multitax_multifactor_", method, ".csv")))



##################################################################################################
##  Random effects   Multiple or single high taxonomic level     relative abundance      ##############
##################################################################################################

formula <- "Cropping + Fertilization + (1|Compartment)"
taxlevel <- "Genus"

method <- "lme"
t1 <- trans_diff$new(dataset = amplicon_16S_microtable, method = method, formula = formula, taxa_level = taxlevel, transformation = "AST", filter_thres = 0.005)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_mixedeff_Genus_", method, ".csv")))

method <- "glmm_beta"
t1 <- trans_diff$new(dataset = amplicon_16S_microtable, method = method, formula = formula, taxa_level = taxlevel, filter_thres = 0.005)
write.csv(t1$res_diff, file.path(output_dir, paste0("Diff_abund_test_mixedeff_Genus_", method, ".csv")))





