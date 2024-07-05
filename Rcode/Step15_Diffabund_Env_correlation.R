


######################################################

formula <- "Cropping+Fertilization"
taxlevel <- "Genus"

# use beta regression method as an example to show the correlation between differential taxa and environmental variables
tmp_transdiff_betareg <- trans_diff$new(dataset = tmp_microtable_rhizo, method = "betareg", formula = formula, taxa_level = taxlevel, filter_thres = 0.001)

# select the feature with extremely significance in rotational cropping and fertilization treatments
select_taxa <- tmp_transdiff_betareg$res_diff %>% .[!.$Factors %in% c("(Intercept)", "(phi)"), ] %>% .[.$Significance %in% c("***"), ] %>% .$Taxa %>% unique
# further select the feature with standard taxonomic names; \\d is a regular expression usage, meaning any number from 0 to 9
select_taxa %<>% .[!grepl("\\d", .)]

# use trans_env class to perform correlation between taxa and environmental variables
tmp_transenv_rhizo <- trans_env$new(dataset = tmp_microtable_rhizo, env_cols = 7:19, standardize = TRUE)
# use_data = "other" means provide customized taxa names via other_taxa parameter
tmp_transenv_rhizo$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = select_taxa)
write.csv(tmp_transenv_rhizo$res_cor, file.path(output_dir, "Diff_abund_env_cor_betareg.csv"))

# visualization
g1 <- tmp_transenv_rhizo$plot_cor(cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_cluster_betareg.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)

# When the cluster plot is added, the output is aplot class, not ggplot class. 
# If you want to change the plot like ggplot2 usage, you should access the object like managing a list object.
# access the element
g1[[1]]
# optimize the x axis text name
g1[[1]] <- g1[[1]] + ggplot2::scale_x_discrete(labels = c(NH4 = expression(NH[4]^'+'-N), NO3 = expression(NO[3]^'-'-N)))
g1
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_cluster_betareg_newxaxistext.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)


# also visualize the differential test result
# select factors and taxa that will be visualized in the figure
tmp <- tmp_transdiff_betareg$res_diff %>% .[.$Taxa %in% select_taxa, ] %>% .[!.$Factors %in% c("(Intercept)", "(phi)"), ]
# assign back to object
tmp_transdiff_betareg$res_diff <- tmp
write.csv(tmp_transdiff_betareg$res_diff, file.path(output_dir, "Diff_abund_diff_betareg.csv"))

# when the formula is found, the plot_diff_bar function can automatically employ the heatmap instead of bar plot.
g1 <- tmp_transdiff_betareg$plot_diff_bar(heatmap_cell = "Estimate", heatmap_lab_fill = "Betareg\nCoef", cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_diff_betareg.png"), g1, base_aspect_ratio = 1, dpi = 300, base_height = 8)















