

#########################################################################################
## Correlation between alpha diversity and environmental variables
#########################################################################################

# select bulk soil
tmp_microtable_bulksoil <- clone(tmp_microtable)
tmp_microtable_bulksoil$sample_table %<>% .[.$Compartment == "Bulk soil", ]
tmp_microtable_bulksoil$tidy_dataset()

# create a trans_env object for the following correlation analysis
# env_cols parameter to select the columns in sample_table of microtable object
tmp <- trans_env$new(dataset = tmp_microtable_bulksoil, env_cols = 9:21)

tmp$cal_cor(add_abund_table = tmp_microtable_bulksoil$alpha_diversity)
write.csv(tmp$res_cor, file.path(output_dir, "AlphaDiv_Env_bulksoil_correlation.csv"))

g1 <- tmp$plot_cor()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Env_bulksoil_correlation.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)








