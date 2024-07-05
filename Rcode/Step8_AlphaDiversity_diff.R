




# Compare two cropping treatments across compartments
measure <- "Shannon"

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "Cropping", by_group = "Compartment")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_bycompart_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_bycompart_Shannon_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_alpha(measure = measure, use_boxplot = FALSE, add_line = TRUE, line_type = 2)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_bycompart_Shannon_errorbar_line.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


# Compare three fertilization treatments across compartments
measure <- "Shannon"

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "Fertilization", by_group = "Compartment")

tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure, y_increase = 0.05)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_Shannon_boxplot_wilcox.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

tmp$cal_diff(method = "anova", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_anova.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_Shannon_boxplot_anova.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

tmp$cal_diff(method = "KW_dunn", measure = measure, KW_dunn_letter = TRUE)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Fertilization_bycompart_KW_dunn.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Fertilization_bycompart_Shannon_boxplot_KW_dunn.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


# two-way anova
tmp <- trans_alpha$new(dataset = tmp_microtable, by_group = "Compartment")

# for single measure
measure <- "Shannon"
tmp$cal_diff(method = "anova", measure = measure, formula = "Cropping*Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova_Shannon.csv"))

## for all measures
tmp <- trans_alpha$new(dataset = tmp_microtable, by_group = "Compartment")
tmp$cal_diff(method = "anova", formula = "Cropping*Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova_allmeasures.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_bycompart_twowayanova_allmeasures_heatmap.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 7)


# linear mixed-effects model

tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "lme", formula = "Cropping+Fertilization+(1|Compartment)")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lme.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lme.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


# linear regression for rhizosphere soil
tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

tmp <- trans_alpha$new(dataset = tmp_microtable_rhizo)
tmp$cal_diff(method = "lm", formula = "Cropping+Fertilization")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lm_Rhizo.csv"))
g1 <- tmp$plot_alpha(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef")
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Cropping_Fertilization_lm_Rhizo.png"), g1, base_aspect_ratio = 1, dpi = 300, base_height = 7)


# one-way anova for rhizosphere soil

measure <- "Shannon"
tmp <- trans_alpha$new(dataset = tmp_microtable_rhizo, group = "Group")
tmp$cal_diff(method = "anova")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_Group_Rhizo_anova.csv"))
g1 <- tmp$plot_alpha(measure = measure, add_sig_text_size = 5)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_Group_Rhizo_anova.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 5)







