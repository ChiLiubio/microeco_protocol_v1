

######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
######################################################
# check output directory
output_dir <- "./Output/1.Amplicon/StageⅢ_coremicrobiome"
if(! dir.exists(output_dir)){
	stop("Please first run the scrip in Step5 !")
}

######################################################

# generate a total core-taxa microtable object according to names
tmp <- clone(tmp_microtable)
tmp$otu_table %<>% .[rownames(.) %in% c(S$taxa_names(), RS$taxa_names(), R$taxa_names()), ]
tmp$tidy_dataset()
# merge samples into compartments
tmp_merge <- tmp$merge_samples(use_group = "Compartment")

# analyze intersection set with trans_venn class
trans_vennobj <- trans_venn$new(tmp_merge, ratio = "numratio", name_joint = "-")
# remove the sets with 0
trans_vennobj$data_summary %<>% .[.$Counts != 0, ]
# output data_summary to file
write.csv(trans_vennobj$data_summary, file.path(output_dir, "Coretaxa_intersection.csv"))

# bar plot for intersections
g1 <- trans_vennobj$plot_bar(sort_samples = FALSE, left_background_fill = "white")
cowplot::save_plot(file.path(output_dir, "Coretaxa_intersection.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 7)


# transform intersection data to microtable object for other visualization
trans_vennobj_mt <- trans_vennobj$trans_comm(use_frequency = TRUE)
trans_vennobj_mt$cal_abund()
write.csv(trans_vennobj_mt$taxa_abund$Genus, file.path(output_dir, "Coretaxa_intersection_Genusratio.csv"))

# Ordered elements in x-axis according to the number
tmp_x <- trans_vennobj$data_summary %>% tibble::rownames_to_column(.) %>% .[order(.$Counts, decreasing = TRUE), ] %>% .$rowname

trans_vennobj_mt_abund <- trans_abund$new(dataset = trans_vennobj_mt, taxrank = "Genus", ntaxa = 10)
g2 <- trans_vennobj_mt_abund$plot_bar(bar_full = TRUE, legend_text_italic = T, xtext_angle = 40, order_x = tmp_x) + ylab("Ratio (%)") + theme(plot.margin = unit(c(1, 0, 0, 2), "cm"))

cowplot::save_plot(file.path(output_dir, "Coretaxa_intersection_Genusratio.png"), g2, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)









