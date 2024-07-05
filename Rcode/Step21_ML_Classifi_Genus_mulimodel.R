

######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/1.Amplicon/StageⅦ_Machine_learning"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "Output/1.Amplicon/StageⅡ_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in StageⅡ !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable)
# calculate relative abundance at each taxonomic level
tmp_microtable$cal_abund(rel = TRUE)


######################################################
# Genus level: relative abundance
taxa_level <- "Genus"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = "Cropping", x.predictors = taxa_level)
# standardize data; "nzv" means removing features with near zero variance
t1$cal_preProcess(method = c("center", "scale", "nzv"))

# feature selection or add other customized or manual selected data into the object
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)

# analyze different models; 'rf': random forest; 'svmRadial': support vector machine; 'rda': regularized discriminant analysis model
# If the function cal_split is not performed, there is no data_train in the object. Thus, the function cal_caretList will use all the samples for the training.
t1$cal_caretList(methodList = c('rf', 'svmRadial', 'rda'))
# analyze a set of resampling results and reshape the metric values
t1$cal_caretList_resamples()
write.csv(t1$res_caretList_resamples$values, file.path(output_dir, paste0("Classification_", taxa_level, "_multiple_models_metrics_raw.csv")))
write.csv(t1$res_caretList_resamples_reshaped, file.path(output_dir, paste0("Classification_", taxa_level, "_multiple_models_metrics_reshaped.csv")))
# visualize the performance of models
g1 <- t1$plot_caretList_resamples() + theme_bw()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_multiple_models_metrics.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)








