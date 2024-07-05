


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

# a continuous variable in tmp_microtable$sample_table
y_response <- "pH"

# create trans_classifier object
t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)
# standardize data; "nzv" means removing features with near zero variance
t1$cal_preProcess(method = c("center", "scale", "nzv"))
write.csv(t1$data_feature, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_preProcess_featuredata.csv")))

# feature selection or add other customized or manual selected data into the object
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_feature, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_featuresel_boruta_featuredata.csv")))

# generate train and test set
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_split_data_test.csv")))

# add set_trainControl to conveniently use trainControl function to pass customized parameters
t1$set_trainControl()
# train the random forest model
t1$cal_train(method = "rf")

# prediction
t1$cal_predict()

t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("IncNodePurity")
cowplot::save_plot(file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# save the object
save(t1, file = file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))





# analyze different models; lm: Linear Regression
t1$cal_caretList(methodList = c('rf', 'lm'))
# analyze a set of resampling results and reshape the metric values
t1$cal_caretList_resamples()
# visualize the performance of models
g1 <- t1$plot_caretList_resamples() + theme_bw()
cowplot::save_plot(file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_multiple_models_metrics.png")), g1, 
	base_aspect_ratio = 1.1, dpi = 300, base_height = 5)







