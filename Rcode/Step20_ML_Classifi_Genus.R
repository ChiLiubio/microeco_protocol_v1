## 
## Classification at Genus level
## 


######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/1.Amplicon/Stage7_Machine_learning"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable)
# calculate relative abundance (i.e. total sum scaling) at each taxonomic level
tmp_microtable$cal_abund(rel = TRUE)


######################################################
# Genus level: relative abundance
taxa_level <- "Genus"

# Cropping treatments: two groups
y_response <- "Cropping"

# create trans_classifier object
t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)


# generate train and test set
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

# ##### Optional: preprocessing. "nzv" means removing features with near zero variance; "center" subtracts the mean; "scale" divides by the standard deviation
# t1$cal_preProcess(method = c("center", "scale", "nzv"))

# Optional: feature selection
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_traindata.csv")))

# add set_trainControl to conveniently use trainControl function to pass customized parameters
t1$set_trainControl()
# train the random forest model
t1$cal_train(method = "rf")

# prediction
t1$cal_predict()
# save the object
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# perform ROC and PR analysis with training data
t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.csv")))
# ROC curve
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
# PR curve
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# perform ROC and PR analysis with prediction data
t1$cal_ROC(input = "pred")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.csv")))
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.csv")))
# Figure 5a
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
# Figure 5b
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# calculate feature importance using varImp function in caret package
t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# calculate feature importance with significance using rfPermute package
t1$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
# Figure 5c
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)




# Fertilization treatments: three groups
y_response <- "Fertilization"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

# generate train and test set
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

# Optional
# feature selection or add other customized or manual selected data into the object
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_feature, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_featuredata.csv")))

# add set_trainControl to conveniently use trainControl function to pass customized parameters
t1$set_trainControl()
# train the random forest model
t1$cal_train(method = "rf")

# prediction
t1$cal_predict()
# save the object
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# perform ROC and PR analysis with training data
t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# perform ROC and PR analysis with prediction data
t1$cal_ROC(input = "pred")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.csv")))
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# calculate feature importance using varImp function in caret package
t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# calculate feature importance with significance using rfPermute package
t1$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)







