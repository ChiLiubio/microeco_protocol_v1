## 
## Classification at ASV level
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
# load a normalized microtable object
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_CLR.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable_CLR)

# filter ASVs based on the customized standards
# As an example, we first perform low-abundance filtering based on the unnormalized raw data, and then conduct screening on the CLR-normalized data
load("Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData")
tmp <- clone(amplicon_16S_microtable)
tmp$filter_taxa(rel_abund = 0.001)
tmp_microtable$tax_table %<>% .[rownames(.) %in% rownames(tmp$tax_table), ]
tmp_microtable$tidy_dataset()

# add ASV to tax_table to generate ASV abundance in the taxa_abund list, which is the input data in the trans_classifier class
tmp_microtable$add_rownames2taxonomy(use_name = "ASV")
# regenerate taxa_abund list with original abundance (i.e., CLR-normalized data); now there is a table called 'ASV' in the taxa_abund list
tmp_microtable$cal_abund(rel = FALSE)


######################################################
# ASV level; Cropping treatments
taxa_level <- "ASV"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = "Cropping", x.predictors = taxa_level)

# Optional
# standardize data; "nzv" means removing features with near zero variance
# t1$cal_preProcess(method = c("center", "scale", "nzv"))

# generate training and testing data set
t1$cal_split(prop.train = 3/4)
# save generated data
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_split_data_test.csv")))

# Optional
# feature selection or add other customized operation to select features
# t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
# save selected feature data to the directory
# write.csv(t1$data_feature, file.path(output_dir, paste0("Classification_", taxa_level, "_featuresel_boruta_featuredata.csv")))


# add set_trainControl to conveniently use trainControl function to pass customized parameters
t1$set_trainControl()
# train the random forest model; For more models, please refer to: https://topepo.github.io/caret/available-models.html
t1$cal_train(method = "rf")

# prediction
t1$cal_predict()
# save the trans_classifier object
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_train_predict.RData")))

# confusion matrix
g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# perform ROC and PR analysis with training data
t1$cal_ROC(input = "train")
# save the data for ROC curve of training
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_ROC_train.csv")))
# save the data for PR curve
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_PR_train.csv")))
# ROC curve
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
# PR curve
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# perform ROC and PR analysis with prediction data
t1$cal_ROC(input = "pred")
# save the data for ROC curve of testing
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_ROC_pred.csv")))
# save the data for PR curve of testing
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_PR_pred.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# calculate feature importance using varImp function in caret package
t1$cal_feature_imp(scale = TRUE)
# In the result, 'Overall' represents MeanDecreaseGini.
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_feature_importance.csv")))
# feature importance plot
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_feature_importance.png")), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# calculate feature importance with significance using rfPermute package
t1$cal_feature_imp(rf_feature_sig = TRUE)
# MeanDecreaseAccuracy.pval represents the p value of MeanDecreaseAccuracy based on the permutations; Similarly, MeanDecreaseGini.pval represents the p value of MeanDecreaseGini
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_feature_importance_withsig_gini.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_model_rf_feature_importance_withsig_accuracy.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)







