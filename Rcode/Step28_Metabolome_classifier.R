## 
## Machine learning for metabolome data
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)
library(ggplot2)
theme_set(theme_bw())
######################################################
output_dir <- "Output/3.Metabolome/Stage10_Metabolome"
# load data
input_path <- file.path(output_dir, "metab_microtable.RData")
if(! file.exists(input_path)){
	stop("Please first run the script in the step27!")
}
load(input_path)
######################################################
set.seed(123)

tmp_microtable <- clone(metab_microtable)
number_raw <- nrow(tmp_microtable$otu_table)

tmp_microtable$filter_taxa(freq = 0.2)
number_filter <- nrow(tmp_microtable$otu_table)

# generate taxa_abund list for further analysis; use raw abundance
tmp_microtable$cal_abund(rel = FALSE)

######################################################
taxa_level <- "class"

# Cropping treatments: two groups
y_response <- "Cropping"

# create trans_classifier object
t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)
# standardize data; "nzv" means removing features with near zero variance; "center" subtracts the mean; "scale" divides by the standard deviation
t1$cal_preProcess(method = c("center", "scale", "nzv"))

# feature selection or add other customized or manual selected data into the object
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)

number_select_Cropping <- ncol(t1$data_feature)

# generate train and test set
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", y_response, "_split_data_test.csv")))

# add set_trainControl to conveniently use trainControl function to pass customized parameters
t1$set_trainControl()
# train the random forest model
t1$cal_train(method = "rf")

# prediction
t1$cal_predict()
# save the object
save(t1, file = file.path(output_dir, paste0("Classification_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# perform ROC and PR analysis with training data
t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_PR_train.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# calculate feature importance with significance using rfPermute package
t1$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)




# Fertilization treatments: three groups
y_response <- "Fertilization"

t2 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)
# standardize data; "nzv" means removing features with near zero variance
t2$cal_preProcess(method = c("center", "scale", "nzv"))

# feature selection or add other customized or manual selected data into the object
t2$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)

number_select_Fertilization <- ncol(t2$data_feature)

# generate train and test set
t2$cal_split(prop.train = 3/4)
write.csv(t2$data_train, file.path(output_dir, paste0("Classification_", y_response, "_split_data_train.csv")))
write.csv(t2$data_test, file.path(output_dir, paste0("Classification_", y_response, "_split_data_test.csv")))

# add set_trainControl to conveniently use trainControl function to pass customized parameters
t2$set_trainControl()
# train the random forest model
t2$cal_train(method = "rf")

# prediction
t2$cal_predict()
# save the object
save(t2, file = file.path(output_dir, paste0("Classification_", y_response, "_model_rf_train_predict.RData")))

g1 <- t2$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# perform ROC and PR analysis with training data
t2$cal_ROC(input = "train")
write.csv(t2$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t2$res_ROC$res_pr, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_PR_train.csv")))
g1 <- t2$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t2$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# calculate feature importance with significance using rfPermute package
t2$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t2$res_feature_imp, file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t2$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t2$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1, 
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)




# select and save the important/informative features for the following analysis
tmp_sel_features <- c(colnames(t1$data_feature), t2$res_feature_imp %>% .[.$MeanDecreaseAccuracy.pval < 0.01, ] %>% rownames)
tmp_sel_features %<>% unique

save(tmp_sel_features, file = file.path(output_dir, "Classification_rf_select_features.RData"))




# PLS-DA for rhizosphere soil

tmp_microtable_rhizo <- clone(tmp_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

t1 <- trans_beta$new(dataset = tmp_microtable_rhizo, group = "Group")
t1$cal_ordination(method = "PLS-DA", scale_species = TRUE)
t1$plot_ordination(plot_color = "Group")
g1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"), loading_arrow = FALSE) +
	theme(panel.grid=element_blank()) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
	geom_hline(yintercept = 0, linetype = "dashed", color = "grey80")

cowplot::save_plot(file.path(output_dir, "Classification_PLSDA_Rhizo.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
imp <- as.data.frame(sort(t1$res_ordination$model@vipVn, decreasing = TRUE))
# Variable Importance
colnames(imp)[1] <- "Importance"
write.csv(imp, file.path(output_dir, "Classification_PLSDA_Rhizo_vip.csv"))













