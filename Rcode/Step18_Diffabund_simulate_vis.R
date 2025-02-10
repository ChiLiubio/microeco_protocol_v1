## 
## Visualization for differential test results of simulated data
##


######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory to save differential test results from each method
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"

output_dir_test <- file.path(output_dir, "2.sim_test")
if(! dir.exists(output_dir)){
	stop("Please first run the last step!")
}



#########################################
library(ggplot2)
library(ggpubr)

# obtain all the method according to the directory names
all_methods <- list.files(output_dir_test)
res_table <- data.frame()

for(i in all_methods){
	# first read data into list, then convert it to data.frame object
	tmp_list <- list()
	tmp_path <- file.path(output_dir_test, i)
	tmp_files <- list.files(tmp_path)
	for(j in tmp_files){
		load(file.path(tmp_path, j))
		tmp_list[[j]] <- c(i, gsub(".RData$", "", j), res[1], res[2])
	}
	tmp <- do.call(rbind, tmp_list)
	res_table <- rbind(res_table, unname(tmp))
}

colnames(res_table) <- c("method", "param", "power", "fdr")

# generate parameters from the file names
res_table <- tidyr::separate_wider_delim(res_table, cols = "param", names = c("param", "iter"), delim = "_run_", cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param", names = c("param", "signum"), delim = "_singlenf_", cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param", names = c("param", "feature_num"), delim = "_nfeature_", cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param", names = c("param", "nsample"), delim = "_nsample_", cols_remove = TRUE)

# try to convert character type to numeric for each column
res_table %<>% dropallfactors(unfac2num = TRUE)



# a customized theme
tmp_theme <- theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 14))

p1 <- res_table %>%
	ggboxplot(x = "method", y = "power", color = "method", palette = "Dark2", add = "mean", xlab = "", ylab = "Power", size = 0.6, width = 0.6)
p1 <- ggpar(p1, legend = "none", x.text.angle = 30) + tmp_theme

p2 <- res_table %>%
	ggboxplot(x = "method", y = "fdr", color = "method", palette = "Dark2", add = "mean", xlab = "", ylab = "FDR", size = 0.6, width = 0.6)
p2 <- ggpar(p2, legend = "none", x.text.angle = 30) + tmp_theme

# merge two plots into one
# Figure 4c and 4d
g1 <- ggarrange(p1, p2, ncol = 1)
cowplot::save_plot(file.path(output_dir, "diff_methods_simulation_powerfdr.png"), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 8)











