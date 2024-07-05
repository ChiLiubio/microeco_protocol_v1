

######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory to save differential test results from each method
output_dir <- "./Output/1.Amplicon/StageⅥ_Diff_abund"
if(! dir.exists(output_dir)){
	stop("Please first run the last step!")
}

dir_simdata <- file.path(output_dir, "1.sim_generate_SparseDOSSA2_data")

output_dir_test <- file.path(output_dir, "2.sim_test")
dir.create(output_dir_test, recursive = TRUE)

#########################################
# differential test on the simulated data

# create a function to calculate power (sensitivity) and fdr
powerfdr <- function(H, test){
	tp <- sum(H != 0 & test != 0)
	fp <- sum(H == 0 & test != 0)
	fn <- sum(H != 0 & test == 0)
	power1 <- tp/(tp + fn)
	fdr1 <- fp/(tp + fp)
	c(power1, fdr1)
}


all_input <- list.files(dir_simdata)
group <- "Group"
taxa_level <- "ASV"

for(method in c("metagenomeSeq", "linda", "ALDEx2_kw", "ancombc2", "DESeq2", "edgeR", "GMPR", "Wrench")){
	dir_difftest_method <- file.path(output_dir_test, method)
	dir.create(dir_difftest_method)

	for(i in all_input){
		# load microtable object data
		load(file.path(dir_simdata, i))
		# make sure 'Group' is a categorical variable
		d1$sample_table$Group %<>% as.character
		
		if(method %in% c("GMPR", "Wrench")){
			# normalize data with trans_norm class
			tmp <- trans_norm$new(d1)
			d1 <- tmp$norm(method = method, condition = group)
			# generate taxa_abund list with a table named ASV; rel = FALSE generates absolute abundance
			d1$cal_abund(rel = FALSE)
			t2 <- trans_diff$new(dataset = d1, method = "wilcox", group = group, taxa_level = taxa_level)
		}else{
			t2 <- trans_diff$new(dataset = d1, method = method, group = group, taxa_level = taxa_level)
		}
		
		# extract significant features
		tmp <- t2$res_diff %>% .[grepl("\\*", .$Significance), ] %>% .[, "Taxa"]
		rej <- as.numeric(gsub("Feature", "", tmp))
		
		# extract parameters from saved file name; regular expression (\\d+) specify the extracted number; \\1 means extracting the thing within first paired brackets
		n_feature <- gsub(".*nfeature_(\\d+)_.*", "\\1", i) %>% as.numeric
		H <- test <- rep(0, n_feature)
		j <- gsub(".*singlenf_(\\d+)_.*", "\\1", i) %>% as.numeric
		H[1:j] <- 1
		test[rej] <- 1
		
		# calculate power and fdr
		res <- powerfdr(H, test)
		save(res, file = paste0(dir_difftest_method, "/", i), compress = TRUE)
	}
}




