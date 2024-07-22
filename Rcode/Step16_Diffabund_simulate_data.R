## 
## Simulate communities for customized benchmarking
## 


######################################################
# load packages
library(SparseDOSSA2)
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################

tmp_microtable_rhizo <- clone(amplicon_16S_microtable)
tmp_microtable_rhizo$sample_table %<>% .[.$Compartment == "Rhizosphere", ]
tmp_microtable_rhizo$tidy_dataset()

# delete the genera without clear taxonomic information and those with low abundance to reduce running time
tmp_microtable_rhizo$tax_table %<>% .[!grepl("\\d+", .$Genus), ]
tmp_microtable_rhizo$tax_table %<>% .[.$Genus != "g__", ]
tmp_microtable_rhizo$tidy_dataset()
tmp_microtable_rhizo$filter_taxa(rel_abund = 0.0002)

#########################################
# fit parameters with SparseDOSSA2 package
fitted_param <- fit_SparseDOSSA2(data = as.matrix(tmp_microtable_rhizo$otu_table))
save(fitted_param, file = file.path(output_dir, "tmp_microtable_rhizo_SparseDOSSA2_fitted_param.RData"))


#########################################
# create a directory and generate simulated data
dir_simdata <- file.path(output_dir, "1.sim_data")
dir.create(dir_simdata)

# simulate single binary metadata; sample number: 30; feature number: 600; significant feature ratio: 0.1
n_sample = 30
n_feature = 600
n_sigfeature_ratio = 0.1
n_sigfeature <- n_sigfeature_ratio * n_feature

# simulation times for each parameter set
n_sim = 10

# load fitted parameters
load(file.path(output_dir, "tmp_microtable_rhizo_SparseDOSSA2_fitted_param.RData"))

for(i in n_sample){
	for(j in n_sigfeature){
		group_sample_num <- i/2
		# binary group
		Group <- c(rep(1, group_sample_num), rep(0, group_sample_num))
		metadata_matrix <- as.matrix(data.frame(Group))
		# fixed effect size: 2
		spike_metadata <- data.frame(metadata_datum = 1, feature_spiked = paste0("Feature", 1:j), associated_property = "abundance", effect_size = 2)

		for(k in seq_len(n_sim)){
			set.seed(k)
			sim1 <- SparseDOSSA2(template = fitted_param, n_sample = i, n_feature = n_feature, new_features = TRUE, 
				spike_metadata = spike_metadata, metadata_matrix = metadata_matrix)

			tmp_sample <- sim1$spike_metadata$metadata_matrix %>% as.data.frame
			tmp_abund <- sim1$simulated_data %>% as.data.frame
			rownames(tmp_sample) <- colnames(tmp_abund)
			tmp_tax <- data.frame(ASV = rownames(tmp_abund))
			rownames(tmp_tax) <- tmp_tax[, 1]
			# prepare microtable object
			tmp_mtobj <- microtable$new(sample_table = tmp_sample, otu_table = tmp_abund, tax_table = tmp_tax)
			tmp_mtobj$tidy_dataset()
			tmp_mtobj$raw_sim <- sim1
			# save each microtable object to RData
			save(tmp_mtobj, file = file.path(dir_simdata, paste0("microtable_nsample_", i, "_nfeature_", n_feature, "_singlenf_", j, "_run_", k, ".RData")), compress = TRUE)
		}
	}
}







