


######################################################
# load packages
library(SparseDOSSA2)
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/StageⅥ_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/StageⅡ_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in StageⅡ !")
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
# fit parameters
fitted_param <- fit_SparseDOSSA2(data = as.matrix(tmp_microtable_rhizo$otu_table))
save(fitted_param, file = file.path(output_dir, "tmp_microtable_rhizo_SparseDOSSA2_fitted_param.RData"))


#########################################
# create a directory and generate simulated data
dir_simdata <- file.path(output_dir, "1.sim_generate_SparseDOSSA2_data")
dir.create(dir_simdata)
# simulate single binary metadata; sample number: 30; feature number: 600; significant feature ratio: 0.1
n_sample = 30
n_feature = 600
n_sigfeature_ratio = c(0.1)
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

			s1 <- sim1$spike_metadata$metadata_matrix %>% as.data.frame
			a1 <- sim1$simulated_data %>% as.data.frame
			rownames(s1) <- colnames(a1)
			t1 <- data.frame(ASV = rownames(a1))
			rownames(t1) <- t1[, 1]
			# prepare microtable object
			d1 <- microtable$new(sample_table = s1, otu_table = a1, tax_table = t1)
			d1$tidy_dataset()
			d1$raw_sim <- sim1

			save(d1, file = file.path(dir_simdata, paste0("microtable_nsample_", i, "_nfeature_", n_feature, "_singlenf_", j, "_run_", k, ".RData")), compress = TRUE)
		}
	}
}







