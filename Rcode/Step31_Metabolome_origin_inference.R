## 
## Metabolite origin inference from V2.1.0
## 


######################################################
# load packages
library(microeco)
library(magrittr)
library(readxl)

# check the microeco package version again
if(packageVersion("microeco") < '2.1.0'){stop("Minimum version of microeco package should be 2.1.0! Current version is ", packageVersion("microeco"), " ! Please reinstall it!")}

if(!require("stringdist")){install.packages("stringdist")}

######################################################
output_dir <- "Output/3.Metabolome/Stage10_Metabolome"
# load metabolome microtable object
input_path <- file.path(output_dir, "metab_microtable.RData")
if(!file.exists(input_path)){
	stop("Please first run the script in step27 !")
}
load(input_path)

# load amplicon 16S microtable object
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)

######################################################


tmp <- trans_metab$new(metab = metab_microtable, microb = amplicon_16S_microtable)

# please download the database from https://zenodo.org/records/18618912
# uncompress the zip file
# Fuzzy matching of metabolite names against those in the database
tmp$cal_match(database_path = "./metorigindb_split_202602")
# View(test$res_match)

# origin inference
# The parameter match_col works by determining metabolite assignment based on the provided "names" (row names) or other columns (e.g., "KEGG_ID"). 
# For "names", the function first checks whether res_match (from the last step) exists; if it does, metabolite identification is performed according to the distance threshold.
tmp$cal_origin(database_path = "./metorigindb_split_202602", match_col = c("names", "HMDB_ID", "KEGG_ID"))
# tmp$res_origin_rawtable and tmp$res_origin_list

# create the metabolites-bacteria network according to the result from the last step
res_network <- tmp$cal_origin_network()

# use trans_network class to visualize the result
test <- trans_network$new()
test$res_network <- res_network
test$cal_module()
test$get_node_table()
head(test$res_node_table)

# test$plot_network(method = "networkD3", node_color = "type")












