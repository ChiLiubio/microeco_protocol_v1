## 
## Multiple normalization methods demonstration
## 


###########################
# load packages
library(microeco)
library(magrittr)
###########################
# output directory
output_dir <- "./Output/1.Amplicon/Stage2_amplicon_microtable"
# load data
input_path <- file.path(output_dir, "amplicon_16S_microtable.RData")
# first check whether saved data path exists
if(! file.exists(input_path)){
	stop("Please first run the script in the last step !")
}
load(input_path)
###########################
# fix random seed to ensure the reproducibility of the analysis
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable)

# use trans_norm class to perform all the normalization methods
tmp <- trans_norm$new(dataset = tmp_microtable)

# Rarefaction
# the minimum of sample sequences is over 20000
tmp_microtable$sample_sums() %>% range
amplicon_16S_microtable_rarefy <- tmp$norm(method = "rarefy", sample.size = 20000)
# save each normalized microtable object to output directory
save(amplicon_16S_microtable_rarefy, file = file.path(output_dir, "amplicon_16S_microtable_rarefy.RData"), compress = TRUE)

# CLR: Centered log-ratio normalization <ISBN:978-0-412-28060-3> <doi: 10.3389/fmicb.2017.02224>
amplicon_16S_microtable_CLR <- tmp$norm(method = "CLR")
save(amplicon_16S_microtable_CLR, file = file.path(output_dir, "amplicon_16S_microtable_CLR.RData"), compress = TRUE)

# rclr: Robust centered log-ratio normalization <doi: doi:10.1128/msystems.00016-19>
amplicon_16S_microtable_RCLR <- tmp$norm(method = "rclr")
save(amplicon_16S_microtable_RCLR, file = file.path(output_dir, "amplicon_16S_microtable_RCLR.RData"), compress = TRUE)

# GMPR: Geometric mean of pairwise ratios <doi: 10.7717/peerj.4600>
amplicon_16S_microtable_GMPR <- tmp$norm(method = "GMPR")
save(amplicon_16S_microtable_GMPR, file = file.path(output_dir, "amplicon_16S_microtable_GMPR.RData"), compress = TRUE)

# CSS: Cumulative sum scaling normalization based on the metagenomeSeq package <doi:10.1038/nmeth.2658>
amplicon_16S_microtable_CSS <- tmp$norm(method = "CSS")
save(amplicon_16S_microtable_CSS, file = file.path(output_dir, "amplicon_16S_microtable_CSS.RData"), compress = TRUE)

# TMM: Trimmed mean of M-values method based on the normLibSizes function of edgeR package <doi: 10.1186/gb-2010-11-3-r25>
amplicon_16S_microtable_TMM <- tmp$norm(method = "TMM")
save(amplicon_16S_microtable_TMM, file = file.path(output_dir, "amplicon_16S_microtable_TMM.RData"), compress = TRUE)

# RLE: Relative log expression.
amplicon_16S_microtable_RLE <- tmp$norm(method = "RLE")
save(amplicon_16S_microtable_RLE, file = file.path(output_dir, "amplicon_16S_microtable_RLE.RData"), compress = TRUE)

# TSS: Total sum scaling. Abundance is divided by the sequencing depth.
amplicon_16S_microtable_TSS <- tmp$norm(method = "TSS")
save(amplicon_16S_microtable_TSS, file = file.path(output_dir, "amplicon_16S_microtable_TSS.RData"), compress = TRUE)

# DESeq2: Median ratio of counts relative to geometric mean per feature based on the DESeq function of DESeq2 package <doi: 10.1186/s13059-014-0550-8>
# Either group or formula parameter should be provided
formula <- "Compartment+Cropping+Fertilization"
amplicon_16S_microtable_DESeq2 <- tmp$norm(method = "DESeq2", group = formula)
save(amplicon_16S_microtable_DESeq2, file = file.path(output_dir, "amplicon_16S_microtable_DESeq2.RData"), compress = TRUE)

# Wrench <doi: 10.1186/s12864-018-5160-5>
# condition parameter is necessary to be provided
# use condition = "Compartment" as the example; it is suggested to be provided according to the analysis context
amplicon_16S_microtable_Wrench <- tmp$norm(method = "Wrench", condition = "Compartment")
save(amplicon_16S_microtable_Wrench, file = file.path(output_dir, "amplicon_16S_microtable_Wrench.RData"), compress = TRUE)






