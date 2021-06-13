#####################################################################################
## Author: Jenny Sjaarda
## Project: proxyMR
##
## Time-stamp: <[settings.R] by JS 2021-06-01 10:42:03 CEST>
##
## Description:
##
##
## History:
##
#####################################################################################

### Set data and project directories

project_dir <- "/data/sgg2/jenny/projects/proxyMR/"
SGG_generic <- "/data/sgg2/jenny/SGG_generic/"


## register clustermq and future plans
#options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")
#future::plan(batchtools_slurm, template = "slurm_batchtools.tmpl")

#source(paste0(SGG_generic,"/scripts/settings.r"))

traits <- c()
MR_method_list <- c("mr_wald_ratio","mr_ivw","mr_ivw_fe","mr_egger_regression","mr_weighted_median")

### define variables
IV_threshold <- 5e-08
prune_threshold <- 0.001
GRS_thresholds <- c(0.1,0.01,0.001)
num_household_bins <- 5

num_IVs_threshold <- 5
household_correlation_threshold <-0.1
irnt=TRUE


#phesant_directory <- read.table(paste0(UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt"), header=T)
#phesant_file_list <- unique(phesant_directory$File)

# Data-Field 6141
# Description:	How are people in household related to participant
# limit to only individuals who responded YES to husband, wife, or partner

relatedness_field <- "6141_1"
#relations_file <- as.character(phesant_directory[which(phesant_directory[,2]==relatedness_field),"File"])

# name of first phesant file
#first_phesant_file <- as.character(phesant_directory[1,"File"])


#time_at_address_file <- "/data/sgg2/jenny/data/UKBB_processed/PHESANT/ukb31459/bin1/out_bin1..tsv"
#time_at_address_raw_file <- "/data/sgg2/jenny/data/UKBB_raw/pheno/ukb31459.csv"
time_at_address_field <- "699"
time_at_address_raw_field <- "699-0.0"

# -10 represents "Less than a year"
# -1 represents "Do not know"
# -3 represents "Prefer not to answer"
time_together_min <- 0
time_together_max <- 50
time_together_interval <- 10

age_min <- 40
age_max <- 70
age_interval <- 10
## functions to fit:

Neale_output_path <- "/data/sgg3/data/neale_files"
Neale_manifest_file_name <- "UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.tsv"
