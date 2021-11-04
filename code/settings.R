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

library(tidyverse)

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
z_prune_threshold <- 0.8
reverse_MR_threshold <- 0.001

num_IVs_threshold <- 5
household_correlation_threshold <-0.1
irnt=TRUE


# Data-Field 6141
# Description:	How are people in household related to participant
# limit to only individuals who responded YES to husband, wife, or partner

relatedness_field <- "6141_1"

# name of first phesant file

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

## directories
UKBB_dir <- "/data/sgg2/jenny/data/UKBB_raw"
UKBB_processed_dir <- "/data/sgg2/jenny/data/UKBB_processed"
Neale_summary_dir <- "/data/sgg2/jenny/data/Neale_UKBB_GWAS"
Neale_output_dir <- "/data/sgg3/data/neale_files"

## files

Neale_SGG_dir_file_cp <- "data/Neale_SGG_directory_15_07_2021.csv"
PHESANT_file_dir_cp <- paste0("data/PHESANT_file_directory_05_10_2021.txt")
UKBB_pheno_dir_cp <- paste0("data/UKBB_pheno_directory_05_10_2021.csv")

Neale_manifest_file <- "UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.tsv"
Neale_variant_file <- paste0(Neale_summary_dir, "/variants.tsv")


## UKBB fields
household_relationships_field <- "6141_1"
time_at_address_field <- "699"
time_at_address_raw_field <- "699-0.0"

## Input DATA
files_custom_names <- c("household_info", "phesant_directory", "relatives", "fam", "sqc", "time_at_address",
                 "time_at_address_raw", "UKBB_directory", "Neale_SGG_dir", "Neale_manifest", "code_process_Neale", "Neale_variants",
                 "UKBB_sample", "main_UKBB_raw")

files <- c(paste0(UKBB_dir,"/pheno/ukb6881.csv"), PHESANT_file_dir_cp,
          paste0(UKBB_dir,"/geno/","ukb1638_rel_s488366.dat"),  paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam"),
          paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt"),
          paste0(UKBB_processed_dir, "/PHESANT/ukb31459/bin1/out_bin1..tsv"), paste0(UKBB_dir, "/pheno/ukb31459.csv"),
          UKBB_pheno_dir_cp,
          Neale_SGG_dir_file_cp, paste0(Neale_output_dir,"/", Neale_manifest_file), "code/process_Neale.sh",
          Neale_variant_file, paste0(UKBB_dir, "/imp/ukb1638_imp_chr1_v2_s487398.sample"),
          paste0(UKBB_dir,"/pheno/ukb21067.csv"))

files_ref <- tibble(custom_names = files_custom_names, files = files)

custom_col <- c("#a6cee3",
                "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                "#6a3d9a", "#ffff99", "#b15928")

## UKBB main file for showing how data was filtered (not actually used in pipeline)

main_UKBB_raw_file <- paste0(UKBB_dir,"/pheno/ukb21067.csv")


