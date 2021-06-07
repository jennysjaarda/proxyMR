library(targets)
library(tarchetypes)
library(data.table)
source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "tidyverse", "data.table"))

UKBB_dir <- "/data/sgg2/jenny/data/UKBB_raw/"

household_relationships_field <- "6141_1"
#time_at_address_file <- "/data/sgg2/jenny/data/UKBB_processed/PHESANT/ukb31459/bin1/out_bin1..tsv"
#time_at_address_raw_file <- "/data/sgg2/jenny/data/UKBB_raw/pheno/ukb31459.csv"
time_at_address_field <- "699"
time_at_address_raw_field <- "699-0.0"

source("code/settings.R")

list(

  ##################################
  ### MAKE HOUSEHOLD PAIRS FILE ####
  ##################################
  tar_target(
    household_info_file,
    paste0(UKBB_dir,"/pheno/ukb6881.csv"),
    format = "file"
  ),
  tar_target(
    household_info,
    read.table(household_info_file, header=T),
  ),
  tar_target(
    relatives_file,
    paste0(UKBB_dir,"/geno/","ukb1638_rel_s488366.dat", header=T),
    format = "file"
  ),
  tar_target(
    relatives,
    read.table(relatives_file),
  ),
  tar_target(
    phesant_directory_file,
    paste0(UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt"),
    format = "file"
  ),
  tar_target(
    phesant_directory,
    read.table(phesant_directory_file, header=T),
  ),
  tar_target(
    hh_pairs,
    pairs_only(household_info),
  ),
  tar_target(
    hh_pairs_kin,
    find_kinship(hh_pairs, relatives),
  ),
  tar_target(
    household_relationships_file,
    as.character(phesant_directory[which(phesant_directory[,2]==household_relationships_field),"File"]),
    format = "file"
  ),
  tar_target(
    household_relationships,
    fread(household_relationships_file, header=T, select=c("userId","sex",household_relationships_field), data.table=F),
  ),
  tar_target(
    hh_pairs_filter,
    filter_pairs(hh_pairs_kin, household_relationships, household_relationships_field),
  ),

  # write_hh_pairs = write.csv(hh_pairs_filter, file_out("analysis/data_setup/household_pairs.csv"),row.names=F, quote=T ),


  ################################
  ######## CREATE PHENO FILE #####
  ################################

  tar_target(
    first_phesant_file,
    as.character(phesant_directory[1,"File"]),
    format = "file"
  ),
  tar_target(
    pheno,
    fread(first_phesant_file, header=T, select=c("userId","age"), data.table=F),
  ),
  tar_target(
    sqc_file,
    paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt"),
    format = "file"
  ),
  tar_target(
    sqc,
    fread(sqc_file, header=F, data.table=F),
  ),
  tar_target(
    fam_file,
    paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam"),
    format = "file"
  ),
  tar_target(
    fam,
    fread(fam_file, header=F, data.table=F),
  ),
  tar_target(
    sqc_munge,
    munge_sqc(sqc, fam),
  ),

  tar_target(
    model_adjustments,
    add_pcs(hh_pairs_filter, pheno, sqc_munge),
  ),
  tar_target(
    joint_model_adjustments,
    merge(model_adjustments[[1]], model_adjustments[[2]], by=c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSE_ID","kinship","HOUSEHOLD_MEMBER1_sex","HOUSEHOLD_MEMBER2_sex")),
  ),




  # write_houshold1 = write.csv(model_adjustments[[1]], file_out("analysis/data_setup/HOUSEHOLD_MEMBER1_pheno_model_adjustments.csv"),row.names=F, quote=T ),
  # write_houshold2 = write.csv(model_adjustments[[2]], file_out("analysis/data_setup/HOUSEHOLD_MEMBER2_pheno_model_adjustments.csv"),row.names=F, quote=T ),
  # write_joint = write.csv(joint_model_adjustments, file_out("analysis/data_setup/JOINT_pheno_model_adjustments.csv"), row.names=F, quote=T ),

  #####################################
  #### ADD TIME AT HOUSEHOLD ##########
  #####################################

  tar_target(
    time_at_address_file,
    "/data/sgg2/jenny/data/UKBB_processed/PHESANT/ukb31459/bin1/out_bin1..tsv",
    format = "file"
  ),
  tar_target(
    time_at_address_raw_file,
    "/data/sgg2/jenny/data/UKBB_raw/pheno/ukb31459.csv",
    format = "file"
  ),
  tar_target(
    time_at_address,
    fread(time_at_address_file,  select=c("userId",time_at_address_field), data.table=F)
  ),
  tar_target(
    time_at_address_raw,
    fread(time_at_address_raw_file, data.table=F)
  ),
  tar_target(
    household_time,
    calc_time_together(joint_model_adjustments,time_at_address,time_at_address_raw)
  ),
  tar_target(
    household_intervals,
    interval_labels(c(seq(time_together_min,time_together_max, by=time_together_interval)))
  ),
  tar_target(
    household_time_munge,
    munge_household_time(household_time, joint_model_adjustments, time_together_min,time_together_max, time_together_interval,
                                                 age_min, age_max, age_interval, num_household_bins)
  )

  # write_household_time = write.csv(household_time_munge, file_out("analysis/data_setup/household_time.csv"), row.names=F, quote=T )

  #####################################
  #### FITLER TRAITS ##################
  #####################################

)


