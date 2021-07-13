library(targets)
library(tarchetypes)

source("R/functions.R")
source("code/settings.R")

options(tidyverse.quiet = TRUE)
options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")


tar_option_set(
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(num_cores = 1, account = "sgg",
                                                        ntasks = 1, partition = "sgg",
                                                        log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))
  ),
  packages = c("tidyverse", "data.table", "cutr", "ukbtools", "rbgen", "bigsnpr", "TwoSampleMR", "ggplot2", "purrr"),
  error = "workspace",
  memory = "transient",
  garbage_collection = TRUE

)

path_UKBB_imp_data <- paste0(UKBB_dir, "/imp") #too large to track

##################################
### INPUT DATA ####
##################################

input_data <- tar_map(
  values = list(
    custom_names = c("household_info", "phesant_directory", "relatives", "fam", "sqc", "time_at_address",
                     "time_at_address_raw", "UKBB_directory", "Neale_SGG_dir", "Neale_manifest", "code_process_Neale", "Neale_variants",
                     "UKBB_sample"),
    files = c(paste0(UKBB_dir,"/pheno/ukb6881.csv"), paste0(UKBB_processed_dir,"/PHESANT/","PHESANT_file_directory.txt"),
              paste0(UKBB_dir,"/geno/","ukb1638_rel_s488366.dat"),  paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam"),
              paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt"),
              paste0(UKBB_processed_dir, "/PHESANT/ukb31459/bin1/out_bin1..tsv"), paste0(UKBB_dir, "/pheno/ukb31459.csv"),
              paste0(UKBB_processed_dir,"/UKBB_pheno_directory.csv"),
              Neale_SGG_dir_file_cp, paste0(Neale_output_dir,"/",Neale_manifest_file), "code/process_Neale.sh",
              Neale_variant_file, paste0(UKBB_dir, "/imp/ukb1638_imp_chr1_v2_s487398.sample"))
  ),
  names = custom_names,
  unlist = FALSE,
  tar_target(path, files, format = "file")
)


list(

  input_data,

  tar_target(
    path_household_relationships,
    as.character(data_phesant_directory[which(data_phesant_directory[,2]==household_relationships_field),"File"]),
    format = "file"
  ),
  tar_target(
    path_first_phesant_file,
    as.character(data_phesant_directory[1,"File"]),
    format = "file"
  ),
  tar_files(path_phesant, as.character(unique(data_phesant_directory$File))),

  ##################################
  ### READ DATA ####################
  ##################################

  tar_target(
    data_household_info,
    read.csv(path_household_info, header=T),
  ),
  tar_target(
    data_relatives,
    read.table(path_relatives, header = T),
  ),
  tar_target(
    data_phesant_directory,
    read.table(path_phesant_directory, header=T),
  ),
  tar_target(
    data_household_relationships,
    fread(path_household_relationships, header=T, select=c("userId","sex",household_relationships_field), data.table=F),
  ),
  tar_target(
    data_id_age,
    fread(path_first_phesant_file, header=T, select=c("userId","age"), data.table=F),
  ),
  tar_target(
    data_id_sex,
    fread(path_first_phesant_file, header=T, select=c("userId","sex"), data.table=F),
  ),
  tar_target(
    data_sqc,
    fread(path_sqc, header=F, data.table=F),
  ),
  tar_target(
    data_fam,
    fread(path_fam, header=F, data.table=F),
  ),
  tar_target(
    data_time_at_address,
    fread(path_time_at_address,  select=c("userId",time_at_address_field), data.table=F)
  ),
  tar_target(
    data_time_at_address_raw,
    fread(path_time_at_address_raw, data.table=F)
  ),
  tar_target(
    data_UKBB_directory,
    read.csv(path_UKBB_directory, header=T)
  ),
  tar_target(
    data_Neale_SGG_dir,
    read.csv(path_Neale_SGG_dir, header=T)
  ),
  tar_target(
    data_Neale_manifest,
    read_tsv(path_Neale_manifest, col_types = cols())
  ),
  tar_target(
    data_UKBB_sample,
    fread(path_UKBB_sample, skip=2, header=F,data.table=F)
  ),

  ##################################
  ### MAKE HOUSEHOLD PAIRS FILE ####
  ##################################

  tar_group_count(
    hh_pairs,
    pairs_only(data_household_info),
    count=20
  ),

  tar_target(
    hh_pairs_kin,
    find_kinship(hh_pairs, data_relatives), pattern = map(hh_pairs)
  ),

  tar_target(
    hh_pairs_filter,
    filter_pairs(hh_pairs_kin, data_household_relationships, household_relationships_field),
  ),


  # write_hh_pairs = write.csv(hh_pairs_filter, file_out("analysis/data_setup/household_pairs.csv"),row.names=F, quote=T ),


  ################################
  ######## CREATE PHENO FILE #####
  ################################


  tar_target(
    sqc_munge,
    munge_sqc(data_sqc, data_fam),
  ),

  tar_target(
    model_adjustments,
    add_pcs(hh_pairs_filter, data_id_age, sqc_munge),
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
    household_time,
    calc_time_together(joint_model_adjustments,data_time_at_address,data_time_at_address_raw, time_at_address_raw_field, time_at_address_field)
  ),
  tar_target(
    household_intervals,
    interval_labels(c(seq(time_together_min,time_together_max, by=time_together_interval)))
  ),
  tar_target(
    household_time_munge,
    munge_household_time(household_time, joint_model_adjustments, time_together_min,time_together_max, time_together_interval,
                                                 age_min, age_max, age_interval, num_household_bins)
  ),

  # write_household_time = write.csv(household_time_munge, file_out("analysis/data_setup/household_time.csv"), row.names=F, quote=T )

  #####################################
  #### FITLER TRAITS ##################
  #####################################

  tar_target(
    traits_corr,
    compute_trait_corr(data_phesant_directory,data_UKBB_directory,hh_pairs_filter)
  ),

  tar_target(
    PCs_corr,
    compute_pc_corr(sqc_munge, hh_pairs_filter)
  ),
  # write_traits_corr = write.csv(trait_corrs, file_out( "output/tables/1.household_correlations.csv"), row.names=F),

  tar_target(
    Neale_SGG_dir_filt,
    SGG_link_with_Neale(data_Neale_SGG_dir)
  ),
  tar_target(
    traits_corr2,
    filter_by_corr(traits_corr,Neale_SGG_dir_filt,household_correlation_threshold)
  ),
  tar_target(
    Neale_to_process,
    organize_Neale(traits_corr2)
  ),
  tar_target(
    path_define_cats,
    write_define_cats(Neale_to_process),
    format = "file"
  ),
  tar_target(
    path_download_list,
    write_download_list(Neale_to_process),
    format = "file"
  ),
  tar_target(
    path_define_cats_filled, ## THIS IS CREATED MANUALLY
    {
      path_define_cats
      "output/tables/define_Neale_categories_filled.csv"
    },
    format = "file"
  ),
  tar_target(
    data_define_cats_filled,
    read.csv(path_define_cats_filled, check.names=F)
  ),
  tar_target(
    traits_corr2_filled,
    download_Neale(data_define_cats_filled,Neale_to_process$download_rest,traits_corr2,
                   path_Neale_manifest, Neale_output_dir)
  ),
  tar_target(
    run_process_Neale,
    {
      traits_corr2_filled
      processx::run(command = "sbatch", c(path_code_process_Neale))
    }
  ),

  tar_target(
    traits_corr2_update,
    {
      run_process_Neale
      stats1 <- 2
      stats2 <- 2
      while ( (stats1) > 1 | (stats2) > 1){
        stats1 <- length(suppressWarnings(system(paste("squeue -n", "process_Neale"), intern = TRUE)))
        stats2 <- length(suppressWarnings(system(paste("squeue -n", "clump_Neale_IVs"), intern = TRUE)))
        print("Still running...")
        Sys.sleep(2)
      }

      update_download_info(traits_corr2_filled, data_Neale_SGG_dir)
    }, deployment = "main"
  ),
  tar_target(
    traits_to_count_IVs,
    pull_traits_to_count_IVs(traits_corr2_update)
  ),
  tar_target(
    IV_list,
    get_IV_list(traits_corr2_update,traits_to_count_IVs$Neale_pheno_ID, data_Neale_manifest,IV_threshold, Neale_output_dir, Neale_summary_dir), pattern = map(traits_to_count_IVs)
  ),
  tar_target(
    path_IV_list,
    write_IV_list(traits_corr2_update, Neale_pheno_ID = traits_to_count_IVs$Neale_pheno_ID, IV_list,
                  IV_threshold, "analysis/data_setup/IV_lists/"), pattern = map(traits_to_count_IVs, IV_list),
    format = "file"
  ),
  tar_target(
    count_IVs,
    tibble(Neale_pheno_ID = traits_to_count_IVs$Neale_pheno_ID, num_IVs = length(IV_list)), pattern = map(traits_to_count_IVs, IV_list)
  ),
  tar_target(
    traits_corr3,
    IV_filter(traits_corr2_update, count_IVs, num_IVs_threshold)
  ),
  # write_traits_corr2 = write.csv(traits_corr3$non_filtered,file_out("output/tables/2.household_correlations.corr_filter.csv"), row.names=F),

  tar_target(
    variant_IV_data,
    reduce_Neale_variant_data(path_Neale_variants, IV_list), deployment = "main"
  ),
  tar_target(
    traits_to_calc_het,
    pull_traits_to_count_IVs(traits_corr3$to_run)
  ),

  tar_target(
    IV_data_summary, # this target used to be called: sex_het_summary
                     # could try running it as it used to be and make sure you get the same result
    {
      ## This function gets info on all IVs for male and females, calcs het between and provides a summary line with number of SNPs that pass filter
      path_IV_list
      summarize_IV_data(traits_corr3$to_run, traits_to_calc_het$Neale_pheno_ID, variant_IV_data,
                        data_Neale_manifest, Neale_summary_dir, Neale_output_dir, IV_threshold)
    }, pattern = map(traits_to_calc_het), iteration = "list"
  ),
  tar_target(
    path_IV_info,
    write_IV_info(IV_data_summary, traits_to_calc_het$Neale_pheno_ID),
    format = "file",
    pattern = map(IV_data_summary, traits_to_calc_het)
  ),

  tar_target(
    path_sex_het,
    write_sex_het(IV_data_summary, traits_to_calc_het$Neale_pheno_ID),
    format = "file",
    pattern = map(IV_data_summary, traits_to_calc_het)
  ),

  tar_target(
    traits_corr4,
    sex_het_filter(traits_corr3$to_run, IV_data_summary, num_IVs_threshold)
  ),
  # filter for continuous traits only
  tar_target(
    traits_corr5,
    continuous_filter(traits_corr4$to_run)
  ),

  # final list of traits to run in pipeline, just a copy of above target
  tar_target(
    traits_final,
    traits_corr5
  ),

  tar_target(
    path_correlations_final_filter,
    write_final_filter(traits_final, "output/tables/household_correlations.final_filter.csv"),
    format = "file"
  ),

  tar_target(
    exposures_to_run,
    pull_traits_to_run(traits_final)
  ),

  tar_target(
    outcomes_to_run,
    pull_traits_to_run(traits_final) #this could be changed to traits_corr2_update, add Neale file name column
  ),

  tar_target(
    IV_indices_to_run,
    pull_IV_indices_to_run(traits_final, traits_to_calc_het), iteration = "list"
  ),

  ## data prep
  tar_target(
    path_outcome_dirs,
    create_trait_dirs(outcomes_to_run$Neale_pheno_ID), pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    pheno_data, # pheno data list could expand beyond just those traits with relevant IVs by changing `outcomes_to_run`
    {
      path_phesant
      path_outcome_dirs
      prep_pheno_data(traits_corr2_update, outcomes_to_run$Neale_pheno_ID,
                data_sqc, data_fam, data_relatives)
    }, pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    path_pheno_data,
    write_pheno_data(pheno_data, outcomes_to_run$Neale_pheno_ID),
    format = "file",
    pattern = map(pheno_data, outcomes_to_run)
  ),

  tar_target(
    exposure_info,
    get_trait_info(traits_final, exposures_to_run$Neale_pheno_ID,
                   data_Neale_manifest, Neale_summary_dir, Neale_output_dir),
    pattern = map(exposures_to_run), iteration = "list"
  ),

  tar_target(
    IV_data_summary_run,
    IV_data_summary[[IV_indices_to_run]],
    pattern = map(IV_indices_to_run), iteration= "list"
  ),

  tar_target(grouping_var, c("time_together_even_bins", "age_even_bins")),

  tar_target(
    summ_stats,
    create_summary_stats(exposures_to_run$Neale_pheno_ID, exposure_info, IV_data_summary_run),
    pattern = map(exposures_to_run, exposure_info, IV_data_summary_run),
    iteration = "list"
  ),

  tar_target(
    IV_genetic_data, ## loads genetic data for each set of IVs, same IVs are used for both males and females
    load_geno(summ_stats[["male_IV_data"]], data_UKBB_sample, path_UKBB_imp_data), iteration = "list",
    pattern = map(summ_stats)
  ),

  tar_target(
    path_household_GWAS,
    {
      path_pheno_data  ### map over all phenos
      path_outcome_dirs
      household_GWAS_all_outcomes(exposure_info, summ_stats, outcomes_to_run, traits_corr2_update,
                         IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge)
    },
    pattern = map(exposure_info, summ_stats, IV_genetic_data),  format = "file"

  ),

  ## MR

  tar_target(
    path_MR_dirs,
    {
      path_outcome_dirs
      create_MR_dirs(outcomes_to_run$Neale_pheno_ID)
    },
    pattern = map(outcomes_to_run), iteration = "list"
  ),

  ## Binned results -> could change the name
  tar_target(
    household_MR_binned,
    household_MR_all_outcomes(exposure_info, summ_stats, outcomes_to_run, gwas_files = path_household_GWAS,
                       traits_corr2_update, grouping_var, MR_method_list),
    pattern = map(exposure_info, summ_stats, IV_genetic_data, path_household_GWAS), iteration = "list"
  ),

  tar_target(
    path_household_MR_binned,
    {
      path_MR_dirs
      write_household_MR(exposure_info, outcomes_to_run, household_MR_binned)
    },

    pattern = map(exposure_info, household_MR_binned), format = "file"
  ),

  tar_target(
    household_harmonised_data,
    harmonise_household_data_all_outcomes(exposure_info, summ_stats, outcomes_to_run, gwas_files = path_household_GWAS, traits_corr2_update),
    pattern = map(exposure_info, summ_stats, IV_genetic_data, path_household_GWAS), iteration = "list" #could remove IV_genetic data
  ),

  ## RUN in full sample only, not binned
  tar_target(
    household_MR_exhaustive,
    household_MR_complete_all_outcomes(exposure_info, household_harmonised_data, outcomes_to_run, MR_method_list),
    pattern = map(exposure_info, household_harmonised_data), iteration = "list"
  ),

  ##summarize into one table, ignore leave-1-out analyses for now
  ## add column for outcome_ID==exposure_ID
  ## outcome_description and exposure description

  tar_target(
    household_MR_exhaustive_summary,
    household_MR_complete_summary(household_MR_exhaustive),
    pattern = map(household_MR_exhaustive)
  ),

  tar_target(
    all_IVs,
    tibble(exposure = exposures_to_run$Neale_pheno_ID, rsid = summ_stats[[1]]$rsid), # it doesn't matter if you take male or female summary stats
    pattern = map(exposures_to_run, summ_stats)

  ),

  ## THIS IS THE EXACT SAME AS VARIANT DATA
  # tar_target(
  #   IV_variant_data,
  #   extract_relevant_variant_rows(Neale_variant_file, snp_list = all_IVs$rsid)
  # ),

  tar_target(
    outcome_stats,
    {
      extract_Neale_outcome(Neale_file, IV_variant_data, outcomes_to_run, traits_corr2_update)
    },

    pattern = map(outcomes_to_run, summ_stats)

  ),

  tar_target(
    standard_MR,
    {
      standard_MR(exposures_to_run, outcomes_to_run)
    },

    pattern = map(exposures_to_run)

  ),

  tar_target(
    PC_gwas_input,
    prep_PC_GWAS(data_id_age, data_id_sex, sqc_munge, data_UKBB_sample)
  ),

  tar_target(
    path_PC_gwas_input,
    write_PC_GWAS_input(PC_gwas_input),
    format = "file"
  ),

  tar_target(
    path_v2_snp_list,
    create_UKBB_v2_snp_list(UKBB_processed_dir),
    format = "file"
  ),

  tar_target(
    bgenie_ukbb_chunks,
    make_ukbb_chunks(path_v2_snp_list, chunk_size=1e6), pattern = map(path_v2_snp_list)
  )






#   traits_corr4 = sex_het_filter(traits_corr3$to_run, sex_het_summary, traits_to_calc_het, !!num_IVs_threshold),
#   write_traits_corr3 = write.csv(traits_corr4$non_filtered, file_out("output/tables/3.household_correlations.numIVs_filter.csv"), row.names=F),
#   write_traits_corr4 = write.csv(traits_corr4$to_run, file_out("output/tables/4.household_correlations.sexhet_filter.csv"), row.names=F),
#   write_traits_corr5 = write.csv(traits_corr5, file_out("output/tables/6.household_correlations.nonbinary_filter.csv"), row.names=F)




)


