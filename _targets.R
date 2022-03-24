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
  packages = c("tidyverse", "data.table", "cutr", "ukbtools", "rbgen", "bigsnpr", "TwoSampleMR",
               "ggplot2", "purrr", "rmeta", "PASWR2", "cowplot", "meta", "strex", "RColorBrewer",
               "forestplot", "R.utils", "MendelianRandomization", "rstatix"),
  error = "workspace",
  memory = "transient",
  storage = "worker",
  garbage_collection = TRUE

)

path_UKBB_imp_data <- paste0(UKBB_dir, "/imp") #too large to track

##################################
### INPUT DATA ####
##################################

input_data <- tar_map(
  values = list(
    custom_names = files_custom_names,
    files = files
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

  ## COUNT RAW DATA

  tar_target(
    data_main_UKBB_raw_count,
    countLines(path_main_UKBB_raw)[[1]]-1 # subtract 1 to account for header
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
  # write_traits_corr = write.csv(trait_corrs, file_out( "output/tables/1.household_correlations.csv"), row.names=F),

  tar_target(
    traits_corr_spearman,
    compute_trait_corr_spearman(data_phesant_directory,data_UKBB_directory,hh_pairs_filter)
  ),

  tar_target(
    PCs_corr,
    compute_pc_corr(sqc_munge, hh_pairs_filter, data_id_sex)
  ),

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
    path_define_cats_filled,
    {
      path_define_cats
      "output/tables/define_Neale_categories_filled.csv" ## THIS FILE IS CREATED MANUALLY
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
      stats1 <- 2
      stats2 <- 2
      while ( (stats1) > 1 | (stats2) > 1){
        stats1 <- length(suppressWarnings(system(paste("squeue -n", "process_Neale"), intern = TRUE)))
        stats2 <- length(suppressWarnings(system(paste("squeue -n", "clump_Neale_IVs"), intern = TRUE)))
        print("Still processing Neale files (extracting IVs and clumping)...")
        Sys.sleep(2)
      }

    }
  ),
  tar_target(
    traits_to_count_IVs,
    {
      run_process_Neale
      pull_traits_to_count_IVs(traits_corr2_filled)
    }
  ),
  tar_target(
    IV_list,
    get_IV_list(traits_corr2_filled,traits_to_count_IVs$Neale_pheno_ID,IV_threshold, Neale_summary_dir), pattern = map(traits_to_count_IVs)
  ),
  tar_target(
    path_IV_list,
    write_IV_list(traits_corr2_filled, Neale_pheno_ID = traits_to_count_IVs$Neale_pheno_ID, IV_list,
                  IV_threshold, "analysis/data_setup/IV_lists/"), pattern = map(traits_to_count_IVs, IV_list),
    format = "file"
  ),
  tar_target(
    count_IVs,
    tibble(Neale_pheno_ID = traits_to_count_IVs$Neale_pheno_ID, num_IVs = length(IV_list)), pattern = map(traits_to_count_IVs, IV_list)
  ),
  tar_target(
    traits_corr3,
    IV_filter(traits_corr2_filled, count_IVs, num_IVs_threshold)
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
    IV_data_summary,
    {
      ## This function gets info on all IVs for male and females, calcs het between and provides a summary line with number of SNPs that pass filter
      path_IV_list
      summarize_IV_data(traits_corr3$to_run, traits_to_calc_het$Neale_pheno_ID, variant_IV_data,
                        Neale_output_dir, IV_threshold)
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

  # filter for non-diet traits only
  tar_target(
    traits_corr5,
    non_diet_filter(traits_corr4$to_run)
  ),

  # filter out left-side traits (highly correlated with right-side)
  # filter out categorical 'Qualifications:', because data is roughly captured with continuous variables
  tar_target(
    traits_corr6,
    non_redundant_filter(traits_corr5)
  ),

  # final list of traits to run in pipeline, just a copy of above target
  tar_target(
    traits_final,
    traits_corr6
  ),

  tar_target(
    path_correlations_final_filter,
    write_final_filter(traits_final, "output/tables/household_correlations.final_filter.csv"),
    format = "file"
  ),

  tar_target(
    outcomes_to_run,
    pull_traits_to_run(traits_final) #this could be changed to traits_corr2_filled
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
    path_rev_filter_dirs,
    create_trait_rev_filt_dirs(outcomes_to_run$Neale_pheno_ID), pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    pheno_data, # pheno data list could expand beyond just those traits with relevant IVs by changing `outcomes_to_run`
    {
      path_phesant
      path_outcome_dirs
      prep_pheno_data(traits_corr2_filled, outcomes_to_run$Neale_pheno_ID,
                      data_sqc, data_fam, data_relatives)
    }, pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    path_pheno_data,
    {
      path_outcome_dirs
      write_pheno_data(pheno_data, outcomes_to_run$Neale_pheno_ID)
    },
    format = "file",
    pattern = map(pheno_data, outcomes_to_run)
  ),

  tar_target(
    exposure_info,
    get_trait_info(traits_final,outcomes_to_run$Neale_pheno_ID,
                   data_Neale_manifest, Neale_summary_dir, Neale_output_dir),
    pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    corr_mat_traits,
    calc_corr_mat_traits(outcomes_to_run, path_pheno_data)
  ),

  tar_target(
    corr_mat_traits_all,
    calc_corr_mat_traits_all(traits_corr2_filled, path_pheno_data, data_sqc, data_fam, data_relatives)
  ),

  tar_target(
    PC_trait_corr,
    calc_pc_trait_corr(outcomes_to_run$Neale_pheno_ID, pheno_data),
    pattern = map(outcomes_to_run, pheno_data), iteration = "list"
  ),

  tar_target(
    PC_traits,
    calc_PC_traits(outcomes_to_run, path_pheno_data)
  ),

  tar_target(
    num_tests_by_PCs,
    calc_num_tests_by_PCs(PC_traits, 0.995)
  ),

  tar_target(
    IV_data_summary_run,
    IV_data_summary[[IV_indices_to_run]],
    pattern = map(IV_indices_to_run), iteration= "list"
  ),

  tar_target(grouping_var, c("time_together_even_bins", "age_even_bins")),

  tar_target(
    summ_stats,
    create_summary_stats(outcomes_to_run$Neale_pheno_ID, exposure_info, IV_data_summary_run),
    pattern = map(outcomes_to_run, exposure_info, IV_data_summary_run),
    iteration = "list"
  ),

  tar_target(
    path_summ_stats,
    {
      create_trait_dirs
      write_summ_stats(outcomes_to_run$Neale_pheno_ID, summ_stats)
    },
    pattern = map(outcomes_to_run, summ_stats), format = "file"
  ),

  tar_target(
    IV_genetic_data, ## loads genetic data for each set of IVs, same IVs are used for both males and females
    load_geno(summ_stats[["male_IV_data"]], data_UKBB_sample, path_UKBB_imp_data), iteration = "list",
    pattern = map(summ_stats)
  ),

  ## raw phenotypic correlation in different bins

  tar_target(
    binned_pheno_corrs,
    {
      path_pheno_data
      calc_binned_pheno_corrs(exposure_info, grouping_var, household_time_munge)
    },
    iteration = "list",
    pattern = map(exposure_info)
  ),

  tar_target(
    binned_pheno_figs,
    {
      create_binned_pheno_figs(binned_pheno_corrs, exposure_info, grouping_var)
    },
    iteration = "list",
    pattern = map(binned_pheno_corrs, exposure_info)
  ),

  ## household MR prep

  tar_target(
    path_household_GWAS,
    {
      path_pheno_data
      path_outcome_dirs
      run_household_GWAS(exposure_info, summ_stats, outcomes_to_run, traits_corr2_filled,
                         IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge)
    },
    pattern = map(exposure_info, summ_stats, IV_genetic_data), format = "file"
  ),

  tar_target(
    household_GWAS,
    read_household_GWAS(path_household_GWAS),
    pattern = map(path_household_GWAS), iteration = "list"
  ),

  tar_target(
    household_harmonised_data,
    harmonise_household_data(exposure_info, summ_stats, outcomes_to_run, gwas_results = household_GWAS, traits_corr2_filled),
    pattern = map(exposure_info, summ_stats, household_GWAS), iteration = "list"
  ),

  tar_target(
    household_harmonised_data_meta, # explore meta-analyzing the SNP/trait associations before the MR
    meta_harmonised_household_data(exposure_info, outcomes_to_run, household_harmonised_data),
    pattern = map(exposure_info, household_harmonised_data), iteration = "list"
  ),

  ## standard MR prep

  tar_target(
    outcome_stats,
    extract_Neale_outcome(outcome_ID = outcomes_to_run$Neale_pheno_ID,
                          both_sexes_file = outcomes_to_run$both_sexes_original_Neale_file, male_file = outcomes_to_run$male_original_Neale_file,
                          female_file = outcomes_to_run$female_original_Neale_file, variant_IV_data),
    pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    path_outcome_stats,
    {
      path_outcome_dirs
      write_outcome_stats(exposure_info, outcome_stats,outcomes_to_run, summ_stats)
    },
    pattern = map(exposure_info, outcome_stats), format = "file"
  ),

  tar_target(
    standard_harmonised_data,
    {
      path_outcome_stats
      harmonise_standard_data(exposure_info, summ_stats, outcomes_to_run, traits_corr2_filled)
    },
    pattern = map(exposure_info, summ_stats), iteration = "list"
  ),

  tar_target(
    standard_harmonised_data_meta,
    meta_harmonised_standard_data(exposure_info, outcomes_to_run, standard_harmonised_data),
    pattern = map(exposure_info, standard_harmonised_data), iteration = "list"
  ),

  tar_target(
    path_outcome_stats_meta,
    {
      path_outcome_dirs
      write_outcome_stats_meta(exposure_info, outcomes_to_run, standard_harmonised_data_meta, summ_stats)
    },
    pattern = map(exposure_info, standard_harmonised_data_meta, summ_stats), format = "file"
  ),


  ## Filter SNPs for evidence of reverse-causality based on same-person two-trait MR (`standard_MR`)

  tar_target(
    standard_harmonised_data_meta_reverse_filter, # filter SNPs for evidence of reverse causation using Ninon's filtering method.
    filter_reverse_SNPs_standard_data(exposure_info, outcomes_to_run, standard_harmonised_data_meta, reverse_MR_threshold),
    pattern = map(exposure_info, standard_harmonised_data_meta), iteration = "list"
  ),

  tar_target(
    num_SNPs_lost_reverse_filter,
    check_num_SNPS_removed_reverse_filter(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter),
    pattern = map(exposure_info, standard_harmonised_data_meta_reverse_filter)
  ),

  tar_target(
    household_harmonised_data_meta_reverse_filter, # using same set of filtered SNPs from `standard_harmonised_data_meta_reverse_filter`, filter meta household data.
    filter_reverse_SNPs_household_data(exposure_info, outcomes_to_run, household_harmonised_data_meta, standard_harmonised_data_meta_reverse_filter),
    pattern = map(exposure_info, household_harmonised_data_meta, standard_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  tar_target(
    standard_harmonised_data_reverse_filter, # using same set of filtered SNPs from `standard_harmonised_data_meta_reverse_filter`, filter sex-specific standard data.
    filter_reverse_SNPs_standard_data_sex_spec(exposure_info, outcomes_to_run, standard_harmonised_data, standard_harmonised_data_meta_reverse_filter),
    pattern = map(exposure_info, standard_harmonised_data, standard_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  tar_target(
    household_harmonised_data_reverse_filter, # using same set of filtered SNPs from `standard_harmonised_data_meta_reverse_filter`, filter sex-specific household data.
    filter_reverse_SNPs_household_data_sex_spec(exposure_info, outcomes_to_run, household_harmonised_data, standard_harmonised_data_meta_reverse_filter),
    pattern = map(exposure_info, household_harmonised_data, standard_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  ## Write filtered and meta-analyzed SNP results

  tar_target(
    path_outcome_stats_filter,
    {
      path_rev_filter_dirs
      write_outcome_stats_filter(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter, standard_harmonised_data_reverse_filter, summ_stats)
    },
    pattern = map(exposure_info, standard_harmonised_data_meta_reverse_filter, standard_harmonised_data_reverse_filter, summ_stats), format = "file"
  ),

  tar_target(
    path_household_GWAS_filter,
    {
      path_rev_filter_dirs
      write_household_GWAS_filter(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter, summ_stats)
    },
    pattern = map(exposure_info, household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter, summ_stats), format = "file"
  ),

  ## Household MR

  ## Household MR is run in the following groups:
  ## (1) binned, sex-specific, (2) binned, joint - MR results meta-analyzed across sexes [from (1)], (3) binned, joint - SNP-effects meta-analyzed across sexes,
  ## (4) non-binned, sex-specific, (5) non-binned, joint - MR results meta-analyed acrsos sexes [from (4)], (6) non-binned, joint - SNP-effect meta-analyzed across sexes.

  ## Notes:
  ## Results 4-6 also have MR figures
  tar_target(
    path_MR_dirs,
    {
      path_outcome_dirs
      create_MR_dirs(outcomes_to_run$Neale_pheno_ID)
    },
    pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    household_MR_binned_sex_specific, # MR results are given binned in full sample and binned by time-together and mean age
    run_binned_household_MR(exposure_info, outcomes_to_run, household_harmonised_data_reverse_filter, grouping_var, MR_method_list = MR_method_list),
    pattern = map(exposure_info, household_harmonised_data_reverse_filter), iteration = "list"
  ),

  tar_target(
    household_MR_binned_MRmeta, # meta-analyze binned results by sex, meta-analyze at MR-level
    meta_binned_household_MR(exposure_info, outcomes_to_run, household_MR_binned_sex_specific),
    pattern = map(exposure_info, household_MR_binned_sex_specific), iteration = "list"
  ),

  tar_target(
    household_MR_binned_SNPmeta, # MR results run in data meta-analyzed at SNP-level. Compare with meta-analyzed MR above to see impact of meta-analyzing MRs vs. meta-analyzing SNPs.
    run_binned_household_MR_SNPmeta(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, grouping_var, MR_method_list = MR_method_list),
    pattern = map(exposure_info, household_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  # tar_target(
  ## No longer run because all effects are standardized first.
  #   household_MR_binned_SNPmeta_std, # MR results run in data meta-analyzed at SNP-level. SNP effects were standardized first.
  #   run_binned_household_MR_SNPmeta_std(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, grouping_var, MR_method_list = MR_method_list),
  #   pattern = map(exposure_info, household_harmonised_data_meta_reverse_filter), iteration = "list"
  # ),

  tar_target(
    couple_MR_vs_trait_corr,
    compare_mr_raw_corr(exposure_info, household_MR_binned_SNPmeta, traits_corr),
    pattern = map(exposure_info, household_MR_binned_SNPmeta)
  ),

  tar_target(
    couple_MR_vs_trait_corr_sig,
    couple_MR_vs_trait_corr %>%
      filter(diff_p < 0.05/num_tests_by_PCs)
    ),

  tar_target(
    couple_MR_vs_trait_corr_sig_neg,
    couple_MR_vs_trait_corr_sig %>% filter(!correlation_larger)
  ),

  tar_target(
    couple_MR_vs_trait_corr_sig_pos,
    couple_MR_vs_trait_corr_sig %>% filter(correlation_larger)
  ),

  tar_target(
    corr_potential_trait_confounders_neg,
    find_potential_trait_confounders_neg(Neale_pheno_ID=couple_MR_vs_trait_corr_sig_neg$exposure_ID, Neale_pheno_ID_corr=couple_MR_vs_trait_corr_sig_neg$couple_r,
                                     standard_MR_summary_SNPmeta, household_MR_summary_SNPmeta, traits_corr, num_tests_by_PCs, corr_mat_traits, corr_trait_threshold),
    pattern = map(couple_MR_vs_trait_corr_sig_neg)
  ),

  tar_target(
    corr_potential_trait_confounders_pos,
    find_potential_trait_confounders_pos(Neale_pheno_ID=couple_MR_vs_trait_corr_sig_pos$exposure_ID, Neale_pheno_ID_corr=couple_MR_vs_trait_corr_sig_pos$couple_r,
                                     standard_MR_summary_SNPmeta, household_MR_summary_SNPmeta, traits_corr, num_tests_by_PCs, corr_mat_traits, corr_trait_threshold),
    pattern = map(couple_MR_vs_trait_corr_sig_pos)
  ),

  tar_target(
    corr_potential_PC_confounders,
    find_potential_PC_confounders(Neale_pheno_ID=couple_MR_vs_trait_corr_sig$exposure_ID, Neale_pheno_ID_corr=couple_MR_vs_trait_corr_sig$couple_r,
                                  PC_trait_corr, PCs_corr, path_pheno_data),
    pattern = map(couple_MR_vs_trait_corr_sig)
  ),

  # decided to scrap the MVMR model with C-sum:

  # tar_target(
  #   potential_trait_confounders_MVMR,
  #   {
  #     path_outcome_stats
  #     path_outcome_stats_filter
  #     path_outcome_stats_meta
  #     run_MVMR_potential_trait_confounders(corr_potential_trait_confounders, household_MR_summary_AM, prune_threshold)
  #   },
  #   pattern = map(corr_potential_trait_confounders), iteration= "list"
  # ),


  # tar_target(
  #   proxyMR_yiyp_adj_SNPmeta,
  #   {
  #     path_household_GWAS_filter
  #     path_outcome_stats_filter
  #     adj_yiyp_xIVs_SNPmeta(exposure_info, household_harmonised_data_meta_reverse_filter, household_MR_summary_BF_sig)
  #   },
  #   map(exposure_info, household_harmonised_data_meta_reverse_filter)
  # ),
  #
  #


  tar_target(
    path_household_MR_binned,
    {
      path_MR_dirs
      write_household_MR(exposure_info, outcomes_to_run, household_MR_binned_SNPmeta)
    },
    pattern = map(exposure_info, household_MR_binned_SNPmeta), format = "file"
  ),

  tar_target(
    household_MR_binned_het_sex_spec, # calc sex-het, and slope and Q-stat among bins, in sex-specific results.
    calc_binned_household_MR_het(exposure_info, outcomes_to_run, household_MR_binned_MRmeta),
    pattern = map(exposure_info, household_MR_binned_MRmeta)
  ),

  tar_target(
    household_MR_binned_het_SNPmeta, # calc slope and Q-stat among bins in result meta-analyzed at SNP-level (joint).
    calc_binned_household_MR_het_SNPmeta(exposure_info, outcomes_to_run, household_MR_binned_SNPmeta),
    pattern = map(exposure_info, household_MR_binned_SNPmeta)
  ),

  ####################################
  ## join the result above together ##
  ####################################

  tar_target(
    household_MR_sex_specific, ## `household_MR_sex_specific` is run in full sample only, not binned by age / time-together categories.
    run_household_MR_comprehensive(exposure_info, outcomes_to_run, household_harmonised_data_reverse_filter, MR_method_list),
    pattern = map(exposure_info, household_harmonised_data_reverse_filter), iteration = "list"
  ),

  tar_target(
    household_MR_SNPmeta, ## `household_MR_SNPmeta` is run in full sample only, not binned by age / time-together categories (meta-analyzed at SNP-level).
    run_household_MR_comprehensive_SNPmeta(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, MR_method_list),
    pattern = map(exposure_info, household_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  tar_target(
    household_MR_SNPmeta_sensitivity, ## `household_MR_SNPmeta_sensitivity` is run in full sample only, not binned by age / time-together categories (meta-analyzed at SNP-level).
    run_household_MR_SNPmeta_sensitivity(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter),
    pattern = map(exposure_info, household_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now, meta-analyze across sexes and calculate heterogeneity statistic between sexes.
    household_MR_summary_MRmeta,
    summarize_household_MR_comprehensive(household_MR_sex_specific, corr_mat_traits),
    pattern = map(household_MR_sex_specific)
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now.
    household_MR_summary_SNPmeta,
    summarize_household_MR_comprehensive_SNPmeta(household_MR_SNPmeta, corr_mat_traits),
    pattern = map(household_MR_SNPmeta)
  ),

  ## Filter Household MR for significance and correlation

  tar_target(
    household_MR_summary_BF_sig,
    find_sig_household_MR_summary(household_MR_summary_SNPmeta, num_tests_by_PCs)
  ),


  tar_target(
    household_MR_summary_corr_filter,
    find_non_corr_household_MR_summary(household_MR_summary_BF_sig, corr_trait_threshold)
  ),



  ## Household MR analyses for only same trait (AM = assortative mating)

  tar_target(
    household_MR_binned_AM_figs,
    create_MR_binned_AM_figs(household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter, household_MR_binned_SNPmeta, household_MR_binned_MRmeta, custom_col),
    pattern = map(household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter, household_MR_binned_SNPmeta, household_MR_binned_MRmeta), iteration = "list"
  ),

  tar_target(
    household_MR_AM_FvsM_fig,
    create_household_MR_AM_FvsM_fig(household_MR_binned_het_sex_spec, custom_col),
  ),

  tar_target(
    ## filter to only MR between same traits, these results meta-analyzed at SNP-level.
    household_MR_summary_AM,
    pull_AM_MRs_household_MR_summary(household_MR_summary_SNPmeta)
  ),

  tar_target(
    ## filter to only MR between same traits, these results meta-analyzed at MR-level.
    ## These results will show the heterogeneity statistics.
    household_MR_summary_AM_MRmeta,
    pull_AM_MRs_household_MR_summary(household_MR_summary_MRmeta)
  ),

  tar_target(
    ## after AM MR results is reduced to only those significant based on `num_tests_by_PCs`, reduce further for testing heterogeneity between sexes and amongst bins to reduce number of tests.
    outcomes_to_run_AM_sig,
    find_AM_sig_exposure_info(household_MR_summary_AM, outcomes_to_run, num_tests_by_PCs)
  ),

  tar_target(
    PC_traits_AM_sig,
    calc_PC_traits(outcomes_to_run_AM_sig, path_pheno_data)
  ),

  tar_target(
    num_tests_by_PCs_AM_sig,
    calc_num_tests_by_PCs(PC_traits_AM_sig, 0.995)
  ),

  ## Household MVMR:  $Y_p \sim X_i + Y_i + X_p$
  ## Decided that this won't work afterall since we don't have additional instruments for X_p, beyond those that are already instruments for X_i. The only way around would be to run a GWAS on X_p, but it's not worth the effort.
  ## If we change our minds the `run_household_MVMR_SNPmeta` needs some minor attention before running.

  # tar_target(
  #   household_MVMR_sex_specific,
  #   {
  #     path_household_GWAS
  #     path_outcome_stats
  #     run_household_MVMR(exposure_info, outcomes_to_run)
  #   },
  #   pattern = map(exposure_info), iteration = "list"
  # ),


  # tar_target(
  #   household_MVMR_SNPmeta,
  #   {
  #     path_household_GWAS_filter
  #     path_outcome_stats_filter
  #     run_household_MVMR_SNPmeta(exposure_info, outcomes_to_run)
  #   },
  #   pattern = map(exposure_info), iteration = "list"
  # ),

  # tar_target(
  #   household_MVMR_summary_MRmeta,
  #   summarize_household_MVMR(household_MVMR_sex_specific, traits_corr2_filled, corr_mat_traits),
  #   pattern = map(household_MVMR_sex_specific)
  # ),

  ## Standard MR

  tar_target(
    standard_MR_sex_specific,
    run_standard_MR_comprehensive(exposure_info, outcomes_to_run, standard_harmonised_data, MR_method_list),
    pattern = map(exposure_info, standard_harmonised_data), iteration = "list"
  ),

  tar_target(
    standard_MR_SNPmeta,
    run_standard_MR_comprehensive_SNPmeta(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter, MR_method_list),
    pattern = map(exposure_info, standard_harmonised_data_meta_reverse_filter), iteration = "list"
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now, meta-analyze across sexes and calculate heterogeneity statistic between sexes.
    standard_MR_summary_MRmeta,
    summarize_standard_MR_comprehensive(standard_MR_sex_specific),
    pattern = map(standard_MR_sex_specific)
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now, meta-analyze across sexes and calculate heterogeneity statistic between sexes.
    standard_MR_summary_SNPmeta,
    summarize_standard_MR_comprehensive_SNPmeta(standard_MR_SNPmeta),
    pattern = map(standard_MR_SNPmeta)
  ),

  tar_target(
    standard_MR_summary_BF_sig,
    find_sig_standard_MR_summary(standard_MR_summary_SNPmeta)
  ),

  ## Compute rho, omega and gamma (without adjusting for $X_i$ SNPs)
  ## [We don't use these results because we only consider those from the adjusted model.]

  # tar_target(
  #   proxyMR_comparison_sex_specific,
  #   run_proxyMR_comparison(exposure_info, household_MR_summary_corr_filter, household_MR_summary_MRmeta, standard_MR_summary_MRmeta, household_MR_summary_AM_MRmeta),
  #   map(exposure_info, household_MR_summary_MRmeta, standard_MR_summary_MRmeta), iteration = "list"
  # ),
  #
  # tar_target(
  #   proxyMR_comparison_SNPmeta,
  #   run_proxyMR_comparison_SNPmeta(exposure_info, household_MR_summary_corr_filter, household_MR_summary_SNPmeta, standard_MR_summary_SNPmeta, household_MR_summary_AM),
  #   map(exposure_info, household_MR_summary_SNPmeta, standard_MR_summary_SNPmeta), iteration = "list"
  # ),
  #
  # tar_target(
  #   proxyMR_paths_summary, ## summarize the different MR paths in each Xi -> Yp
  #   summarize_proxyMR_paths(proxyMR_comparison_SNPmeta)
  # ),
  #
  # tar_target(
  #   proxyMR_comparison_summary, ## would need to reformat the function slightly for the sex-specific results
  #   summarize_proxyMR_comparison(proxyMR_comparison_SNPmeta, traits_corr2_filled)
  # ),
  #
  # tar_target(
  #   proxyyMR_IV_overlap,
  #   {
  #     path_summ_stats
  #     find_proxyMR_IV_overlap(exposure_info, proxyMR_paths_summary)
  #   },
  #   map(exposure_info)
  #
  # ),

  ## Adj $Y_p \sim Y_i$ associations for $X_i$ and reperform proxyMR

  tar_target(
    proxyMR_yiyp_adj_sex_specific,
    {
      path_outcome_stats_filter
      path_household_GWAS_filter
      adj_yiyp_xIVs_sex_specific(exposure_info, household_harmonised_data_reverse_filter, household_MR_summary_corr_filter)
    },
    map(exposure_info, household_harmonised_data_reverse_filter)
  ),

  tar_target(
    proxyMR_yiyp_adj_SNPmeta,
    {
      path_household_GWAS_filter
      path_outcome_stats_filter
      adj_yiyp_xIVs_SNPmeta(exposure_info, household_harmonised_data_meta_reverse_filter, household_MR_summary_corr_filter)
    },
    map(exposure_info, household_harmonised_data_meta_reverse_filter)
  ),

  ## Run the formal comparisons between gamma, rho and omega.
  ## Might need to add target(s) for sex-specific results if reviewers request.

  tar_target(
    proxyMR_comparison_yiyp_adj_SNPmeta,
    run_proxyMR_comparison_adj_yiyp_SNPmeta(exposure_info, household_MR_summary_corr_filter, household_MR_summary_SNPmeta, standard_MR_summary_SNPmeta, household_MR_summary_AM, proxyMR_yiyp_adj_SNPmeta),
    map(exposure_info, household_MR_summary_SNPmeta, standard_MR_summary_SNPmeta), iteration = "list"
  ),

  tar_target(
    proxyMR_MR_paths_summary_yiyp_adj_SNPmeta, ## summarize the different MR paths in each Xi -> Yp
    summarize_proxyMR_paths(proxyMR_comparison_yiyp_adj_SNPmeta)
  ),

  tar_target(
    proxyMR_comparison_summary_yiyp_adj_SNPmeta, ## compare the effects
    summarize_proxyMR_comparison_SNPmeta(proxyMR_comparison_yiyp_adj_SNPmeta, traits_corr2_filled)
  ),

  tar_target(
    proxyMR_comparison_summary_yiyp_adj_SNPmeta_prune, ##  if there is a pair A-B and C-D and if max(corr(A,C)*corr(B,D),corr(A,D)*corr(B,C)) > `corr_trait_threshold` [0.8]
    prune_proxyMR_comparison_SNPmeta(proxyMR_comparison_summary_yiyp_adj_SNPmeta, household_MR_summary_corr_filter, corr_mat_traits, corr_trait_threshold)
  ),

  ## Add z's to proxyMR
  ## Did not end up performing this analysis because so much of the variance was explained by rho and gamma (i.e. direct effects)

  # tar_group_count(
  #   MV_z, ## z's are based on standard_MR: x -> z -> y
  #   find_MV_z(household_MR_summary_BF_sig, standard_MR_summary),
  #   count=350
  # ),
  #
  # tar_target(
  #   MV_z_corr_filter, # filter to only z's that have correlation < 0.8 (play with this a bit).
  #   corr_filter_MV_z(MV_z, corr_mat_traits, z_prune_threshold),
  #   map(MV_z)
  # ),
  #
  # tar_target(
  #   z_summ_stats,
  #   {
  #     path_outcome_stats
  #     pull_z_summ_stats(MV_z_corr_filter)
  #   },
  #   map(MV_z_corr_filter)
  # ),
  #
  # tar_target(
  #   z_summ_stats_pruned,
  #   prune_z_summ_stats(MV_z_corr_filter, z_summ_stats, prune_threshold),
  #   map(MV_z_corr_filter, z_summ_stats)
  # ),

  ## Create some summary figures - input for the below figures needs to be a sex-specific result.

  # tar_target(
  #   proxyMR_prod_comparison_fig_yiyp_adj_all,
  #   # creates a 3x3 grid of results
  #   # gam vs rho, rho vs omega, gam vs. omega in M, F and joint
  #   create_proxy_prod_comparison_fig_all(proxyMR_comparison_summary_yiyp_adj_SNPmeta)
  # ),

  # tar_target(
  #   proxyMR_sex_comparison_fig_yiyp_adj,
  #   create_proxy_sex_comparison_fig(proxyMR_comparison_summary_yiyp_adj_SNPmeta)
  # ),

  tar_target(
    proxyMR_prod_comparison_fig_yiyp_adj,
    create_proxy_prod_comparison_fig(proxyMR_comparison_summary_yiyp_adj_SNPmeta)
  ),

  # Run a GWAS of PCs
  # Based on results of PC correlation decided not to pursue this analysis.
  # If we come back to it, it will need to be rerun in the proper UKB subgroup.

  # tar_target(
  #   PC_gwas_input,
  #   prep_PC_GWAS(data_id_age, data_id_sex, sqc_munge, data_UKBB_sample)
  # ),
  #
  # tar_target(
  #   path_PC_gwas_input,
  #   write_PC_gwas_input(PC_gwas_input),
  #   format = "file"
  # ),
  #
  # tar_target(
  #   path_v2_snps,
  #   create_UKBB_v2_snp_file_list(UKBB_processed_dir),
  #   format = "file"
  # ),
  #
  # tar_target(
  #   file_v2_snps,
  #   path_v2_snps
  # ),
  #
  # tar_target(
  #   bgenie_ukbb_chunks,
  #     make_ukbb_chunks(file_v2_snps, chunk_size=5e4), pattern = map(file_v2_snps)
  # ),
  #
  # tar_target(
  #   bgenie_ukbb_chunks_run,
  #   as_tibble(bgenie_ukbb_chunks)
  # ),
  #
  # tar_target(
  #   path_bgenie_GWAS_dir,
  #   create_bgenie_GWAS_dir()
  # ),
  #
  # tar_target(
  #   path_bgenie_pcs,
  #   launch_bgenie(bgenie_ukbb_chunks_run$chr, phenofile = path_PC_gwas_input, UKBB_dir,
  #                 bgenie_ukbb_chunks_run$chr_char, bgenie_ukbb_chunks_run$start, bgenie_ukbb_chunks_run$end, bgenie_ukbb_chunks_run$chunk_num,
  #                 output_dir = path_bgenie_GWAS_dir, output_prefix = "PC_GWAS"),
  #   pattern = map(bgenie_ukbb_chunks_run), format = "file"
  # ),
  #
  # tar_target(
  #   path_bgenie_pcs_unzip,
  #   unzip_bgenie(path_bgenie_pcs),
  #   pattern = map(path_bgenie_pcs), format = "file"
  # ),

  ## Assessing the impact of some potential / important confounders braodly

  ## We did this by calculating the correlation due to confounding with two approachs (1) using the correlation as input and (2) causal effects as input.

  tar_target(
    corr_impact_by_PCs,
    calc_corr_impact_by_PCs(traits_corr, PCs_corr, PC_trait_corr,
                            path_pheno_data),
    pattern = map(PC_trait_corr)
  ),

  tar_target(
    corr_impact_by_coords,
    calc_corr_impact_by_coords(outcomes_to_run, traits_corr, corr_mat_traits,
                               traits_corr2_filled, path_pheno_data, data_sqc, data_fam, data_relatives)
  ),

  tar_target(
    confounder_traits, # Neale IDs of confounder traits to test. Running this differently than coor/PC analysis because there is no summing across traits.
    c("738", "845", "189_irnt", "20016_irnt", "6160_1", "6147_1", "1180", "129_irnt", "130_irnt")
    # household income, age completed education, townsend deprivation index, Fluid intelligence score
    # Leisure/social activities: Sports club or gym, Reason for glasses/contact lenses: For short-sightedness, i.e. only or mainly for distance viewing such as driving, cinema etc (called 'myopia')
    # Morning/evening person (chronotype)

  ),

  tar_target(
    corr_impact_by_traits,
    calc_corr_impact_by_traits(outcomes_to_run, traits_corr, corr_mat_traits_all, confounder_traits,
                               traits_corr2_filled, path_pheno_data, data_sqc, data_fam, data_relatives),
    pattern = map(confounder_traits)
  ),

  tar_target(
    corr_impact_by_traits_causal,
    calc_corr_impact_by_traits_causal(confounder_traits, standard_MR_summary_SNPmeta,
                                      household_MR_summary_SNPmeta, traits_corr),
    pattern = map(confounder_traits)
  ),


  tar_target(
    corr_impact_by_coords_causal,
    calc_corr_impact_by_coords_causal(corr_impact_by_traits_causal),

  ),

)


