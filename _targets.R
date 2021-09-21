library(targets)
library(tarchetypes)

source("R/functions.R")
source("code/settings.R")

options(tidyverse.quiet = TRUE)
options(clustermq.scheduler = "slurm", clustermq.template = "slurm_clustermq.tmpl")


tar_option_set(
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(num_cores = 1, account = "sgg",
                                                        ntasks = 4, partition = "sgg",
                                                        log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))
  ),
  packages = c("tidyverse", "data.table", "cutr", "ukbtools", "rbgen", "bigsnpr", "TwoSampleMR",
               "ggplot2", "purrr", "rmeta", "PASWR2", "cowplot", "meta", "strex", "RColorBrewer",
               "forestplot", "R.utils", "MendelianRandomization"),
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

  tar_target(
    PCs_corr,
    compute_pc_corr(sqc_munge, hh_pairs_filter, data_id_sex)
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
    pull_traits_to_count_IVs(traits_corr2_filled)
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
    get_trait_info(traits_final, exposures_to_run$Neale_pheno_ID,
                   data_Neale_manifest, Neale_summary_dir, Neale_output_dir),
    pattern = map(exposures_to_run), iteration = "list"
  ),

  tar_target(
    PC_traits,
    calc_PC_traits(exposure_info, data_sqc, data_fam, data_relatives)
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
    create_summary_stats(exposures_to_run$Neale_pheno_ID, exposure_info, IV_data_summary_run),
    pattern = map(exposures_to_run, exposure_info, IV_data_summary_run),
    iteration = "list"
  ),

  tar_target(
    path_summ_stats,
    {
      create_trait_dirs
      write_summ_stats(exposures_to_run$Neale_pheno_ID, summ_stats)
    },
    pattern = map(exposures_to_run, summ_stats), format = "file"
  ),

  tar_target(
    IV_genetic_data, ## loads genetic data for each set of IVs, same IVs are used for both males and females
    load_geno(summ_stats[["male_IV_data"]], data_UKBB_sample, path_UKBB_imp_data), iteration = "list",
    pattern = map(summ_stats)
  ),

  tar_target(
    path_household_GWAS,
    {
      path_pheno_data
      path_outcome_dirs # delete all GWAS files and this target if you need to rerun below because it is not saved as `format = "file"`
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

  ## MR

  tar_target(
    path_MR_dirs,
    {
      path_outcome_dirs
      create_MR_dirs(outcomes_to_run$Neale_pheno_ID)
    },
    pattern = map(outcomes_to_run), iteration = "list"
  ),

  tar_target(
    household_MR_binned, # MR results are given binned in full sample and binned by time-together and mean age
    run_binned_household_MR(exposure_info, outcomes_to_run, household_harmonised_data, grouping_var, MR_method_list = MR_method_list),
    pattern = map(exposure_info, household_harmonised_data), iteration = "list"
  ),

  tar_target(
    household_MR_binned_meta, # meta-analyze binned results by sex
    meta_binned_household_MR(exposure_info, outcomes_to_run, household_MR_binned),
    pattern = map(exposure_info, household_MR_binned), iteration = "list"
  ),

  tar_target(
    household_MR_binned_het, # calc sex-het, and slope and Q-stat among bins, sex-specific and binned
    calc_binned_household_MR_het(exposure_info, outcomes_to_run, household_MR_binned_meta),
    pattern = map(exposure_info, household_MR_binned_meta)
  ),

  tar_target(
    household_MR_binned_AM_figs,
    create_MR_binned_AM_figs(household_harmonised_data, household_MR_binned_meta, custom_col),
    pattern = map(household_harmonised_data, household_MR_binned_meta), iteration = "list"
  ),

  tar_target(
    household_MR_AM_FvsM_fig,
    create_household_MR_AM_FvsM_fig(household_MR_binned_het, custom_col),
  ),


  tar_target(
    path_household_MR_binned,
    {
      path_MR_dirs
      write_household_MR(exposure_info, outcomes_to_run, household_MR_binned_meta)
    },
    pattern = map(exposure_info, household_MR_binned_meta), format = "file"
  ),


  tar_target(
    household_MR, ## `household_MR` is run in full sample only, not binned by age / time-together categories.
    run_household_MR_comprehensive(exposure_info, outcomes_to_run, household_harmonised_data, MR_method_list),
    pattern = map(exposure_info, household_harmonised_data), iteration = "list"
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now, meta-analyze across sexes and calculate heterogeneity statistic between sexes.
    household_MR_summary,
    summarize_household_MR_comprehensive(household_MR),
    pattern = map(household_MR)
  ),

  tar_target(
    household_MR_summary_BF_sig,
    find_sig_household_MR_summary(household_MR_summary)
  ),

  tar_target(
    ## filter to only MR between same traits.
    household_MR_summary_AM,
    pull_AM_MRs_household_MR_summary(household_MR_summary)
  ),

  tar_target(
    ## after AM MR results is reduced to only those significant based on `num_tests_by_PCs`, reduce further for testing heterogeneity between sexes and amongst bins to reduce number of tests.
    exposure_info_AM_sig,
    find_AM_sig_exposure_info(household_MR_summary_AM, exposure_info, num_tests_by_PCs)
  ),

  tar_target(
    PC_traits_AM_sig,
    calc_PC_traits(exposure_info_AM_sig, data_sqc, data_fam, data_relatives)
  ),

  tar_target(
    num_tests_by_PCs_AM_sig,
    calc_num_tests_by_PCs(PC_traits_AM_sig, 0.995)
  ),

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
      write_outcome_stats(outcomes_to_run$Neale_pheno_ID, outcome_stats, exposures_to_run, summ_stats)
    },
    pattern = map(outcomes_to_run, outcome_stats), format = "file"
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
    standard_MR,
    run_standard_MR_comprehensive(exposure_info, outcomes_to_run, standard_harmonised_data, MR_method_list),
    pattern = map(exposure_info, standard_harmonised_data), iteration = "list"
  ),

  tar_target(
    ## summarize into one table, ignore leave-1-out analyses for now, meta-analyze across sexes and calculate heterogeneity statistic between sexes.
    standard_MR_summary,
    summarize_standard_MR_comprehensive(standard_MR),
    pattern = map(standard_MR)
  ),

  tar_target(
    standard_MR_summary_BF_sig,
    find_sig_standard_MR_summary(standard_MR_summary)
  ),

  tar_target(
    proxyMR_comparison,
    run_proxyMR_comparison(exposure_info, household_MR_summary_BF_sig, household_MR_summary, standard_MR_summary, household_MR_summary_AM),
    map(exposure_info, household_MR_summary, standard_MR_summary), iteration = "list"
  ),

  tar_target(
    proxyMR_yiyp_adj,
    {
      path_household_GWAS
      adj_yiyp_xIVs(exposure_info, household_harmonised_data, household_MR_summary_BF_sig)
    },
    map(exposure_info, household_harmonised_data)
  ),


  tar_target(
    proxyMR_comparison_yiyp_adj,
    run_proxyMR_comparison(exposure_info, household_MR_summary_BF_sig, household_MR_summary, standard_MR_summary, household_MR_summary_AM, proxyMR_yiyp_adj),
    map(exposure_info, household_MR_summary, standard_MR_summary), iteration = "list"
  ),

  tar_group_count(
    MV_z, ## z's are based on standard_MR: x -> z -> y
    find_MV_z(household_MR_summary_BF_sig, standard_MR_summary),
    count=350
  ),

  tar_target(
    z_summ_stats,
    {
      path_outcome_stats
      pull_z_summ_stats(MV_z)
    },
    map(MV_z)
  ),

  tar_target(
    z_summ_stats_pruned,
    {
      path_outcome_stats
      prune_z_summ_stats(MV_z, z_summ_stats, prune_threshold)
    },
    map(MV_z, z_summ_stats)
  ),

  tar_target(
    proxyMR_MR_paths_summary, ## summarize the different MR paths in each Xi -> Yp
    summarize_proxyMR_paths(proxyMR_comparison)
  ),

  tar_target(
    proxyyMR_IV_overlap,
    {
      path_summ_stats
      find_proxyMR_IV_overlap(proxyMR_MR_paths_summary)
    }

  ),

  tar_target(
    proxyMR_comparison_summary, ## used to be `proxyMR_figure_data`
    summarize_proxyMR_comparison(proxyMR_comparison, traits_corr2_filled)
  ),


  tar_target(
    proxyMR_prod_comparison_fig,
    create_proxy_prod_comparison_fig(proxyMR_comparison_summary)
  ),

  tar_target(
    proxyMR_sex_comparison_fig,
    create_proxy_sex_comparison_fig(proxyMR_comparison_summary)
  ),

  tar_target(
    PC_gwas_input,
    prep_PC_GWAS(data_id_age, data_id_sex, sqc_munge, data_UKBB_sample)
  ),

  tar_target(
    path_PC_gwas_input,
    write_PC_gwas_input(PC_gwas_input),
    format = "file"
  ),

  tar_target(
    path_v2_snps,
    create_UKBB_v2_snp_file_list(UKBB_processed_dir),
    format = "file"
  ),

  tar_target(
    file_v2_snps,
    path_v2_snps
  ),

  tar_target(
    bgenie_ukbb_chunks,
      make_ukbb_chunks(file_v2_snps, chunk_size=5e4), pattern = map(file_v2_snps)
  ),

  tar_target(
    bgenie_ukbb_chunks_run,
    as_tibble(bgenie_ukbb_chunks)
  ),

  tar_target(
    path_bgenie_GWAS_dir,
    create_bgenie_GWAS_dir()
  ),

  tar_target(
    path_bgenie_pcs,
    launch_bgenie(bgenie_ukbb_chunks_run$chr, phenofile = path_PC_gwas_input, UKBB_dir,
                  bgenie_ukbb_chunks_run$chr_char, bgenie_ukbb_chunks_run$start, bgenie_ukbb_chunks_run$end, bgenie_ukbb_chunks_run$chunk_num,
                  output_dir = path_bgenie_GWAS_dir, output_prefix = "PC_GWAS"),
    pattern = map(bgenie_ukbb_chunks_run), format = "file"
  ),

  tar_target(
    path_bgenie_pcs_unzip,
    unzip_bgenie(path_bgenie_pcs),
    pattern = map(path_bgenie_pcs), format = "file"
  )





#   traits_corr4 = sex_het_filter(traits_corr3$to_run, sex_het_summary, traits_to_calc_het, !!num_IVs_threshold),
#   write_traits_corr3 = write.csv(traits_corr4$non_filtered, file_out("output/tables/3.household_correlations.numIVs_filter.csv"), row.names=F),
#   write_traits_corr4 = write.csv(traits_corr4$to_run, file_out("output/tables/4.household_correlations.sexhet_filter.csv"), row.names=F),
#   write_traits_corr5 = write.csv(traits_corr5, file_out("output/tables/6.household_correlations.nonbinary_filter.csv"), row.names=F)




)


