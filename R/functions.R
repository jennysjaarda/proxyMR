

download_neale_files  <-  function(phenotype_ids,
                                   irnt = TRUE,
                                   max_regex = 200,
                                   chunk_size = 10,
                                   max_downloads = 500,
                                   output_path = '/data/sgg3/data/neale_files',
                                   reference_file = list.files(path = output_path,
                                                               pattern = 'UKBB.*Manifest.*[.]tsv')[2],
                                   sex = c('both_sexes', 'female', 'male'),
                                   category = 'unsorted') {
  if (length(phenotype_ids) > max_regex) {
    n_downloads  =  download_neale_files(
      phenotype_ids[1:max_regex],
      irnt = irnt,
      max_regex = max_regex,
      chunk_size = chunk_size,
      max_downloads = max_downloads,
      output_path = output_path,
      reference_file = reference_file,
      sex = sex,
      category = category
    )

    download_neale_files(
      phenotype_ids[(max_regex + 1):length(phenotype_ids)],
      irnt = irnt,
      max_regex = max_regex,
      chunk_size = chunk_size,
      max_downloads = max_downloads - n_downloads,
      output_path = output_path,
      reference_file = reference_file,
      sex = sex,
      category = category
    )
  }
  if (length(sex) > 1)
    sex = paste(sex, collapse = '|') %>%
      paste0('(', ., ')')

  phenotype_ids  =  paste0('^', phenotype_ids, ifelse(irnt, '(_irnt|)$', '(_raw|)$')) %>%
    paste(collapse = '|')

  existing_files  =  list.files(output_path,
                                recursive = TRUE,
                                pattern = '[.]gz') %>%
    str_match('[^/]+$') %>%
    c

  to_download  =  read_tsv(reference_file) %>%
    filter(str_detect(.$'Phenotype Code', phenotype_ids)) %>%
    dplyr::select(wget = 'wget command') %>%
    mutate(wget = str_replace(wget, '[.]bgz$', '.gz')) %>%
    filter(str_detect(.$wget, paste0('-O .*[.]', sex, '([.]v2)?[.]tsv[.]gz$'))) %>%
    transmute(
      url = str_match(wget, '^wget (http.*) -O')[, 2],
      folder = paste(
        output_path,
        str_match(
          wget,
          '-O .*[.](both_sexes|female|male)([.]v2)?[.]tsv[.]gz$'
        )[, 2],
        category,
        sep = '/'
      ),
      filename = str_match(wget, '-O (.*[.]gz$)')[, 2]
    ) %>%
    mutate(url = str_replace(url, '[.]bgz.*$', '.bgz?raw=1')) %>%

    filter(!(filename %in% existing_files) |
             file.info(paste0(folder, "/", filename))$size < 5e+8)  %>%


    mutate(filename = paste0(folder,
                             '/',
                             filename))

  to_download$folder %>%
    unique %>%
    walk(dir.create, showWarnings = FALSE, recursive = TRUE)

  if (nrow(to_download) != 0) {
    for (chunk in 1:ceiling(nrow(to_download) / chunk_size)) {
      if (chunk * chunk_size > max_downloads)
        stop('Maximum number of downloads reached.')
      download.file(to_download[((chunk - 1) * chunk_size + 1):min(nrow(to_download), (chunk * chunk_size)),]$url,
                    to_download[((chunk - 1) * chunk_size + 1):min(nrow(to_download), (chunk * chunk_size)),]$filename,
                    method = 'libcurl')
    }
  }
  nrow(to_download)
}


pairs_only <- function(household_info){

  household_info <- household_info %>% mutate(group_count = ave(household_info$group, household_info$group,  FUN = length))
  #table(household_info$group_count) ## number of individuals in each group household sizes
  household_info <- household_info[which(household_info$group_count==2),]

  household_pairs <- household_info

  household_list <- unique(household_pairs$group)
  pairs <- numeric()
  for(house in household_list)
  {

    house_ids <- household_pairs[which(household_pairs$group==house),]
    for(i in 1:(dim(house_ids)[1]-1))
    {
      for(j in i:dim(house_ids)[1])
      {
        if (i==j) next
        else     (pairs <- rbind(pairs,cbind(house_ids[i,1],house_ids[j,1], house)))
      }
    }
  }

  colnames(pairs)[1] <- "HOUSEHOLD_MEMBER1"  ##Index
  colnames(pairs)[2] <- "HOUSEHOLD_MEMBER2"  ##Household member
  colnames(pairs)[3] <- "HOUSE_ID"  ##Household member

  pairs <- as.data.frame(pairs)
  return(pairs)
}

find_kinship <- function(household_pairs, relatives){
  pairs <- household_pairs
  pairs$kinship <- 0
  for (i in 1:dim(relatives)[1])
  {
    ind1 <- relatives[i,1]
    ind2 <- relatives[i,2]
    if(any((pairs[,1]==ind1 & pairs[,2]==ind2) | (pairs[,2]==ind1 & pairs[,1]==ind2)))
    {
      index <- which((pairs[,1]==ind1 & pairs[,2]==ind2) | (pairs[,2]==ind1 & pairs[,1]==ind2))
      pairs$kinship[index] <- relatives[i,5]
      print(i)
    }
  }
  return(pairs)
}

filter_pairs <- function(pairs,household_relationships,field){
  pairs$sex1 <- household_relationships[["sex"]][match(pairs[["HOUSEHOLD_MEMBER1"]], household_relationships$userId)]
  pairs$sex2 <- household_relationships[["sex"]][match(pairs[["HOUSEHOLD_MEMBER2"]], household_relationships$userId)]

  #merge(pairs, household_relationships, by.x=, by.y=)
  pairs_full_sex <- pairs[which(!is.na(pairs$sex1) & !is.na(pairs$sex2)),]

  pairs_unrelated <- subset(pairs_full_sex, pairs_full_sex$kinship<0.05)

  pairs_unrelated_hetero <- pairs_unrelated[-which(pairs_unrelated$sex1==pairs_unrelated$sex2),]
  pairs_unrelated_hetero2 <- pairs_unrelated_hetero[-which(is.na(pairs_unrelated_hetero$sex1) | is.na(pairs_unrelated_hetero$sex2)),]

  temp1 <- pairs_unrelated_hetero[which(pairs_unrelated_hetero$sex1==0),c("HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1", "HOUSE_ID", "kinship", "sex2", "sex1", "tar_group")]

  colnames(temp1) <- colnames(pairs_unrelated_hetero)

  temp2 <- pairs_unrelated_hetero[which(pairs_unrelated_hetero$sex1==1),c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSE_ID", "kinship", "sex1", "sex2", "tar_group")]
  pairs_filter <- rbind(temp2, temp1)
  ## 79150      6

  pairs_filter$house_rel_member1 <- household_relationships[[field]][match(pairs_filter[["HOUSEHOLD_MEMBER1"]], household_relationships$userId)]
  pairs_filter$house_rel_member2 <- household_relationships[[field]][match(pairs_filter[["HOUSEHOLD_MEMBER2"]], household_relationships$userId)]
  pairs_filter2 <- pairs_unrelated_hetero[which(pairs_filter$house_rel_member1==TRUE & pairs_filter$house_rel_member2==TRUE ) ,]
  ## 75048     6


  colnames(pairs_filter2)[which(colnames(pairs_filter2)=="sex1")] <- "HOUSEHOLD_MEMBER1_sex"
  colnames(pairs_filter2)[which(colnames(pairs_filter2)=="sex2")] <- "HOUSEHOLD_MEMBER2_sex"
  return(pairs_filter2)

}

munge_sqc <- function(sqc,fam){
  sqc_sub <- sqc[,c(26:65)]
  colnames(sqc_sub) <- c(sapply(1:40, function(x) {paste0("PC_",x)}))
  sqc_sub$ID <- fam[,1]
  return(sqc_sub)
}

add_pcs <- function(pairs,pheno,sqc_sub){
  ## get rid of tar_grouping -> don't need it anymore
  pairs <- pairs %>% dplyr::select(-tar_group)
  out <- list()
  for(i in 1:dim(pairs)[1])
  {
    if(pairs[["HOUSEHOLD_MEMBER1_sex"]][i]==0)
    {
      male_ID <-pairs[["HOUSEHOLD_MEMBER2"]] [i]
      female_ID <-pairs[["HOUSEHOLD_MEMBER1"]] [i]
      pairs[["HOUSEHOLD_MEMBER1"]][i] <-male_ID
      pairs[["HOUSEHOLD_MEMBER2"]][i] <-female_ID
    }

  }
  pairs[["HOUSEHOLD_MEMBER1_sex"]] <- 1
  pairs[["HOUSEHOLD_MEMBER2_sex"]] <- 0

  for(merge_by in c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2"))
  {

    merge1_temp <- merge(pairs, pheno, by.x=merge_by, by.y="userId")
    #length(which(!pairs[,1] %in% pheno1_sub[,1]))
    ### when merge_by="HOUSEHOLD_MEMBER1" -> 3 people were in household file and not pheno
    #length(which(!pairs[,2] %in% pheno1_sub[,1]))
    ### when merge_by="HOUSEHOLD_MEMBER2" -> 1 person were in household file and not pheno

    merge2_temp <- merge(merge1_temp, sqc_sub, by.x=merge_by, by.y="ID")
    #length(which(!merge1_temp[,"HOUSEHOLD_MEMBER1"] %in% sqc_sub[,"ID"]))
    ### when merge_by="INDEX" -> 2009 people were in household file and not geno PC file
    #length(which(!merge1_temp[,"HOUSEHOLD_MEMBER2"] %in% sqc_sub[,"ID"]))
    ### when merge_by="HOUSEHOLD_MEMBER" -> 2799 people were in household file and not geno PC file

    colnames(merge2_temp) <- c(colnames(merge2_temp)[1:6], paste0(merge_by, "_", colnames(merge2_temp)[-c(1:6)]))

    out[[paste0(merge_by, "_pheno_data")]] <- merge2_temp
    #write.csv(merge2_temp, paste0("analysis/data_setup/", merge_by,"_pheno_model_adjustments", ".csv"),row.names=F, quote=T )

  }

  return(out)
}

calc_time_together <- function(pheno_cov,time_at_household,time_at_household_raw, time_at_address_raw_field, time_at_address_field){

  households_temp <- pheno_cov[,1:3]

  time_at_household_raw[[time_at_address_raw_field]][time_at_household_raw[[time_at_address_raw_field]] == -10 ] <- NA
  time_at_household_raw[[time_at_address_raw_field]][time_at_household_raw[[time_at_address_raw_field]] == -3 ] <- NA
  time_at_household_raw[[time_at_address_raw_field]][time_at_household_raw[[time_at_address_raw_field]] == -1 ] <- 0 #less than 1 year

  households_temp$house_time_member1 <- time_at_household[[time_at_address_field]][match(households_temp[["HOUSEHOLD_MEMBER1"]], time_at_household$userId)]
  households_temp$house_time_member2 <- time_at_household[[time_at_address_field]][match(households_temp[["HOUSEHOLD_MEMBER2"]], time_at_household$userId)]
  households_temp$house_time_raw_member1 <- time_at_household_raw[[time_at_address_raw_field]][match(households_temp[["HOUSEHOLD_MEMBER1"]], time_at_household_raw$eid)]
  households_temp$house_time_raw_member2 <- time_at_household_raw[[time_at_address_raw_field]][match(households_temp[["HOUSEHOLD_MEMBER2"]], time_at_household_raw$eid)]

  households_temp <- transform(households_temp, time_together = pmin(house_time_member1, house_time_member2))
  households_temp <- transform(households_temp, time_together_raw = pmin(house_time_raw_member1, house_time_raw_member2))
  return(households_temp)
}

munge_household_time <- function(household_time, pheno_cov,
                                 time_together_min,time_together_max,time_together_interval,
                                 age_min, age_max, age_interval, num_household_bins){



  household_time[["time_together_even_width"]] <- cutr::smart_cut(household_time[["time_together_raw"]],seq(time_together_min,time_together_max, by=time_together_interval))
  household_time[["time_together_even_bins"]] <- cutr::smart_cut(household_time[["time_together_raw"]],list(num_household_bins,"balanced"),"groups")

  household_time$age_member1 <- pheno_cov[["HOUSEHOLD_MEMBER1_age"]][match(household_time[["HOUSEHOLD_MEMBER1"]], pheno_cov[["HOUSEHOLD_MEMBER1"]])]
  household_time$age_member2 <- pheno_cov[["HOUSEHOLD_MEMBER2_age"]][match(household_time[["HOUSEHOLD_MEMBER2"]], pheno_cov[["HOUSEHOLD_MEMBER2"]])]

  household_time[["mean_age"]] <- rowMeans(household_time[c('age_member1', 'age_member2')], na.rm=TRUE)

  household_time[["age_even_width"]] <- cutr::smart_cut(household_time[["mean_age"]],seq(age_min,age_max, by=age_interval))
  household_time[["age_even_bins"]] <- cutr::smart_cut(household_time[["mean_age"]],list(num_household_bins,"balanced"),"groups")

  return(household_time)

}


### File A1 ----
## Description: tests phenotypic correlations between household pairs that have been
## filtered for kinship < 0.05, and to be opposite-sex pairs
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

compute_trait_corr <- function(phesant_directory,UKBB_directory,pairs_filter){
  unique_phes_ids <- unique(phesant_directory[,2])
  trait_corr <- numeric()
  prev_file <- ""
  for (id in unique_phes_ids)
  {
    id <- as.character(id)
    phenotype_sub <- unlist(strsplit (id,"_"))[1]

    UKBB_directory[,"field.showcase"] <- as.character(UKBB_directory[,"field.showcase"])
    file_with_pheno <- which(UKBB_directory[,"field.showcase"]==id)
    file_name <- names(which.max(table(as.character(UKBB_directory[file_with_pheno,"file"]))))
    if(length(file_with_pheno)==0)
    {
      file_with_pheno <- which(UKBB_directory[,"field.showcase"]==phenotype_sub)
      file_name <- names(which.max(table(as.character(UKBB_directory[file_with_pheno,"file"]))))
    }
    description <- unlist(strsplit (as.character(UKBB_directory[file_with_pheno[1],"col.name"]),"_datacoding"))[1]
    current_file <- as.character(phesant_directory[which(phesant_directory[,2]==id),"File"])[1]

    #sgg_file <- sgg_paths[grepl(file_name,current_file)][1]
    #flow <- gsub(basename(sgg_file),"variable-flow-counts-all.txt", sgg_file)
    #file_dir <- gsub("/([^/]*)$" , "", sgg_file)
    #bin_num <- substrRight(file_dir,1)
    #log <- read.table(paste0(file_dir,"/out_bin",bin_num, "..log"))

    if(prev_file!=current_file)
    {tsv_data <- fread(current_file, header=TRUE, sep='\t',data.table=F)}
    r2 <- NA
    pairs_filter_copy <- pairs_filter
    pairs_filter_copy$trait1 <- tsv_data[[id]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER1"]], tsv_data$userId)]
    pairs_filter_copy$trait2 <- tsv_data[[id]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER2"]], tsv_data$userId)]
    complete_pairs <- pairs_filter_copy[which(complete.cases(pairs_filter_copy$trait1 ,pairs_filter_copy$trait2 )),]
    n_completed_pairs <- length(pairs_filter_copy[which(complete.cases(pairs_filter_copy$trait2, pairs_filter_copy$trait1)),"trait2"])
    if(n_completed_pairs>1)
    {
      sd1 <- sd(complete_pairs$trait1,na.rm=TRUE)
      sd2 <- sd(complete_pairs$trait2,na.rm=TRUE)
      if(sd1!=0 & sd2!=0)
      {
        r <- cor(pairs_filter_copy$trait1, pairs_filter_copy$trait2, method = c("pearson"),use = "pairwise.complete.obs")
        r2 <- r^2
      }
    }

    prev_file <- current_file
    trait_row <- cbind(id, phenotype_sub,description, n_completed_pairs,r2)
    trait_corr <- rbind(trait_corr, trait_row)
  }
  colnames(trait_corr) <- c("ID", "ID_sub", "description", "N_pairs","r2")
  cat(paste0("Household phenotypic correlations successfully computed, and saved to:\n",
             "'output/tables/1.household_correlations.csv'.\n"))
  return(trait_corr)
}

## File A2 ----

## Description: filter phenotype correlations for > specified threshold (in settings).
## Check which phenotypes have been downloaed, and create list of need to be downloaded
## and defined here: "output/tables/define_Neale_categories.csv".
## Manually categorize this list and save as: "output/tables/define_Neale_categories_filled.csv"
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

SGG_link_with_Neale <- function(Neale_SGG_dir){

  Neale_SGG_dir_filt <- subset(Neale_SGG_dir, Neale_SGG_dir$phesant_processed=="YES")
  unique_Neale_ids <- unique((Neale_SGG_dir_filt[,"phenotype"]))

  both_sex_avail <- numeric()
  for(id in unique_Neale_ids)
  {
    id <- as.character(id)
    sex <- c(as.character(Neale_SGG_dir_filt[which(Neale_SGG_dir_filt[,"phenotype"]==id),"phenotype_sex"]))
    if("male" %in% sex & "female" %in% sex & "both" %in% sex)
    {
      both_sex_avail <- c(both_sex_avail,id)
    }
  }
  cat(paste0("There are: ", length(both_sex_avail), " phenotypes with Neale summary stats with joint
  and sex-specific data that have also been processed in the SGG database.\n\n")) ###1243 traits
  Neale_SGG_dir_filt2 <- subset(Neale_SGG_dir_filt, Neale_SGG_dir_filt$phenotype %in% both_sex_avail)
  return(Neale_SGG_dir_filt2)
}

filter_by_corr <- function(trait_corr,Neale_SGG_dir_filt2,household_correlation_threshold){

  merge_temp <- merge(trait_corr, Neale_SGG_dir_filt2, by.x="ID", by.y="sgg_phesant_name", fill=T)
  merge_temp$r2 <- as.numeric(as.character(merge_temp$r2))
  corr_traits <- merge_temp[which(sqrt(merge_temp$r2) > household_correlation_threshold),]
  corr_traits <- corr_traits[,-which(colnames(corr_traits) %in% c("description.x"))]
  colnames(corr_traits) <- c("SGG_PHESANT_ID", "SGG_PHESANT_ID_sub","N_pairs","r2",
                             "Neale_pheno_ID","Neale_pheno_ID_sub", "Neale_file_sex","description",
                             "variable_type", "SGG_request", "PHESANT_processed", "SGG_location",
                             "Neale_downloaded","category","define_category:T/F",
                             "v2_exists", "v2_downloaded")
  return(corr_traits)
}

organize_Neale <- function(traits_corr_filter){
  corr_traits_sub <- subset(traits_corr_filter, traits_corr_filter$Neale_file_sex=="both")
  define_cats <- corr_traits_sub[,c("Neale_pheno_ID","Neale_pheno_ID_sub","description","variable_type","Neale_downloaded","category","define_category:T/F")]
  download_rest <- subset(define_cats,define_cats[["Neale_downloaded"]]=="PARTIALLY" & define_cats[["define_category:T/F"]]==FALSE)
  #write.csv(corr_traits, "pull_Neale_phenotypes.csv", row.names=F)
  date <- Sys.Date()
  define_cats2 <- subset(define_cats,define_cats[["define_category:T/F"]]==TRUE)

  if(dim(define_cats2)[1]==0 & dim(download_rest)[1]==0)
  {

    cat("All necessary Neale files have been downloaded, proceed in analysis.\n\n")
  }
  if(dim(define_cats2)[1]!=0)
  {

    cat(paste0("Some Neale files need to be categorized, please categorize appropriately within file:
      'output/tables/define_Neale_categories.csv'.\n\n"))
  }

  if(dim(download_rest)[1]!=0)
  {

    cat(paste0("Some Neale files need to be downloaded, but have been categorized previously,
      likely they were updated based on Neale modifications to the PHESANT analysis.
      See here: https://www.dropbox.com/s/ayto55r0z1eigwa/phenotype_update_note.txt?dl=0, for details.\n\n"))
  }

  cat(paste0("Finished linking phenotypic correlations file with Neale files.\n\n"))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have valid Neale summary stats in both sexes seperately and together,
  and have household correlations greater than specified r2 threshold (see 'scripts/settings.r' for chosen threshold),\n",
             "filtered phenotypic correlations have been saved to:
  'output/tables/2.household_correlations.corr_filter.csv'.\n\n"))

  return(list("download_rest"=download_rest, "define_cats"=define_cats2))
}


write_define_cats <- function(Neale_to_process){
  write.csv(Neale_to_process$define_cats, "output/tables/define_Neale_categories.csv", row.names=F)
  return("output/tables/define_Neale_categories.csv")
}

write_download_list <- function(Neale_to_process){
  write.csv(Neale_to_process$download_rest, "analysis/download_Neale_list.csv", row.names=F)
  return("analysis/download_Neale_list.csv")
}

### File A3 ----

## Description: download missing neale files using categories defined here:
## "output/tables/define_Neale_categories_filled.csv".
## Update categories in phenotypic correlations file: output/tables/pheno_household_filter_corr.csv
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

download_Neale <- function(filled_cats,download_rest,traits_corr_filter, reference_file_name){
  full_dl_list <- rbind(filled_cats, download_rest)
  select <- dplyr::select

  for(category
      in c("body", "brain", "diet", "disease", "disease_proxy", "lifestyle", "parental_pheno"))
  {
    IDs <- as.character(full_dl_list[which(full_dl_list[,"category"]==category),"Neale_pheno_ID"])
    download_neale_files( IDs,
                          category = category , reference_file=reference_file_name)
  }

  corr_traits <- traits_corr_filter
  levels(corr_traits$category) <- union (levels(corr_traits$category), levels(filled_cats$category))
  index <- which(corr_traits[["define_category:T/F"]]==TRUE)
  pheno_missing_cat <- corr_traits[index,"Neale_pheno_ID"]
  corr_traits[index, "category"] <- as.factor(filled_cats[["category"]][match(pheno_missing_cat, filled_cats$Neale_pheno_ID)])

  cat(paste0("Missing Neale files have been successfully downloaded and categorized,\n",
             "now they need to be processed (clumped for LD and extracted for p<0.1).\n"))

  return(corr_traits)
}

update_download_info <- function(corr_traits,Neale_SGG_dir){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]
  irnt=TRUE

  index <- which(corr_traits[["v2_exists"]]==TRUE)
  index_df <- corr_traits[index,c("Neale_pheno_ID","Neale_file_sex")]
  index_df2 <- merge(index_df, Neale_SGG_dir, by.x=c("Neale_pheno_ID","Neale_file_sex"), by.y=c("phenotype","phenotype_sex"),all.x=TRUE)
  now_downloaded <- index_df2$v2_downloaded
  corr_traits[index, "v2_downloaded"] <- now_downloaded
  return(corr_traits)

}

file_in_out <- function(traits_corr2_update,i,reference_file,IV_threshold, Neale_summary_dir){

  corr_traits_both <- traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),]

  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))
  irnt <- TRUE
  ID <- corr_traits_both[i,"Neale_pheno_ID"]
  phenotype_ids  =  paste0( '^', ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
    paste( collapse = '|' )

  file_name_temp  =  reference_file %>%
    filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")

  file_name <- file_name_temp[["Phenotype Code"]][1]

  folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]

  corr_traits_both <- traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"]
  out_file <- paste0("analysis/data_setup/IV_lists/", Neale_id, "_IVs_", IV_threshold,"_both_sexes.txt")

  return(list(in_file = folder, out_file = out_file))

}

pull_traits_to_count_IVs <- function(traits_corr2_update){
  output <- tibble(Neale_pheno_ID =traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),"Neale_pheno_ID"])
  return(output)
}

pull_traits_to_calc_het <- function(traits_corr3_run){
  output <- tibble(Neale_pheno_ID =traits_corr3_run[,"Neale_pheno_ID"])
  return(output)
}

pull_traits_to_run <- function(traits_corr5){
  output <- tibble(Neale_pheno_ID =traits_corr5[,"Neale_pheno_ID"])
  return(output)
}

pull_IV_indices_to_run <- function(traits_final, traits_to_calc_het){

  indicies <- which(traits_to_calc_het$Neale_pheno_ID %in% traits_final$Neale_pheno_ID)
  return(indicies)
}

get_IV_list <- function(corr_traits, Neale_pheno_ID, reference_file, IV_threshold, Neale_output_path, Neale_summary_dir){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]

  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                      recursive = TRUE,
                                      pattern = '[.]gz' )
  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))
  irnt=TRUE

  i <- which(corr_traits_both[["Neale_pheno_ID"]]==Neale_pheno_ID)
  category <- corr_traits_both[i,"category"]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"] #same as Neale_pheno_ID
  ID <- corr_traits_both[i,"Neale_pheno_ID"]
  phenotype_ids  =  paste0( '^', ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
    paste( collapse = '|' )

  file_name_temp  =  reference_file %>%
    filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")

  file_name <- file_name_temp[["Phenotype Code"]][1]
  folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
  #folder <- IVs_full[which(grepl(paste0(file_name ), IVs_full) & grepl("both_sexes", IVs))]
  if(length(folder)==2){folder <- folder[grep("v2",folder)]} ## use version 2 if it exists
  name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]
  if(length(name)==2){name <- name[grep("v2",name)]} ## use version 2 if it exists

  threshold <- IV_threshold ## only extract IVs
  build_data <- numeric()
  for(chr in 1:22)

  {
    original_data <- read.table(paste0(folder,"/", name,"_unpruned_chr", chr, ".txt"), header=T)

    if(file.exists(paste0(folder,"/", name,"_unpruned_chr", chr, ".clumped")))
    {
      clump_data <- read.table(paste0(folder,"/", name,"_unpruned_chr", chr, ".clumped"),header = T)
      clumped_SNPs <- as.character(clump_data[,"SNP"])
      IV_temp <- as.character(original_data[which(original_data[["P"]]<threshold & original_data[["SNP"]] %in% clumped_SNPs),"SNP"])
      build_data <- c(build_data, IV_temp) #
      #write.table(IV_temp, file_out, append=T, row.names=F, col.names=F, quote=F)
    }

  }

  return(IV_snp_list = build_data)
}


summarize_IVs <- function(corr_traits, i, reference_file, Neale_summary_dir){

  irnt=TRUE
  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]
  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))
  category <- corr_traits_both[i,"category"]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"]
  ID <- corr_traits_both[i,"Neale_pheno_ID"]
  phenotype_ids  =  paste0( '^', ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
    paste( collapse = '|' )

  file_name_temp  =  reference_file %>%
    filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")

  file_name <- file_name_temp[["Phenotype Code"]][1]
  folder <- IVs_full[which(grepl(file_name, IVs_full) & grepl("both_sexes", IVs))]
  name <- IVs[which(grepl(file_name, IVs) & grepl("both_sexes", IVs))]
  IV_snp_list <- numeric()
  threshold <- IV_threshold
  #if(file.exists(paste0( "analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt")))
  if(file.info(paste0( "analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"))$size!=0)
  {
    GW_ivs <- fread( paste0( "analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"),data.table=F, header=F)
    out <- cbind(i,dim(GW_ivs)[1])
  } else out <- cbind(i,0)

  return(as.data.frame(out))

}

IV_filter <- function(corr_traits, count_IVs, num_IVs_threshold){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]

  corr_traits_both <- merge(corr_traits_both, count_IVs, by="Neale_pheno_ID")
  to_run <- corr_traits_both[which(corr_traits_both$num_IVs >=num_IVs_threshold),]
  return(list(to_run = to_run, non_filtered = corr_traits_both))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have at least minimum number of IVs according to specified threshold
  (see 'scripts/settings.r' for chosen threshold),\n",
             "filtered phenotypic correlations have been saved to:
  'output/tables/3.household_correlations.numIVs_filter.csv'.\n"))

}

reduce_Neale_variant_data <- function(path_Neale_variants, variants_to_extract){

  variant_data <- fread(path_Neale_variants,data.table=F)


  SNP_rows <- which(variant_data[,"rsid"] %in% variants_to_extract)
  reduced_data <- variant_data[SNP_rows,]

  return(reduced_data)

}

# calc_sex_het
summarize_IV_data <- function(traits, Neale_pheno_ID, variant_data, reference_file, Neale_summary_dir, Neale_output_path, IV_threshold){

  irnt=TRUE
  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                      recursive = TRUE,
                                      pattern = '[.]gz' )

  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))

  #traits$num_IVs_pass_het <- NA

  result <- NA

  trait_ID <- Neale_pheno_ID
  i <- which(traits[["Neale_pheno_ID"]]==trait_ID)

  phenotype_ids  =  paste0( '^', trait_ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
    paste( collapse = '|' )
  file_name_temp  =  reference_file %>%
    filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")
  file_name <- file_name_temp[["Phenotype Code"]][1]

  #trait_info <- read.table(paste0(pheno_dir,"/trait_info.txt"), header=F, row.names=1)

  IV_list_both_sexes <- fread(paste0( "analysis/data_setup/IV_lists/", trait_ID, "_IVs_", IV_threshold,"_both_sexes.txt"), data.table=F, header=F)
  SNP_rows <- which(variant_data[,"rsid"] %in% IV_list_both_sexes[,1])
  IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
  IV_file_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]

  if(length(IV_folder)==2){IV_folder <- IV_folder[grep("v2",IV_folder)]} ## use version 2 if it exists
  if(length(IV_file_name)==2){IV_file_name <- IV_file_name[grep("v2",IV_file_name)]} ## use version 2 if it exists

  for(exposure_sex in c("both_sexes", "male", "female"))
  {

    specific_IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl(paste0("[.]",exposure_sex), IVs))]
    specific_IV_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl(paste0("[.]",exposure_sex), IVs))]
    if(length(specific_IV_folder)==2){specific_IV_folder <- specific_IV_folder[grep("v2",specific_IV_folder)]} ## use version 2 if it exists
    if(length(specific_IV_name)==2){specific_IV_name <- specific_IV_name[grep("v2",specific_IV_name)]} ## use version 2 if it exists

    original_neale_temp <- existing_files_full[which(grepl(paste0("/",gsub(".IVs", "",specific_IV_name)), existing_files_full))]
    original_neale_file <- paste0(Neale_output_path, "/", original_neale_temp)
    assign(paste0(exposure_sex, "_IV_folder"), specific_IV_folder)
    assign(paste0(exposure_sex, "_original_Neale_file"), original_neale_file)

  }

  for(file in c("male", "female"))
  {

    temp <- fread(get(paste0(file,"_original_Neale_file")), data.table=F)
    temp <- temp %>% rename(minor_allele_Neale = minor_allele) %>% rename(AC_Neale  = AC) %>%
      rename(minor_AF_Neale = minor_AF)


    variants <- variant_data[SNP_rows,"variant"]
    neale_data_temp <- temp[which(temp$variant %in% variants),]
    var_data_temp <- variant_data[SNP_rows,]
    exp_var_merge <- merge(neale_data_temp,var_data_temp,by="variant")
    #IV_cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "AF","n_complete_samples")
    #exp_var_merge <- exp_var_merge[,IV_cols]

    assign(paste0(file, "_IV_data"), exp_var_merge)

  }

  #write.table(male_IV_data, male_out_file,row.names=F, col.names=T, quote=F)
  #write.table(female_IV_data, female_out_file,row.names=F, col.names=T, quote=F)


  IV_list_both_sexes$p_het <- NA
  colnames(IV_list_both_sexes)[1] <- "SNP"
  for(k in 1:length(SNP_rows))
  {
    m_snp <- which(male_IV_data[["rsid"]]==IV_list_both_sexes[["SNP"]][k])
    f_snp <- which(female_IV_data[["rsid"]]==IV_list_both_sexes[["SNP"]][k])

    beta_F <-  female_IV_data[f_snp,"beta"]
    beta_M <-  male_IV_data[m_snp,"beta"]
    SE_F <- female_IV_data[f_snp,"se"]
    SE_M <- male_IV_data[m_snp,"se"]
    n_M <- male_IV_data[m_snp,"n_complete_samples"]
    n_F <- female_IV_data[f_snp,"n_complete_samples"]
    #  se <- sqrt( (SE_F^2/n_F) + (SE_M^2/n_M) )
    se <- sqrt( (SE_F^2) + (SE_M^2) )
    t <- (beta_F-beta_M)/se
    p_het <- 2*pnorm(-abs(t))
    IV_list_both_sexes$p_het[k] <- p_het
  }



  IV_list_filter <- IV_list_both_sexes[which(IV_list_both_sexes$p_het > 0.05/length(SNP_rows)),]

  num_pass_filter <- length(which(IV_list_both_sexes$p_het > 0.05/length(SNP_rows)))

  colnames(IV_list_both_sexes) <- c("SNP", "P-het")


  #write.table(IV_list_both_sexes, het_output_file, row.names=F, col.names=T, quote=F)

  char_row <- data.frame(lapply(traits[i,], as.character), stringsAsFactors=FALSE)
  sex_het_summary <- as.data.frame(t(unlist(c(char_row, num_pass_filter))))
  output <- list(male_IV_data = male_IV_data, female_IV_data = female_IV_data, IV_list_both_sexes = IV_list_both_sexes, sex_het_summary = sex_het_summary)
  return(output)

}

write_IV_list <- function(traits_corr2_update, Neale_pheno_ID, IV_list, IV_threshold, dir) {

  corr_traits_both <- traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),]

  #for(i in 1:dim(traits_to_count_IVs)[1]){

    Neale_id <- Neale_pheno_ID
    out_file <- paste0("analysis/data_setup/IV_lists/", Neale_id, "_IVs_", IV_threshold,"_both_sexes.txt")
    if(file.exists(out_file)) {file.remove(out_file)}
    write.table(IV_list, out_file, append=F, row.names=F, col.names=F, quote=F)
    return(out_file)
  #}

}

write_IV_info <- function(IV_data_summary, Neale_pheno_ID) {

    trait_ID <- Neale_pheno_ID
    male_out_file <- paste0( "analysis/data_setup/IV_info/", trait_ID, "_IVs_5e-08_male.txt")
    female_out_file <- paste0( "analysis/data_setup/IV_info/", trait_ID, "_IVs_5e-08_female.txt")
    write.table(IV_data_summary$male_IV_data, male_out_file,row.names=F, col.names=T, quote=F)
    write.table(IV_data_summary$female_IV_data, female_out_file,row.names=F, col.names=T, quote=F)
    out_names <- c(male_out_file, female_out_file)


}

write_sex_het <- function(IV_data_summary, Neale_pheno_ID){


    trait_ID <- Neale_pheno_ID
    output_file <- paste0( "analysis/data_setup/sex_heterogeneity/", trait_ID, "_sex_het.txt")
    write.table(IV_data_summary$IV_list_both_sexes, output_file, row.names=F, col.names=T, quote=F)
    #IV_list_both_sexes is in spot 3
    return(output_file)

}

write_pheno_data <- function(pheno_data, Neale_pheno_ID){

  trait_ID <- Neale_pheno_ID
  pheno_dir <- paste0("analysis/traitMR")

  male_file <- paste0(pheno_dir,"/pheno_files/phesant/", trait_ID, "_male.txt")
  female_file <- paste0(pheno_dir,"/pheno_files/phesant/", trait_ID, "_female.txt")

  write.table(pheno_data$unrelated_male_data, male_file,row.names=F, quote=F)
  write.table(pheno_data$unrelated_female_data, female_file, row.names=F, quote=F)

  return(c(male_file, female_file))

}

sex_het_filter <- function(corr_traits, sex_het_summary, num_IVs_threshold){
  list_length <- 4
  build_df <- numeric()
  for(i in 1:dim(corr_traits)[1]){
    data_out <- sex_het_summary[[(i)]][[4]] #sex_het_summary is in spot 4
    build_df <- rbind(build_df, data_out)
  }
  colnames(build_df) <- c(colnames(corr_traits),"num_IVs_pass_het")
  result_df <- as.data.frame(build_df)
  result_df$num_IVs_pass_het <- as.numeric(as.character(result_df$num_IVs_pass_het))
  to_run <- result_df[which(result_df$num_IVs_pass_het >=num_IVs_threshold),]
  removed <- result_df[-which(result_df$num_IVs_pass_het >=num_IVs_threshold),]
  return(list(to_run = to_run, non_filtered = result_df))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have at least minimum number of IVs according to specified threshold after removing IVs with
  with significant evidence of heterogeneity between sexes (p<0.05/[number of GW-sig IVs]),
  filtered phenotypic correlations have been saved to:
  'output/tables/4.household_correlations.sexhet_filter.csv'.\n"))
}

continuous_filter <- function(traits){
  output <- traits %>% filter(variable_type=="ordinal" | variable_type=="continuous_irnt")
  return(output)
}

write_final_filter <- function(traits_final, output_name){
  output_name <- ("output/tables/household_correlations.final_filter.csv")
  write.csv(traits_final, "output/tables/household_correlations.final_filter.csv", row.names=F)
  return(output_name)
}


check_valid_GRS_input <- function(traits,reference_file, Neale_output_path, Neale_summary_dir){

  irnt=TRUE
  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                      recursive = TRUE,
                                      pattern = '[.]gz' )
  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))

  result <- numeric()
  for(i in 1:dim(traits)[1])
  {

    trait_ID <- as.character(traits[i,"Neale_pheno_ID"])
    phenotype_ids  =  paste0( '^', trait_ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
      paste( collapse = '|' )
    file_name_temp  =  reference_file %>%
      filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")
    file_name <- file_name_temp[["Phenotype Code"]][1]

    #trait_info <- read.table(paste0(pheno_dir,"/trait_info.txt"), header=F, row.names=1)

    IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
    IV_file_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]

    if(length(IV_folder)==2){IV_folder <- IV_folder[grep("v2",IV_folder)]} ## use version 2 if it exists
    if(length(IV_file_name)==2){IV_file_name <- IV_file_name[grep("v2",IV_file_name)]} ## use version 2 if it exists

    for(exposure_sex in c("both_sexes", "male", "female"))
    {

      specific_IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl(paste0("[.]",exposure_sex), IVs))]
      specific_IV_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl(paste0("[.]",exposure_sex), IVs))]
      if(length(specific_IV_folder)==2){specific_IV_folder <- specific_IV_folder[grep("v2",specific_IV_folder)]} ## use version 2 if it exists
      if(length(specific_IV_name)==2){specific_IV_name <- specific_IV_name[grep("v2",specific_IV_name)]} ## use version 2 if it exists

      original_neale_temp <- existing_files_full[which(grepl(paste0("/",gsub(".IVs", "",specific_IV_name)), existing_files_full))]
      original_neale_file <- paste0(Neale_output_path, "/", original_neale_temp)
      assign(paste0(exposure_sex, "_IV_folder"), specific_IV_folder)
      assign(paste0(exposure_sex, "_original_Neale_file"), original_neale_file)

    }

    no_GRS_snps <- FALSE

    GRS_length <- numeric()
    for(file in c("male", "female"))
    {
      specific_IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl(paste0("[.]",file), IVs))]
      specific_IV_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl(paste0("[.]",file), IVs))]
      if(length(specific_IV_folder)==2){specific_IV_folder <- specific_IV_folder[grep("v2",specific_IV_folder)]} ## use version 2 if it exists
      if(length(specific_IV_name)==2){specific_IV_name <- specific_IV_name[grep("v2",specific_IV_name)]} ## use version 2 if it exists

      temp <- fread(paste0(get(paste0(file,"_IV_folder")), "/",specific_IV_name, "_unpruned.txt"), data.table=F, nrow=5)
      GRS_length <- c(GRS_length, dim(temp)[1])

    }
    if(any(GRS_length==0))
    {
      no_GRS_snps <- TRUE
    }

    out <- cbind(i,no_GRS_snps)
    result <- rbind(result,out)
  }

  return(result)
}

valid_GRS_filter <- function(traits,valid_GRS_summary){


  result <- valid_GRS_summary[order(valid_GRS_summary[,1]),]
  traits$GRS_list <- result[,2]

  to_run <- traits[-which(result[,2]==1),]


  cat(paste0("The following trait ID(s): ", paste(as.character(traits[which(result[,2]==1),"Neale_pheno_ID"]), collapse="; "), ",
  with pheno description(s): ", paste(as.character(traits[which(result[,2]==1),"description"]), collapse="; "), ",
  do not have any SNPs that pass the clump/GRS filters in at least one sex (p<0.1 and 'low_quality_SNP'==FALSE),
  and have been removed from all subsequent analyses.\n\n"))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have a valid sex-specific input for the GRS calculation (p<0.1 and 'low_quality_SNP'==FALSE),
  filtered phenotypic correlations have been saved to:
  'output/tables/5.household_correlations.baseGRS_filter.csv'.\n"))

  return(to_run)

}

create_trait_dirs <- function(Neale_pheno_ID){

  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description

  print(trait_ID)
  dir.create(paste0("output/tables/traitMR/", trait_ID), showWarnings = FALSE)
  dir.create(paste0("output/figures/traitMR/", trait_ID), showWarnings = FALSE)

  pheno_dir <- paste0("analysis/traitMR")

  ### create_relevant directories
  dir.create(paste0(pheno_dir, "/household_MR/", trait_ID), showWarnings = FALSE)
  dir.create(paste0(pheno_dir, "/household_GWAS/", trait_ID), showWarnings = FALSE)


  dir.create(paste0(pheno_dir, "/IVs/Neale/", trait_ID), showWarnings = FALSE)




}

prep_pheno_data <- function(traits, Neale_pheno_ID, sqc, fam, relatives){

  i <- which(traits[["Neale_pheno_ID"]]==Neale_pheno_ID & traits[["Neale_file_sex"]]=="both")
  category <- as.character(traits[i,"category"])
  trait_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description
  phes_ID <- as.character(traits[i,"SGG_PHESANT_ID"])
  variable_type <- as.character(traits[i,"variable_type"])

   ## process pheno data to have IID column and no quotes.
  phesant_file <- as.character(traits[i,"SGG_location"]) #as.character(phesant_directory[which(phesant_directory[,2]==phes_ID),"File"])[1]
  tsv_data <- fread(phesant_file, header=TRUE, sep='\t',data.table=F,select=c("userId","sex","age", phes_ID))

  if(class(tsv_data[,phes_ID])=="logical")
  {
    tsv_data[,phes_ID] <- as.numeric(tsv_data[,phes_ID]) ## 1 turns to true
  }

  pheno_male <- tsv_data[which(tsv_data[["sex"]]==1),]
  pheno_female <- tsv_data[which(tsv_data[["sex"]]==0),]

  pheno_male_comp <- pheno_male[complete.cases(pheno_male), ]
  pheno_female_comp <- pheno_female[complete.cases(pheno_female), ]

  colnames(pheno_male_comp)[1] <- "IID"
  colnames(pheno_female_comp)[1] <- "IID"

  sqc_sub <- sqc[,c(26:65)]
  colnames(sqc_sub) <- c(sapply(1:40, function(x) {paste0("PC_",x)}))
  ## PCs are column 26-65

  ##fam file in same order as sample QC (sqc) file
  sqc_sub$ID <- fam[,1]

  pheno_male_full <- merge(pheno_male_comp, sqc_sub, by.x="IID", by.y="ID")
  pheno_female_full <- merge(pheno_female_comp, sqc_sub, by.x="IID", by.y="ID")

  maleIDs <- as.integer(pheno_male_full$IID)
  femaleIDs <- as.integer(pheno_female_full$IID)


  related_males <- ukb_gen_samples_to_remove(relatives, ukb_with_data = maleIDs)
  related_females <- ukb_gen_samples_to_remove(relatives, ukb_with_data = femaleIDs)

  unrelated_male_data <- pheno_male_full[-which(pheno_male_full$IID %in% related_males ),]
  unrelated_female_data <- pheno_female_full[-which(pheno_female_full$IID %in% related_females ),]

  cat(paste0("Successfully prepared phenotype for processing (created directories, phenotype file inputs, trait summary, etc).\n"))

  pheno_data = list(unrelated_male_data = unrelated_male_data, unrelated_female_data = unrelated_female_data)
  return(pheno_data)
}

get_trait_info <- function(traits, Neale_pheno_ID, data_Neale_manifest, Neale_summary_dir, Neale_output_path){

  reference_file <- data_Neale_manifest

  ## find Neale file names and IV folder:
  i <- which(traits[["Neale_pheno_ID"]]==Neale_pheno_ID)
  category <- as.character(traits[i,"category"])
  trait_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description
  phes_ID <- as.character(traits[i,"SGG_PHESANT_ID"])
  variable_type <- as.character(traits[i,"variable_type"])
  phesant_file <- as.character(traits[i,"SGG_location"]) #as.character(phesant_directory[which(phesant_directory[,2]==phes_ID),"File"])[1]


  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))

  irnt=TRUE
  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                      recursive = TRUE,
                                      pattern = '[.]gz' )

  phenotype_ids  =  paste0( '^', trait_ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
    paste( collapse = '|' )
  file_name_temp  =  reference_file %>%
    filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")
  file_name <- file_name_temp[["Phenotype Code"]][1]
  IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
  IV_file_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]

  if(length(IV_folder)==2){IV_folder <- IV_folder[grep("v2",IV_folder)]} ## use version 2 if it exists
  if(length(IV_file_name)==2){IV_file_name <- IV_file_name[grep("v2",IV_file_name)]} ## use version 2 if it exists

  for(exposure_sex in c("both_sexes", "male", "female"))
  {

    specific_IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl(paste0("[.]",exposure_sex), IVs))]
    specific_IV_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl(paste0("[.]",exposure_sex), IVs))]
    if(length(specific_IV_folder)==2){specific_IV_folder <- specific_IV_folder[grep("v2",specific_IV_folder)]} ## use version 2 if it exists
    if(length(specific_IV_name)==2){specific_IV_name <- specific_IV_name[grep("v2",specific_IV_name)]} ## use version 2 if it exists

    original_neale_temp <- existing_files_full[which(grepl(paste0("/",gsub(".IVs", "",specific_IV_name)), existing_files_full))]
    original_neale_file <- paste0(Neale_output_path, "/", original_neale_temp)
    assign(paste0(exposure_sex, "_IV_folder"), specific_IV_folder)
    assign(paste0(exposure_sex, "_original_Neale_file"), original_neale_file)


  }

  num_IVs <- traits[i,"num_IVs_pass_het"]
  description <- as.character(traits[i,"description"])
  trait_info <- rbind(i,trait_ID,phes_ID,phesant_file,category,description, variable_type, num_IVs,
                      both_sexes_original_Neale_file, male_original_Neale_file, female_original_Neale_file,
                      both_sexes_IV_folder, male_IV_folder, female_IV_folder)


  num_IVs <- traits[i,"num_IVs_pass_het"]
  description <- as.character(traits[i,"description"])


  trait_info <- rbind(i,trait_ID,phes_ID,phesant_file,category,description, variable_type, num_IVs,
                      both_sexes_original_Neale_file, male_original_Neale_file, female_original_Neale_file,
                      both_sexes_IV_folder, male_IV_folder, female_IV_folder)


  return(trait_info)

}

extract_trait_info <- function(pheno_data){

  output <- pheno_data$trait_info
  return(output)

}


write_data_prep <- function(traits, traits_to_run, out1, out2){

  for(i in 1:dim(traits_to_run)[1]){
    print(i)
    data_out <- readd(data_prep, subtargets=i)
    category <- as.character(traits[i,"category"])
    trait_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description
    pheno_dir <- paste0("analysis/traitMR")
    #pheno_dir <- paste0(project_dir, "/analysis/traitMR")
    write.table(data_out[[1]]$unrelated_male_data, paste0(pheno_dir,"/pheno_files/phesant/", trait_ID, "_male.txt"),row.names=F, quote=F)
    write.table(data_out[[1]]$unrelated_female_data, paste0(pheno_dir,"/pheno_files/phesant/", trait_ID, "_female.txt"), row.names=F, quote=F)
    write.table(data_out[[1]]$trait_info, paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), row.names=T, col.names=F, quote=T)
  }
}

create_summary_stats <- function(Neale_pheno_ID, trait_info, IV_data_summary){

  #i <- which(traits[["Neale_pheno_ID"]]==Neale_pheno_ID)
  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description
  pheno_dir <- paste0("analysis/traitMR/")
  #trait_info <-  read.table(paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), header=F, row.names=1,check.names=F)

  #variant_data <- fread(variant_file_full_name,data.table=F)

  het_stats <- IV_data_summary$IV_list_both_sexes #fread(paste0( "analysis/data_setup/sex_heterogeneity/", trait_ID, "_sex_het.txt"), header=T, data.table=F)
  IV_list_filter <- het_stats[which(het_stats[["P-het"]] > 0.05/dim(het_stats)[1]),]

  for(sex in c("male", "female"))
  {
    #file_name <- paste0( "analysis/data_setup/IV_info/", trait_ID, "_IVs_5e-08_", file,".txt")
    temp <- IV_data_summary[[paste0(sex, "_IV_data")]] #read.table(file_name,header=T)
    SNP_rows <- which(temp[,"rsid"] %in% IV_list_filter[,1])
    temp <- temp[SNP_rows,]
    IV_cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "AF","n_complete_samples")
    exp_var_merge <- temp[,IV_cols]
    assign(paste0(sex, "_IV_data"), exp_var_merge)
  }


  #write.table(male_IV_data,paste0(pheno_dir,"/IVs/", trait_ID, "/male_IVs.txt"), row.names=F, col.names=T, quote=F)
  #write.table(female_IV_data,paste0(pheno_dir,"/IVs/", trait_ID, "/female_IVs.txt"), row.names=F, col.names=T, quote=F)

  ###############################



  cat(paste0("Successfully prepared MR inputs from Neale summary stats.
  MR inputs are sex-specific but based on GW-significance in both-sex file
  (only including SNPs which passed filters - see pipeline part A)\n"))
  output = list(male_IV_data = male_IV_data, female_IV_data = female_IV_data)
  return(output)

}

write_summ_stats <- function(summ_stats_create, traits, traits_to_run, GRS_thresholds, out1, out2){

  pheno_dir <- paste0("analysis/traitMR/")
  for(i in 1:dim(traits)[1]){
    print(i)
    data_out <- readd(summ_stats_create, subtargets=i)
    trait_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description


    write.table(data_out[[1]][["male_IV_data"]],paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/male_IVs.txt"), row.names=F, col.names=T, quote=F)
    write.table(data_out[[1]][["female_IV_data"]],paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/female_IVs.txt"), row.names=F, col.names=T, quote=F)


    for(sex in c("male", "female")){

      for(GRS_threshold in GRS_thresholds){

        if(file.exists(paste0( pheno_dir,"/GRS/",trait_ID, "/", sex,"/", GRS_threshold,"/", sex,"_GRS_",GRS_threshold,".txt")))
        {file.remove(paste0( pheno_dir,"/GRS/",trait_ID, "/", sex, "/", GRS_threshold,"/", sex,"_GRS_",GRS_threshold,".txt"))}
      }

      for (chr in 1:22){

        for(GRS_threshold in GRS_thresholds){


          list_name <- paste0(sex, "_chr", chr, "_GRSthresh_", GRS_threshold)

          write.table(data_out[[1]][["GRS_list"]][[list_name]], paste0( pheno_dir,"/GRS/", trait_ID, "/", sex, "/", GRS_threshold,"/", sex,"_GRS_",GRS_threshold,"_chr",chr,".txt"), append=F, row.names=F, col.names=T, quote=F)
          write.table(data_out[[1]][["GRS_list"]][[list_name]], paste0( pheno_dir,"/GRS/", trait_ID, "/", sex, "/", GRS_threshold,"/", sex,"_GRS_",GRS_threshold,".txt"), append=T, row.names=F, col.names=F, quote=F)
        }

      }
    }

    ## add colnames
    colnames_GRS <- colnames(summ_stats_create[[i]][["GRS_list"]][[list_name]])
    for(sex in c("male", "female"))
    {
      for(GRS_threshold in GRS_thresholds)
      {
        file <- read.table(paste0( pheno_dir, "/GRS/", trait_ID, "/", sex, "/", GRS_threshold, "/", sex, "_GRS_", GRS_threshold, ".txt"), header=F)
        colnames(file) <- colnames_GRS
        write.table(file, paste0( pheno_dir, "/GRS/", trait_ID, "/", sex, "/", GRS_threshold, "/", sex, "_GRS_", GRS_threshold, ".txt"), row.names=F, col.names=T, quote=F)
      }

    }

    #write.table(male_IV_data,paste0(pheno_dir,"/IVs/", trait_ID, "/male_IVs.txt"), row.names=F, col.names=T, quote=F)
    #write.table(female_IV_data,paste0(pheno_dir,"/IVs/", trait_ID, "/female_IVs.txt"), row.names=F, col.names=T, quote=F)

  }


}


read_shiny_data <- function(traits, i){

  trait_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description
  pheno_dir <- paste0("analysis/traitMR/")
  trait_info <-  read.table(paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), header=F, row.names=1,check.names=F)
  male_IVs <- fread(paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/male_IVs.txt"),header=T, data.table=F)
  female_IVs <- fread(paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/female_IVs.txt"),header=T, data.table=F)
  output <- list(trait_info = trait_info, male_IVs = male_IVs, female_IVs = female_IVs)
  return(list(output))

}


find_ntile_groups <- function(ntile_in, ntile_out, n){

  join <- cbind(ntile_in, ntile_out)
  groups <- unique(ntile_out)
  out <- numeric()
  for(i in 1:n){
    min <- min(join[which(join[,2]==i),1])
    max <- max(join[which(join[,2]==i),1])
    temp <- cbind(i, paste0(min, "-", max))
    out <- rbind(out, temp)
  }

  full <- out[match(ntile_out, out[,1]),2]
  return(full)

}



### FUNCTION LIST:

#1 extract_Neale_IV
## Description: extract IV's from Neale files and clump

#2 extract_Neale_outcome
## Description: extract statistics for a given set of SNPs from Neale files

#3 extract_bgen
## Description: extract genotypes from bgen files

#4 proper_sex_pairs
## Description: returns list of rows to be included given a sex pairing (eg. male/female)

#5 full_model_summary
## Description: returns summary of model with removed rows previously removed because of NA

#6 extract_model_stats
## Description: returns single row of beta, se, pval for given set coefficients

#7 pretty_round
## Description: rounds to 2 decimal places if >1, keeps 2 significant figures if <1,
## and uses scientific notation if <0.001

#8 make_beta_95ci
## Description: returns character in the form: "beta (95ci)"

#9 include_MR_NA
## Description: include NA rows in MR output

#10 summarize_mr_result
## Description: summarizes mr_results

#11 plink_clump
## Description: runs plink clumping function in shell

#12 load_geno
## Description: loads genotype data for a list of IVs with input (i) which corresponds to row in traits file.

#13 meta
## Description: meta analyzes a list of effects/se's that are organized into a df and removes any
## "problem rows", uses 'rmeta' pacakge, outputs a list containing: effect, se, ci, pval
#####################################################################################


### extract independent set of IVs from Neale_files
### first remove any duplicated SNP and SNP not in HRC SNP list
extract_Neale_IV <- function(Neale_file, variant_file, filter_in,pval_threshold,prune_threshold){
  exposure_raw <- fread(paste0(Neale_file),data.table=F)
  exposure_cols  <- c("variant","beta","se","pval","n_complete_samples")
  exposure_raw <- exposure_raw[, exposure_cols]

  variant_data <- fread(paste0(variant_file),data.table=F)
  variant_cols <- c("rsid", "ref", "alt", "AF","chr")
  variant_data <- variant_data[,variant_cols]

  exp_var_merge <- cbind(exposure_raw,variant_data)

  data_IV_no_duplicates <- exp_var_merge[!(duplicated(exp_var_merge$rsid) | duplicated(exp_var_merge$rsid, fromLast=TRUE)),]
  data_IV_filter <- data_IV_no_duplicates[which(data_IV_no_duplicates$rsid %in% filter_in), ]

  data_IV <- data_IV_filter[which(data_IV_filter[["pval"]] < pval_threshold),]
  IV_cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "AF","n_complete_samples")
  data_IV <- data_IV[,IV_cols]

  colnames(data_IV) <- c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf","samplesize")
  data_IV_format <- format_data(data_IV, type="exposure")

  data_prune <- numeric()
  for(chr in 1:22)
  {
    data_IV_temp <- data_IV_format[which(data_IV_format[["chr.exposure"]]==chr),]
    fn <- tempfile(tmpdir = tempdir())
    snps <- data_IV_temp$SNP
    if(length(snps)==0) next
    pvals <- data_IV_temp$pval.exposure
    write.table(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F)
    refdat=paste0("/data/sgg2/jenny/data/1000G/chr",chr,"/1000G_EUR_chr",chr,"_filt")

    snp_clump <- plink_clump(refdat, fn, prune_threshold)

    data_prune_temp <- data_IV_temp[which(data_IV_temp$SNP %in% snp_clump),]
    data_prune <- rbind(data_prune, data_prune_temp)
  }
  return(data_prune)
}

extract_Neale_outcome <- function(Neale_file, variant_file, IVs){
  outcome_raw <- fread(paste0(Neale_file),data.table=F)
  outcome_cols  <- c("variant","beta","se","pval","n_complete_samples")
  outcome_raw <- outcome_raw[, outcome_cols]

  variant_data <- fread(paste0(variant_file),data.table=F)
  variant_cols <- c("rsid", "ref", "alt", "AF","chr")
  variant_data <- variant_data[,variant_cols]

  outcome_var_merge <- cbind(outcome_raw,variant_data)

  outcome_filter <- outcome_var_merge[which(outcome_var_merge$rsid %in% IVs), ]

  cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "AF", "n_complete_samples")
  outcome_filter <- outcome_filter[,cols]
  colnames(outcome_filter) <- c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf", "samplesize")
  outcome_format <- format_data(outcome_filter, type="outcome")

  return(outcome_format)
}

extract_bgen <- function(snps, file){

  num_snps <- length(snps)
  bgen_data <- bgen.load(file=file, rsids = snps)
  snp_map <- bgen_data$variants
  snp_map$AF <- NA
  geno <- numeric()

  for(i in 1:dim(snp_map)[1])
  {
    snp <- as.character(snp_map[["rsid"]][i])
    snp_temp <- bgen_data$data[snp,,]
    snp_temp <- as.data.frame(snp_temp)
    snp_temp[,4] <- snp_temp[,2]+2*snp_temp[,3]
    snp_map$AF[i] <- mean(snp_temp[,4], rm.na=T)/2
    snp_out <- snp_temp[,4]
    geno <- cbind(geno, snp_out)
  }
  colnames(geno) <- bgen_data$variants$rsid
  return(list(geno_data=geno,snp_map=snp_map))
}

proper_sex_pairs <- function(pheno,exposure_sex,outcome_sex,index,opp_index){
  sex_exp_index <- rep(FALSE,dim(pheno)[1])
  if(exposure_sex=="female"){sex_exp_index[which(pheno[[paste0(opp_index, "_sex")]]==0)] <- TRUE}
  if(exposure_sex=="male"){sex_exp_index[which(pheno[[paste0(opp_index, "_sex")]]==1)] <- TRUE}
  if(exposure_sex=="both_sexes"){sex_exp_index[] <- TRUE}

  sex_outcome_index <- rep(FALSE,dim(pheno)[1])
  if(outcome_sex=="female"){sex_outcome_index[which(pheno[[paste0(index, "_sex")]]==0)] <- TRUE}
  if(outcome_sex=="male"){sex_outcome_index[which(pheno[[paste0(index, "_sex")]]==1)] <- TRUE}
  if(outcome_sex=="both_sexes"){sex_outcome_index[] <- TRUE}

  sex_joint_index <- sex_exp_index & sex_outcome_index
  return(sex_joint_index)
}

full_model_summary <- function(model){

  nr <- length(which(is.na(coef(model))))
  nc <- 4
  rnames <- names(which(summary(model)$aliased))
  cnames <- colnames(summary(model)$coefficients)
  mat_na <- matrix(data = NA,nrow = nr,ncol = nc,
                   dimnames = list(rnames,cnames))
  mat_coef <- rbind(summary(model)$coefficients,mat_na)
  return(mat_coef)

}

extract_model_stats <- function(model_summary, coef_list){
  summary <- numeric()
  for(coef in coef_list)
  {
    beta <- model_summary[coef,1]
    se <- model_summary[coef,2]
    pval <- model_summary[coef,4]
    temp <- cbind(beta, se, pval)
    colnames(temp) <- paste0(coef, c("_beta", "_se","_pval"))
    summary <- cbind (summary,  temp)
  }
  return(summary)

}

pretty_round <- function(x){
  x <- as.numeric(x)
  if(!is.na(x))
  {
    if(abs(x)>1) {out<- sprintf("%.2f", round(x,2))}
    if(abs(x)<=1 & abs(x)>=0.001) {
      om = floor(log10(abs(x)))
      dp = 2-om-1
      out<- sprintf(paste("%.",dp,"f", sep=""), signif(x,2))}
    if(abs(x)<0.001) {out <- sprintf("%.1e", signif(x, 2))}
  } else out <- NA
  return(out)
}

make_beta_95ci <- function(beta, se){
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  beta_round <- pretty_round(beta)
  lower_ci <- pretty_round(beta-(1.96*se))
  upper_ci <- pretty_round(beta+(1.96*se))
  out <- paste0(beta_round," (",lower_ci,", ", upper_ci, ")")
  return(out)
}

include_MR_NA <- function(MR_res){
  out <- MR_res
  used_methods <- MR_res[,"method"]
  method_col <- which(colnames(MR_res)=="method")
  total_col <- dim(MR_res)[2]
  attempted_methods <- levels(used_methods)[which(!levels(used_methods) %in%  used_methods)]
  for(i in attempted_methods)
  {
    empty_row <- c(rep(NA,method_col-1), i, rep(NA,total_col-method_col))
    out <- rbind(out, empty_row)
  }
  return(out)
}

summarize_mr_result <- function(exposure, outcome, exposure_sex, outcome_sex, nsnps, MR_res, het_test, egger_test,n_exposure,n_outcome){
  out <- cbind(exposure, outcome, exposure_sex, outcome_sex, nsnps, n_exposure,n_outcome)
  if(nsnps==1)
  {
    egger_pval <- NA
    het_IVW_pval <- NA
    for(method in c("Wald ratio","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      if(any(is.na(MR_res[row,c("b","se","pval")])))
      {
        out <- cbind(out, NA, NA)
      } else {
        b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
        pval <- pretty_round(MR_res[row,"pval"])
        out <- cbind(out, b.ci, pval)
      }
    }
    out <- cbind(out, egger_pval,het_IVW_pval)
  }
  if(nsnps>1)
  {

    for(method in c("Inverse variance weighted","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
      pval <- pretty_round(MR_res[row,"pval"])
      out <- cbind(out, b.ci, pval)
    }
    het_IVW_pval <- pretty_round(het_test[which(het_test[,"method"]=="Inverse variance weighted"),"Q_pval"])
    egger_pval <- pretty_round(egger_test[1,"pval"])
    out <- cbind(out, egger_pval,het_IVW_pval)
  }
  colnames(out) <- c("exposure", "outcome", "exposure_sex","outcome_sex","N_snps", "N_exposure_GWAS", "N_outcome_GWAS",
                     "IVW/Wald_summary", "IVW/Wald_pval","Egger_summary", "Egger_pval", "Egger_int_pval", "Het_IVW_Qpval")
  return(out)
}

plink_clump <- function(bfile, fn, prune_threshold){

  fun2 <- paste0(
    "plink",
    " --bfile ", bfile,
    " --clump ", fn,
    " --clump-p1 ", clump_p1,
    " --clump-p2 ", clump_p2,
    " --clump-r2 ", prune_threshold,
    " --clump-kb ", clump_kb,
    " --out ", fn
  )
  system(fun2)
  a <- read.table(paste(fn, ".clumped", sep=""), header =T)
  unlink(paste(fn, "*", sep=""))
  a_out <- a[,"SNP"]
  return(a_out)
}

load_geno <- function(snp_data, sample_file, UKBB_imp_dir){ # UKBB_dir should be added to variable list

  ####### load IV list
  # IV_data <- fread(IV_data_file, header=T, data.table=F) ## IV list is the same for males and females, but the effects will be different

  ####### Extract SNP froms bgen files
  snp_map <- numeric()
  geno_data <- numeric()
  for(chr in 1:22)
  {
    snps <- snp_data[which(snp_data[,"chr"]==chr),"rsid"]
    if(length(snps)==0) next
    bgen_file=paste(UKBB_imp_dir, "/_001_ukb_imp_chr", chr, "_v2.bgen", sep="")
    bgen_chr_extract <- extract_bgen(snps,bgen_file)
    geno_data <- cbind(geno_data, bgen_chr_extract$geno_data)
    snp_map <- rbind(snp_map,bgen_chr_extract$snp_map)
    cat(paste0("Finished loading chr: ", chr, ".\n"))

  }
  row.names(geno_data) <- sample_file[,1]
  return(list(geno_data = geno_data, snp_map = snp_map))
}



meta <- function(result_data, remove_rows, effect_col, se_col,n){

  if(length(remove_rows)!=0)
  {
    result_data <- result_data[-remove_rows,]
  }

  effects <- result_data[,effect_col]
  ses <- result_data[,se_col]
  meta.result=meta.summaries(d=effects, se=ses,method=c("fixed"), conf.level=0.95)
  b_meta=round(meta.result[3]$summary,digits=10)
  b_meta_se=round(meta.result[4]$se.summary,digits=10)
  lowerbound=b_meta-b_meta_se*1.96
  upperbound=b_meta+b_meta_se*1.96
  meta_p=round(2*(pt(abs(b_meta/b_meta_se),((n)-meta.result$het[2]),lower.tail=FALSE)),digits=10)

  return(list(effect=b_meta, se=b_meta_se,
              lower_ci=lowerbound, upper_ci=upperbound,
              pval=meta_p))

}

interval_labels <- function(intervals){

  labels <- numeric()
  for(j in 1:length(intervals))
  {
    if(j!=length(intervals))
    {
      label_j <- paste0(intervals[j],"-",(intervals[j+1]-1))
    }
    if(j==length(intervals))
    {
      label_j <- paste0(intervals[j],"+")
    }
    labels <- c(labels, label_j)

  }
  return(labels)
}

define_models <- function(traits){

  i_col <- numeric()
  exposure_sex_col <- numeric()
  trait_ID_col <- numeric()
  phenotype_file <- numeric()
  phenotype <- numeric()
  phenotype_description <- numeric()
  IV_file <- numeric()
  output_file <- numeric()
  for(i in 1:dim(traits)[1]){

    trait_ID <- as.character(traits[i,"Neale_pheno_ID"])
    pheno_dir <- paste0("analysis/traitMR/" )
    trait_info <-  read.table(paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), header=F, row.names=1,check.names=F)
    phesant_ID <- as.character(trait_info["phes_ID",1])

    for(exposure_sex in c("male", "female")){
      i_col <- c(i_col,i)
      trait_ID_col <- c(trait_ID_col, trait_ID)
      if(exposure_sex=="male"){outcome_sex="female"}
      if(exposure_sex=="female"){outcome_sex="male"}
      exposure_sex_col <- c(exposure_sex_col, exposure_sex)
      IV_file <- c(IV_file, paste0(pheno_dir,"/IVs/Neale/",trait_ID,"/", exposure_sex,"_IVs.txt"))

      phenotype_file <- c(phenotype_file,paste0(pheno_dir,"/pheno_files/phesant/", trait_ID, "_", outcome_sex,".txt"))
      phenotype = c(phenotype,phesant_ID)
      phenotype_description = c(phenotype_description,"phesant")
      output_file = c(output_file, paste0(pheno_dir, "/household_GWAS/", trait_ID, "/outcome_",outcome_sex,"/", "phesant","_" ,exposure_sex, "-",outcome_sex, "_GWAS.csv"))

    }
  }
  output <- tibble(
    i = i_col,
    trait_ID = trait_ID_col,
    exposure_sex = exposure_sex_col,
    phenotype_file = phenotype_file,
    phenotype_col = phenotype,
    phenotype_description = phenotype_description,
    IV_file = IV_file,
    gwas_outcome_file = output_file)
  return(output)

}

summarize_gwas <- function(geno, outcome, covar){


  gwas <- big_univLinReg(as_FBM(geno), y.train = outcome,
                         covar.train = covar_from_df(covar))

  gwas$pval <- predict(gwas, log10 = FALSE)
  gwas$n <- dim(geno)[1]
  gwas$group_AF <- colMeans(geno)/2

  return(gwas)

}

run_household_GWAS <- function(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_update,
                               IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge){

  exposure_ID <- as.character(exposure_info["trait_ID",1])
  cat(paste0("Calculating household GWAS for phenotype `", exposure_ID, "` as exposure and phenotype `", outcome_ID, "` as outcome...\n"))
  pheno_dir <- paste0("analysis/traitMR/")

  IV_data <- summ_stats # IVs are same for males and females
  household_time <- household_time_munge
  pheno_cov <- joint_model_adjustments

  IV_geno <- IV_genetic_data[[1]]
  snp_map <- IV_genetic_data[[2]]

  phesant_ID <- as.character(exposure_info["phes_ID",1])

  outcome_traits <- traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),]
  outcome_phes_ID <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==outcome_ID), "SGG_PHESANT_ID"])

  for(exposure_sex in c("male", "female")){

    cat(paste0("Processing ", exposure_sex, "s...\n"))
    if(exposure_sex=="male"){outcome_sex="female"}
    if(exposure_sex=="female"){outcome_sex="male"}
    if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
    if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
    opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")
    IV_data_sex <- IV_data[[paste0(exposure_sex, "_IV_data")]]


    pheno_data_sex <- pheno_data[[paste0("unrelated_", outcome_sex, "_data")]] %>% dplyr::select(IID, !!outcome_phes_ID)


    household_intervals <- levels(household_time[[grouping_var]])

    genetic_IDs <- tibble(IID = as.character(rownames(IV_geno)))

    # to reduce to only genetic IDs with good genetic data
    temp1 <- merge(pheno_cov,genetic_IDs, by.x=index, by.y="IID")
    # to reudce to only IDs with phenotype data
    temp2 <- merge(temp1, pheno_data_sex, by.x=opp_index, by.y="IID")
    # to reduce to only those in a household pair
    temp3 <- merge(temp2, household_time[,c("HOUSEHOLD_MEMBER1",grouping_var)], by="HOUSEHOLD_MEMBER1")
    final_data <- temp3
    colnames(final_data) <- c(colnames(pheno_cov), "outcome", grouping_var)

    pheno_run <- final_data[,c(grep("_age", names(final_data)),grep("_PC_", names(final_data)), which(names(final_data)=="outcome"))]
    pheno_run$outcome <- scale(pheno_run$outcome)


    IIDs_keep <- as.character(format(final_data[[index]], scientific = F))
    geno_data_sub <- IV_geno[as.character(IIDs_keep),]

    template <- cbind(exposure_ID, outcome_ID, SNP = colnames(geno_data_sub), grouping_var, bin = "all", exposure_sex, outcome_sex)

    gwas_result <- summarize_gwas(geno_data_sub, pheno_run$outcome, pheno_run[,!names(pheno_run) %in% c("outcome")])

    outcome_GWAS <- cbind(template, gwas_result)

    for(bin in household_intervals)
    {

      bin_sub <- final_data[which(final_data[[grouping_var]]==bin),]
      bin_pheno_run <- bin_sub[,c(grep("_age", names(bin_sub)),grep("_PC_", names(bin_sub)), which(names(bin_sub)=="outcome"))]

      bin_pheno_run$outcome <- scale(bin_pheno_run$outcome)
      bin_pheno_run <- bin_pheno_run[complete.cases(bin_pheno_run),] #if all values are the same then NA's will be produced

      bin_IIDs_keep <- bin_sub[[index]]
      bin_geno_data_sub <- IV_geno[as.character(bin_IIDs_keep),]

      template <- cbind(exposure_ID, outcome_ID, SNP = colnames(geno_data_sub), grouping_var, bin, exposure_sex, outcome_sex)

      if(dim(bin_pheno_run)[1]!=0){

        bin_gwas_result <- summarize_gwas(bin_geno_data_sub, bin_pheno_run$outcome, bin_pheno_run[,!names(bin_pheno_run) %in% c("outcome")])

      } else bin_gwas_result <- as.data.frame(sapply(c(NA, NA, NA, NA, 0, NA), rep, 9))

      outcome_GWAS <- rbind(outcome_GWAS, cbind(template, bin_gwas_result))
    }

    colnames(outcome_GWAS)[match(c("estim", "std.err", "score", "pval", "n"), colnames(outcome_GWAS))] <-
      c("geno_index_beta", "geno_index_se", "geno_index_score", "geno_index_pval", "bin_n")

    outcome_GWAS_snp_info <- full_join(outcome_GWAS, snp_map,by = c("SNP"="rsid"))

    colnames(outcome_GWAS_snp_info)[match(c("AF"), colnames(outcome_GWAS_snp_info))] <-
      c("UKBB_AF")

    assign(paste0(exposure_sex, "_", outcome_sex, "_GWAS"), outcome_GWAS_snp_info)
  }

  output <- rbind(male_female_GWAS, female_male_GWAS)
  return(output)
}

household_GWAS_bin <- function(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_update,
                                           IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge){

  for(grouping_var in grouping_var_list){

    cat(paste0("Running regression for each GWAS in all participants and bins of household pairs for group `", grouping_var, "`...\n"))

    group_result <- run_household_GWAS(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_update,
                       IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge)

    assign(paste0(grouping_var, "_GWAS"), group_result)
  }
  output <- rbind(time_together_even_bins_GWAS, age_even_bins_GWAS)
  return(output)

}



household_GWAS_across_phenos <- function(exposure_info, summ_stats, outcomes_to_run, traits_corr2_update,
                               IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- as.character(exposure_info["trait_ID",1])
  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("Calculating household GWAS for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    cat(paste0("Loading phenotype data for phenotype `", outcome_ID, "` and performing GWAS...\n"))

    male_file <- paste0(pheno_dir,"/pheno_files/phesant/", outcome_ID, "_male.txt")
    female_file <- paste0(pheno_dir,"/pheno_files/phesant/", outcome_ID, "_female.txt")

    male_pheno_data <- fread( male_file,header=T, data.table=F)
    female_pheno_data <- fread( female_file,header=T, data.table=F)

    pheno_data <- list(unrelated_male_data = male_pheno_data, unrelated_female_data = female_pheno_data)

    outcome_result <- household_GWAS_bin(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_update,
                                       IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge)

    output_file_i <- paste0(pheno_dir, "/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")

    write.csv(outcome_result, output_file_i, row.names = F)

    #output_list[[outcome_ID]] <- outcome_result ## getting to large, slowing down the analysis
    output_files <- c(output_files, output_file_i)

    cat(paste0("Finished GWAS for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }
  #return(output_list)
  return(output_files)

}


run_household_MR <- function(trait_info, summ_stats, household_GWAS_result, grouping_var) {
  #household_GWAS_result <- "gwas_age_bins"
  #grouping_var <- "age_even_bins"

  #number <- household_GWAS_result[[1]]$outcome_gwas_out
  trait_ID <- as.character(trait_info["trait_ID",1])

  pheno_dir <- paste0("analysis/traitMR/")
  print(trait_ID)
  # trait_info <-  read.table(paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), header=F, row.names=1,check.names=F)
  trait_description <- as.character(trait_info["description",1])
  phesant_ID <- as.character(trait_info["phes_ID",1])

  for(exposure_sex in c("male", "female")){
    cat(paste0("Processing ", exposure_sex, "s...\n"))


  if(exposure_sex=="male"){outcome_sex="female"}
  if(exposure_sex=="female"){outcome_sex="male"}
  if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
  if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
  opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")
  ####### load IV list


  IV_data <- summ_stats[[paste0(exposure_sex, "_IV_data")]]
  colnames(IV_data) <- c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf","samplesize")
  data_IV_format <- format_data(IV_data, type="exposure")
  outcome <- "phenotype" #phenotype_description
  #j <- which(models_to_run[["exposure_sex"]]==exposure_sex & models_to_run[["trait_ID"]]==trait_ID) # j represents the model
  #outcome_gwas <- get(household_GWAS_result)[[j]]$outcome_gwas_out
  outcome_gwas <- household_GWAS_result[[1]]$outcome_gwas_out
  household_intervals <- levels(outcome_gwas[[grouping_var]])

  outcome_gwas <- outcome_gwas %>% mutate_if(is.factor,as.character) %>%
    mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)
  outcome_dat_full <- format_data(outcome_gwas, type="outcome",
                                  snp_col = "SNP",
                                  beta_col = paste0("geno_index_beta"),
                                  se_col = paste0("geno_index_se"),
                                  effect_allele_col = "allele1",
                                  other_allele_col = "allele0",
                                  pval_col = paste0("geno_index_pval"),
                                  eaf_col = "AF",
                                  phenotype_col = grouping_var
  )
  harmonise_dat_full <- harmonise_data(
    exposure_dat = data_IV_format,
    outcome_dat = outcome_dat_full, action=1
  )

  full_MR_summary <-  numeric()
  bin_summary <- numeric()

  household_intervals_num <- strex::str_first_number(household_intervals)

  household_intervals_num[which(is.na(household_intervals_num))] <- 1e3
  for(bin in household_intervals[order(household_intervals_num)]){

    gwas_sub <- outcome_gwas[which(outcome_gwas[[grouping_var]]==bin),] %>% mutate_if(is.factor,as.character) %>%
      mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)

    n_bin <- max(gwas_sub$n)
    outcome_dat <- format_data(gwas_sub, type="outcome",
                               snp_col = "SNP",
                               beta_col = paste0("geno_index_beta"),
                               se_col = paste0("geno_index_se"),
                               effect_allele_col = "allele1",
                               other_allele_col = "allele0",
                               pval_col = paste0("geno_index_pval"),
                               eaf_col = "AF"
    )

    harmonise_dat <- harmonise_data(
      exposure_dat = data_IV_format,
      outcome_dat = outcome_dat, action=1
    )


    original_MR <- mr(harmonise_dat, method=MR_method_list)
    if(dim(original_MR)[1]!=0){
      MR_res <- include_MR_NA(original_MR)
      MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
      MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")
      #MR_dir <- paste0(pheno_dir, "/household_MR/exposure_", exposure_sex)
      #write.csv(MR_res, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR_bin_",bin, ".csv"), row.names=F)

      bin_result <- c(bin, outcome, n_bin, exposure_sex, outcome_sex, MR_res[MR_ivw_row,"b"],MR_res[MR_ivw_row,"se"], MR_res[MR_ivw_row,"pval"])
    } else bin_result <- c(bin, outcome, n_bin, exposure_sex, outcome_sex, NA,NA, NA)
    bin_summary <- rbind(bin_summary, bin_result)



    if(bin=="all"){

      het_test <- mr_heterogeneity(harmonise_dat, method_list=c("mr_egger_regression", "mr_ivw"))
      egger_test <- mr_pleiotropy_test(harmonise_dat)
      leave1out_test <- mr_leaveoneout(harmonise_dat)

      #write.csv(het_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-het.csv"), row.names=F)
      #write.csv(egger_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-Egger.csv"), row.names=F)
      #write.csv(leave1out_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-leave_out.csv"), row.names=F)

      gwas_sub <- outcome_gwas[which(outcome_gwas[[grouping_var]]==bin),] %>% mutate_if(is.factor,as.character) %>%
        mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)
      ngwas <- max(gwas_sub$n, na.rm=T)
      nsnps <- dim(harmonise_dat)[1]
      n_neale <- max(data_IV_format$samplesize.exposure, na.rm=T)
      temp_summary <- summarize_mr_result (paste0(trait_ID,"_INDEX"), paste0(trait_ID,"_HOUSEHOLD"), exposure_sex, outcome_sex, nsnps,MR_res, het_test, egger_test,n_neale, ngwas)

      ## Test for reverse causality
      harmonise_dat_sensitivity <- harmonise_dat[which(harmonise_dat[["pval.exposure"]] < harmonise_dat[["pval.outcome"]]),]
      MR_res <- include_MR_NA(mr(harmonise_dat_sensitivity, method=MR_method_list))
      #write.csv(MR_res, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-sensitivity.csv"), row.names=F)

      MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
      MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")
      nsnps_sensitivity <- dim(harmonise_dat_sensitivity)[1]
      correct_row <- ifelse(nsnps_sensitivity==1, MR_wald_row, MR_ivw_row)
      sensitivity_result <- cbind(nsnps_sensitivity, make_beta_95ci(MR_res[correct_row,"b"],MR_res[correct_row,"se"]),pretty_round(MR_res[correct_row,"pval"]))
      temp_summary <- cbind(temp_summary, sensitivity_result)

      ## Make MR plot
      #output_figure_dir <- paste0(project_dir, "/output/figures/traitMR/",trait_ID)
      #pdf(file=paste0(output_figure_dir,"/", trait_ID, "_", exposure_sex, "-",outcome_sex, "_household_MR", ".pdf"))
      mr_title <- bquote(atop(.(paste0("Estimate of the assortative mating effect of ")),
                              italic(.(trait_description)) ~ .(paste0(' for ', exposure_sex, "s on ", outcome_sex,"s"))))
      mr_plot <- mr_scatter_plot_custom(original_MR, harmonise_dat, mr_title, exposure_sex, outcome_sex)
      #dev.off()

      MR_summary <- cbind(outcome, temp_summary)
      full_MR_summary <- rbind(full_MR_summary, MR_summary)
    }
  }

  num_cols <- length(colnames(full_MR_summary))

  colnames(full_MR_summary)[num_cols-1] <- "IVW/Wald_summary_sensitivity"
  colnames(full_MR_summary)[num_cols] <- "IVW/Wald_pval_sensitivity"
  colnames(full_MR_summary)[num_cols-2] <- "N_snps_sensitivity" #test for reverse causation
  #output_table_dir <- paste0(project_dir, "/output/tables/traitMR/", trait_ID)
  #write.csv(full_MR_summary, paste0(output_table_dir, "/", trait_ID, "_household_MR.csv"), row.names=F, quote=T)

  colnames(bin_summary) <- c("bin","outcome","n", "exposure_sex","outcome_sex", "IVW_beta", "IVW_se", "IVW_pval")
  #write.csv(bin_summary, paste0(pheno_dir, "/household_MR/", phenotype_description,"_", exposure_sex, "-",outcome_sex, "_household_MR_bin.csv"), row.names=F, quote=T)
  }

  out <- list(harmonise_dat_full = harmonise_dat_full, bin_summary = bin_summary, full_MR_summary = full_MR_summary, mr_plot = mr_plot, leave1out_test = leave1out_test)
  return(list(out))

}

calc_Q_stat <- function(household_MR_result, trait_ID){

  bin_result_temp <- household_MR_result[[1]]$bin_summary
  full_MR_temp <- household_MR_result[[1]]$full_MR_summary
  full_MR_temp <- full_MR_temp[,-1]
  bin_result_temp <- as.data.frame(bin_result_temp)

  bin_result_temp <- bin_result_temp %>% mutate_if(is.factor,as.character) %>%
    mutate_at(vars(IVW_beta, IVW_se, IVW_pval, n), as.numeric)

  all_row <- which(bin_result_temp$bin=="all")
  all_sum <- c(bin_result_temp[6,"IVW_beta"],bin_result_temp[6,"IVW_se"])
  names(all_sum) <- c("MR_est", "MR_se")

  bin_result_temp <- bin_result_temp %>%
    dplyr::filter(bin != "all")
  meta_bin <- metagen(TE = as.numeric(IVW_beta), seTE = as.numeric(IVW_se), studlab = bin, data = bin_result_temp)

  Q_stat <- meta_bin$Q
  Q_pval <- meta_bin$pval.Q
  Q_df <- meta_bin$df.Q
  Q_sum <- c(Q_stat, Q_pval, Q_df)
  names(Q_sum) <- c("Q_stat", "Q_pval", "Q_df")


  out <- as.data.frame(t(c(trait_ID, full_MR_temp, all_sum, Q_sum))) %>%
    mutate_if(is.factor,as.character) %>%
    as_tibble()

  colnames(out)[1] <- "trait_ID"
  return(out)

}

mr_sex_het <- function(calc_Q_stat_age, calc_Q_stat_tt, traits_summary){
  mr_summary <- calc_Q_stat_age
  mr_summary2 <- calc_Q_stat_tt

  traits <- unique(mr_summary$trait_ID)

  mr_summary$trait_description <- NA

  mr_summary$sex_het <- NA
  mr_summary$MR_meta_est <- NA
  mr_summary$MR_meta_se <- NA
  mr_summary$MR_meta_L95 <- NA
  mr_summary$MR_meta_U95 <- NA
  mr_summary$MR_meta_pval <- NA

  for(trait in traits){

    trait_description <- as.character(traits_summary[which(traits_summary$Neale_pheno_ID==trait), "description"])

    male_row <- which(mr_summary$trait_ID==trait & mr_summary$exposure_sex=="male")
    female_row <- which(mr_summary$trait_ID==trait & mr_summary$exposure_sex=="female")

    beta_F <-  as.numeric(mr_summary[female_row,"MR_est"])
    beta_M <-  as.numeric(mr_summary[male_row,"MR_est"])
    SE_F <- as.numeric(mr_summary[female_row,"MR_se"])
    SE_M <- as.numeric(mr_summary[male_row,"MR_se"])
    n_M <- as.numeric(mr_summary[male_row,"N_outcome_GWAS"])
    n_F <- as.numeric(mr_summary[male_row,"N_outcome_GWAS"])

    se <- sqrt( (SE_F^2) + (SE_M^2) )
    t <- (beta_F-beta_M)/se
    p_het <- 2*pnorm(-abs(t))

    mr_summary$sex_het[male_row] <- p_het
    mr_summary$sex_het[female_row] <- p_het

    effects <- c(beta_F, beta_M)
    ses <- c(SE_F, SE_M)
    n <- sum(n_M, n_F)
    meta.result=meta.summaries(d=effects, se=ses,method=c("fixed"), conf.level=0.95)
    b_meta=round(meta.result[3]$summary,digits=10)
    b_meta_se=round(meta.result[4]$se.summary,digits=10)
    lowerbound=b_meta-b_meta_se*1.96
    upperbound=b_meta+b_meta_se*1.96
    meta_p=round(2*(pt(abs(b_meta/b_meta_se),((n)-meta.result$het[2]),lower.tail=FALSE)),digits=10)

    for(sex_row in c(male_row, female_row)){
      mr_summary$trait_description[sex_row] <- trait_description
      mr_summary$MR_meta_est[sex_row] <- b_meta
      mr_summary$MR_meta_se[sex_row] <- b_meta_se
      mr_summary$MR_meta_L95[sex_row] <- lowerbound
      mr_summary$MR_meta_U95[sex_row] <- upperbound
      mr_summary$MR_meta_pval[sex_row] <- meta_p

    }
  }

  mr_summary <- mr_summary %>%
    rename_at(vars(starts_with("Q_")),function(x) paste0("median_age_", x))

  mr_summary2 <- mr_summary2 %>% dplyr::select(c(trait_ID, exposure_sex, starts_with("Q_"))) %>%
    rename_at(vars(starts_with("Q_")),function(x) paste0("time_together_", x))

  out <- full_join(mr_summary, mr_summary2, by = c("trait_ID", "exposure_sex"))

  return(out)

}


#saveRDS(out, "code/shiny/mr_summary.rds")

mr_scatter_plot_custom <-  function (mr_results, dat, mr_title, exposure_sex, outcome_sex){
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"),
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] *
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] *
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] &
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure,
                                                       d$beta.outcome, d$se.exposure, d$se.outcome,
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure,
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome,
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure,
                                                                y = beta.outcome)) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome -se.outcome, ymax = beta.outcome + se.outcome),
                                                  colour = "grey", width = 0) +
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
                           ggplot2::geom_point() + #ggplot2::aes(text = paste("SNP:", SNP))) +
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
                           theme_minimal() + labs(title=mr_title) + theme(plot.title = element_text(hjust = 0.5)) +
                           ggplot2::scale_colour_manual(values = c("#a6cee3",
                                                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                                                                   "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test",
                                                                                                                     x = paste0("SNP effect on ", exposure_sex, "s (general population)"), y = paste("SNP effect on",
                                                                                                                                                                                                     outcome_sex, "partner")) + ggplot2::theme(legend.position = "bottom",
                                                                                                                                                                                                                                               legend.direction = "horizontal") + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
                           theme(legend.title = element_blank())
                       })
  mrres
}

MR_MLE <- function(sqc_munge, i, trait_ID, exposure_sex, phenotype_file, phenotype_col, phenotype_description,IV_file,sample_file,pheno_cov,grouping_var,household_time,output){

  cat(paste0("Running MLE MR for phenotype: ", trait_ID, "...\n"))
  cat(paste0("Loading genotype data...\n"))
  pheno_dir <- paste0("analysis/traitMR/")
  IV_geno_out <- load_geno(sample_file, pheno_dir, IV_file)
  IV_geno <- IV_geno_out[[1]]
  snp_map <- IV_geno_out[[2]]

  trait_info <- read.table(paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), header=F, row.names=1,check.names=F)
  phesant_ID <- as.character(trait_info["phes_ID",1])

  cat(paste0("Processing ", exposure_sex, "s...\n"))
  if(exposure_sex=="male"){outcome_sex="female"}
  if(exposure_sex=="female"){outcome_sex="male"}
  if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
  if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
  opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")

  ####### load IV list
  IV_data <- fread(IV_file, header=T, data.table=F)

  ####### load phenotypes: GRS and raw
  phenotype <- fread( gsub(paste0("_", outcome_sex, ".txt"), paste0("_", exposure_sex, ".txt"),
                           phenotype_file),header=T, data.table=F, select=c("IID","age",phenotype_col))

  phenotype_partner <- fread( phenotype_file,header=T, data.table=F, select=c("IID","age",phenotype_col))
  colnames(phenotype) <- c("IID", "age", "outcome")
  ukbb_pheno <- merge(phenotype, sqc_munge, by.x="IID", by.y="ID")

  ukbb_pheno[["outcome"]] <- scale(ukbb_pheno[["outcome"]])

  # temp1 <- merge(pheno_cov,ukbb_pheno, by.x=opp_index, by.y="IID")
  #
  # temp1 <- merge(pheno_cov,geno_sub, by.x=index, by.y="IID")
  # temp2 <- merge(temp1, phenotype, by.x=opp_index, by.y="IID")
  # temp3 <- merge(temp2, household_time[,c("HOUSEHOLD_MEMBER1",grouping_var)], by="HOUSEHOLD_MEMBER1")
  #

  ukbb_pheno[["outcome_resid"]] <- scale(resid(lm(outcome ~ I(age^2) + . - IID, data=ukbb_pheno)))

  mle_geno_out <- numeric()

  ##########################################################
  ## Run the likelihood functions ##########################
  ##########################################################

  ############################################################
  # define the likelihood function for use in entire cohorot #

  LL <- function(y0,gam, C_pam, sigma) {
    #-sum(dnorm(model_data$outcome-y0 -gam * model_data$PRS*exp(model_data$age*C_pam), 0,sigma, log=TRUE))

    -sum(dnorm(outcome_resid, y0+ gam * geno*exp(-age*C_pam), sigma, log=TRUE))
  }


  model_data <- ukbb_pheno
  for( k in 1:dim(IV_data)[1]) #mle_geno_out <- foreach(k = 1:dim(IV_data)[1],.combine=rbind) %dopar% #
  {
    snp <- IV_data[["rsid"]][k]
    col_extract <- which(colnames(IV_geno)==snp)
    geno_sub <- as.data.frame(IV_geno[,col_extract])
    colnames(geno_sub)[1] <- "geno"
    geno_sub$IID <- as.numeric(row.names(geno_sub))

    ############################################################
    model_data_geno <- merge(model_data,geno_sub, by.x="IID", by.y="IID")
    #model_data_geno <- merge(model_data,geno_sub, by.x="IID", by.y="IID")
    model_data_geno <- model_data_geno %>%
      mutate_at(c("outcome"), scale) %>%
      mutate_at(c("geno","outcome_resid"), scale)

    cor_geno_init <- cor(model_data_geno[["outcome_resid"]], model_data_geno[["geno"]], method = c("pearson"),use = "pairwise.complete.obs")
    if(cor_geno_init < 0){
      model_data_geno$geno <- 2-model_data_geno$geno
    }

    cor_geno <- cor(model_data_geno[["outcome_resid"]], model_data_geno[["geno"]], method = c("pearson"),use = "pairwise.complete.obs")

    ## null model
    l0 <- suppressWarnings(mle2(LL,start = list(y0=0, gam=cor_geno,sigma=1), fixed=list(C_pam=0),
                                data=model_data_geno, method="Nelder-Mead"))

    ## model fitted with "C-pam"
    ## use the fitted correlation in the model with "C-pam"
    l0_gam <- summary(l0)@coef["gam",c(1,2,4)]

    l1 <- suppressWarnings(mle2(LL,start = list(y0=0, gam=l0_gam[1],C_pam=0,sigma=1),
                                data=model_data_geno, method="Nelder-Mead"))

    ## Fitted C values (se based on Hessian matrix) for l1
    # sometimes LRT returns a negative value for 'diff_log_geno' so the results pval is 1 and the SD cannot be calculated,
    # so in these cases we use the Hessian result instead:

    # extract estimate, se, and pval
    l1_gam <- suppressWarnings(summary(l1)@coef["gam",c(1,2,4)])
    l1_C <- suppressWarnings(summary(l1)@coef["C_pam",c(1,2,4)])

    # calculate the p-value
    diff_C <-  2*(stats4::logLik(l1)[1]- stats4::logLik(l0)[1])
    LRT_C_pval <- pchisq(diff_C, df = 1,  lower.tail=FALSE)

    # calculate the SE on "C-pam"
    deg_free<- dim(model_data_geno)[1] -1
    T_C <- abs(qt(LRT_C_pval/2, deg_free))
    LRT_C_se <- abs(l1_C[1]/T_C)

    cat(paste0("Finished estimating C for ", k, " of ", dim(IV_data)[1], " SNPs.\n"))
    row_temp <- cbind(snp, cor_geno_init, l0_gam[1],
                      l1_gam[1],l1_gam[2],l1_gam[3],
                      l1_C[1],l1_C[2],l1_C[3],
                      diff_C,LRT_C_se,LRT_C_pval)

    mle_geno_out <- rbind(mle_geno_out, row_temp) #comment this line out if using %dopar%
  }


  colnames(mle_geno_out) <- c("SNP", "correlation","gamma0",
                              "gamma", "gamma_se", "gamma_p",
                              "C", "C_se", "C_p",
                              "C_LRT_diff","C_LRT_se", "C_LRT_p")
  mle_geno_out <- as.data.frame(mle_geno_out)

  indx <- (sapply(mle_geno_out, is.factor ) | sapply(mle_geno_out, is.character )) & colnames(mle_geno_out)!="SNP"
  mle_geno_out[indx] <- lapply(mle_geno_out[indx], function(x) as.numeric(as.character(x)))
  #### meta-anlayze C-parameters

  mle_geno_out[,"C_se_for_meta"] <- mle_geno_out[,"C_LRT_se"]
  LRT_probs <- which(mle_geno_out[,"C_LRT_p"]==1 )
  mle_geno_out[LRT_probs,"C_se_for_meta"] <-  mle_geno_out[LRT_probs,"C_se"]
  both_probs <-  which(mle_geno_out[,"C_se"]=="NaN" & mle_geno_out[,"C_LRT_p"]==1)

  l_n=dim(model_data)[1]
  l_meta <- meta(mle_geno_out,both_probs,"C","C_se_for_meta",l_n)

  c_meta <- l_meta$effect
  c_pass_sig <- ifelse(l_meta$p > 0.05/2/dim(traits)[1], FALSE, TRUE)
  if(c_pass_sig){c_meta <- 0}


  #########################################################################################
  # define the likelihood function for use in couple pairs while estimating convergence ###

  LL_am <- function(y0,alpha,sigma , rho, pi_pam) {

    -sum(dnorm(outcome_partner_resid,y0+ gam * exp(-c_meta*index_age)* alpha
               * geno*exp(-(rho*index_age_meeting + pi_pam*time_together_raw))
               , sigma, log=TRUE))
  }

  models_cols <- c(index, opp_index, "outcome_partner_resid","geno", paste0(index,"_age"),"time_together_raw", "index_age_meeting")
  mle_geno_pairs_out <- numeric()

  for(k in 1:dim(IV_data)[1])
  {
    snp <- IV_data[["rsid"]][k]
    col_extract <- which(colnames(IV_geno)==snp)
    geno_sub <- as.data.frame(IV_geno[,col_extract])
    colnames(geno_sub)[1] <- "geno"
    geno_sub$IID <- as.numeric(row.names(geno_sub))
    model_data_geno <- merge(household_pheno,geno_sub, by.x=index, by.y="IID")

    cor_geno <- mle_geno_out[which(mle_geno_out[,"SNP"]==snp),"correlation"]
    if(cor_geno < 0){
      model_data_geno$geno <- 2-model_data_geno$geno
    }

    model_data_geno <- model_data_geno %>%
      filter(!is.na(time_together_raw) & !is.na(index_age_meeting)) %>%
      mutate_at(c("geno","outcome_partner_resid"), scale2)

    gam <- mle_geno_out[which(mle_geno_out[,"SNP"]==snp),"gamma"]
    gam <- ifelse(c_pass_sig, mle_geno_out[which(mle_geno_out[,"SNP"]==snp),"gamma"], mle_geno_out[which(mle_geno_out[,"SNP"]==snp),"gamma0"])

    model_data_geno <- model_data_geno[,models_cols]
    colnames(model_data_geno)[which(colnames(model_data_geno)==paste0(index,"_age"))] <- "index_age"

    ## null model
    m0 <- suppressWarnings(mle2(LL_am,start = list(y0=0,sigma=1), fixed=list(alpha=0,rho=0, pi_pam=0),
                                data=model_data_geno, method="Nelder-Mead"))

    ## model fitted with alpha
    m1 <- suppressWarnings(mle2(LL_am,start = list(y0=0,sigma=1,alpha=0), fixed=list(rho=0, pi_pam=0),
                                data=model_data_geno, method="Nelder-Mead"))
    ## Fitted alpha values (se based on Hessian matrix) for m1
    m1_alpha <- suppressWarnings(summary(m1)@coef["alpha",c(1,2,4)])

    ## model fitted with alpha and rho
    m2 <- suppressWarnings(mle2(LL_am,start = list(y0=0,sigma=1, alpha=m1_alpha[1], rho=0), fixed=list(pi_pam=0),
                                data=model_data_geno, method="Nelder-Mead"))
    ## Fitted alpha values and rho parameters (se based on Hessian matrix) for m2
    m2_rho <- suppressWarnings(summary(m2)@coef["rho",c(1,2,4)])
    m2_alpha <- suppressWarnings(summary(m2)@coef["alpha",c(1,2,4)])

    ## model fitted with all parameters
    m3 <- suppressWarnings(mle2(LL_am,start = list(y0=0,sigma=1,alpha=m2_alpha[1], rho=m2_rho[1],pi_pam=0),
                                data=model_data_geno, method="Nelder-Mead"))
    m3_pi <- suppressWarnings(summary(m3)@coef["pi_pam",c(1,2,4)])

    # calculate the LRT based pvalue's and se's
    # combine individual SNP results
    LRT_result <- numeric()
    row_temp <- cbind(snp)
    count=0
    deg_free <- dim(model_data_geno)[1] -1
    for(pam in c("alpha", "rho", "pi_pam"))
    {
      diff_log <-  2*(stats4::logLik(get(paste0("m",count+1)))[1]- stats4::logLik(get(paste0("m",count)))[1])
      LRT_pval <- pchisq(diff_log, df = 1,  lower.tail=FALSE)
      T_stat <- abs(qt(LRT_pval/2, deg_free))
      estimate <- suppressWarnings(summary(get(paste0("m",count+1)))@coef[pam,"Estimate"])
      hessian_se <- suppressWarnings(summary(get(paste0("m",count+1)))@coef[pam,"Std. Error"])
      hessian_p <- suppressWarnings(summary(get(paste0("m",count+1)))@coef[pam,"Pr(z)"])
      LRT_se <- abs(estimate/T_stat)
      count <- count+1
      LRT_result <- rbind(LRT_result, cbind(pam, estimate, diff_log, LRT_pval, LRT_se))
      row_temp <- cbind(row_temp,estimate,hessian_se,hessian_p,
                        diff_log,LRT_se,LRT_pval)
    }
    row.names(LRT_result) <- c("alpha", "rho", "pi_pam")

    # sometimes LRT returns a negative value for 'diff_log'
    # so the resulting pval is 1 and the SD cannot be calculated,
    # so in these cases we use the Hessian result instead

    cat(paste0("Finished estimating assortative mating and convergence effects for ", k, " of ", dim(IV_data)[1], " SNPs.\n"))

    mle_geno_pairs_out <- rbind(mle_geno_pairs_out, row_temp)
  }





}
