

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


  existing_files <- numeric()
  for(path in paste0( output_path, "/", c("both_sexes", "male", "female")))
  {

      existing_files_temp  =  list.files( path,
                                        recursive = TRUE,
                                        pattern = '[.]gz' ) %>%
      str_match( '[^/]+$' ) %>%
      c

    existing_files <- c(existing_files, existing_files_temp)
  }


  to_download  =  read_tsv(reference_file) %>%
    filter(str_detect(`Phenotype Code`, phenotype_ids)) %>%
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
  return(as.data.frame(trait_corr))
}


## got below function from: https://github.com/MRCIEU/PHESANT/blob/master/WAS/testContinuous.r
irnt_phesant <- function(pheno) {
  set.seed(1234)
  numPhenos = length(which(!is.na(pheno)))
  quantilePheno = (rank(pheno, na.last="keep", ties.method="random")-0.5)/numPhenos
  phenoIRNT = qnorm(quantilePheno)
  return(phenoIRNT);
}



ivt <- function(x){
  out <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(out)
}


compute_pc_corr <- function(sqc_munge, pairs_filter, data_id_sex){


  euro_phes_data <- merge(sqc_munge, data_id_sex, by.x = "ID", by.y = "userId")

  sqc_munge <- euro_phes_data %>% mutate_at(vars(starts_with("PC_")), ivt)
  PCs <- colnames(sqc_munge)[-which(colnames(sqc_munge) %in% c("ID", "sex"))]
  trait_corr <- numeric()
  ID_sub <- NA
  for (PC in PCs)
  {
    id <- as.character(PC)
    r2 <- NA
    pairs_filter_copy <- pairs_filter
    tsv_data <- sqc_munge %>% dplyr::select(ID, !!PC)
    pairs_filter_copy$trait1 <- tsv_data[[PC]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER1"]], tsv_data$ID)]
    pairs_filter_copy$trait2 <- tsv_data[[PC]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER2"]], tsv_data$ID)]
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

    trait_row <- cbind("22009", "22009", PC, n_completed_pairs,r2)
    trait_corr <- rbind(trait_corr, trait_row)
  }
  colnames(trait_corr) <- c("ID", "ID_sub", "description", "N_pairs","r2")
  cat(paste0("Household PC correlations successfully computed.\n"))
  return(as.data.frame(trait_corr))


}


## File A2 ----

## Description: filter phenotype correlations for > specified threshold (in settings).
## Check which phenotypes have been downloaded, and create list of need to be downloaded
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
  cat(paste0("There are: ", length(unique(both_sex_avail)), " phenotypes with Neale summary stats with joint
  and sex-specific data that have also been processed in the SGG database.\n\n")) ###1278 traits
  Neale_SGG_dir_filt2 <- subset(Neale_SGG_dir_filt, Neale_SGG_dir_filt$phenotype %in% both_sex_avail)
  return(Neale_SGG_dir_filt2)
}

filter_by_corr <- function(traits_corr,Neale_SGG_dir_filt2,household_correlation_threshold){

  merge_temp <- merge(traits_corr, Neale_SGG_dir_filt2, by.x="ID", by.y="sgg_phesant_name", fill=T)
  merge_temp$r2 <- as.numeric(as.character(merge_temp$r2))
  corr_traits <- merge_temp[which(sqrt(merge_temp$r2) > household_correlation_threshold),]
  corr_traits <- corr_traits[,-which(colnames(corr_traits) %in% c("description.x"))]
  colnames(corr_traits) <- c("SGG_PHESANT_ID", "SGG_PHESANT_ID_sub","N_pairs","r2",
                             "Neale_pheno_ID","Neale_pheno_ID_sub", "Neale_file_sex","description",
                             "variable_type", "SGG_request", "PHESANT_processed", "SGG_location",
                             "Neale_downloaded","category","Neale_file_location", "define_category:T/F",
                             "v2_exists", "v2_downloaded", "variant_flag_updated")
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

    cat(paste0("Some Neale files need to be categorized and downloaded, please categorize appropriately within file:
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

download_Neale <- function(filled_cats,download_rest,traits_corr_filter, reference_file_name, Neale_output_dir){
  full_dl_list <- rbind(filled_cats, download_rest)

  for(category
      in c("body", "brain", "diet", "disease", "disease_proxy", "lifestyle", "parental_pheno"))
  {
    IDs <- as.character(full_dl_list[which(full_dl_list[,"category"]==category),"Neale_pheno_ID"])
    if(!length(IDs)==0){
      download_neale_files( IDs,
                            category = category , reference_file=reference_file_name)
    }

  }

  corr_traits <- traits_corr_filter
  levels(corr_traits$category) <- union (levels(corr_traits$category), levels(filled_cats$category))
  index <- which(corr_traits[["define_category:T/F"]]==TRUE)
  pheno_missing_cat <- corr_traits[index,"Neale_pheno_ID"]
  corr_traits[index, "category"] <- as.character(filled_cats[["category"]][match(pheno_missing_cat, filled_cats$Neale_pheno_ID)])



  existing_files_full <- numeric()
  existing_files <- numeric()

  for(path in paste0( Neale_output_dir, "/", c("both_sexes", "male", "female")))
  {

    existing_files_full_temp  =  list.files( path,
                                             recursive = TRUE,
                                             pattern = '[.]gz', full.names = TRUE )

    existing_files_temp  =  list.files( path,
                                        recursive = TRUE,
                                        pattern = '[.]gz' ) %>%
      str_match( '[^/]+$' ) %>%
      c

    existing_files_full <- c(existing_files_full, existing_files_full_temp)
    existing_files <- c(existing_files, existing_files_temp)
  }

  for(i in index){

    Neale_ID <- corr_traits[["Neale_pheno_ID"]][i]
    full_path <- existing_files_full[grepl(Neale_ID, existing_files_full )]
    Neale_paths <- paste(full_path, collapse=";")
    corr_traits[i, "Neale_file_location"] <- Neale_paths
    if(!is.null(Neale_paths)){
      corr_traits[i, "Neale_downloaded"] <- "YES"
      corr_traits[i, "define_category:T/F"] <- "FALSE"

      if(grepl("v2", Neale_paths )){
        corr_traits[i, "v2_exists"] <- "TRUE"
        corr_traits[i, "v2_downloaded"] <- "TRUE"

        ## restrict to only v2 files
        for(sex_cat in c("male", "female", "both_sexes")){
          sex_specific_paths <- full_path[grepl(paste0("[.]", sex_cat, "[.]"), full_path)]
          if(length(grep("v2", sex_specific_paths))==1){
            file_to_remove <- sex_specific_paths[!grepl("v2", sex_specific_paths)]
            full_path <- full_path[-which(full_path==file_to_remove)]
          }

        }
        Neale_paths <- paste(full_path, collapse=";")


      }
    }

  }

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

file_in_out <- function(traits_corr2_filled,i,reference_file,IV_threshold, Neale_summary_dir){

  corr_traits_both <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]

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

  corr_traits_both <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"]
  out_file <- paste0("analysis/data_setup/IV_lists/", Neale_id, "_IVs_", IV_threshold,"_both_sexes.txt")

  return(list(in_file = folder, out_file = out_file))

}

pull_traits_to_count_IVs <- function(traits_corr2_filled){
  output <- tibble(Neale_pheno_ID =traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),"Neale_pheno_ID"])
  return(output)
}

pull_traits_to_calc_het <- function(traits_corr3_run){
  output <- tibble(Neale_pheno_ID =traits_corr3_run[,"Neale_pheno_ID"])
  return(output)
}

pull_traits_to_run <- function(traits_final){


  both_sexes_original_Neale_file <- numeric()
  female_original_Neale_file <- numeric()
  male_original_Neale_file <- numeric()


  for(i in 1:dim(traits_final)[1]){

    Neale_files <- traits_final[i,"Neale_file_location"]
    Neale_file_list <- unlist(str_split(Neale_files, ";"))
    male_original_Neale_file <- c(male_original_Neale_file, Neale_file_list[grepl(paste0("\\.", "male", "\\."), Neale_file_list)])
    female_original_Neale_file <- c(female_original_Neale_file, Neale_file_list[grepl(paste0("\\.", "female", "\\."), Neale_file_list)])
    both_sexes_original_Neale_file <- c(both_sexes_original_Neale_file, Neale_file_list[grepl(paste0("\\.", "both_sexes", "\\."), Neale_file_list)])

  }

  output <- tibble(Neale_pheno_ID =traits_final[,"Neale_pheno_ID"],
                   both_sexes_original_Neale_file = both_sexes_original_Neale_file,
                   male_original_Neale_file = male_original_Neale_file,
                   female_original_Neale_file = female_original_Neale_file)

  return(output)
}

pull_IV_indices_to_run <- function(traits_final, traits_to_calc_het){

  indicies <- which(traits_to_calc_het$Neale_pheno_ID %in% traits_final$Neale_pheno_ID)
  return(indicies)
}

get_IV_clump_folder <- function(Neale_file, Neale_summary_dir){

  replace_ending <- str_replace(Neale_file, ".tsv.gz", ".IVs")
  basefolder <- basename(replace_ending)
  IV_clump_path <- file.path(Neale_summary_dir, "IVs", "clump", basefolder)
  return(IV_clump_path)
}


get_IV_list <- function(corr_traits, Neale_pheno_ID, IV_threshold, Neale_summary_dir){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]

  i <- which(corr_traits_both[["Neale_pheno_ID"]]==Neale_pheno_ID)
  category <- corr_traits_both[i,"category"]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"] #same as Neale_pheno_ID
  ID <- corr_traits_both[i,"Neale_pheno_ID"]
  Neale_files <- corr_traits_both[i,"Neale_file_location"]
  Neale_file_list <- unlist(str_split(Neale_files, ";"))
  Neale_file_both_sexes <- Neale_file_list[grepl("both_sexes", Neale_file_list)]

  IV_folder <- get_IV_clump_folder(Neale_file_both_sexes, Neale_summary_dir)
  name <- basename(IV_folder)

  threshold <- IV_threshold ## only extract IVs
  build_data <- numeric()
  for(chr in 1:22)

  {
    original_data <- read.table(paste0(IV_folder,"/", name,"_unpruned_chr", chr, ".txt"), header=T)

    if(file.exists(paste0(IV_folder,"/", name,"_unpruned_chr", chr, ".clumped")))
    {
      clump_data <- read.table(paste0(IV_folder,"/", name,"_unpruned_chr", chr, ".clumped"),header = T)
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
  reduced_data$original_row <- SNP_rows
  return(reduced_data)

}

# calc_sex_het
summarize_IV_data <- function(traits, Neale_pheno_ID, variant_data, Neale_summary_dir, IV_threshold){

  trait_ID <- Neale_pheno_ID

  IV_list_both_sexes <- fread(paste0( "analysis/data_setup/IV_lists/", trait_ID, "_IVs_", IV_threshold,"_both_sexes.txt"), data.table=F, header=F)
  SNP_rows <- which(variant_data[,"rsid"] %in% IV_list_both_sexes[,1])

  result <- NA

  i <- which(traits[["Neale_pheno_ID"]]==trait_ID)
  category <- traits[i,"category"]
  Neale_id <- traits[i,"Neale_pheno_ID"] #same as Neale_pheno_ID
  ID <- traits[i,"Neale_pheno_ID"]
  Neale_files <- traits[i,"Neale_file_location"]
  Neale_file_list <- unlist(str_split(Neale_files, ";"))


  for(exposure_sex in c("both_sexes", "male", "female"))
  {

    Neale_file_exposure_sex <- Neale_file_list[grepl(paste0("\\.", exposure_sex, "\\."), Neale_file_list)]
    IV_folder_exposure_sex <- get_IV_clump_folder(Neale_file_exposure_sex, Neale_summary_dir)
    assign(paste0(exposure_sex, "_IV_folder"), IV_folder_exposure_sex)
    assign(paste0(exposure_sex, "_original_Neale_file"), Neale_file_exposure_sex)

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
  char_row <- data.frame(lapply(traits[i,], as.character), stringsAsFactors=FALSE)
  sex_het_summary <- as.data.frame(t(unlist(c(char_row, num_pass_filter))))
  output <- list(male_IV_data = male_IV_data, female_IV_data = female_IV_data, IV_list_both_sexes = IV_list_both_sexes, sex_het_summary = sex_het_summary)
  return(output)

}

write_IV_list <- function(traits_corr2_filled, Neale_pheno_ID, IV_list, IV_threshold, dir) {

  corr_traits_both <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]

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

non_diet_filter <- function(traits){ # filter out only binary diet traits?
  output <- traits %>% filter(!category=="diet") %>% filter(!grepl("intake", description))
  return(output)
}

non_qual_filter <- function(traits){
  output <- traits %>% filter(!grepl("Qualifications:", description))
  return(output)
}

non_left_side_filter <- function(traits){
  output <- traits %>% filter(!grepl("left", description))
  return(output)
}

non_redundant_filter <- function(traits){
  temp <- non_qual_filter(traits)
  output <- non_left_side_filter(temp)
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

  dir.create(paste0(pheno_dir, "/standard_MR/", trait_ID), showWarnings = FALSE)
  dir.create(paste0(pheno_dir, "/standard_GWAS/", trait_ID), showWarnings = FALSE)

  dir.create(paste0(pheno_dir, "/IVs/Neale/", trait_ID), showWarnings = FALSE)

}

create_trait_rev_filt_dirs <- function(Neale_pheno_ID){

  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description

  print(trait_ID)

  pheno_dir <- paste0("analysis/traitMR")

  dir.create(paste0(pheno_dir, "/household_GWAS_rev_filter/", trait_ID), showWarnings = FALSE)
  dir.create(paste0(pheno_dir, "/standard_GWAS_rev_filter/", trait_ID), showWarnings = FALSE)

}

create_MR_dirs <- function(Neale_pheno_ID){

  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description

  print(trait_ID)
   pheno_dir <- paste0("analysis/traitMR")

  ### create_relevant directories
  dir.create(paste0(pheno_dir, "/household_MR/", trait_ID, "/univariate_MR"), showWarnings = FALSE)
  dir.create(paste0(pheno_dir, "/household_MR/", trait_ID, "/multivariate_MR"), showWarnings = FALSE)

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

  cat(paste0("Successfully prepared phenotype for processing (created directories, phenotype file inputs, trait summary, etc).\n"))

  pheno_data = list(unrelated_male_data = pheno_male_full, unrelated_female_data = pheno_female_full)
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


  Neale_files <- traits[i,"Neale_file_location"]
  Neale_file_list <- unlist(str_split(Neale_files, ";"))


  for(exposure_sex in c("both_sexes", "male", "female"))
  {

    Neale_file_exposure_sex <- Neale_file_list[grepl(paste0("\\.", exposure_sex, "\\."), Neale_file_list)]
    IV_folder_exposure_sex <- get_IV_clump_folder(Neale_file_exposure_sex, Neale_summary_dir)
    assign(paste0(exposure_sex, "_IV_folder"), IV_folder_exposure_sex)
    assign(paste0(exposure_sex, "_original_Neale_file"), Neale_file_exposure_sex)

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


  trait_info <- as.data.frame(trait_info)
  trait_info <- tibble::rownames_to_column(trait_info, "Value") %>% rename(Info = V1)
  return(trait_info)

}

extract_trait_info <- function(pheno_data){

  output <- pheno_data$trait_info
  return(output)

}

calc_corr_mat_traits <- function(outcomes_to_run, path_pheno_data){

  full_df <- numeric()

  for(i in 1:dim(outcomes_to_run)[1]){

    Neale_pheno_ID <- outcomes_to_run$Neale_pheno_ID[i]
    phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

    male_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_male.txt")))]
    female_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_female.txt")))]
    male_pheno_data <- fread(male_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))
    female_pheno_data <- fread(female_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))

    joint_data <- rbind(male_pheno_data, female_pheno_data)

    if(i==1){
      full_df <- joint_data
    } else full_df <- merge(full_df, joint_data[,c(1, 4)], all = TRUE)

    print(i)
  }


  calc_PC_data <- full_df[,-c(1:3)]

  cor_matrix <- cor(calc_PC_data, use = "pairwise.complete.obs")

  return(cor_matrix)


}

calc_corr_mat_traits_all <- function(traits_all, path_pheno_data, sqc, fam, relatives){

  full_df <- numeric()
  traits <- traits_all %>% filter(Neale_file_sex=="both")

  for(i in 1:dim(traits)[1]){

    category <- as.character(traits[i,"category"])
    Neale_pheno_ID <- as.character(traits[i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description

    pheno_data <- prep_pheno_data(traits_all, Neale_pheno_ID, sqc, fam, relatives)

    joint_data <- rbind(pheno_data[[1]], pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

    if(i==1){
      full_df <- joint_data
    } else full_df <- merge(full_df, joint_data[,c(1, 4)], all = TRUE)

    print(i)
  }

  calc_PC_data <- full_df[,-c(1:3)]

  cor_matrix <- cor(calc_PC_data, use = "pairwise.complete.obs")

  return(cor_matrix)

}

calc_PC_traits <- function(outcomes_to_run, path_pheno_data){

  cor_matrix <- calc_corr_mat_traits (outcomes_to_run, path_pheno_data)

  res.pca <- prcomp((cor_matrix), scale = TRUE)
  # to visualize the PC results
  # fviz_eig(res.pca)

  return(res.pca)

}

calc_num_tests_by_PCs <- function(prcomp_result, threshold){

  num_tests <- which(summary(prcomp_result)$importance[3,] > threshold)[1]
  return(num_tests[[1]])
}

calc_pc_trait_corr <- function(Neale_pheno_ID, pheno_data){

  full_data <- rbind(pheno_data[[1]], pheno_data[[2]])
  PCs <- colnames(full_data)[which(startsWith(colnames(full_data), "PC_"))]
  result <- numeric()
  for (PC in PCs)
  {

    phes_ID <- gsub("_irnt", "", Neale_pheno_ID)
    full_data[["PC_ivt"]] <- ivt(full_data[[PC]])
    # cor_joint <- cor(full_data[[phes_ID]], full_data[[PC]])
    cor_joint <- cor(full_data[[phes_ID]], full_data[["PC_ivt"]])

    output_row <- cbind(Neale_pheno_ID, PC, cor_joint)
    for(sex_data in names(pheno_data)){

      dat <- pheno_data[[sex_data]]
      cor_sex <- cor(dat[[phes_ID]], dat[[PC]])
      output_row <- cbind(output_row, cor_sex)
    }

    result <- rbind(result, output_row)
  }

  colnames(result) <- c("Neale_pheno_ID", "PC", "corr_both_sexes", "corr_male", "corr_female")

  result <- result %>% as_tibble() %>% type_convert() %>% mutate_at(c("Neale_pheno_ID"), as.character)
  return(result)
}

calc_corr_impact_by_PCs <- function(traits_corr, PCs_corr, PC_trait_corr,
                                    path_pheno_data){

  result <- numeric()

  Neale_pheno_ID <- PC_trait_corr$Neale_pheno_ID[1]
  phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

  ###########################################################
  ## Get the number of pairs in the correlation between the
  ## trait `Neale_pheno_ID` and PCs
  ## to use in the SE calculation
  ###########################################################

  male_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_male.txt")))]
  female_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_female.txt")))]
  male_pheno_data <- fread(male_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))
  female_pheno_data <- fread(female_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))

  joint_data <- rbind(male_pheno_data, female_pheno_data)

  trait_PC_n_pairs <- dim(joint_data)[1] # represents number of pairs in the correlation between the PCs and trait of interest

  ###########################################################
  ## Couple correlations for trait of interest (and se and sample size)
  ###########################################################

  trait_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"r2"])
  trait_couple_r <- sqrt(trait_couple_r2)

  trait_couple_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"N_pairs"])
  trait_couple_r_se <- sqrt((1-trait_couple_r^2)/(trait_couple_n_pairs-2))


  PCs <- PCs_corr$description

  for(PC in PCs){

    ###########################################################
    ## UKB correlation between for trait of interest and PC
    ## And calculate SE on R^2 term (formula from Zoltan)
    ## We need SE on R^2 term because it is R^2 that is used in the `correlation_due_to_confounding` formula
    ###########################################################

    trait_PC_r <- PC_trait_corr[which(PC_trait_corr$PC==PC),"corr_both_sexes"][[1]]
    trait_PC_r2_se <- sqrt(4*trait_PC_r^2*(1-trait_PC_r^4)/trait_PC_n_pairs)

    ###########################################################
    ## Couple correlations for PCs
    ## (and se and sample size)
    ###########################################################

    PC_couple_r2 <- as.numeric(PCs_corr[which(PCs_corr$description==PC),"r2"])
    PC_couple_r <- sqrt(PC_couple_r2)
    PC_couple_n_pairs <- as.numeric(PCs_corr[which(PCs_corr$description==PC),"N_pairs"])
    PC_couple_r_se <- sqrt((1-PC_couple_r^2)/(PC_couple_n_pairs-2))

    ###########################################################
    ## Calculate the correlation due to confounding
    ## and the corresponding SE
    ###########################################################

    correlation_due_to_confounding <- trait_PC_r^2 *PC_couple_r
    correlation_due_to_confounding_se <- sqrt(variance_of_product(trait_PC_r^2, trait_PC_r2_se, PC_couple_r, PC_couple_r_se))

    ## Bind results together
    temp_row <- cbind(Neale_pheno_ID, PC, trait_couple_r, trait_couple_r_se, PC_couple_r, trait_PC_r, correlation_due_to_confounding, correlation_due_to_confounding_se)
    result <- rbind(result, temp_row)
  }

  colnames(result) <- c("Neale_pheno_ID", "PC", "trait_couple_corr", "trait_couple_corr_se", "PC_couple_corr", "trait_PC_corr",
                        "corr_due_to_confounding", "corr_due_to_confounding_se")

  result <- result %>% as_tibble() %>% type_convert() %>% mutate_at(c("Neale_pheno_ID"), as.character)

  # sum up along PCs
  # remove home location and birth place coordinates

  temp <- result %>% group_by(Neale_pheno_ID) %>% mutate(corr_due_to_confounding_all = sum(corr_due_to_confounding)) %>%
    mutate(corr_due_to_confounding_all_se = sqrt(sum(corr_due_to_confounding_se^2))) %>% ### CHECK THIS
    filter(Neale_pheno_ID != "130_irnt") %>% filter(Neale_pheno_ID != "129_irnt") %>%
    filter(Neale_pheno_ID != "22702_irnt") %>% filter(Neale_pheno_ID != "22704_irnt")

  return(temp)

}

calc_corr_impact_by_coords <- function(outcomes_to_run, traits_corr, corr_mat_traits,
                                       traits_all, path_pheno_data, sqc, fam, relatives){

  result <- numeric()

  NC_pheno_data <- prep_pheno_data(traits_all, "129_irnt", sqc, fam, relatives)
  NC_joint_data <- rbind(NC_pheno_data[[1]], NC_pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

  EC_pheno_data <- prep_pheno_data(traits_all, "130_irnt", sqc, fam, relatives)
  EC_joint_data <- rbind(EC_pheno_data[[1]], EC_pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

  for(i in 1:dim(outcomes_to_run)[1]){

    Neale_pheno_ID <- outcomes_to_run$Neale_pheno_ID[i]
    phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

    ###########################################################
    ## Couple correlations for trait_i (and se and sample size)
    ###########################################################

    trait_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"r2"])
    trait_couple_r <- sqrt(trait_couple_r2)

    trait_couple_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"N_pairs"])
    trait_couple_r_se <- sqrt((1-trait_couple_r^2)/(trait_couple_n_pairs-2))

    ###########################################################
    ## Get the number of pairs in the correlation between the
    ## trait `i` and N/E coordinates
    ## to use in the SE calculation
    ###########################################################

    trait_pheno_data <- prep_pheno_data(traits_all, Neale_pheno_ID, sqc, fam, relatives)
    trait_joint_data <- rbind(trait_pheno_data[[1]], trait_pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

    full_df_NC <- merge(NC_joint_data, trait_joint_data[,c(1, 4)], all = TRUE)
    full_df_EC <- merge(EC_joint_data, trait_joint_data[,c(1, 4)], all = TRUE)

    full_df_complete_NC <- full_df_NC[complete.cases(full_df_NC), ]
    full_df_complete_EC <- full_df_EC[complete.cases(full_df_EC), ]

    trait_NC_n_pairs <- dim(full_df_complete_NC)[1] # represents number of pairs in the correlation between the NC and trait `i`
    trait_EC_n_pairs <- dim(full_df_complete_EC)[1] # represents number of pairs in the correlation between the EC and trait `i`

    ###########################################################
    ## UKB correlation between for trait_i and N/E coordinataes
    ## And calculate SE on R^2 term (formula from Zoltan)
    ## We need SE on R^2 term because it is R^2 that is used in the `correlation_due_to_confounding` formula
    ###########################################################

    trait_NC_r <- corr_mat_traits[phes_ID, "129"]
    trait_EC_r <- corr_mat_traits[phes_ID, "130"]

    trait_NC_r2_se <- sqrt(4*trait_NC_r^2*(1-trait_NC_r^4)/trait_NC_n_pairs)
    trait_EC_r2_se <- sqrt(4*trait_EC_r^2*(1-trait_EC_r^4)/trait_EC_n_pairs)

    ###########################################################
    ## Couple correlations for N/E coord
    ## (and se and sample size)
    ###########################################################

    NC_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID=="129"),"r2"])
    EC_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID=="130"),"r2"])

    NC_couple_r <- sqrt(NC_couple_r2)
    EC_couple_r <- sqrt(EC_couple_r2)

    NC_couple_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID=="129"),"N_pairs"])
    EC_couple_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID=="130"),"N_pairs"])

    NC_couple_r_se <- sqrt((1-NC_couple_r^2)/(NC_couple_n_pairs-2))
    EC_couple_r_se <- sqrt((1-EC_couple_r^2)/(EC_couple_n_pairs-2))

    ###########################################################
    ## Calculate the correlation due to confounding
    ## and the corresponding SE
    ###########################################################

    correlation_due_to_confounding_NC <- trait_NC_r^2 *NC_couple_r
    correlation_due_to_confounding_EC <- trait_EC_r^2 *EC_couple_r

    correlation_due_to_confounding_NC_se <- sqrt(variance_of_product(trait_NC_r^2, trait_NC_r2_se, NC_couple_r, NC_couple_r_se))
    correlation_due_to_confounding_EC_se <- sqrt(variance_of_product(trait_EC_r^2, trait_EC_r2_se, EC_couple_r, EC_couple_r_se))


    ## Bind results together
    temp_row_NC <- cbind(Neale_pheno_ID, "129", trait_couple_r, trait_couple_r_se, NC_couple_r, trait_NC_r, correlation_due_to_confounding_NC, correlation_due_to_confounding_NC_se)
    temp_row_EC <- cbind(Neale_pheno_ID, "130", trait_couple_r, trait_couple_r_se, EC_couple_r, trait_EC_r, correlation_due_to_confounding_EC, correlation_due_to_confounding_EC_se)

    result_temp <- rbind(temp_row_NC, temp_row_EC)
    result <- rbind(result, result_temp)

  }

  colnames(result) <- c("Neale_pheno_ID", "coordinate", "trait_couple_corr", "trait_couple_corr_se", "coordinate_couple_corr", "trait_coordiante_corr", "corr_due_to_confounding", "corr_due_to_confounding_se")
  result <- result %>% as_tibble() %>% type_convert() %>% mutate_at(c("Neale_pheno_ID"), as.character)

  # sum up along coordinates
  # remove home location and birth place coordinates

  temp <- result %>% group_by(Neale_pheno_ID) %>% mutate(corr_due_to_confounding_all = sum(corr_due_to_confounding)) %>%
    mutate(corr_due_to_confounding_all_se = sqrt(sum(corr_due_to_confounding_se^2))) %>%
    filter(Neale_pheno_ID != "130_irnt") %>% filter(Neale_pheno_ID != "129_irnt") %>%
    filter(Neale_pheno_ID != "22702_irnt") %>% filter(Neale_pheno_ID != "22704_irnt")

  return(temp)
}


calc_corr_impact_by_traits <- function(outcomes_to_run, traits_corr, corr_mat_traits, trait_interest,
                                       traits_all, path_pheno_data, sqc, fam, relatives){

  result <- numeric()

  trait_interest_pheno_data <- prep_pheno_data(traits_all, trait_interest, sqc, fam, relatives)
  trait_interest_joint_data <- rbind(trait_interest_pheno_data[[1]], trait_interest_pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

  for(i in 1:dim(outcomes_to_run)[1]){

    Neale_pheno_ID <- outcomes_to_run$Neale_pheno_ID[i]
    phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

    ###########################################################
    ## Couple correlations for trait_i (and se and sample size)
    ###########################################################

    i_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"r2"])
    i_couple_r <- sqrt(i_couple_r2)

    i_couple_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"N_pairs"])
    i_couple_r_se <- sqrt((1-i_couple_r^2)/(i_couple_n_pairs-2))

    ###########################################################
    ## Get the number of pairs in the correlation between the
    ## trait of interest and trait `i`
    ## to use in the SE calculation
    ###########################################################

    i_pheno_data <- prep_pheno_data(traits_all, Neale_pheno_ID, sqc, fam, relatives)
    i_joint_data <- rbind(i_pheno_data[[1]], i_pheno_data[[2]]) %>% dplyr::select(-starts_with("PC_")) # remove PC columns

    full_df <- merge(trait_interest_joint_data, i_joint_data[,c(1, 4)], all = TRUE)
    full_df_complete <- full_df[complete.cases(full_df), ]

    trait_interest_i_n_pairs <- dim(full_df_complete)[1] # represents number of pairs in the correlation between the trait of interest and trait `i`

    ###########################################################
    ## UKB correlation between for trait_i and trait of interest
    ## And calculate SE on R^2 term (formula from Zoltan)
    ## We need SE on R^2 term because it is R^2 that is used in the `correlation_due_to_confounding` formula
    ###########################################################

    trait_interest_phes_ID <- gsub("_irnt", "", trait_interest)
    trait_interest_i_r <- corr_mat_traits[phes_ID, trait_interest_phes_ID]
    trait_interest_i_r2_se <- sqrt(4*trait_interest_i_r^2*(1-trait_interest_i_r^4)/trait_interest_i_n_pairs)

    ###########################################################
    ## Couple correlations for trait of interest
    ## (and se and sample size)
    ###########################################################

    trait_interest_couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID==trait_interest_phes_ID),"r2"])
    trait_interest_couple_r <- sqrt(trait_interest_couple_r2)

    trait_interest_n_pairs <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"N_pairs"])
    trait_interest_couple_r_se <- sqrt((1-trait_interest_couple_r^2)/(trait_interest_n_pairs-2))

    ###########################################################
    ## Calculate the correlation due to confounding
    ## and the corresponding SE
    ###########################################################

    correlation_due_to_confounding <- trait_interest_i_r^2 *trait_interest_couple_r
    correlation_due_to_confounding_se <- sqrt(variance_of_product(trait_interest_i_r^2, trait_interest_i_r2_se, trait_interest_couple_r, trait_interest_couple_r_se))

    temp_row <- cbind(Neale_pheno_ID, trait_interest, i_couple_r, i_couple_r_se, trait_interest_couple_r, trait_interest_i_r, correlation_due_to_confounding, correlation_due_to_confounding_se)
    result <- rbind(result, temp_row)

  }

  colnames(result) <- c("trait_j", "trait_interest", "trait_j_couple_corr", "trait_j_couple_corr_se", "trait_interest_couple_corr", "trait_j_trait_inter_corr", "corr_due_to_confounding", "corr_due_to_confounding_se")
  result <- result %>% as_tibble() %>% type_convert() %>% mutate_at(c("trait_j", "trait_interest"), as.character)

  temp <- result %>%
    filter(trait_j != "130_irnt") %>% filter(trait_j != "129_irnt") %>%
    filter(trait_j != "22702_irnt") %>% filter(trait_j != "22704_irnt") %>%
    filter(trait_interest!=trait_j)


  return(temp)
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
    write.table(data_out[[1]]$trait_info, paste0(pheno_dir,"/trait_info/", trait_ID, "_trait_info.txt"), row.names=F, quote=T)
  }
}

create_summary_stats <- function(Neale_pheno_ID, trait_info, IV_data_summary){

  #i <- which(traits[["Neale_pheno_ID"]]==Neale_pheno_ID)
  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description
  pheno_dir <- paste0("analysis/traitMR/")

  het_stats <- IV_data_summary$IV_list_both_sexes #fread(paste0( "analysis/data_setup/sex_heterogeneity/", trait_ID, "_sex_het.txt"), header=T, data.table=F)
  IV_list_filter <- het_stats[which(het_stats[["P-het"]] > 0.05/dim(het_stats)[1]),]

  for(sex in c("male", "female"))
  {
    #file_name <- paste0( "analysis/data_setup/IV_info/", trait_ID, "_IVs_5e-08_", file,".txt")
    temp <- IV_data_summary[[paste0(sex, "_IV_data")]] #read.table(file_name,header=T)
    SNP_rows <- which(temp[,"rsid"] %in% IV_list_filter[,1])
    temp <- temp[SNP_rows,]

    temp <- temp %>% mutate(alt_AF = case_when(minor_allele_Neale == alt ~ minor_AF_Neale,
                                               TRUE ~ 1 - minor_AF_Neale))

    IV_cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "alt_AF","n_complete_samples")
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

write_summ_stats <- function(Neale_pheno_ID, summ_stats){

  trait_ID <- Neale_pheno_ID ## this is the Neale_id, used to be pheno_description
  pheno_dir <- paste0("analysis/traitMR/")

  pheno_dir <- paste0("analysis/traitMR/")
  male_file <- paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/male_IVs.txt")
  female_file <- paste0(pheno_dir,"/IVs/Neale/", trait_ID, "/female_IVs.txt")

  write.table(summ_stats[["male_IV_data"]],male_file, row.names=F, col.names=T, quote=F)
  write.table(summ_stats[["female_IV_data"]],female_file, row.names=F, col.names=T, quote=F)

  return(c(male_file, female_file))

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

extract_relevant_variant_rows <- function(variant_file, snp_list){
  variant_data <- fread(paste0(variant_file),data.table=F)
  extract_rows <- which(variant_data$rsid %in% snp_list)
  variant_data_sub <- variant_data[extract_rows,]
  variant_data_sub$original_row <- extract_rows

  return(variant_data_sub)

}

extract_ind_Neale_file <- function(Neale_file, variant_data){

  # `variant_data` is data frame with a set of variants to extract (one row/variant) which must have at least one column `original_row`,
  # which is the row from the Neale variant file indicating which SNPs to extract

  outcome_raw <- fread(paste0(Neale_file),data.table=F)
  outcome_cols  <- c("variant","beta","se","pval","n_complete_samples", "minor_AF")


  outcome_raw_sub <- outcome_raw[variant_data$original_row,outcome_cols]

  variant_cols <- c("rsid", "ref", "alt", "chr", "minor_allele")
  variant_data_sub <- variant_data[,variant_cols]

  outcome_var_merge <- cbind(outcome_raw_sub,variant_data_sub)

  temp <- outcome_var_merge %>% mutate(alt_AF = case_when(minor_allele == alt ~ minor_AF,
                                                          TRUE ~ 1 - minor_AF))


  cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "alt_AF", "n_complete_samples")
  outcome_filter <- temp[,cols]
  colnames(outcome_filter) <- c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","alt_AF", "samplesize")
  #outcome_format <- format_data(outcome_filter, type="outcome")

  return(outcome_filter)
}

extract_Neale_outcome <- function(outcome_ID, both_sexes_file, male_file, female_file, variant_data){

  outcome_list <- list()
  count <- 0
  for(sex in c("both_sexes", "male", "female")){

    Neale_file <- get(paste0(sex, "_file"))
    if(!is.null(Neale_file)){

      Neale_stats <- extract_ind_Neale_file(Neale_file, variant_data)
      Neale_stats$outcome_ID <- outcome_ID
      Neale_stats$sex <- sex
      outcome_list[[paste0(outcome_ID, "_", sex, "_summary_stats")]] <- Neale_stats

    }

  }
  return(outcome_list)

}

write_outcome_stats <- function(exposure_info, extract_Neale_outcome_result, outcomes_to_run, summ_stats){

  pheno_dir <- paste0("analysis/traitMR/" )

  outcome_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)


  file_list <- numeric()
  for(i in 1:dim(outcomes_to_run)[1]){
    IV_stats <- summ_stats[[i]]
    exposure_ID <- outcomes_to_run$Neale_pheno_ID[i]
    rsids <- IV_stats[[1]]$rsid

    GWAS_file_i <- paste0(pheno_dir, "/standard_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
    file_list <- c(GWAS_file_i, file_list)
    outcome_result_i <- numeric()

    for(sex in c("both_sexes", "male", "female")){
      outcome_stats_sex <- extract_Neale_outcome_result[[paste0(outcome_ID, "_", sex, "_summary_stats")]]
      outcome_stats_sex_sub <- outcome_stats_sex[which(outcome_stats_sex$SNP %in% rsids),]
      if(!all(rsids %in% outcome_stats_sex_sub$SNP)) stop(paste0("Missing outcome info for some IVs for exposure `", exposure_ID, "` in outcome `", outcome_ID, "`."))
      outcome_stats_sex_sub$sex <- sex
      outcome_stats_sex_sub$exposure_ID <- exposure_ID
      outcome_result_i <- rbind(outcome_result_i, outcome_stats_sex_sub)
    }

    write.csv(outcome_result_i, GWAS_file_i, row.names = F)
    cat(paste0("Finished writing standard GWAS statistics for exposure ", i, " of ", dim(outcomes_to_run)[1], " (from Neale).\n\n" ))
  }

  return(file_list)

}

write_outcome_stats_filter <- function(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter, standard_harmonised_data_reverse_filter, summ_stats){


  pheno_dir <- paste0("analysis/traitMR/" )

  file_list <- numeric()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  snp_chr <- summ_stats[[1]] %>% dplyr::select(rsid, chr)

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[i]
    GWAS_file_i <- paste0(pheno_dir, "/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
    file_list <- c(GWAS_file_i, file_list)

    meta_GWAS <- standard_harmonised_data_meta_reverse_filter[[i]] %>% dplyr::select(SNP, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)
    male_GWAS <- standard_harmonised_data_reverse_filter[[i]][[paste0("exp_male_harmonised_data_filter")]] %>% dplyr::select(SNP, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)
    female_GWAS <- standard_harmonised_data_reverse_filter[[i]][[paste0("exp_female_harmonised_data_filter")]] %>% dplyr::select(SNP, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)

    outcome_result_i <- numeric()

    for(sex in c("meta", "male", "female")){

      dat <- get(paste0(sex, "_GWAS"))
      dat$chr <- snp_chr[match(dat$SNP, snp_chr$rsid), "chr"]
      dat <- dat %>% dplyr::select(SNP, chr, everything()) %>% mutate(outcome_ID=outcome_ID) %>% mutate(sex = !!sex) %>% mutate(exposure_ID=exposure_ID)


      colnames(dat)[match(c("beta.outcome", "se.outcome", "pval.outcome", "other_allele.outcome", "effect_allele.outcome", "eaf.outcome", "samplesize.outcome"), colnames(dat))] <-
        c("beta", "se", "pval", "other_allele", "effect_allele", "eaf", "samplesize")

      outcome_result_i <- rbind(outcome_result_i, dat)


    }


    write.csv(outcome_result_i, GWAS_file_i, row.names = F)
    cat(paste0("Finished writing standard GWAS statistics for exposure ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(file_list)
}

write_household_GWAS_filter <- function(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter, summ_stats){

  pheno_dir <- paste0("analysis/traitMR/" )

  file_list <- numeric()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  snp_chr <- summ_stats[[1]] %>% dplyr::select(rsid, chr)

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[i]
    GWAS_file_i <- paste0(pheno_dir, "/household_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
    file_list <- c(GWAS_file_i, file_list)

    meta_GWAS <- household_harmonised_data_meta_reverse_filter[[i]] %>% dplyr::select(exposure_ID, outcome_ID, SNP, grouping_var, bin, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)
    male_GWAS <- household_harmonised_data_reverse_filter[[i]][[paste0("exp_male_harmonised_data_filter")]] %>% dplyr::select(exposure_ID, outcome_ID, SNP, grouping_var, bin, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)
    female_GWAS <- household_harmonised_data_reverse_filter[[i]][[paste0("exp_female_harmonised_data_filter")]] %>% dplyr::select(exposure_ID, outcome_ID, SNP, grouping_var, bin, beta.outcome, se.outcome, pval.outcome, other_allele.outcome, effect_allele.outcome, eaf.outcome, samplesize.outcome)

    outcome_result_i <- numeric()

    for(sex in c("meta", "male", "female")){

      if(sex=="meta"){outcome_sex <- "meta"}
      if(sex=="male"){outcome_sex <- "female"}
      if(sex=="female"){outcome_sex <- "male"}
      dat <- get(paste0(sex, "_GWAS"))
      dat$chr <- snp_chr[match(dat$SNP, snp_chr$rsid), "chr"]
      dat <- dat %>% mutate(exposure_sex = !!sex) %>% mutate(outcome_sex = !!outcome_sex) %>% dplyr::select(exposure_ID, outcome_ID, SNP, chr, grouping_var, bin, exposure_sex, outcome_sex, everything())


      colnames(dat)[match(c("beta.outcome", "se.outcome", "pval.outcome", "other_allele.outcome", "effect_allele.outcome", "eaf.outcome", "samplesize.outcome"), colnames(dat))] <-
        c("beta", "se", "pval", "other_allele", "effect_allele", "eaf", "samplesize")

      outcome_result_i <- rbind(outcome_result_i, dat)


    }


    write.csv(outcome_result_i, GWAS_file_i, row.names = F)
    cat(paste0("Finished writing household GWAS statistics for exposure ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(file_list)


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
    egger_pval_round <- NA
    het_IVW_pval_round <- NA

    for(method in c("Wald ratio","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      if(any(is.na(MR_res[row,c("b","se","pval")])))
      {
        out <- cbind(out, NA, NA, NA, NA, NA)
      } else {
        b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
        pval <- pretty_round(MR_res[row,"pval"])
        out <- cbind(out, MR_res[row,"b"],MR_res[row,"se"], MR_res[row,"pval"], b.ci, pval)
      }
    }
    out <- cbind(out, egger_pval,het_IVW_pval, egger_pval_round, het_IVW_pval_round)
  }
  if(nsnps>1)
  {

    for(method in c("Inverse variance weighted","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
      pval <- pretty_round(MR_res[row,"pval"])
      out <- cbind(out, MR_res[row,"b"],MR_res[row,"se"], MR_res[row,"pval"], b.ci, pval)
    }
    het_IVW_pval <- het_test[which(het_test[,"method"]=="Inverse variance weighted"),"Q_pval"]
    egger_pval <- egger_test[1,"pval"]

    het_IVW_pval_round <- pretty_round(het_IVW_pval)
    egger_pval_round <- pretty_round(egger_pval)

    out <- cbind(out, egger_pval,het_IVW_pval, egger_pval_round, het_IVW_pval_round)

  }
  colnames(out) <- c("exposure", "outcome", "exposure_sex","outcome_sex","N_snps", "N_exposure_GWAS", "N_outcome_GWAS",
                     "IVW_beta", "IVW_se", "IVW_pval", "IVW_summary", "IVW_pval_round",
                     "Egger_beta", "Egger_se", "Egger_pval", "Egger_summary", "Egger_pval_round",
                     "Egger_int_pval", "Het_IVW_Qpval", "Egger_int_pval_round", "Het_IVW_Qpval_round")
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

summarize_gwas <- function(geno, outcome, covar, variable_type, tempdir = "/data/sgg2/jenny"){

  temp_file <- tempfile(tmpdir = tempdir)
  gwas <- suppressWarnings(big_univLinReg(as_FBM(geno, backingfile = temp_file), y.train = outcome,
                         covar.train = covar_from_df(covar)))

  unlink(paste0(temp_file, ".bk"))
  gwas$pval <- predict(gwas, log10 = FALSE)
  gwas$n <- dim(geno)[1]
  gwas$group_AF <- colMeans(geno)/2

  return(gwas)

}

household_GWAS_bin <- function(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_filled,
                               IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  variable_type <- exposure_info %>% filter(Value=="variable_type") %>% pull(Info)

  cat(paste0("Calculating household GWAS for phenotype `", exposure_ID, "` as exposure and phenotype `", outcome_ID, "` as outcome...\n"))
  pheno_dir <- paste0("analysis/traitMR/")

  IV_data <- summ_stats # IVs are same for males and females
  household_time <- household_time_munge
  pheno_cov <- joint_model_adjustments

  IV_geno <- IV_genetic_data[[1]]
  snp_map <- IV_genetic_data[[2]]

  phesant_ID <- exposure_info %>% filter(Value=="phes_ID") %>% pull(Info)

  outcome_traits <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]
  outcome_phes_ID <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==outcome_ID), "SGG_PHESANT_ID"])

  for(exposure_sex in c("male", "female")){

    cat(paste0("Processing ", exposure_sex, "s in each bin...\n"))
    if(exposure_sex=="male"){outcome_sex="female"}
    if(exposure_sex=="female"){outcome_sex="male"}
    if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
    if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
    opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")
    IV_data_sex <- IV_data[[paste0(exposure_sex, "_IV_data")]]


    pheno_data_sex <- pheno_data[[paste0("unrelated_", outcome_sex, "_data")]] %>% dplyr::select(IID, !!outcome_phes_ID)


    household_intervals <- levels(household_time[[grouping_var]])

    genetic_IDs <- tibble(IID = as.character(rownames(IV_geno)))

    ### SHOULD REALLY USE THIS FORM OF GENETIC_IDs where IIDs are numeric.
    ### male ID 4e6 does not merge with pheno_cov when IID is a character.
    ## see dim(pheno_cov) vs. dim(temp1)
    #genetic_IDs <- tibble(IID = as.numeric(rownames(IV_geno)))

    # to reduce to only genetic IDs with good genetic data
    temp1 <- merge(pheno_cov,genetic_IDs, by.x=index, by.y="IID")
    # to reudce to only IDs with phenotype data
    temp2 <- merge(temp1, pheno_data_sex, by.x=opp_index, by.y="IID")
    # to reduce to only those in a household pair
    temp3 <- merge(temp2, household_time[,c("HOUSEHOLD_MEMBER1",grouping_var)], by="HOUSEHOLD_MEMBER1")
    final_data <- temp3
    colnames(final_data) <- c(colnames(pheno_cov), "outcome", grouping_var)

    pheno_run <- final_data[,c(grep("_age", names(final_data)),grep("_PC_", names(final_data)), which(names(final_data)=="outcome"))]

    if(variable_type!="binary" | variable_type!="ordinal") {pheno_run$outcome <- scale(pheno_run$outcome)}


    IIDs_keep <- as.character(format(final_data[[index]], scientific = F))
    geno_data_sub <- IV_geno[as.character(IIDs_keep),]

    template <- cbind(exposure_ID, outcome_ID, SNP = colnames(geno_data_sub), grouping_var, bin = "all", exposure_sex, outcome_sex)

    gwas_result <- summarize_gwas(geno_data_sub, pheno_run$outcome, pheno_run[,!names(pheno_run) %in% c("outcome")])

    outcome_GWAS <- cbind(template, gwas_result)

    for(bin in household_intervals)
    {

      bin_sub <- final_data[which(final_data[[grouping_var]]==bin),]
      bin_pheno_run <- bin_sub[,c(grep("_age", names(bin_sub)),grep("_PC_", names(bin_sub)), which(names(bin_sub)=="outcome"))]

      if(variable_type!="binary") {pheno_run$outcome <- scale(pheno_run$outcome)}
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

household_GWAS_group <- function(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_filled,
                                           IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge){

  for(grouping_var in grouping_var_list){

    cat(paste0("Running regression for each GWAS in all participants and bins of household pairs for group `", grouping_var, "`...\n"))

    group_result <- household_GWAS_bin(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_filled,
                       IV_genetic_data, joint_model_adjustments, grouping_var, household_time_munge)

    assign(paste0(grouping_var, "_GWAS"), group_result)
  }
  output <- rbind(time_together_even_bins_GWAS, age_even_bins_GWAS)
  return(output)

}


run_household_GWAS <- function(exposure_info, summ_stats, outcomes_to_run, traits_corr2_filled,
                               IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")

  cat(paste0("\nCalculating household GWAS for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    GWAS_file_i <- paste0(pheno_dir, "/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
    output_files <- c(output_files, GWAS_file_i)

    #if(file.exists(GWAS_file_i)) {
    #  cat(paste0("Skipping `", outcome_ID, "` because GWAS results already exist...\n\n"))
    #  next
    #}

    male_file <- paste0(pheno_dir,"/pheno_files/phesant/", outcome_ID, "_male.txt")
    female_file <- paste0(pheno_dir,"/pheno_files/phesant/", outcome_ID, "_female.txt")

    male_pheno_data <- fread( male_file,header=T, data.table=F)
    female_pheno_data <- fread( female_file,header=T, data.table=F)

    pheno_data <- list(unrelated_male_data = male_pheno_data, unrelated_female_data = female_pheno_data)

    outcome_result_i <- household_GWAS_group(exposure_info, summ_stats, pheno_data, outcome_ID, traits_corr2_filled,
                                         IV_genetic_data, joint_model_adjustments, grouping_var_list, household_time_munge)


    write.csv(outcome_result_i, GWAS_file_i, row.names = F)

    cat(paste0("Finished GWAS for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_files)

}

read_household_GWAS <- function(gwas_file_list){ #output from `run_household_GWAS`

  output <- list()

  for(file in gwas_file_list){
    household_GWAS_result <- fread(file)
    exposure_ID <- household_GWAS_result$exposure_ID[1]
    outcome_ID <- household_GWAS_result$outcome_ID[1]
    household_GWAS_result <- as.data.frame(household_GWAS_result)
    output[[paste0(outcome_ID, "_vs_", exposure_ID, "_household_GWAS")]] <- household_GWAS_result
  }

  return(output)

}

###########################################
### harmonise exposure and outcome data

harmonise_household_data_ind <- function(exposure_info, summ_stats, traits_corr2_filled, outcome_ID, household_GWAS_result) {

  output_list <- list()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  outcome_traits <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]
  outcome_description <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==outcome_ID), "description"])
  exposure_description <- exposure_info %>% filter(Value=="description") %>% pull(Info)

  for(exposure_sex in c("male", "female")){

    cat(paste0("Harmonising data for ", exposure_sex, "s as exposure...\n"))

    if(exposure_sex=="male"){outcome_sex="female"}
    if(exposure_sex=="female"){outcome_sex="male"}
    if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
    if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
    opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")

    ## FORMAT IV DATA
    IV_data <- summ_stats[[paste0(exposure_sex, "_IV_data")]]
    colnames(IV_data)[match(c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "alt_AF", "n_complete_samples"), colnames(IV_data))] <-
      c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf","samplesize")
    IV_data$phenotype <- exposure_ID
    data_IV_format <- format_data(IV_data, type="exposure", phenotype_col = "phenotype")

    ## FILTER OUTCOME DATA - it doesn't matter which `grouping_var` we choose because we are selecting "all"
    outcome_gwas <- household_GWAS_result %>% filter(exposure_sex == !!exposure_sex) %>%
      mutate(phenotype = paste0(grouping_var, ".", bin))

    outcome_dat <- format_data(outcome_gwas, type="outcome",
                               snp_col = "SNP",
                               beta_col = paste0("geno_index_beta"),
                               se_col = paste0("geno_index_se"),
                               effect_allele_col = "allele1",
                               other_allele_col = "allele0",
                               pval_col = paste0("geno_index_pval"),
                               eaf_col = "group_AF",
                               phenotype_col = "phenotype",
                               samplesize_col = "bin_n"
    )

    harmonise_dat <- harmonise_data(
      exposure_dat = data_IV_format,
      outcome_dat = outcome_dat, action=1
    )

    harmonise_dat$exposure_ID <- exposure_ID
    harmonise_dat$outcome_ID <- outcome_ID

    harmonise_dat$exposure_sex <- exposure_sex
    harmonise_dat$outcome_sex <- outcome_sex

    harmonise_dat$exposure_description <- exposure_description
    harmonise_dat$outcome_description <- outcome_description

    grouping_var_bin <- harmonise_dat %>% separate(outcome, c("grouping_var", "bin"), "[.]", extra = "merge") %>% dplyr::select("grouping_var", "bin")
    harmonise_dat$grouping_var <- grouping_var_bin$grouping_var
    harmonise_dat$bin <- grouping_var_bin$bin
    #output_table_dir <- paste0(project_dir, "/output/tables/traitMR/", trait_ID)
    #write.csv(full_MR_summary, paste0(output_table_dir, "/", trait_ID, "_household_MR.csv"), row.names=F, quote=T)

    output_list[[paste0("exp_", exposure_sex, "_harmonised_data")]] <- harmonise_dat

  }

  return(output_list)

}

harmonise_household_data <- function(exposure_info, summ_stats, outcomes_to_run, gwas_results, traits_corr2_filled){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nHarmonising household MR input for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    #GWAS_file_i <- gwas_files[[i]] ## should be the same as: `paste0(pheno_dir, "/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")`
    #household_GWAS_result <- fread(GWAS_file_i, sep=",", header = T, data.table=F)

    household_GWAS_result <- gwas_results[[i]]

    harmonise_result <- harmonise_household_data_ind(exposure_info, summ_stats, traits_corr2_filled, outcome_ID, household_GWAS_result)

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data")]] <- harmonise_result

    cat(paste0("Finished harmonising data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)

}

meta_harmonised_household_data <- function(exposure_info, outcomes_to_run, household_harmonised_data){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nMeta-analyzing harmonised household data across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    household_harmonised_data_i <- household_harmonised_data[[i]]
    household_harmonised_data_M <- household_harmonised_data_i[[paste0("exp_", "male", "_harmonised_data")]]
    household_harmonised_data_F <- household_harmonised_data_i[[paste0("exp_", "female", "_harmonised_data")]]
    household_harmonised_data_join <- rbind(household_harmonised_data_M, household_harmonised_data_F)

    exposure_ID <- household_harmonised_data_join$exposure_ID[1]
    outcome_ID <- household_harmonised_data_join$outcome_ID[1]
    exposure_description <- household_harmonised_data_join$exposure_description[1]
    outcome_description <- household_harmonised_data_join$outcome_description[1]


    for(data_type in c("exposure", "outcome")){

      temp <- household_harmonised_data_join %>% rename_all(~stringr::str_replace(., paste0(".", data_type),".data_type"))
      meta_temp <- temp %>% group_by(SNP, bin, grouping_var, effect_allele.data_type, other_allele.data_type, outcome) %>% group_modify(~ summarize_sex_specific_results(.x$beta.data_type, .x$se.data_type))
      meta_temp_n <- temp %>% mutate(samplesize.data_type_adj = case_when(is.na(beta.data_type) | is.na(se.data_type) ~ 0L,
                                                                           TRUE ~ samplesize.data_type)) %>% group_by(SNP, bin, grouping_var) %>%
        summarise(samplesize.data_type_meta = sum(samplesize.data_type_adj))

      meta_temp_eaf <- temp %>% group_by(SNP, bin, grouping_var) %>% mutate(samplesize.data_type_meta = sum(samplesize.data_type)) %>% ungroup() %>%
        mutate(eaf.data_type_weighted = eaf.data_type*(samplesize.data_type/samplesize.data_type_meta)) %>% group_by(SNP, bin, grouping_var) %>% summarise(eaf.data_type_meta = sum(eaf.data_type_weighted))
      meta_temp_join <- list(meta_temp, meta_temp_n, meta_temp_eaf) %>% reduce(left_join)

      if(data_type=="exposure"){
        meta_temp_join <- meta_temp_join %>% ungroup() %>% dplyr::select(-outcome, -bin, -grouping_var) %>% unique()
      }
      format_dat <- format_data(meta_temp_join, type=data_type,
                                 snp_col = "SNP",
                                 beta_col = paste0("meta_beta"),
                                 se_col = paste0("meta_se"),
                                 effect_allele_col = "effect_allele.data_type",
                                 other_allele_col = "other_allele.data_type",
                                 pval_col = paste0("meta_pval"),
                                 eaf_col = "eaf.data_type_meta",
                                 phenotype_col = "outcome",
                                 samplesize_col = "samplesize.data_type_meta"
      )

      assign(paste0(data_type, "_dat"), format_dat)

    }

    harmonise_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat, action=1
    )


    harmonise_dat$exposure_ID <- exposure_ID
    harmonise_dat$outcome_ID <- outcome_ID

    harmonise_dat$exposure_description <- exposure_description
    harmonise_dat$outcome_description <- outcome_description

    grouping_var_bin <- harmonise_dat %>% separate(outcome, c("grouping_var", "bin"), "[.]", extra = "merge") %>% dplyr::select("grouping_var", "bin")
    harmonise_dat$grouping_var <- grouping_var_bin$grouping_var
    harmonise_dat$bin <- grouping_var_bin$bin

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta")]] <- harmonise_dat
    cat(paste0("Finished meta-analyzing harmonised household data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)
}

harmonise_standard_data_ind <- function(exposure_info, summ_stats, traits_corr2_filled, outcome_ID, standard_GWAS_result) {

  output_list <- list()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  outcome_traits <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]
  outcome_description <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==outcome_ID), "description"])
  exposure_description <- exposure_info %>% filter(Value=="description") %>% pull(Info)

  for(exposure_sex in c("male", "female")){

    cat(paste0("Harmonising data for ", exposure_sex, "s as exposure...\n"))

    if(exposure_sex=="male"){outcome_sex="female"}
    if(exposure_sex=="female"){outcome_sex="male"}
    if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
    if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
    opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")

    ## FORMAT IV DATA
    IV_data <- summ_stats[[paste0(exposure_sex, "_IV_data")]]
    colnames(IV_data)[match(c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "alt_AF", "n_complete_samples"), colnames(IV_data))] <-
      c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf","samplesize")
    IV_data$phenotype <- exposure_ID
    data_IV_format <- format_data(IV_data, type="exposure", phenotype_col = "phenotype")

    ## FILTER OUTCOME DATA - it doesn't matter which `grouping_var` we choose because we are selecting "all"
    outcome_gwas <- standard_GWAS_result %>% filter(sex == !!exposure_sex)

    outcome_dat <- format_data(outcome_gwas, type="outcome",
                               snp_col = "SNP",
                               beta_col = paste0("beta"),
                               se_col = paste0("se"),
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               pval_col = paste0("pval"),
                               eaf_col = "alt_AF",
                               phenotype_col = "phenotype",
                               samplesize_col = "samplesize"
    )

    harmonise_dat <- harmonise_data(
      exposure_dat = data_IV_format,
      outcome_dat = outcome_dat, action=1
    )

    harmonise_dat$exposure_ID <- exposure_ID
    harmonise_dat$outcome_ID <- outcome_ID

    harmonise_dat$exposure_sex <- exposure_sex
    harmonise_dat$outcome_sex <- exposure_sex # standard MR is in only one sex

    harmonise_dat$exposure_description <- exposure_description
    harmonise_dat$outcome_description <- outcome_description

    #output_table_dir <- paste0(project_dir, "/output/tables/traitMR/", trait_ID)
    #write.csv(full_MR_summary, paste0(output_table_dir, "/", trait_ID, "_household_MR.csv"), row.names=F, quote=T)

    output_list[[paste0("exp_", exposure_sex, "_harmonised_data")]] <- harmonise_dat

  }

  return(output_list)

}

harmonise_standard_data <- function(exposure_info, summ_stats, outcomes_to_run, traits_corr2_filled){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nHarmonising standard MR input for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){


    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    GWAS_file_i <- paste0(pheno_dir, "/standard_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")

    standard_GWAS_result <- fread(GWAS_file_i, sep=",", header = T, data.table=F)

    harmonise_result <- harmonise_standard_data_ind(exposure_info, summ_stats, traits_corr2_filled, outcome_ID, standard_GWAS_result)

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data")]] <- harmonise_result

    cat(paste0("Finished harmonising data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)

}

meta_harmonised_standard_data <- function(exposure_info, outcomes_to_run, standard_harmonised_data){


  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nMeta-analyzing harmonised standard data across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    standard_harmonised_data_i <- standard_harmonised_data[[i]]
    standard_harmonised_data_M <- standard_harmonised_data_i[[paste0("exp_", "male", "_harmonised_data")]]
    standard_harmonised_data_F <- standard_harmonised_data_i[[paste0("exp_", "female", "_harmonised_data")]]
    standard_harmonised_data_join <- rbind(standard_harmonised_data_M, standard_harmonised_data_F)

    exposure_ID <- standard_harmonised_data_join$exposure_ID[1]
    outcome_ID <- standard_harmonised_data_join$outcome_ID[1]
    exposure_description <- standard_harmonised_data_join$exposure_description[1]
    outcome_description <- standard_harmonised_data_join$outcome_description[1]


    for(data_type in c("exposure", "outcome")){

      temp <- standard_harmonised_data_join %>% rename_all(~stringr::str_replace(., paste0(".", data_type),".data_type"))
      meta_temp <- temp %>% group_by(SNP, effect_allele.data_type, other_allele.data_type, outcome) %>% group_modify(~ summarize_sex_specific_results(.x$beta.data_type, .x$se.data_type))
      meta_temp_n <- temp %>% mutate(samplesize.data_type_adj = case_when(is.na(beta.data_type) | is.na(se.data_type) ~ 0L,
                                                                          TRUE ~ samplesize.data_type)) %>% group_by(SNP) %>%
        summarise(samplesize.data_type_meta = sum(samplesize.data_type_adj))

      meta_temp_eaf <- temp %>% group_by(SNP) %>% mutate(samplesize.data_type_meta = sum(samplesize.data_type)) %>% ungroup() %>%
        mutate(eaf.data_type_weighted = eaf.data_type*(samplesize.data_type/samplesize.data_type_meta)) %>% group_by(SNP) %>% summarise(eaf.data_type_meta = sum(eaf.data_type_weighted))
      meta_temp_join <- list(meta_temp, meta_temp_n, meta_temp_eaf) %>% reduce(left_join)

      if(data_type=="exposure"){
        meta_temp_join <- meta_temp_join %>% ungroup() %>% dplyr::select(-outcome) %>% unique()
      }
      format_dat <- format_data(meta_temp_join, type=data_type,
                                snp_col = "SNP",
                                beta_col = paste0("meta_beta"),
                                se_col = paste0("meta_se"),
                                effect_allele_col = "effect_allele.data_type",
                                other_allele_col = "other_allele.data_type",
                                pval_col = paste0("meta_pval"),
                                eaf_col = "eaf.data_type_meta",
                                phenotype_col = "outcome",
                                samplesize_col = "samplesize.data_type_meta"
      )

      assign(paste0(data_type, "_dat"), format_dat)

    }

    harmonise_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat, action=1
    )

    harmonise_dat <- steiger_filtering(harmonise_dat)

    harmonise_dat$exposure_ID <- exposure_ID
    harmonise_dat$outcome_ID <- outcome_ID

    harmonise_dat$exposure_description <- exposure_description
    harmonise_dat$outcome_description <- outcome_description

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta")]] <- harmonise_dat
    cat(paste0("Finished meta-analyzing harmonised standard data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}

filter_reverse_SNPs_standard_data <- function(exposure_info, outcomes_to_run, standard_harmonised_data_meta, reverse_MR_threshold){

  output_list <- list()
  output_files <- numeric()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nFiltering meta-anlayzed, harmonised standard SNP data for evidence of reverse causality across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  reverse_t_threshold  =  stats::qnorm(reverse_MR_threshold)


  for(i in 1:dim(outcomes_to_run)[1]){

    dat <- standard_harmonised_data_meta[[i]]
    snps_before_filter <- dim(dat)[1]
    outcome_ID <- dat$outcome_ID[1]
    dat_filter <- dat %>%
      mutate(z_score.exposure = beta.exposure / se.exposure) %>%
      mutate(z_score.outcome = beta.outcome / se.outcome) %>%
      mutate(std_beta.exposure = z_score.exposure / sqrt(samplesize.exposure)) %>%
      mutate(std_beta.outcome = z_score.outcome / sqrt(samplesize.outcome)) %>%
      mutate(std_se.exposure = 1 / sqrt(samplesize.exposure)) %>%
      mutate(std_se.outcome = 1 / sqrt(samplesize.outcome)) %>%
      mutate(std_t_stat = (abs(std_beta.exposure) - abs(std_beta.outcome)) /
               sqrt(std_se.exposure^2 + std_se.outcome^2)) %>%
      # mutate(t_stat = (abs(beta.exposure) - abs(beta.outcome)) /
      #        sqrt(se.exposure^2 + se.outcome^2)) %>%
      filter(std_t_stat > reverse_t_threshold)

    snps_after_filter <- dim(dat_filter)[1]
    dat_filter$snps_before_filter <- snps_before_filter
    dat_filter$snps_after_filter <- snps_after_filter

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta_filter")]] <- dat_filter
    cat(paste0("Finished filtering meta-analyzed harmonised standard data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)
}

check_num_SNPS_removed_reverse_filter <- function(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  snps_before_filter <- numeric()
  snps_after_filter <- numeric()
  outcome_ID <- numeric()

  for(i in 1:dim(outcomes_to_run)[1]){

    dat <- standard_harmonised_data_meta_reverse_filter[[i]]
    snps_before_filter <- c(snps_before_filter, dat$snps_before_filter[1])
    snps_after_filter <- c(snps_after_filter, dat$snps_after_filter[1])
    outcome_ID <- c(outcome_ID, dat$outcome_ID[1])

  }

  output <- tibble(
    exposure_ID = exposure_ID,
    outcome_ID = outcome_ID,
    num_snps_before_rev_MR_filter = snps_before_filter,
    num_snps_after_rev_MR_filter = snps_after_filter
  )

  return(output)

}


filter_reverse_SNPs_household_data_sex_spec <- function(exposure_info, outcomes_to_run, household_harmonised_data, standard_harmonised_data_meta_reverse_filter){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  output_list <- list()
  cat(paste0("\nFiltering sex-specific, harmonised household SNP data for evidence of reverse causality across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    sex_filter <- list()
    for(exposure_sex in c("male", "female")){

      dat <- household_harmonised_data[[i]][[paste0("exp_", exposure_sex, "_harmonised_data")]]
      outcome_ID <- dat$outcome_ID[1]

      snps_to_keep <- standard_harmonised_data_meta_reverse_filter[[i]]$SNP
      dat_filt <- dat[which(dat$SNP %in% snps_to_keep),]

      sex_filter[[paste0("exp_", exposure_sex, "_harmonised_data_filter")]] <- dat_filt



    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_filter")]] <- sex_filter

    cat(paste0("Finished filtering sex-specific harmonised household data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)
}

filter_reverse_SNPs_standard_data_sex_spec <- function(exposure_info, outcomes_to_run, standard_harmonised_data, standard_harmonised_data_meta_reverse_filter){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  output_list <- list()
  cat(paste0("\nFiltering sex-specific, harmonised standard SNP data for evidence of reverse causality across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    sex_filter <- list()
    for(exposure_sex in c("male", "female")){

      dat <- standard_harmonised_data[[i]][[paste0("exp_", exposure_sex, "_harmonised_data")]]
      outcome_ID <- dat$outcome_ID[1]

      snps_to_keep <- standard_harmonised_data_meta_reverse_filter[[i]]$SNP
      dat_filt <- dat[which(dat$SNP %in% snps_to_keep),]

      sex_filter[[paste0("exp_", exposure_sex, "_harmonised_data_filter")]] <- dat_filt

    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_filter")]] <- sex_filter

    cat(paste0("Finished filtering sex-specific standard household data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)

}


filter_reverse_SNPs_household_data <- function(exposure_info, outcomes_to_run, household_harmonised_data_meta, standard_harmonised_data_meta_reverse_filter){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  output_list <- list()
  cat(paste0("\nFiltering meta-anlayzed, harmonised household SNP data for evidence of reverse causality across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    dat <- household_harmonised_data_meta[[i]]
    outcome_ID <- dat$outcome_ID[1]

    snps_to_keep <- standard_harmonised_data_meta_reverse_filter[[i]]$SNP
    dat_filt <- dat[which(dat$SNP %in% snps_to_keep),]

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta_filter")]] <- dat_filt
    cat(paste0("Finished filtering meta-analyzed harmonised household data for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)
}


binned_household_MR_ind <- function(exposure_info, outcome_ID, household_harmonised_data, grouping_var, MR_method_list) {

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  exposure_description <- exposure_info %>% filter(Value=="description") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR/")

  for(exposure_sex in c("male", "female")){

    cat(paste0("\nProcessing ", exposure_sex, "s in each bin...\n"))

    if(exposure_sex=="male"){outcome_sex="female"}
    if(exposure_sex=="female"){outcome_sex="male"}
    if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
    if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
    opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")

    household_harmonised_data_sex <- household_harmonised_data[[paste0("exp_", exposure_sex, "_harmonised_data_filter")]] %>% filter(grouping_var == !!grouping_var)

    household_intervals <- levels(as.factor(household_harmonised_data_sex$bin))
    household_intervals_num <- strex::str_first_number(household_intervals)
    household_intervals_num[which(is.na(household_intervals_num))] <- 1e3

    full_MR_summary <-  numeric()
    bin_summary <- numeric()

    for(bin in household_intervals[order(household_intervals_num)]){

      harmonise_sub <- household_harmonised_data_sex[which(household_harmonised_data_sex[["bin"]]==bin),]

      harmonise_sub_std <- harmonise_sub %>%
        mutate(z_score.exposure = beta.exposure / se.exposure) %>%
        mutate(z_score.outcome = beta.outcome / se.outcome) %>%
        mutate(std_beta.exposure = z_score.exposure / sqrt(samplesize.exposure)) %>%
        mutate(std_beta.outcome = z_score.outcome / sqrt(samplesize.outcome)) %>%
        mutate(std_se.exposure = 1 / sqrt(samplesize.exposure)) %>%
        mutate(std_se.outcome = 1 / sqrt(samplesize.outcome)) %>%
        mutate(std_t_stat = (abs(std_beta.exposure) - abs(std_beta.outcome)) /
                 sqrt(std_se.exposure^2 + std_se.outcome^2))

      harmonise_sub_std <- harmonise_sub_std %>% dplyr::select(-beta.exposure, -se.exposure, -beta.outcome, -se.outcome) %>%
        rename(beta.exposure = std_beta.exposure) %>% rename(se.exposure = std_se.exposure) %>%
        rename(beta.outcome = std_beta.outcome) %>% rename(se.outcome = std_se.outcome)

      harmonise_sub <- harmonise_sub_std

      n_bin <- max(harmonise_sub$samplesize.outcome )
      n_exposure <- max(harmonise_sub$samplesize.exposure )

      original_MR <- mr(harmonise_sub, method=MR_method_list)

      if(dim(original_MR)[1]!=0){
        MR_res <- include_MR_NA(original_MR)
        MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
        MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")

        bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, exposure_sex, outcome_sex, MR_res[MR_ivw_row,"nsnp"], MR_res[MR_ivw_row,"b"],MR_res[MR_ivw_row,"se"], MR_res[MR_ivw_row,"pval"])
      } else bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, exposure_sex, outcome_sex, NA, NA,NA, NA)

      bin_summary <- rbind(bin_summary, bin_result)


    }


    colnames(bin_summary) <- c("grouping_var", "bin","exposure_ID", "outcome_ID","n_exposure", "n_outcome", "exposure_sex","outcome_sex", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")
    assign(paste0(exposure_sex, "_", outcome_sex, "_MR"), bin_summary)
  }

  output <- rbind(male_female_MR, female_male_MR)

  numeric_cols <- c("n_exposure", "n_outcome", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")
  output <- as_tibble(output) %>% mutate_at(numeric_cols, as.numeric )
  return(output)


}

binned_household_MR_ind_SNPmeta <- function(exposure_info, outcome_ID, household_harmonised_data, grouping_var, MR_method_list) {

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  exposure_description <- exposure_info %>% filter(Value=="description") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR/")

  household_harmonised_data_group <- household_harmonised_data %>% filter(grouping_var == !!grouping_var)

  household_intervals <- levels(as.factor(household_harmonised_data_group$bin))
  household_intervals_num <- strex::str_first_number(household_intervals)
  household_intervals_num[which(is.na(household_intervals_num))] <- 1e3

  full_MR_summary <-  numeric()
  bin_summary <- numeric()

  for(bin in household_intervals[order(household_intervals_num)]){

    harmonise_sub <- household_harmonised_data_group[which(household_harmonised_data_group[["bin"]]==bin),]

    n_bin <- max(harmonise_sub$samplesize.outcome )
    n_exposure <- max(harmonise_sub$samplesize.exposure )

    original_MR <- mr(harmonise_sub, method=MR_method_list)

    if(dim(original_MR)[1]!=0){
      MR_res <- include_MR_NA(original_MR)
      MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
      MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")

      bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, MR_res[MR_ivw_row,"nsnp"], MR_res[MR_ivw_row,"b"],MR_res[MR_ivw_row,"se"], MR_res[MR_ivw_row,"pval"])
    } else bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, NA, NA,NA, NA)

    bin_summary <- rbind(bin_summary, bin_result)


  }


    colnames(bin_summary) <- c("grouping_var", "bin","exposure_ID", "outcome_ID","n_exposure", "n_outcome", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")


  output <- bin_summary

  numeric_cols <- c("n_exposure", "n_outcome", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")
  output <- as_tibble(output) %>% mutate_at(numeric_cols, as.numeric )
  return(output)


}

binned_household_MR_ind_SNPmeta_std <- function(exposure_info, outcome_ID, household_harmonised_data, grouping_var, MR_method_list) {

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  exposure_description <- exposure_info %>% filter(Value=="description") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR/")

  household_harmonised_data_group <- household_harmonised_data %>% filter(grouping_var == !!grouping_var)

  household_intervals <- levels(as.factor(household_harmonised_data_group$bin))
  household_intervals_num <- strex::str_first_number(household_intervals)
  household_intervals_num[which(is.na(household_intervals_num))] <- 1e3

  full_MR_summary <-  numeric()
  bin_summary <- numeric()

  for(bin in household_intervals[order(household_intervals_num)]){

    harmonise_sub <- household_harmonised_data_group[which(household_harmonised_data_group[["bin"]]==bin),]

    harmonise_sub_std <- harmonise_sub %>%
      mutate(z_score.exposure = beta.exposure / se.exposure) %>%
      mutate(z_score.outcome = beta.outcome / se.outcome) %>%
      mutate(std_beta.exposure = z_score.exposure / sqrt(samplesize.exposure)) %>%
      mutate(std_beta.outcome = z_score.outcome / sqrt(samplesize.outcome)) %>%
      mutate(std_se.exposure = 1 / sqrt(samplesize.exposure)) %>%
      mutate(std_se.outcome = 1 / sqrt(samplesize.outcome)) %>%
      mutate(std_t_stat = (abs(std_beta.exposure) - abs(std_beta.outcome)) /
               sqrt(std_se.exposure^2 + std_se.outcome^2))

    harmonise_sub_std <- harmonise_sub_std %>% dplyr::select(-beta.exposure, -se.exposure, -beta.outcome, -se.outcome) %>%
      rename(beta.exposure = std_beta.exposure) %>% rename(se.exposure = std_se.exposure) %>%
      rename(beta.outcome = std_beta.outcome) %>% rename(se.outcome = std_se.outcome)

    n_bin <- max(harmonise_sub$samplesize.outcome )
    n_exposure <- max(harmonise_sub$samplesize.exposure )

    original_MR <- mr(harmonise_sub, method=MR_method_list)
    original_MR_std <- mr(harmonise_sub_std, method=MR_method_list)

    if(dim(original_MR_std)[1]!=0){
      MR_res <- include_MR_NA(original_MR_std)
      MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
      MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")

      bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, MR_res[MR_ivw_row,"nsnp"], MR_res[MR_ivw_row,"b"],MR_res[MR_ivw_row,"se"], MR_res[MR_ivw_row,"pval"])
    } else bin_result <- c(grouping_var, bin, exposure_ID, outcome_ID, n_exposure, n_bin, NA, NA,NA, NA)

    bin_summary <- rbind(bin_summary, bin_result)


  }


  colnames(bin_summary) <- c("grouping_var", "bin","exposure_ID", "outcome_ID","n_exposure", "n_outcome", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")


  output <- bin_summary

  numeric_cols <- c("n_exposure", "n_outcome", "nsnp", "IVW_beta", "IVW_se", "IVW_pval")
  output <- as_tibble(output) %>% mutate_at(numeric_cols, as.numeric )
  return(output)


}

# `run_household_MR_binned` runs the MR in individual bins (by age or time_together) and returns a table with results for each bin.
run_binned_household_MR <- function(exposure_info, outcomes_to_run, household_harmonised_data, grouping_var, MR_method_list = MR_method_list){

  output_list <- list()
  output_files <- numeric()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nCalculating household MR in age and time-together bins for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  household_MR_result <- numeric()

  for(i in 1:dim(outcomes_to_run)[1]){

    household_MR_result <- numeric()
    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    household_harmonised_data_i <- household_harmonised_data[[i]]
    for(group in grouping_var){

      cat(paste0("\nRunning MR in bins of household pairs for group `", group, "`...\n"))

      # run for each group

      group_MR_result <- binned_household_MR_ind(exposure_info, outcome_ID, household_harmonised_data_i, group, MR_method_list)
      household_MR_result <- rbind(household_MR_result, group_MR_result)
    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR")]] <- household_MR_result

    cat(paste0("Finished MR for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)

}

run_binned_household_MR_SNPmeta <- function(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, grouping_var, MR_method_list = MR_method_list){

  output_list <- list()
  output_files <- numeric()

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  pheno_dir <- paste0("analysis/traitMR")
  cat(paste0("\nCalculating household MR in meta-analyzed SNPs in age and time-together bins for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  household_MR_result <- numeric()

  for(i in 1:dim(outcomes_to_run)[1]){

    household_MR_result <- numeric()
    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    household_harmonised_data_i <- household_harmonised_data_meta_reverse_filter[[i]]
    for(group in grouping_var){

      cat(paste0("\nRunning MR in bins of household pairs for group `", group, "`...\n"))

      # run for each group

      # change function to `binned_household_MR_ind_SNPmeta` if you dont want standardized MR.

      group_MR_result <- binned_household_MR_ind_SNPmeta_std(exposure_info, outcome_ID, household_harmonised_data_i, group, MR_method_list)
      household_MR_result <- rbind(household_MR_result, group_MR_result)
    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR")]] <- household_MR_result

    cat(paste0("Finished MR for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))
  }

  return(output_list)

}

compare_mr_raw_corr <- function(exposure_info, household_MR_binned_joint_std, traits_corr){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  result_name <- paste0(exposure_ID, "_", "vs_", exposure_ID, "_MR")

  MR_result <- household_MR_binned_joint_std[[result_name]] %>% filter(grouping_var == "age_even_bins") %>% filter(bin == "all")

  phes_ID <- gsub("_irnt", "", exposure_ID)

  couple_r2 <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"r2"])
  couple_r <- sqrt(couple_r2)
  n_pairs <- as.numeric(traits_corr[which(traits_corr$ID==phes_ID),"N_pairs"])
  r_se <- sqrt((1-couple_r^2)/(n_pairs-2))
  MR_result$n_pairs <- n_pairs
  MR_result$couple_r2 <- couple_r2
  MR_result$couple_r <- couple_r
  MR_result$couple_r_se <- r_se

  MR_result$diff_z <- z.test(abs(MR_result[["IVW_beta"]]), MR_result[["IVW_se"]], abs(MR_result[["couple_r"]]), MR_result[["couple_r_se"]], alternative = "less")$statistic
  MR_result$diff_p <- z.test(abs(MR_result[["IVW_beta"]]), MR_result[["IVW_se"]], abs(MR_result[["couple_r"]]), MR_result[["couple_r_se"]], alternative = "less")$p

  MR_result$correlation_larger <- ifelse(MR_result$couple_r > MR_result$IVW_beta, TRUE, FALSE)

  MR_result <- MR_result %>% dplyr::select(-grouping_var, -bin)
  return(MR_result)

}

find_potential_trait_confounders <- function(Neale_pheno_ID, Neale_pheno_ID_corr, standard_MR_summary_SNPmeta,
                                             household_MR_summary_SNPmeta, traits_corr, num_tests_by_PCs, corr_mat_traits, corr_trait_threshold){

  # here `Neale_pheno_ID` is used as the outcome_ID, because we are trying to find potential confounders that have an impact on this `Neale_pheno_ID`.
  # `Neale_pheno_ID_corr` refers to the couple correlation for `Neale_pheno_ID`.

  phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

  pheno_i_MR <- standard_MR_summary_SNPmeta %>% filter(outcome_ID==Neale_pheno_ID)

  AM_MR <- household_MR_summary_SNPmeta %>% filter(same_trait) %>% dplyr::select(exposure_ID, IVW_beta, IVW_se) %>%
    rename(exposure_ID_AM_IVW_beta = IVW_beta) %>%
    rename(exposure_ID_AM_IVW_se = IVW_se)

  corr_mat_traits_sub <- as.data.frame(corr_mat_traits) %>% dplyr::select(!!phes_ID) %>% rownames_to_column(var = "exposure_ID") %>% as_tibble() %>%
    dplyr::rename(between_trait_correlation = !!phes_ID)

  ## calculate correlation due to confounding
  ## filter to only confounders with correlation with index trait < `corr_trait_threshold`
  ## filter to only confounders with
  output <- pheno_i_MR %>% dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, IVW_beta, IVW_se, IVW_pval) %>%
    left_join(AM_MR, by = c("exposure_ID" = "exposure_ID")) %>%
    mutate(corr_due_to_confounding = IVW_beta^2*exposure_ID_AM_IVW_beta) %>%
    mutate(corr_due_to_confounding_ratio = corr_due_to_confounding/Neale_pheno_ID_corr) %>%
    mutate(sig_confounder = ifelse(IVW_pval < 0.05/num_tests_by_PCs, TRUE, FALSE)) %>%
    mutate(exposure_ID_phes = gsub("_irnt", "", exposure_ID)) %>%
    left_join(corr_mat_traits_sub, by = c("exposure_ID_phes" = "exposure_ID")) %>%
    filter(abs(between_trait_correlation) < corr_trait_threshold) %>% arrange(-corr_due_to_confounding_ratio) %>%
    filter(sig_confounder)

  # prune remaining by prioritizing those with highest `corr_due_to_confounding_ratio`

  output_pruned <- output
  counter <- 1
  while (counter < dim(output_pruned)[1]){

    index_trait <- output_pruned$exposure_ID[counter]
    index_trait_phes_ID <- gsub("_irnt", "", index_trait)

    corr_mat_traits_index_trait <- as.data.frame(corr_mat_traits) %>% dplyr::select(!!index_trait_phes_ID) %>% rownames_to_column(var = "exposure_ID") %>% as_tibble() %>%
      dplyr::rename(index_vs_trait_correlation = !!index_trait_phes_ID)

    temp <- output_pruned %>% left_join(corr_mat_traits_index_trait, by = c("exposure_ID_phes" = "exposure_ID"))

    remove_rows <- which(temp[["index_vs_trait_correlation"]] >= corr_trait_threshold)[-1]

    if(length(remove_rows)!=0){
      output_pruned <- temp[-remove_rows,]
    }

    counter <- 1 + counter

  }

  output_pruned <- output_pruned %>% dplyr::select(-index_vs_trait_correlation)

  return(output_pruned)

}

find_potential_PC_confounders <- function(Neale_pheno_ID, Neale_pheno_ID_corr, PC_trait_corr, PCs_corr, path_pheno_data){


  ###########################################################
  ## Get the number of pairs in the correlation between the
  ## trait `Neale_pheno_ID` and PCs
  ## to use in the SE calculation
  ###########################################################

  phes_ID <- gsub("_irnt", "", Neale_pheno_ID)

  male_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_male.txt")))]
  female_pheno_file <- path_pheno_data[which(endsWith(path_pheno_data, paste0("/phesant/", Neale_pheno_ID, "_female.txt")))]
  male_pheno_data <- fread(male_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))
  female_pheno_data <- fread(female_pheno_file, header=TRUE, select=c("IID","sex","age", phes_ID))

  joint_data <- rbind(male_pheno_data, female_pheno_data)

  trait_PC_n_pairs <- dim(joint_data)[1] # represents number of pairs in the correlation between the PCs and trait of interest

  output <- bind_rows(PC_trait_corr) %>% filter(Neale_pheno_ID == !!Neale_pheno_ID) %>%
    dplyr::select(Neale_pheno_ID, PC, corr_both_sexes) %>%
    rename(trait_PC_corr = corr_both_sexes) %>%
    mutate(trait_PC_n_pairs = trait_PC_n_pairs) %>%
    left_join(PCs_corr, by = c("PC" = "description")) %>% mutate_at("r2", as.numeric) %>%
    mutate(PC_couple_corr = sqrt(r2)) %>%
    dplyr::select(-ID, -ID_sub, -r2) %>%
    rename(PC_couple_n_pairs = N_pairs) %>%
    mutate(corr_due_to_confounding = trait_PC_corr^2*PC_couple_corr) %>%
    mutate(corr_due_to_confounding_ratio = corr_due_to_confounding/Neale_pheno_ID_corr)

  return(output)

}


meta_binned_household_MR <- function(exposure_info, outcomes_to_run, household_MR_binned){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nMeta-analyzing binned household MR result across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))


  output_list <- list()
  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    household_MR_binned_i <- household_MR_binned[[i]]
    household_MR_binned_i_all <- household_MR_binned_i %>% dplyr::filter(bin=="all" & grouping_var=="time_together_even_bins")

    male_beta <- household_MR_binned_i_all[which(household_MR_binned_i_all$exposure_sex=="male"),"IVW_beta"][[1]]
    female_beta <- household_MR_binned_i_all[which(household_MR_binned_i_all$exposure_sex=="female"),"IVW_beta"][[1]]
    male_se <- household_MR_binned_i_all[which(household_MR_binned_i_all$exposure_sex=="male"),"IVW_se"][[1]]
    female_se <- household_MR_binned_i_all[which(household_MR_binned_i_all$exposure_sex=="female"),"IVW_se"][[1]]

    sex_het_p <- z.test_p(male_beta, male_se, female_beta, female_se)

    meta_result <- household_MR_binned_i %>% group_by(grouping_var, bin) %>% group_modify(~ summarize_sex_specific_results(.x$IVW_beta, .x$IVW_se))
    meta_counts <- household_MR_binned_i %>% group_by(grouping_var, bin) %>% mutate(meta_n_outcome = sum(n_outcome)) %>% mutate(meta_n_exposure = sum(n_exposure)) %>% dplyr::select(grouping_var, bin, meta_n_outcome, meta_n_exposure)

    colnames(meta_result)[-c(1:2)] <- paste0("IVW", "_", colnames(meta_result)[-c(1:2)])

    meta_result <- full_join(meta_result, meta_counts) %>% unique()
    summarized_result <- full_join(household_MR_binned_i, meta_result)

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR_meta")]] <- summarized_result

    cat(paste0("Finished meta-anlayzing across bins for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}


calc_binned_household_MR_het <- function(exposure_info, outcomes_to_run, household_MR_binned_meta){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nCalculating difference in household MR between sexes and among age and time-together bins for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))


  exposure_i_result <- numeric()
  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    household_MR_binned_meta_i <- household_MR_binned_meta[[i]]
    household_MR_binned_meta_i_all <- household_MR_binned_meta_i %>% dplyr::filter(bin=="all" & grouping_var=="time_together_even_bins")


    sex_het_p <- household_MR_binned_meta_i_all$IVW_sex_het_pval[[1]]

    for(group in c("age_even_bins", "time_together_even_bins")){

      group_i_result <- numeric()
      for(exposure_sex in c("male", "female")){

        if(exposure_sex=="male"){outcome_sex="female"}
        if(exposure_sex=="female"){outcome_sex="male"}

        bin_result_temp <- household_MR_binned_meta_i %>% filter(grouping_var==group) %>% filter(exposure_sex==!!exposure_sex)

        all_row <- which(bin_result_temp$bin=="all")
        all_sum <- c(bin_result_temp[all_row,"IVW_beta"][[1]],bin_result_temp[all_row,"IVW_se"][[1]], bin_result_temp[all_row,"IVW_pval"][[1]],
                     bin_result_temp[all_row,"IVW_meta_beta"][[1]],bin_result_temp[all_row,"IVW_meta_se"][[1]], bin_result_temp[all_row,"IVW_meta_pval"][[1]])
        names(all_sum) <- c("IVW_beta", "IVW_se", "IVW_pval", "IVW_meta_beta", "IVW_meta_se", "IVW_meta_pval")

        result_to_analyze <- bin_result_temp %>% filter(bin != "all")

        meta_bin <- metagen(TE = as.numeric(IVW_beta), seTE = as.numeric(IVW_se), studlab = bin, data = result_to_analyze)
        Q_stat <- meta_bin$Q
        Q_pval <- meta_bin$pval.Q

        result_to_analyze <- result_to_analyze %>% separate(bin, c("bin_start_temp", "bin_stop_temp"), ",", remove = F) %>% mutate(bin_start = substring(bin_start_temp, 2)) %>%
          mutate(bin_stop = str_sub(bin_stop_temp,1,nchar(bin_stop_temp)-1)) %>% rowwise() %>%  mutate(bin_median = median(c(as.numeric(bin_start), as.numeric(bin_stop))))

        bin_lm_weight <- lm(IVW_beta ~ bin_median, data = result_to_analyze, weights = 1/(IVW_se^2))
        bin_lm <- lm(IVW_beta ~ bin_median, data = result_to_analyze)


        lm_summary_weight <- summary(bin_lm_weight)$coefficients["bin_median",c(1,2,4)]
        lm_summary <- summary(bin_lm)$coefficients["bin_median",c(1,2,4)]

        diff_sum <- c(sex_het_p, Q_stat, Q_pval, lm_summary, lm_summary_weight)

        names(diff_sum) <- c("sex_het_pval", "Q_stat", "Q_pval",
                             "bin_slope_beta", "bin_slope_se", "bin_slope_pval",
                             "bin_slope_beta_wt", "bin_slope_se_wt", "bin_slope_pval_wt")

        same_trait <- ifelse(exposure_ID==outcome_ID, TRUE, FALSE)
        description_sum <- c(exposure_ID, outcome_ID, exposure_sex, outcome_sex, group, same_trait)
        names(description_sum) <- c("exposure_ID", "outcome_ID", "exposure_sex", "outcome_sex", "grouping_var", "same_trait")

        ## generate the Q-stat and slope for the meta-analyzed results

        meta_bin_meta <- metagen(TE = as.numeric(IVW_meta_beta), seTE = as.numeric(IVW_meta_se), studlab = bin, data = result_to_analyze)
        Q_stat_meta <- meta_bin_meta$Q
        Q_pval_meta <- meta_bin_meta$pval.Q
        bin_lm_meta <- lm(IVW_meta_beta ~ bin_median, data = result_to_analyze)
        bin_lm_meta_weight <- lm(IVW_meta_beta ~ bin_median, data = result_to_analyze, weights = 1/(IVW_meta_se^2))

        lm_summary_meta <- summary(bin_lm_meta)$coefficients["bin_median",c(1,2,4)]
        lm_summary_meta_weight <-  summary(bin_lm_meta_weight)$coefficients["bin_median",c(1,2,4)]


        meta_summ <- c(Q_stat_meta, Q_pval_meta, lm_summary_meta, lm_summary_meta_weight)
        names(meta_summ) <- c("Q_stat_meta", "Q_pval_meta",
                              "bin_slope_meta_beta", "bin_slope_meta_se", "bin_slope_meta_pval",
                              "bin_slope_meta_beta_wt", "bin_slope_meta_se_wt", "bin_slope_meta_pval_wt")

        out_temp <- as.data.frame(t(c(description_sum, all_sum, diff_sum, meta_summ))) %>%
          mutate_if(is.factor,as.character) %>%
          as_tibble()

        group_i_result <- rbind( group_i_result, out_temp)

      }

      ## get the heterogeneity p-value for the slope

      for(slope in c("unweighted", "weighted")){

        appendage <- ifelse(slope=="unweighted", "", "_wt")

        betas <- as.numeric(group_i_result[[paste0("bin_slope_beta", appendage)]])
        ses <-  as.numeric(group_i_result[[paste0("bin_slope_se", appendage)]])
        bin_slope_sex_het_p = z.test_p(betas[1], ses[1], betas[2], ses[2])

        assign(paste0("bin_slope_sex_het_p", appendage), bin_slope_sex_het_p) #group_i_result$bin_slope_sex_het_pval <- bin_slope_sex_het_p


      }


      group_i_result$bin_slope_sex_het_pval <- bin_slope_sex_het_p
      group_i_result$bin_slope_sex_het_pval_wt <- bin_slope_sex_het_p_wt

      exposure_i_result <- rbind(exposure_i_result, group_i_result)

    }
    cat(paste0("Finished calculating heterogeneity statistics for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  #meta_result <- exposure_i_result %>% group_by(outcome_ID, grouping_var) %>% ## DO THE GUESS THING
  # group_modify(~ summarize_sex_specific_results(as.numeric(.x$bin_slope_beta), as.numeric(.x$bin_slope_se)))

  #colnames(meta_result)[-c(1:2)] <- paste0("bin_slope_", colnames(meta_result)[-c(1:2)])

  #output <- left_join(exposure_i_result, meta_result)

  output <- exposure_i_result %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  return(output)


}


calc_binned_household_MR_het_SNPmeta <- function(exposure_info, outcomes_to_run, household_MR_binned_joint){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nCalculating difference in household MR between sexes and among age and time-together bins for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))


  exposure_i_result <- numeric()
  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]

    household_MR_binned_joint_i <- household_MR_binned_joint[[i]]
    household_MR_binned_joint_i_all <- household_MR_binned_joint_i %>% dplyr::filter(bin=="all" & grouping_var=="time_together_even_bins")


    for(group in c("age_even_bins", "time_together_even_bins")){

      group_i_result <- numeric()


      bin_result_temp <- household_MR_binned_joint_i %>% filter(grouping_var==group)

      all_row <- which(bin_result_temp$bin=="all")
      all_sum <- c(bin_result_temp[all_row,"IVW_beta"][[1]],bin_result_temp[all_row,"IVW_se"][[1]], bin_result_temp[all_row,"IVW_pval"][[1]])
      names(all_sum) <- c("IVW_beta", "IVW_se", "IVW_pval")

      result_to_analyze <- bin_result_temp %>% filter(bin != "all")

      meta_bin <- metagen(TE = as.numeric(IVW_beta), seTE = as.numeric(IVW_se), studlab = bin, data = result_to_analyze)
      Q_stat <- meta_bin$Q
      Q_pval <- meta_bin$pval.Q

      result_to_analyze <- result_to_analyze %>% separate(bin, c("bin_start_temp", "bin_stop_temp"), ",", remove = F) %>% mutate(bin_start = substring(bin_start_temp, 2)) %>%
        mutate(bin_stop = str_sub(bin_stop_temp,1,nchar(bin_stop_temp)-1)) %>% rowwise() %>%  mutate(bin_median = median(c(as.numeric(bin_start), as.numeric(bin_stop))))

      bin_lm_weight <- lm(IVW_beta ~ bin_median, data = result_to_analyze, weights = 1/(IVW_se^2))
      bin_lm <- lm(IVW_beta ~ bin_median, data = result_to_analyze)


      lm_summary_weight <- summary(bin_lm_weight)$coefficients["bin_median",c(1,2,4)]
      lm_summary <- summary(bin_lm)$coefficients["bin_median",c(1,2,4)]

      diff_sum <- c(Q_stat, Q_pval, lm_summary, lm_summary_weight)

      names(diff_sum) <- c("Q_stat", "Q_pval",
                           "bin_slope_beta", "bin_slope_se", "bin_slope_pval",
                           "bin_slope_beta_wt", "bin_slope_se_wt", "bin_slope_pval_wt")

      same_trait <- ifelse(exposure_ID==outcome_ID, TRUE, FALSE)
      exposure_sex <- "joint"
      outcome_sex <- "joint"
      description_sum <- c(exposure_ID, outcome_ID, exposure_sex, outcome_sex, group, same_trait)
      names(description_sum) <- c("exposure_ID", "outcome_ID", "exposure_sex", "outcome_sex", "grouping_var", "same_trait")

      out_temp <- as.data.frame(t(c(description_sum, all_sum, diff_sum))) %>%
        mutate_if(is.factor,as.character) %>%
        as_tibble()

      group_i_result <- rbind( group_i_result, out_temp)
      exposure_i_result <- rbind(exposure_i_result, group_i_result)
    }


    cat(paste0("Finished calculating heterogeneity statistics for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  output <- exposure_i_result %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  return(output)


}


household_MR_comprehensive_ind <- function(harmonise_dat, MR_method_list){


  exposure_ID <- harmonise_dat$exposure_ID[1]
  outcome_ID <- harmonise_dat$outcome_ID[1]
  exposure_sex <- harmonise_dat$exposure_sex[1]
  outcome_sex <- harmonise_dat$outcome_sex[1]

  if(is.null(exposure_sex)){
    exposure_sex <- NA
    outcome_sex <- NA
  }


  harmonise_dat <- harmonise_dat %>%
    mutate(z_score.exposure = beta.exposure / se.exposure) %>%
    mutate(z_score.outcome = beta.outcome / se.outcome) %>%
    mutate(std_beta.exposure = z_score.exposure / sqrt(samplesize.exposure)) %>%
    mutate(std_beta.outcome = z_score.outcome / sqrt(samplesize.outcome)) %>%
    mutate(std_se.exposure = 1 / sqrt(samplesize.exposure)) %>%
    mutate(std_se.outcome = 1 / sqrt(samplesize.outcome)) %>%
    mutate(std_t_stat = (abs(std_beta.exposure) - abs(std_beta.outcome)) /
             sqrt(std_se.exposure^2 + std_se.outcome^2))

  harmonise_dat <- harmonise_dat %>% dplyr::select(-beta.exposure, -se.exposure, -beta.outcome, -se.outcome) %>%
    rename(beta.exposure = std_beta.exposure) %>% rename(se.exposure = std_se.exposure) %>%
    rename(beta.outcome = std_beta.outcome) %>% rename(se.outcome = std_se.outcome)


  exposure_description <- harmonise_dat$exposure_description[1]
  outcome_description <- harmonise_dat$outcome_description[1]

  het_test <- mr_heterogeneity(harmonise_dat, method_list=c("mr_egger_regression", "mr_ivw"))
  egger_test <- mr_pleiotropy_test(harmonise_dat)
  leave1out_test <- mr_leaveoneout(harmonise_dat)

  ngwas <- max(harmonise_dat$samplesize.outcome, na.rm=T)
  nsnps <- dim(harmonise_dat)[1]
  n_neale <- max(harmonise_dat$samplesize.exposure, na.rm=T)

  # Standard MR
  original_MR <- mr(harmonise_dat, method=MR_method_list)
  MR_res <- include_MR_NA(original_MR)

  MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
  MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")

  ## Test for reverse causality
  harmonise_dat_sensitivity <- harmonise_dat[which(harmonise_dat[["pval.exposure"]] < harmonise_dat[["pval.outcome"]]),]
  if(dim(harmonise_dat_sensitivity)[1]!=0){
    MR_res_sensitivity <- include_MR_NA(mr(harmonise_dat_sensitivity, method=MR_method_list))
    nsnps_sensitivity <- dim(harmonise_dat_sensitivity)[1]
    correct_row <- ifelse(nsnps_sensitivity==1, MR_wald_row, MR_ivw_row)
    sensitivity_result <- cbind(nsnps_sensitivity, MR_res_sensitivity[correct_row,"b"],MR_res_sensitivity[correct_row,"se"], MR_res_sensitivity[correct_row,"pval"], make_beta_95ci(MR_res_sensitivity[correct_row,"b"],MR_res_sensitivity[correct_row,"se"]),pretty_round(MR_res_sensitivity[correct_row,"pval"]))

  }

  if(dim(harmonise_dat_sensitivity)[1]==0){
    nsnps_sensitivity <- dim(harmonise_dat_sensitivity)[1]
    sensitivity_result <- cbind(nsnps_sensitivity, NA, NA, NA, NA, NA)

  }

  temp_summary <- summarize_mr_result (paste0(exposure_ID,"_INDEX"), paste0(outcome_ID,"_HOUSEHOLD"), exposure_sex, outcome_sex, nsnps, MR_res, het_test, egger_test,n_neale, ngwas)


  temp_summary <- cbind(temp_summary, sensitivity_result)

  same_trait <- exposure_ID == outcome_ID
  MR_summary <- cbind(exposure_ID, outcome_ID,exposure_description, outcome_description, same_trait, temp_summary)


  num_cols <- length(colnames(MR_summary))

  colnames(MR_summary)[num_cols-1] <- "IVW_sensitivity_summary"
  colnames(MR_summary)[num_cols] <- "IVW_sensitivity_pval_round"

  colnames(MR_summary)[num_cols-2] <- "IVW_sensitivity_pval"
  colnames(MR_summary)[num_cols-3] <- "IVW_sensitivity_beta"
  colnames(MR_summary)[num_cols-4] <- "IVW_sensitivity_se"

  colnames(MR_summary)[num_cols-5] <- "N_snps_sensitivity" #test for reverse causation

  mr_plot <- household_MR_plot(harmonise_dat, original_MR)

  leave1out_test$exposure_sex <- exposure_sex
  leave1out_test$outcome_sex <- outcome_sex


  MR_summary <- MR_summary %>% as_tibble() %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)
  colnames(MR_summary) <- gsub("_Wald", "", colnames(MR_summary)) #remove the Wald from the name because Wald results are only returned if there is only 1 SNP in the MR and we know there are more than this for each MR.

  return(list(MR_summary = MR_summary, leave1out_test = leave1out_test, MR_plot = mr_plot))

}

run_household_MR_comprehensive <- function(exposure_info, outcomes_to_run, household_harmonised_data, MR_method_list){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning complete MR analyses for all outcomes with phenotype `", exposure_ID, "`\nas exposure (i.e. Egger, leave-one-out, sensitivity, MR plot. \n[This is only being run in full sample, not in individual bins.]\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    for(exposure_sex in c("male", "female")){

      harmonised_dat_i_sex <- household_harmonised_data[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_filter")]][[paste0("exp_", exposure_sex, "_harmonised_data_filter")]]
      harmonised_dat_sub <- harmonised_dat_i_sex %>% filter(grouping_var == "age_even_bins") %>% filter(bin == "all")
      MR_complete_i_sex <- household_MR_comprehensive_ind(harmonised_dat_sub, MR_method_list)
      MR_complete_i[[paste0("exp_", exposure_sex, "_MR_complete")]] <- MR_complete_i_sex
    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR_complete")]] <- MR_complete_i

    cat(paste0("Finished computing full MR results for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}

run_household_MR_comprehensive_SNPmeta <- function(exposure_info, outcomes_to_run, household_harmonised_data_meta_reverse_filter, MR_method_list){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning complete MR analyses for all outcomes with phenotype `", exposure_ID, "`\nas exposure (i.e. Egger, leave-one-out, sensitivity, MR plot. \n[This is only being run in full sample, not in individual bins.]\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    harmonised_dat_i <- household_harmonised_data_meta_reverse_filter[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta_filter")]]
    harmonised_dat_sub <- harmonised_dat_i %>% filter(grouping_var == "age_even_bins") %>% filter(bin == "all")
    MR_complete_i <- household_MR_comprehensive_ind(harmonised_dat_sub, MR_method_list)


    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR_complete")]] <- MR_complete_i

    cat(paste0("Finished computing full MR results for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}


run_household_MVMR <- function(exposure_info, outcomes_to_run){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning MVMR analyses for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()

  exposure_ID_prev <- ""
  outcome_ID_prev <- ""

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    for(exposure_sex in c("male", "female")){


      #######################################
      ## Standard GWAS results with Y IVs ###
      #######################################

      GWAS_file_yIVs_1 <- paste0("analysis/traitMR/standard_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
      xi_YIV_GWAS_results <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == !!exposure_sex)
      xi_YIV_GWAS_results_both <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == "both_sexes")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xi for only Y IVs.

      GWAS_file_yIVs_2 <- paste0("analysis/traitMR/standard_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
      yi_YIV_GWAS_results <- fread(GWAS_file_yIVs_2, data.table = F) %>% filter(sex == !!exposure_sex)
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yi for only Y IVs.

      #######################################
      ## Standard GWAS results with X IVs ###
      #######################################

      GWAS_file_xIVs_1 <- paste0("analysis/traitMR/standard_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
      xi_XIV_GWAS_results <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == !!exposure_sex)
      xi_XIV_GWAS_results_both <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == "both_sexes")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xi for only X IVs.

      GWAS_file_XIVs_2 <- paste0("analysis/traitMR/standard_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yi_XIV_GWAS_results <- fread(GWAS_file_XIVs_2, data.table = F) %>% filter(sex == !!exposure_sex)
      ## this will pull the standard GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yi for only X IVs.


      ## Prune the XIV and YIVs based on associations with Xi

      if(exposure_ID!=exposure_ID_prev | outcome_ID!=outcome_ID_prev){

        snps <- c(xi_YIV_GWAS_results_both$SNP, xi_XIV_GWAS_results_both$SNP)
        chr <- c(xi_YIV_GWAS_results_both$chr, xi_XIV_GWAS_results_both$chr)
        pvals <- c(xi_YIV_GWAS_results_both$pval, xi_XIV_GWAS_results_both$pval)


        to_prune_mat <- tibble(SNP = snps,
                               chr = chr,
                               pval = pvals) %>% unique()

        pruned_mat <- IV_clump(to_prune_mat, prune_threshold)

      }

      #####################################################
      ## Household GWAS results with X and Y IVs for Xp ###
      #####################################################

      hh_GWAS_file_yIVs_1 <- paste0("analysis/traitMR/household_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
      xp_YIV_GWAS_results <- fread(hh_GWAS_file_yIVs_1, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xp for only Y IVs.

      hh_GWAS_file_xIVs_1 <- paste0("analysis/traitMR/household_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
      xp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs_1, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xp for only X IVs.

      yi_GWAS <- rbind(yi_YIV_GWAS_results, yi_XIV_GWAS_results)
      xi_GWAS <- rbind(xi_YIV_GWAS_results, xi_XIV_GWAS_results)
      xp_GWAS <- rbind(xp_YIV_GWAS_results, xp_XIV_GWAS_results)

      ## join the xi and yi data together
      standard_GWAS_data <- full_join(yi_GWAS, xi_GWAS, by = "SNP", suffix = c("_yi","_xi")) %>% dplyr::select(-exposure_ID_xi, -exposure_ID_yi) %>% unique()

      ## join the xi and yi data with the xp data
      xp_GWAS_sub <- xp_GWAS %>% dplyr::select(SNP, geno_index_beta, geno_index_se, allele1) %>% rename(beta_xp = geno_index_beta,  se_xp = geno_index_se, effect_allele_xp = allele1)
      mv_X_data <- full_join(standard_GWAS_data, xp_GWAS_sub, by = "SNP") %>% unique()

      # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
      # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
      mv_X_data <- mv_X_data %>% mutate(beta_xp_align = case_when(effect_allele_xp == effect_allele_yi ~ beta_xp,
                                                                  TRUE ~ -1*beta_xp))


      #####################################################
      ## Household GWAS results with X and Y IVs for Yp ###
      #####################################################

      hh_GWAS_file_yIVs_2 <- paste0("analysis/traitMR/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
      yp_YIV_GWAS_results <- fread(hh_GWAS_file_yIVs_2, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yp for only Y IVs.

      hh_GWAS_file_xIVs_2 <- paste0("analysis/traitMR/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs_2, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yp for only X IVs.

      yp_GWAS <- rbind(yp_YIV_GWAS_results, yp_XIV_GWAS_results)
      yp_GWAS_sub <- yp_GWAS %>% dplyr::select(SNP, geno_index_beta, geno_index_se, allele1) %>% rename(beta_yp = geno_index_beta,  se_yp = geno_index_se, effect_allele_yp = allele1)

      mv_data <- full_join(mv_X_data, yp_GWAS_sub, by = "SNP") %>% unique()

      # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
      # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
      mv_data <- mv_data %>% mutate(beta_yp_align = case_when(effect_allele_yp == effect_allele_yi ~ beta_yp,
                                                              TRUE ~ -1*beta_yp))



      mv_data_pruned <- mv_data %>% filter(SNP %in% pruned_mat$SNP)
      mv_data_format <- mr_mvinput(bx = cbind(mv_data_pruned[["beta_yi"]], mv_data_pruned[["beta_xi"]], mv_data_pruned[["beta_xp_align"]]), bxse = cbind(mv_data_pruned[["se_yi"]], mv_data_pruned[["se_xi"]], mv_data_pruned[["se_xp"]]),
                                   by = mv_data_pruned[["beta_yp_align"]], byse = mv_data_pruned[["se_yp"]])



      mv_result <- mr_mvivw(mv_data_format) # see `str(mv_result)` to know how to access the results
      betas <- mv_result@Estimate
      ses <- mv_result@StdError
      pvalues <- mv_result@Pvalue
      mv_yi <- c(betas[1], ses[1], pvalues[1])
      names(mv_yi) <- paste0("yi", c("_beta", "_se", "_pval"))
      mv_xi <- c(betas[2], ses[2], pvalues[2])
      names(mv_xi) <- paste0("xi", c("_beta", "_se", "_pval"))
      mv_xp <- c(betas[3], ses[3], pvalues[3])
      names(mv_xp) <- paste0("xp", c("_beta", "_se", "_pval"))


      temp_row <- c(mv_yi, mv_xi, mv_xp)
      description <- c(exposure_ID, outcome_ID, exposure_sex)
      names(description) <- c("exposure_ID", "outcome_ID", "exposure_sex")

      out_temp <- as.data.frame(t(c(description, temp_row)))

      exposure_ID_prev <- exposure_ID
      outcome_ID_prev <- outcome_ID

      result <- rbind(result, out_temp)

    }

  }

  result_out <- result %>%
    as_tibble() %>% type_convert() %>%
    mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  return(result_out)

}


run_household_MVMR_SNPmeta <- function(exposure_info, outcomes_to_run){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning MVMR analyses for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()

  exposure_ID_prev <- ""
  outcome_ID_prev <- ""

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    #for(exposure_sex in c("male", "female")){


      #######################################
      ## Standard GWAS results with Y IVs ###
      #######################################

      GWAS_file_yIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
      xi_YIV_GWAS_results <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == "meta")
      #xi_YIV_GWAS_results_both <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == "both_sexes")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xi for only Y IVs.

      GWAS_file_yIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
      yi_YIV_GWAS_results <- fread(GWAS_file_yIVs_2, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yi for only Y IVs.

      #######################################
      ## Standard GWAS results with X IVs ###
      #######################################

      GWAS_file_xIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
      xi_XIV_GWAS_results <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == "meta")
      #xi_XIV_GWAS_results_both <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == "both_sexes")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xi for only X IVs.

      GWAS_file_XIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yi_XIV_GWAS_results <- fread(GWAS_file_XIVs_2, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yi for only X IVs.


      ## Prune the XIV and YIVs based on associations with Xi

      if(exposure_ID!=exposure_ID_prev | outcome_ID!=outcome_ID_prev){

        snps <- c(xi_YIV_GWAS_results$SNP, xi_XIV_GWAS_results$SNP)
        chr <- c(xi_YIV_GWAS_results$chr, xi_XIV_GWAS_results$chr)
        pvals <- c(xi_YIV_GWAS_results$pval, xi_XIV_GWAS_results$pval)


        to_prune_mat <- tibble(SNP = snps,
                               chr = chr,
                               pval = pvals) %>% unique()

        pruned_mat <- IV_clump(to_prune_mat, prune_threshold)

      }

      #####################################################
      ## Household GWAS results with X and Y IVs for Xp ###
      #####################################################

      hh_GWAS_file_yIVs_1 <- paste0("analysis/traitMR/household_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
      xp_YIV_GWAS_results <- fread(hh_GWAS_file_yIVs_1, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xp for only Y IVs.

      hh_GWAS_file_xIVs_1 <- paste0("analysis/traitMR/household_GWAS/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
      xp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs_1, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xp for only X IVs.

      yi_GWAS <- rbind(yi_YIV_GWAS_results, yi_XIV_GWAS_results)
      xi_GWAS <- rbind(xi_YIV_GWAS_results, xi_XIV_GWAS_results)
      xp_GWAS <- rbind(xp_YIV_GWAS_results, xp_XIV_GWAS_results)

      ## join the xi and yi data together
      standard_GWAS_data <- full_join(yi_GWAS, xi_GWAS, by = "SNP", suffix = c("_yi","_xi")) %>% dplyr::select(-exposure_ID_xi, -exposure_ID_yi) %>% unique()

      ## join the xi and yi data with the xp data
      xp_GWAS_sub <- xp_GWAS %>% dplyr::select(SNP, geno_index_beta, geno_index_se, allele1) %>% rename(beta_xp = geno_index_beta,  se_xp = geno_index_se, effect_allele_xp = allele1)
      mv_X_data <- full_join(standard_GWAS_data, xp_GWAS_sub, by = "SNP") %>% unique()

      # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
      # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
      mv_X_data <- mv_X_data %>% mutate(beta_xp_align = case_when(effect_allele_xp == effect_allele_yi ~ beta_xp,
                                                                  TRUE ~ -1*beta_xp))


      #####################################################
      ## Household GWAS results with X and Y IVs for Yp ###
      #####################################################

      hh_GWAS_file_yIVs_2 <- paste0("analysis/traitMR/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
      yp_YIV_GWAS_results <- fread(hh_GWAS_file_yIVs_2, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yp for only Y IVs.

      hh_GWAS_file_xIVs_2 <- paste0("analysis/traitMR/household_GWAS/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs_2, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yp for only X IVs.

      yp_GWAS <- rbind(yp_YIV_GWAS_results, yp_XIV_GWAS_results)
      yp_GWAS_sub <- yp_GWAS %>% dplyr::select(SNP, geno_index_beta, geno_index_se, allele1) %>% rename(beta_yp = geno_index_beta,  se_yp = geno_index_se, effect_allele_yp = allele1)

      mv_data <- full_join(mv_X_data, yp_GWAS_sub, by = "SNP") %>% unique()

      # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
      # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
      mv_data <- mv_data %>% mutate(beta_yp_align = case_when(effect_allele_yp == effect_allele_yi ~ beta_yp,
                                                              TRUE ~ -1*beta_yp))



      mv_data_pruned <- mv_data %>% filter(SNP %in% pruned_mat$SNP)
      mv_data_format <- mr_mvinput(bx = cbind(mv_data_pruned[["beta_yi"]], mv_data_pruned[["beta_xi"]], mv_data_pruned[["beta_xp_align"]]), bxse = cbind(mv_data_pruned[["se_yi"]], mv_data_pruned[["se_xi"]], mv_data_pruned[["se_xp"]]),
                                   by = mv_data_pruned[["beta_yp_align"]], byse = mv_data_pruned[["se_yp"]])



      mv_result <- mr_mvivw(mv_data_format) # see `str(mv_result)` to know how to access the results
      betas <- mv_result@Estimate
      ses <- mv_result@StdError
      pvalues <- mv_result@Pvalue
      mv_yi <- c(betas[1], ses[1], pvalues[1])
      names(mv_yi) <- paste0("yi", c("_beta", "_se", "_pval"))
      mv_xi <- c(betas[2], ses[2], pvalues[2])
      names(mv_xi) <- paste0("xi", c("_beta", "_se", "_pval"))
      mv_xp <- c(betas[3], ses[3], pvalues[3])
      names(mv_xp) <- paste0("xp", c("_beta", "_se", "_pval"))


      temp_row <- c(mv_yi, mv_xi, mv_xp)
      description <- c(exposure_ID, outcome_ID, exposure_sex)
      names(description) <- c("exposure_ID", "outcome_ID", "exposure_sex")

      out_temp <- as.data.frame(t(c(description, temp_row)))

      exposure_ID_prev <- exposure_ID
      outcome_ID_prev <- outcome_ID

      result <- rbind(result, out_temp)

    #}

  }

  result_out <- result %>%
    as_tibble() %>% type_convert() %>%
    mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  return(result_out)

}



summarize_sex_specific_results  <- function(d,se){

  if(any(is.na(d)) | any(is.na(se))){

    to_rm <- union(which(is.na(d)), which(is.na(se)))
    d <- d[-to_rm]
    se <- se[-to_rm]
  }

  if(length(d)!=0){
    meta.result <- meta.summaries(d, se, method=c("fixed"), conf.level=0.95)

    b_meta=round(meta.result[3]$summary,digits=10)
    b_meta_se=round(meta.result[4]$se.summary,digits=10)
    lowerbound=b_meta-b_meta_se*1.96
    upperbound=b_meta+b_meta_se*1.96

    # n_total <- sum(n)
    ## meta_p=round(2*(pt(abs(b_meta/b_meta_se),((n_total)-meta.result$het[2]),lower.tail=FALSE)),digits=10)
    # The p-value above is based on t-distribution and the one below is based on z. Since we are working with such large samples it shouldn't make any difference.

    meta_p <- 2*pnorm(-abs(b_meta/b_meta_se))

    se <- sqrt( (se[1]^2) + (se[2]^2) )
    t <- (d[1]-d[2])/se
    p_het <- 2*pnorm(-abs(t))

    # same as above
    # p_het2 <- z.test_p(d[1], d[2], se[1], se[2])
  } else {

    b_meta <- NA
    b_meta_se <- NA
    lowerbound <- NA
    upperbound <- NA
    meta_p <- NA
    p_het <- NA
  }

  tibble(meta_beta = b_meta,
         meta_se = b_meta_se,
         meta_l95 = lowerbound,
         meta_u95 = upperbound,
         meta_pval = meta_p,
         sex_het_pval = p_het
  )

}



summarize_household_MR_comprehensive <- function(household_MR, corr_mat_traits){

  exposure_ID <- household_MR[[1]][["exp_male_MR_complete"]][["MR_summary"]][,"exposure_ID"][[1]]
  cat(paste0("Summarizing and meta-analyzing household MR results across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()
  for(i in 1:length(household_MR)){
    outcome_ID <- household_MR[[i]][["exp_male_MR_complete"]][["MR_summary"]][,"outcome_ID"][[1]]
    exposure_IDX <- exposure_ID
    outcome_IDX <- outcome_ID
    if(endsWith(exposure_ID, "_irnt")) {
      exposure_IDX <- gsub("_irnt", "", exposure_ID)
    }
    if(endsWith(outcome_ID, "_irnt")) {
      outcome_IDX <- gsub("_irnt", "", outcome_ID)
    }
    corr_traits <- corr_mat_traits[exposure_IDX, outcome_IDX]
    male_result <- household_MR[[i]][["exp_male_MR_complete"]][["MR_summary"]]
    female_result <- household_MR[[i]][["exp_female_MR_complete"]][["MR_summary"]]
    result_i <- rbind(male_result, female_result)
    result_i$corr_traits <- corr_traits
    result <- rbind(result, result_i)

  }

  result <- as_tibble(result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  meta_result <- result %>% group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$IVW_beta, .x$IVW_se))

  colnames(meta_result)[-c(1:2)] <- paste0("IVW", "_", colnames(meta_result)[-c(1:2)])
  summarized_result <- full_join(result, meta_result)

  return(summarized_result)

}

summarize_household_MR_comprehensive_SNPmeta <- function(household_MR_joint, corr_mat_traits){

  exposure_ID <- household_MR_joint[[1]][["MR_summary"]][,"exposure_ID"][[1]]
  cat(paste0("Summarizing and meta-analyzing household MR results across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()
  for(i in 1:length(household_MR_joint)){
    outcome_ID <- household_MR_joint[[i]][["MR_summary"]][,"outcome_ID"][[1]]
    exposure_IDX <- exposure_ID
    outcome_IDX <- outcome_ID
    if(endsWith(exposure_ID, "_irnt")) {
      exposure_IDX <- gsub("_irnt", "", exposure_ID)
    }
    if(endsWith(outcome_ID, "_irnt")) {
      outcome_IDX <- gsub("_irnt", "", outcome_ID)
    }
    corr_traits <- corr_mat_traits[exposure_IDX, outcome_IDX]
    result_i <- household_MR_joint[[i]][["MR_summary"]]
    result_i$corr_traits <- corr_traits
    result <- rbind(result, result_i)

  }

  result <- as_tibble(result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID", "exposure_sex", "outcome_sex"), as.character)

  return(result)

}

summarize_household_MVMR <- function(household_MVMR, traits_corr2_filled, corr_mat_traits){

  outcome_traits <- traits_corr2_filled[which(traits_corr2_filled[["Neale_file_sex"]]=="both"),]
  exposure_ID <- household_MVMR[["exposure_ID"]][1]
  exposure_description <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==exposure_ID), "description"])

  cat(paste0("Summarizing and meta-analyzing household MMVR results across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  meta_list <- list()
  for(var in c("yi", "xi", "xp")){

    meta_temp <- household_MVMR %>% dplyr::select(exposure_ID, outcome_ID, starts_with(var)) %>% rename_all(~stringr::str_replace(., paste0("^", var, "_"),"")) %>%
      dplyr::group_by(exposure_ID, outcome_ID) %>%
      # note can't easily use the `z.test_p` function because it isn't set up for vector of beta and se's, but rather one argument for each
      dplyr::group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$beta, .x$se))

    colnames(meta_temp)[-c(1:2)] <- paste0(var, "_", colnames(meta_temp)[-c(1:2)])
    meta_list[[var]] <- meta_temp

  }

  meta_result <- reduce(meta_list, full_join)

  summarized_result <- full_join(household_MVMR, meta_result)
  summarized_result$exposure_description <- NA
  summarized_result$outcome_description <- NA
  summarized_result$corr_traits <- NA

  for(i in 1:dim(summarized_result)[1]){
    outcome_ID <- household_MVMR[["outcome_ID"]][i]
    exposure_IDX <- exposure_ID
    outcome_IDX <- outcome_ID
    if(endsWith(exposure_ID, "_irnt")) {
      exposure_IDX <- gsub("_irnt", "", exposure_ID)
    }
    if(endsWith(outcome_ID, "_irnt")) {
      outcome_IDX <- gsub("_irnt", "", outcome_ID)
    }
    corr_traits <- corr_mat_traits[exposure_IDX, outcome_IDX]
    summarized_result$corr_traits[i] <- corr_traits

    outcome_description <- as.character(outcome_traits[which(outcome_traits[["Neale_pheno_ID"]]==outcome_ID), "description"])
    summarized_result$exposure_description[i] <- exposure_description
    summarized_result$outcome_description[i] <- outcome_description

  }

  summarized_result <- summarized_result  %>% mutate(same_trait = ifelse(exposure_ID==outcome_ID, TRUE, FALSE)) %>% select(exposure_ID, outcome_ID, exposure_description, outcome_description, same_trait, everything())
  return(summarized_result)

}

run_standard_MR_comprehensive <- function(exposure_info, outcomes_to_run, standard_harmonised_data, MR_method_list){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning complete MR analyses for all outcomes with phenotype `", exposure_ID, "as exposure\n(i.e. Egger, leave-one-out, sensitivity, MR plot in same individual (i.e. standard MR). \n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    if(outcome_ID != exposure_ID){
      for(exposure_sex in c("male", "female")){

        standard_dat_i_sex <- standard_harmonised_data[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data")]][[paste0("exp_", exposure_sex, "_harmonised_data")]]
        MR_complete_i_sex <- household_MR_comprehensive_ind(standard_dat_i_sex, MR_method_list) ##the inner function for running the MR is the same for household and standard MR.
        MR_complete_i[[paste0("exp_", exposure_sex, "_MR_complete")]] <- MR_complete_i_sex
      }
    }


    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR_complete")]] <- MR_complete_i

    cat(paste0("Finished computing full MR results for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}


run_standard_MR_comprehensive_SNPmeta <- function(exposure_info, outcomes_to_run, standard_harmonised_data_meta_reverse_filter, MR_method_list){

  output_list <- list()
  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)

  cat(paste0("\nRunning complete MR analyses for all outcomes with phenotype `", exposure_ID, "as exposure\n(i.e. Egger, leave-one-out, sensitivity, MR plot in same individual (i.e. standard MR). \n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_complete_i <- list()

    if(outcome_ID != exposure_ID){

      standard_dat_i <- standard_harmonised_data_meta_reverse_filter[[paste0(outcome_ID, "_vs_", exposure_ID, "_harmonised_data_meta_filter")]]
      MR_complete_i <- household_MR_comprehensive_ind(standard_dat_i, MR_method_list) ##the inner function for running the MR is the same for household and standard MR.

    }

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR_complete")]] <- MR_complete_i

    cat(paste0("Finished computing full MR results for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_list)

}

write_household_MR <- function(exposure_info, outcomes_to_run, household_MR){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  output_files <- numeric()
  pheno_dir <- paste0("analysis/traitMR")

  cat(paste0("\nWriting household MR results for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  for(i in 1:dim(outcomes_to_run)[1]){

    outcome_ID <- outcomes_to_run$Neale_pheno_ID[[i]]
    MR_file_i <- paste0(pheno_dir, "/household_MR/", outcome_ID, "/univariate_MR/", outcome_ID, "_vs_", exposure_ID, "_MR.csv")

    output_files <- c(output_files, MR_file_i)

    MR_result_i <- household_MR[[paste0(outcome_ID, "_vs_", exposure_ID, "_MR")]]
    write.csv(MR_result_i, MR_file_i, row.names = F)

    cat(paste0("Finished writing MR results for outcome ", i, " of ", dim(outcomes_to_run)[1], ".\n\n" ))

  }

  return(output_files)

}

summarize_standard_MR_comprehensive <- function(standard_MR){

  exposure_ID <- standard_MR[[1]][["exp_male_MR_complete"]][["MR_summary"]][,"exposure_ID"][[1]]
  cat(paste0("Summarizing and meta-analyzing standard MR results across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()
  for(i in 1:length(standard_MR)){
    if(!length(standard_MR[[i]])==0){
      male_result <- standard_MR[[i]][["exp_male_MR_complete"]][["MR_summary"]]
      female_result <- standard_MR[[i]][["exp_female_MR_complete"]][["MR_summary"]]
      result_i <- rbind(male_result, female_result)
      result <- rbind(result, result_i)
    }
  }

  result <- as_tibble(result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  meta_result <- result %>% group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$IVW_beta, .x$IVW_se))

  colnames(meta_result)[-c(1:2)] <- paste0("IVW", "_", colnames(meta_result)[-c(1:2)])
  summarized_result <- full_join(result, meta_result)

  return(summarized_result)

}


summarize_standard_MR_comprehensive_SNPmeta <- function(standard_MR_joint){

  exposure_ID <- standard_MR_joint[[1]][["MR_summary"]][,"exposure_ID"][[1]]
  cat(paste0("Summarizing and meta-analyzing standard MR results across sexes for all outcomes with phenotype `", exposure_ID, "` as exposure.\n\n"))

  result <- numeric()
  for(i in 1:length(standard_MR_joint)){
    if(!length(standard_MR_joint[[i]])==0){
       result_i <- standard_MR_joint[[i]][["MR_summary"]]
      result <- rbind(result, result_i)
    }
  }

  result <- as_tibble(result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID", "exposure_sex", "outcome_sex"), as.character)


  return(result)

}

## data -> data_set with a column that has Neale_ID that we want to replace with description (eg. "exposure_ID")
## traits_corr2_filled -> reference data set which has both Neale_ID and descrition
## data_Neale_ID_col -> the column name of the Neale_ID (eg. "exposure_ID")
## output_column -> the output column name of the description (eg. "description")

replace_Neale_ID <- function(data, traits_corr2_filled, data_Neale_ID_col, output_column){

  descriptions <- traits_corr2_filled[match(data[[data_Neale_ID_col]], traits_corr2_filled$Neale_pheno_ID), "description"]
  data[[data_Neale_ID_col]] <- descriptions
  colnames(data)[which(colnames(data)==data_Neale_ID_col)] <- output_column
  return(data)

}

find_sig_standard_MR_summary <- function(standard_MR_summary_joint){

  denom <- dim(standard_MR_summary_joint)[1] #divide by 2 because each meta result is there twice (one row/sex)
  sig_only <- standard_MR_summary_joint %>% filter(IVW_pval < 0.05/denom)
  return(sig_only)
}

find_non_corr_household_MR_summary <- function(household_MR_summary_joint, corr_trait_threshold){

  # filter to only pairs that have correlation < 0.8
  output <- household_MR_summary_joint[which(abs(household_MR_summary_joint$corr_traits)<corr_trait_threshold),]

  return(output)
}

find_sig_household_MR_summary <- function(household_MR_summary_corr_filter){

  denom <- dim(household_MR_summary_corr_filter)[1]
  sig_only <- household_MR_summary_corr_filter %>% filter(IVW_pval < 0.05/denom)
  return(sig_only)
}

pull_AM_MRs_household_MR_summary <- function(household_MR_summary){

  output <- household_MR_summary %>% filter(same_trait=="TRUE")
  return(output)
}

find_AM_sig_exposure_info <- function(household_MR_summary_AM, outcomes_to_run, num_tests_by_PCs){


  AM_sig_traits <- household_MR_summary_AM %>% filter(IVW_pval < 0.05/num_tests_by_PCs) %>% pull(exposure_ID) %>% unique()

  outcomes_to_run_sub <- outcomes_to_run[which(outcomes_to_run$Neale_pheno_ID %in% AM_sig_traits),]

  return(outcomes_to_run_sub)

}

z.test_p <- Vectorize(function(x, sigma.x, y, sigma.y) {z.test(x, sigma.x, y, sigma.y)$p},
                        vectorize.args = c("x", "sigma.x", "y", "sigma.y"))

z.test_z <- Vectorize(function(x, sigma.x, y, sigma.y) {z.test(x, sigma.x, y, sigma.y)$statistic},
                      vectorize.args = c("x", "sigma.x", "y", "sigma.y"))

variance_of_product <- function(x, x_se, y, y_se){

  x^2*y_se^2 + y^2*x_se^2 + x_se^2*y_se^2
}

variance_of_sum <- function(x_se, y_se){
  x_se^2 + y_se^2
}

run_proxyMR_comparison <- function(exposure_info, household_MR_summary_BF_sig, household_MR_summary, standard_MR_summary, household_MR_summary_AM){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(exposure_ID==!!exposure_ID) %>% filter(exposure_ID!=outcome_ID)

  summarized_result <- as_tibble(numeric())

  cat(paste0("Comparing gamma, rho and omega estimates for phenotype `", exposure_ID, "` as exposure.\n\n"))


  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){
      outcome_ID <- MR_sub$outcome_ID[[i]]
      for(exposure_sex in c("male", "female")){
        if(exposure_sex=="male"){outcome_sex="female"}
        if(exposure_sex=="female"){outcome_sex="male"}

        # p = partner / outcome

        hh_MR_sub <- household_MR_summary %>% filter(exposure_sex==!!exposure_sex) %>%
          filter(outcome_ID==!!outcome_ID)

        cols_interst <- c("IVW_beta", "IVW_se", "IVW_pval", "N_exposure_GWAS", "N_outcome_GWAS", "N_snps")


        ## Sex-specific proxy MR
        xiyp_summary <- hh_MR_sub %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyp_', names(.)))

        ## Sex-specific AM MR
        xixp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!exposure_ID) %>% filter(outcome_ID==!!exposure_ID) %>%
          filter(exposure_sex==!!exposure_sex) %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xixp_', names(.)))
        yiyp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!outcome_ID) %>% filter(outcome_ID==!!outcome_ID) %>%
          filter(exposure_sex==!!exposure_sex) %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('yiyp_', names(.)))

        ## Sex-specific standard MR
        ## (y is the outcome, x is the exposure)
        ## exposure and outcome sex are the same, i.e. exposure/outcome sex are irrelevant

        xiyi_summary <- standard_MR_summary %>% filter(exposure_sex==!!exposure_sex) %>% filter(outcome_ID==!!outcome_ID) %>%
          dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyi_', names(.)))
        xpyp_summary <- standard_MR_summary %>% filter(outcome_sex==!!outcome_sex) %>% filter(outcome_ID==!!outcome_ID) %>%
          dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xpyp_', names(.)))



        summary_cols <- MR_sub %>% slice(i) %>% dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex)
        row_i <- cbind(summary_cols, xiyp_summary, xixp_summary, yiyp_summary, xiyi_summary, xpyp_summary)

        summarized_result <- rbind(summarized_result, row_i)
      }


    }


    summarized_result <- as_tibble(summarized_result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

    expsoure_sex_temp <- summarized_result$exposure_sex
    outcome_sex_temp <- ifelse(expsoure_sex_temp=="male", "female", "male")

    ## outcome sex previously was the same as exposure sex coming from the standard MR (where exposure and outcome sex would be the same)
    summarized_result$outcome_sex <- outcome_sex_temp


    # Calc products and their SEs

    summarized_result <- summarized_result %>% #dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, contains(var), contains("N_outcome_GWAS"), contains("N_snps")) %>% rename_all(~stringr::str_replace(., paste0("_", var, "_"),"_")) %>%
      mutate(xixp_xpyp_beta = xixp_IVW_beta*xpyp_IVW_beta) %>%
      mutate(xixp_xpyp_se = sqrt(variance_of_product(xixp_IVW_beta, xixp_IVW_se, xpyp_IVW_beta, xpyp_IVW_se))) %>%

      mutate(xiyi_yiyp_beta = xiyi_IVW_beta*yiyp_IVW_beta) %>%
      mutate(xiyi_yiyp_se = sqrt(variance_of_product(xiyi_IVW_beta, xiyi_IVW_se, yiyp_IVW_beta, yiyp_IVW_se))) %>%

      mutate(gam_beta = xixp_xpyp_beta) %>% mutate(gam_se = xixp_xpyp_se) %>%
      mutate(rho_beta = xiyi_yiyp_beta) %>% mutate(rho_se = xiyi_yiyp_se) %>%
      mutate(omega_beta = xiyp_IVW_beta) %>% mutate(omega_se = xiyp_IVW_se) %>%
      mutate(gam_rho_beta = gam_beta + rho_beta) %>%
      mutate(gam_rho_se = sqrt(variance_of_sum(gam_se, rho_se)))

    # Meta analyzed MRs estimates and gam, rho and omega across sexes and calculate heterogeneity between them

    meta_list <- list()
    for(var in c("xixp_IVW", "xpyp_IVW", "xiyi_IVW", "yiyp_IVW", "xiyp_IVW", "gam", "rho", "omega", "gam_rho")){

      meta_temp <- summarized_result %>% dplyr::select(exposure_ID, outcome_ID, starts_with(var)) %>% rename_all(~stringr::str_replace(., paste0("^", var, "_"),"")) %>%
        dplyr::group_by(exposure_ID, outcome_ID) %>%
        # note can't easily use the `z.test_p` function because it isn't set up for vector of beta and se's, but rather one argument for each
        dplyr::group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$beta, .x$se))

      colnames(meta_temp)[-c(1:2)] <- paste0(var, "_", colnames(meta_temp)[-c(1:2)])
      meta_list[[var]] <- meta_temp

    }

    meta_result <- reduce(meta_list, full_join)

    summarized_result_prev <- summarized_result
    summarized_result <- full_join(summarized_result_prev, meta_result)

    # Calc difference between estimates (gamma, rho and omega)
    prod_diff_list <- list()

    prod_diff_result <- summarized_result %>%
      mutate(omega_vs_gam_pval = z.test_p(gam_beta, gam_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_rho_pval = z.test_p(rho_beta, rho_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_gam_rho_pval = z.test_p(gam_rho_beta, gam_rho_se, omega_beta, omega_se)) %>%
      mutate(gam_vs_rho_pval = z.test_p(rho_beta, rho_se, gam_beta, gam_se)) %>%


      mutate(omega_vs_gam_meta_pval = z.test_p(gam_meta_beta, gam_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(omega_vs_rho_meta_pval = z.test_p(rho_meta_beta, rho_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(omega_vs_gam_rho_meta_pval = z.test_p(gam_rho_meta_beta, gam_rho_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(gam_vs_rho_meta_pval = z.test_p(rho_meta_beta, rho_meta_se, gam_meta_beta, gam_meta_se)) %>%

      dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, starts_with("gam"), starts_with("rho"), starts_with("omega"))


    summarized_result_sub <- summarized_result %>% dplyr::select(-starts_with("gam")) %>% dplyr::select(-starts_with("omega")) %>% dplyr::select(-starts_with("rho")) %>% dplyr::select(-starts_with("xiyi_yiyp")) %>% dplyr::select(-starts_with("xixp_xpyp"))
    output <- list(MR_results = summarized_result_sub, proxy_MR_comparison = prod_diff_result)

  } else output <- NULL

  return(output)

}

run_proxyMR_comparison_SNPmeta <- function(exposure_info, household_MR_summary_BF_sig, household_MR_summary_joint, standard_MR_summary_joint, household_MR_summary_AM){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(exposure_ID==!!exposure_ID) %>% filter(exposure_ID!=outcome_ID)

  summarized_result <- as_tibble(numeric())

  cat(paste0("Comparing gamma, rho and omega estimates for phenotype `", exposure_ID, "` as exposure.\n\n"))


  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){
      outcome_ID <- MR_sub$outcome_ID[[i]]
      #exposure_sex <- MR_sub$exposure_sex[[i]]
      #if(exposure_sex=="male"){outcome_sex="female"}
      #if(exposure_sex=="female"){outcome_sex="male"}

      # p = partner / outcome

      hh_MR_sub <- household_MR_summary_joint %>% #filter(exposure_sex==!!exposure_sex) %>%
        filter(outcome_ID==!!outcome_ID)

      cols_interst <- c("IVW_beta", "IVW_se", "IVW_pval", "N_exposure_GWAS", "N_outcome_GWAS", "N_snps")


      ## Sex-specific proxy MR
      xiyp_summary <- hh_MR_sub %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyp_', names(.)))

      ## Sex-specific AM MR
      xixp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!exposure_ID) %>% filter(outcome_ID==!!exposure_ID) %>% #filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xixp_', names(.)))
      yiyp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!outcome_ID) %>% filter(outcome_ID==!!outcome_ID) %>% #filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('yiyp_', names(.)))

      ## Sex-specific standard MR
      ## (y is the outcome, x is the exposure)
      ## exposure and outcome sex are the same, i.e. exposure/outcome sex are irrelevant

      xiyi_summary <- standard_MR_summary_joint %>% filter(outcome_ID==!!outcome_ID) %>% # %>% filter(exposure_sex==!!exposure_sex)
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyi_', names(.)))
      xpyp_summary <- standard_MR_summary_joint %>% filter(outcome_ID==!!outcome_ID) %>% # %>% filter(outcome_sex==!!outcome_sex)
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xpyp_', names(.)))


      summary_cols <- MR_sub %>% slice(i) %>% dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description)
      row_i <- cbind(summary_cols, xiyp_summary, xixp_summary, yiyp_summary, xiyi_summary, xpyp_summary)

      summarized_result <- rbind(summarized_result, row_i)
    }


    summarized_result <- as_tibble(summarized_result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

    #expsoure_sex_temp <- summarized_result$exposure_sex
    #outcome_sex_temp <- ifelse(expsoure_sex_temp=="male", "female", "male")

    ## outcome sex previously was the same as exposure sex coming from the standard MR (where exposure and outcome sex would be the same)
    #summarized_result$outcome_sex <- outcome_sex_temp


    # Calc products and their SEs

    summarized_result <- summarized_result %>% #dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, contains(var), contains("N_outcome_GWAS"), contains("N_snps")) %>% rename_all(~stringr::str_replace(., paste0("_", var, "_"),"_")) %>%
      mutate(xixp_xpyp_beta = xixp_IVW_beta*xpyp_IVW_beta) %>%
      mutate(xixp_xpyp_se = sqrt(variance_of_product(xixp_IVW_beta, xixp_IVW_se, xpyp_IVW_beta, xpyp_IVW_se))) %>%

      mutate(xiyi_yiyp_beta = xiyi_IVW_beta*yiyp_IVW_beta) %>%
      mutate(xiyi_yiyp_se = sqrt(variance_of_product(xiyi_IVW_beta, xiyi_IVW_se, yiyp_IVW_beta, yiyp_IVW_se))) %>%

      mutate(gam_beta = xixp_xpyp_beta) %>% mutate(gam_se = xixp_xpyp_se) %>%
      mutate(rho_beta = xiyi_yiyp_beta) %>% mutate(rho_se = xiyi_yiyp_se) %>%
      mutate(omega_beta = xiyp_IVW_beta) %>% mutate(omega_se = xiyp_IVW_se) %>%
      mutate(gam_rho_beta = gam_beta + rho_beta) %>%
      mutate(gam_rho_se = sqrt(variance_of_sum(gam_se, rho_se)))

    # Meta analyzed MRs estimates and gam, rho and omega across sexes and calculate heterogeneity between them

    summarized_result_prev <- summarized_result

    # Calc difference between estimates (gamma, rho and omega)
    prod_diff_list <- list()

    prod_diff_result <- summarized_result %>%
      mutate(omega_vs_gam_pval = z.test_p(gam_beta, gam_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_rho_pval = z.test_p(rho_beta, rho_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_gam_rho_pval = z.test_p(gam_rho_beta, gam_rho_se, omega_beta, omega_se)) %>%
      mutate(gam_vs_rho_pval = z.test_p(rho_beta, rho_se, gam_beta, gam_se)) %>%

      dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, starts_with("gam"), starts_with("rho"), starts_with("omega"))


    summarized_result_sub <- summarized_result %>% dplyr::select(-starts_with("gam")) %>% dplyr::select(-starts_with("omega")) %>% dplyr::select(-starts_with("rho")) %>% dplyr::select(-starts_with("xiyi_yiyp")) %>% dplyr::select(-starts_with("xixp_xpyp"))
    output <- list(MR_results = summarized_result_sub, proxy_MR_comparison = prod_diff_result)

  } else output <- NULL

  return(output)

}


adj_yiyp_xIVs_sex_specific <- function(exposure_info, household_harmonised_data, household_MR_summary_BF_sig){

  outcome_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID!=outcome_ID)
  summarized_result <- as_tibble(numeric())

  cat(paste0("Adjusting all relevant YiYp MRs of `", outcome_ID, "` for effects from Xi to Yp (for all relevant exposures).\n\n"))
  result <- numeric()
  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){

      exposure_ID <- MR_sub$exposure_ID[[i]]

      for(exposure_sex in c("male", "female")){

        if(exposure_sex=="male"){outcome_sex="female"}
        if(exposure_sex=="female"){outcome_sex="male"}
        harmonised_dat_name <- paste0(outcome_ID, "_vs_", outcome_ID, "_harmonised_data_filter")
        expsosure_sex_name <- paste0("exp_", exposure_sex, "_harmonised_data_filter")
        univar_harmonised_dat <- household_harmonised_data[[harmonised_dat_name]][[expsosure_sex_name]] %>% filter(bin == "all" & grouping_var == "time_together_even_bins")

        hh_GWAS_file_xIVs <- paste0("analysis/traitMR/household_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
        yp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs, data.table = F) %>% filter(exposure_sex==!!exposure_sex & bin == "all" & grouping_var == "time_together_even_bins")
        ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yp for only X IVs.


        GWAS_file_yIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
        xi_YIV_GWAS_results <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == !!exposure_sex)
        ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xi for only Y IVs.

        GWAS_file_yIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
        yi_YIV_GWAS_results <- fread(GWAS_file_yIVs_2, data.table = F) %>% filter(sex == !!exposure_sex)
        yi_YIV_GWAS_results_both <- fread(GWAS_file_yIVs_2, data.table = F) %>% filter(sex == "meta")

        ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yi for only Y IVs.

        ## need to also pull the effect of G on Xi for X IVs
        ## and the effect of G on Yi on for X IVs
        ## i.e. we want to include X IVs in the MV MR as well: Y_p ~ Y_i + X_i

        GWAS_file_xIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
        xi_XIV_GWAS_results <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == !!exposure_sex)
        ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xi for only X IVs.

        GWAS_file_xIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
        yi_XIV_GWAS_results <- fread(GWAS_file_xIVs_2, data.table = F) %>% filter(sex == !!exposure_sex)
        yi_XIV_GWAS_results_both <- fread(GWAS_file_xIVs_2, data.table = F) %>% filter(sex == "meta")

        ## this will pull the standard GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yi for only X IVs.

        ## Pruning should be based on meta-analyzed pvals from effect on `outcome_ID`
        ## First subset `yi_YIV_GWAS_results` to be only those set of SNPs that show no evidence of reverse causality with X

        yi_YIV_GWAS_results_both_sub <- yi_YIV_GWAS_results_both %>% filter(SNP %in% xi_YIV_GWAS_results$SNP)
        yi_YIV_GWAS_results_sub <- yi_YIV_GWAS_results %>% filter(SNP %in% xi_YIV_GWAS_results$SNP)


        snps <- c(yi_YIV_GWAS_results_both_sub$SNP, yi_XIV_GWAS_results_both$SNP)
        chr <- c(yi_YIV_GWAS_results_both_sub$chr, yi_XIV_GWAS_results_both$chr)
        pvals <- c(yi_YIV_GWAS_results_both_sub$pval, yi_XIV_GWAS_results_both$pval)

        # 865
        to_prune_mat <- tibble(SNP = snps,
                               chr = chr,
                               pval = pvals) %>% unique()

        # 487
        pruned_mat <- IV_clump(to_prune_mat, prune_threshold)

        ## Subset `xi_XIV_GWAS_results` to be only those set of SNPs that show no evidence of reverse causality with Y
        xi_XIV_GWAS_results_sub <- xi_XIV_GWAS_results %>% filter(SNP %in% yi_XIV_GWAS_results$SNP)

        yi_GWAS <- rbind(yi_YIV_GWAS_results_sub, yi_XIV_GWAS_results)
        xi_GWAS <- rbind(xi_YIV_GWAS_results, xi_XIV_GWAS_results_sub)

        mv_X_data <- full_join(yi_GWAS, xi_GWAS, by = "SNP", suffix = c("_yi","_xi")) %>% dplyr::select(-exposure_ID_xi, -exposure_ID_yi) %>% unique()

        mv_Y_data_YIV <- univar_harmonised_dat[,c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")] %>% filter(SNP %in% xi_YIV_GWAS_results$SNP)
        mv_Y_data_XIV <- yp_XIV_GWAS_results[,c("SNP", "beta", "se", "effect_allele")] # allele1 is effect allele

        colnames_MV_y <- c("SNP", "beta_yp", "se_yp", "effect_allele_yp")
        colnames(mv_Y_data_YIV) <- colnames_MV_y
        colnames(mv_Y_data_XIV) <- colnames_MV_y

        mv_Y_data <- rbind(mv_Y_data_YIV, mv_Y_data_XIV) %>% unique()

        mv_data <-  full_join(mv_X_data, mv_Y_data, by = "SNP")

        # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
        # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
        mv_data <- mv_data %>% mutate(beta_yp_align = case_when(effect_allele_yp == effect_allele_yi ~ beta_yp,
                                                                TRUE ~ -1*beta_yp))

        mv_data_pruned <- mv_data %>% filter(SNP %in% pruned_mat$SNP)
        mv_data_format <- mr_mvinput(bx = cbind(mv_data_pruned[["beta_yi"]], mv_data_pruned[["beta_xi"]]), bxse = cbind(mv_data_pruned[["se_yi"]], mv_data_pruned[["se_xi"]]),
                                     by = mv_data_pruned[["beta_yp_align"]], byse = mv_data_pruned[["se_yp"]])



        mv_result <- mr_mvivw(mv_data_format) # see `str(mv_result)` to know how to access the results
        betas <- mv_result@Estimate
        ses <- mv_result@StdError
        pvalues <- mv_result@Pvalue
        mv_yi <- c(betas[1], ses[1], pvalues[1])
        names(mv_yi) <- paste0("yi", c("_beta", "_se", "_pval"))
        mv_xi <- c(betas[2], ses[2], pvalues[2])
        names(mv_xi) <- paste0("xi", c("_beta", "_se", "_pval"))

        temp_row <- c(mv_yi, mv_xi)
        description <- c(exposure_ID, outcome_ID, exposure_sex)
        names(description) <- c("exposure_ID", "outcome_ID", "exposure_sex")

        out_temp <- as.data.frame(t(c(description, temp_row)))

        result <- rbind(result, out_temp)

      }

    }

    result_out <- result %>%
      as_tibble() %>% type_convert() %>%
      mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  } else result_out <- NULL

  return(result_out)

}



adj_yiyp_xIVs_SNPmeta <- function(exposure_info, household_harmonised_data_meta_reverse_filter, household_MR_summary_BF_sig){

  outcome_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID!=outcome_ID)
  summarized_result <- as_tibble(numeric())

  cat(paste0("Adjusting all relevant YiYp MRs of `", outcome_ID, "` for effects from Xi to Yp (for all relevant exposures).\n\n"))
  result <- numeric()
  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){

      ## The numbers below represent using exposure_ID #17 (i.e. `exposure_info <- tar_read(exposure_info, branches = 17)[[1]]` and `i=1`)
      exposure_ID <- MR_sub$exposure_ID[[i]]
      harmonised_dat_name <- paste0(outcome_ID, "_vs_", outcome_ID, "_harmonised_data_meta_filter")
      univar_harmonised_dat <- household_harmonised_data_meta_reverse_filter[[harmonised_dat_name]] %>% filter(bin == "all" & grouping_var == "time_together_even_bins")

      #1 - 516 IVs
      hh_GWAS_file_xIVs <- paste0("analysis/traitMR/household_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yp_XIV_GWAS_results <- fread(hh_GWAS_file_xIVs, data.table = F) %>% filter(exposure_sex=="meta" & bin == "all" & grouping_var == "time_together_even_bins")
      ## this will pull the household GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yp for only X IVs.

      #2 - 384 IVs
      GWAS_file_yIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", outcome_ID, "_GWAS.csv")
      xi_YIV_GWAS_results <- fread(GWAS_file_yIVs_1, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `outcome_ID`, i.e. the effect of G on Xi for only Y IVs.

      #3 - 409 IVs
      GWAS_file_yIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", outcome_ID, "_GWAS.csv")
      yi_YIV_GWAS_results <- fread(GWAS_file_yIVs_2, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `outcome_ID` for IVs from `outcome_ID`, i.e. the effect of G on Yi for only Y IVs.

      ## need to also pull the effect of G on Xi for X IVs
      ## and the effect of G on Yi on for X IVs
      ## i.e. we want to include X IVs in the MV MR as well: Y_p ~ Y_i + X_i

      #4 - 521 IVs
      GWAS_file_xIVs_1 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", exposure_ID, "/", exposure_ID, "_vs_", exposure_ID, "_GWAS.csv")
      xi_XIV_GWAS_results <- fread(GWAS_file_xIVs_1, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `exposure_ID` for IVs from `exposure_ID`, i.e. the effect of G on Xi for only X IVs.

      #5 - 516 IVs
      GWAS_file_xIVs_2 <- paste0("analysis/traitMR/standard_GWAS_rev_filter/", outcome_ID, "/", outcome_ID, "_vs_", exposure_ID, "_GWAS.csv")
      yi_XIV_GWAS_results <- fread(GWAS_file_xIVs_2, data.table = F) %>% filter(sex == "meta")
      ## this will pull the standard GWAS results for phenotype `outcome_ID` for IVs from `exposure_ID`, i.e. the effect of G on Yi for only X IVs.

      ## Pruning should be based on meta-analyzed pvals from effect on `outcome_ID`
      ## First subset `yi_YIV_GWAS_results` to be only those set of SNPs that show no evidence of reverse causality with X

      yi_YIV_GWAS_results_sub <- yi_YIV_GWAS_results %>% filter(SNP %in% xi_YIV_GWAS_results$SNP)

      snps <- c(yi_YIV_GWAS_results_sub$SNP, yi_XIV_GWAS_results$SNP)
      chr <- c(yi_YIV_GWAS_results_sub$chr, yi_XIV_GWAS_results$chr)
      pvals <- c(yi_YIV_GWAS_results_sub$pval, yi_XIV_GWAS_results$pval)


      ## 865
      to_prune_mat <- tibble(SNP = snps,
                             chr = chr,
                             pval = pvals) %>% unique()

      ## 487
      pruned_mat <- IV_clump(to_prune_mat, prune_threshold)

      ## Subset `xi_XIV_GWAS_results` to be only those set of SNPs that show no evidence of reverse causality with Y
      xi_XIV_GWAS_results_sub <- xi_XIV_GWAS_results %>% filter(SNP %in% yi_XIV_GWAS_results$SNP)


      yi_GWAS <- rbind(yi_YIV_GWAS_results_sub, yi_XIV_GWAS_results)
      xi_GWAS <- rbind(xi_YIV_GWAS_results, xi_XIV_GWAS_results_sub)

      mv_X_data <- full_join(yi_GWAS, xi_GWAS, by = "SNP", suffix = c("_yi","_xi")) %>% dplyr::select(-exposure_ID_xi, -exposure_ID_yi) %>% unique()

      mv_Y_data_YIV <- univar_harmonised_dat[,c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")] %>% filter(SNP %in% xi_YIV_GWAS_results$SNP)
      mv_Y_data_XIV <- yp_XIV_GWAS_results[,c("SNP", "beta", "se", "effect_allele")]

      colnames_MV_y <- c("SNP", "beta_yp", "se_yp", "effect_allele_yp")
      colnames(mv_Y_data_YIV) <- colnames_MV_y
      colnames(mv_Y_data_XIV) <- colnames_MV_y

      mv_Y_data <- rbind(mv_Y_data_YIV, mv_Y_data_XIV) %>% unique()

      mv_data <-  full_join(mv_X_data, mv_Y_data, by = "SNP")

      # the GWAS data with effects on Xi and Yi should all be aligned because they come from Neale database and same effect alleles were used across all GWAS.
      # double check that alleles are aligned for x and y in mv_data. Y data is effect on Yp (run in house) and is possible that alleles aren't aligned
      mv_data <- mv_data %>% mutate(beta_yp_align = case_when(effect_allele_yp == effect_allele_yi ~ beta_yp,
                                                              TRUE ~ -1*beta_yp))

      mv_data_pruned <- mv_data %>% filter(SNP %in% pruned_mat$SNP)
      mv_data_format <- mr_mvinput(bx = cbind(mv_data_pruned[["beta_yi"]], mv_data_pruned[["beta_xi"]]), bxse = cbind(mv_data_pruned[["se_yi"]], mv_data_pruned[["se_xi"]]),
                                   by = mv_data_pruned[["beta_yp_align"]], byse = mv_data_pruned[["se_yp"]])



      mv_result <- mr_mvivw(mv_data_format) # see `str(mv_result)` to know how to access the results
      betas <- mv_result@Estimate
      ses <- mv_result@StdError
      pvalues <- mv_result@Pvalue
      mv_yi <- c(betas[1], ses[1], pvalues[1])
      names(mv_yi) <- paste0("yi", c("_beta", "_se", "_pval"))
      mv_xi <- c(betas[2], ses[2], pvalues[2])
      names(mv_xi) <- paste0("xi", c("_beta", "_se", "_pval"))

      temp_row <- c(mv_yi, mv_xi)

      exposure_sex <- "meta"
      description <- c(exposure_ID, outcome_ID, exposure_sex)
      names(description) <- c("exposure_ID", "outcome_ID", "exposure_sex")

      out_temp <- as.data.frame(t(c(description, temp_row)))

      result <- rbind(result, out_temp)

    }

    result_out <- result %>%
      as_tibble() %>% type_convert() %>%
      mutate_at(c("exposure_ID", "outcome_ID"), as.character)

  } else result_out <- NULL

  return(result_out)

}


run_proxyMR_comparison_adj_yiyp <- function(exposure_info, household_MR_summary_BF_sig, household_MR_summary, standard_MR_summary, household_MR_summary_AM, proxyMR_yiyp_adj){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(exposure_ID==!!exposure_ID) %>% filter(exposure_ID!=outcome_ID)

  summarized_result <- as_tibble(numeric())

  cat(paste0("Comparing gamma, rho and omega estimates for phenotype `", exposure_ID, "` as exposure.\n\n"))


  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){
      outcome_ID <- MR_sub$outcome_ID[[i]]
      exposure_sex <- MR_sub$exposure_sex[[i]]
      if(exposure_sex=="male"){outcome_sex="female"}
      if(exposure_sex=="female"){outcome_sex="male"}

      # p = partner / outcome

      hh_MR_sub <- household_MR_summary %>% filter(exposure_sex==!!exposure_sex) %>%
        filter(outcome_ID==!!outcome_ID)

      cols_interst <- c("IVW_beta", "IVW_se", "IVW_pval", "N_exposure_GWAS", "N_outcome_GWAS", "N_snps")


      ## Sex-specific proxy MR
      xiyp_summary <- hh_MR_sub %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyp_', names(.)))

      ## Sex-specific AM MR
      xixp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!exposure_ID) %>% filter(outcome_ID==!!exposure_ID) %>%
        filter(exposure_sex==!!exposure_sex) %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xixp_', names(.)))
      yiyp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!outcome_ID) %>% filter(outcome_ID==!!outcome_ID) %>%
        filter(exposure_sex==!!exposure_sex) %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('yiyp_', names(.)))

      yiyp_summary_adj <- proxyMR_yiyp_adj %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID==!!exposure_ID) %>%
        filter(exposure_sex==!!exposure_sex) %>% dplyr::select(yi_beta, yi_se, yi_pval) %>% setNames(c('yiyp_IVW_beta_adj', 'yiyp_IVW_se_adj', 'yiyp_IVW_pval_adj'))


      ## Sex-specific standard MR
      ## (y is the outcome, x is the exposure)
      ## exposure and outcome sex are the same, i.e. exposure/outcome sex are irrelevant

      xiyi_summary <- standard_MR_summary %>% filter(exposure_sex==!!exposure_sex) %>% filter(outcome_ID==!!outcome_ID) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyi_', names(.)))
      xpyp_summary <- standard_MR_summary %>% filter(outcome_sex==!!outcome_sex) %>% filter(outcome_ID==!!outcome_ID) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xpyp_', names(.)))



      summary_cols <- MR_sub %>% slice(i) %>% dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex)
      row_i <- cbind(summary_cols, xiyp_summary, xixp_summary, yiyp_summary, yiyp_summary_adj, xiyi_summary, xpyp_summary)

      summarized_result <- rbind(summarized_result, row_i)
    }


    summarized_result <- as_tibble(summarized_result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

    expsoure_sex_temp <- summarized_result$exposure_sex
    outcome_sex_temp <- ifelse(expsoure_sex_temp=="male", "female", "male")

    ## outcome sex previously was the same as exposure sex coming from the standard MR (where exposure and outcome sex would be the same)
    summarized_result$outcome_sex <- outcome_sex_temp

    # Calc products and their SEs

    summarized_result <- summarized_result %>% #dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, contains(var), contains("N_outcome_GWAS"), contains("N_snps")) %>% rename_all(~stringr::str_replace(., paste0("_", var, "_"),"_")) %>%
      mutate(xixp_xpyp_beta = xixp_IVW_beta*xpyp_IVW_beta) %>%
      mutate(xixp_xpyp_se = sqrt(variance_of_product(xixp_IVW_beta, xixp_IVW_se, xpyp_IVW_beta, xpyp_IVW_se))) %>%

      mutate(xiyi_yiyp_beta = xiyi_IVW_beta*yiyp_IVW_beta_adj) %>%
      mutate(xiyi_yiyp_se = sqrt(variance_of_product(xiyi_IVW_beta, xiyi_IVW_se, yiyp_IVW_beta_adj, yiyp_IVW_se_adj))) %>%

      mutate(gam_beta = xixp_xpyp_beta) %>% mutate(gam_se = xixp_xpyp_se) %>%
      mutate(rho_beta = xiyi_yiyp_beta) %>% mutate(rho_se = xiyi_yiyp_se) %>%
      mutate(omega_beta = xiyp_IVW_beta) %>% mutate(omega_se = xiyp_IVW_se) %>%
      mutate(gam_rho_beta = gam_beta + rho_beta) %>%
      mutate(gam_rho_se = sqrt(variance_of_sum(gam_se, rho_se)))

    # Meta analyzed MRs estimates and gam, rho and omega across sexes and calculate heterogeneity between them

    meta_list <- list()
    for(var in c("xixp_IVW", "xpyp_IVW", "xiyi_IVW", "yiyp_IVW", "xiyp_IVW", "gam", "rho", "omega", "gam_rho")){

      meta_temp <- summarized_result %>% dplyr::select(exposure_ID, outcome_ID, starts_with(var)) %>% rename_all(~stringr::str_replace(., paste0("^", var, "_"),"")) %>%
        dplyr::group_by(exposure_ID, outcome_ID) %>%
        # note can't easily use the `z.test_p` function because it isn't set up for vector of beta and se's, but rather one argument for each
        dplyr::group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$beta, .x$se))

      if(var == "yiyp_IVW"){
        meta_temp_adj <- summarized_result %>% dplyr::select(exposure_ID, outcome_ID, starts_with(var)) %>% rename_all(~stringr::str_replace(., paste0("^", var, "_"),"")) %>%
          dplyr::group_by(exposure_ID, outcome_ID) %>%
          # note can't easily use the `z.test_p` function because it isn't set up for vector of beta and se's, but rather one argument for each
          dplyr::group_by(exposure_ID, outcome_ID) %>% group_modify(~ summarize_sex_specific_results(.x$beta_adj, .x$se_adj))

        colnames(meta_temp_adj)[-c(1:2)] <- paste0(var, "_", colnames(meta_temp)[-c(1:2)], "_adj")
        meta_list[[paste0(var, "_adj")]] <- meta_temp_adj
      }

      colnames(meta_temp)[-c(1:2)] <- paste0(var, "_", colnames(meta_temp)[-c(1:2)])
      meta_list[[var]] <- meta_temp

    }

    meta_result <- reduce(meta_list, full_join)

    summarized_result_prev <- summarized_result
    summarized_result <- full_join(summarized_result_prev, meta_result)

    # Calc difference between estimates (gamma, rho and omega)
    prod_diff_list <- list()

    prod_diff_result <- summarized_result %>%
      mutate(omega_vs_gam_pval = z.test_p(gam_beta, gam_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_rho_pval = z.test_p(rho_beta, rho_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_gam_rho_pval = z.test_p(gam_rho_beta, gam_rho_se, omega_beta, omega_se)) %>%
      mutate(gam_vs_rho_pval = z.test_p(rho_beta, rho_se, gam_beta, gam_se)) %>%


      mutate(omega_vs_gam_meta_pval = z.test_p(gam_meta_beta, gam_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(omega_vs_rho_meta_pval = z.test_p(rho_meta_beta, rho_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(omega_vs_gam_rho_meta_pval = z.test_p(gam_rho_meta_beta, gam_rho_meta_se, omega_meta_beta, omega_meta_se)) %>%
      mutate(gam_vs_rho_meta_pval = z.test_p(rho_meta_beta, rho_meta_se, gam_meta_beta, gam_meta_se)) %>%

      dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, starts_with("gam"), starts_with("rho"), starts_with("omega"))


    summarized_result_sub <- summarized_result %>% dplyr::select(-starts_with("gam")) %>% dplyr::select(-starts_with("omega")) %>% dplyr::select(-starts_with("rho")) %>% dplyr::select(-starts_with("xiyi_yiyp")) %>% dplyr::select(-starts_with("xixp_xpyp"))
    output <- list(MR_results = summarized_result_sub, proxy_MR_comparison = prod_diff_result)

  } else output <- NULL

  return(output)


}

run_proxyMR_comparison_adj_yiyp_SNPmeta <- function(exposure_info, household_MR_summary_BF_sig, household_MR_summary_SNPmeta, standard_MR_summary_SNPmeta, household_MR_summary_AM, proxyMR_yiyp_adj_SNPmeta){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(exposure_ID==!!exposure_ID) %>% filter(exposure_ID!=outcome_ID)

  summarized_result <- as_tibble(numeric())

  cat(paste0("Comparing gamma, rho and omega estimates for phenotype `", exposure_ID, "` as exposure.\n\n"))


  if(dim(MR_sub)[1]!=0){

    for(i in 1:dim(MR_sub)[1]){
      outcome_ID <- MR_sub$outcome_ID[[i]]
      #exposure_sex <- MR_sub$exposure_sex[[i]]
      #outcome_sex <- MR_sub$outcome_sex[[i]]

      #if(exposure_sex=="male"){outcome_sex="female"}
      #if(exposure_sex=="female"){outcome_sex="male"}

      # p = partner / outcome

      hh_MR_sub <- household_MR_summary_SNPmeta %>% #filter(exposure_sex==!!exposure_sex) %>%
        filter(outcome_ID==!!outcome_ID)

      cols_interst <- c("IVW_beta", "IVW_se", "IVW_pval", "N_exposure_GWAS", "N_outcome_GWAS", "N_snps")


      ## Sex-specific proxy MR
      xiyp_summary <- hh_MR_sub %>% dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyp_', names(.)))

      ## Sex-specific AM MR
      xixp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!exposure_ID) %>% filter(outcome_ID==!!exposure_ID) %>% # filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xixp_', names(.)))
      yiyp_summary <- household_MR_summary_AM %>% filter(exposure_ID==!!outcome_ID) %>% filter(outcome_ID==!!outcome_ID) %>% # filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('yiyp_', names(.)))

      yiyp_summary_adj <- proxyMR_yiyp_adj_SNPmeta %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID==!!exposure_ID) %>% # filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(yi_beta, yi_se, yi_pval) %>% setNames(c('yiyp_IVW_beta_adj', 'yiyp_IVW_se_adj', 'yiyp_IVW_pval_adj'))


      ## standard MR results for Xi -> Yi,
      ## since results are meta-analyzed at SNP-level, these two MRs are the same
      ## (y is the outcome, x is the exposure)

      xiyi_summary <- standard_MR_summary_SNPmeta %>% filter(outcome_ID==!!outcome_ID) %>% # filter(exposure_sex==!!exposure_sex) %>%
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xiyi_', names(.)))
      xpyp_summary <- standard_MR_summary_SNPmeta %>% filter(outcome_ID==!!outcome_ID) %>% # %>% filter(outcome_sex==!!outcome_sex)
        dplyr::select(all_of(cols_interst)) %>% setNames(paste0('xpyp_', names(.)))



      summary_cols <- MR_sub %>% slice(i) %>% dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex)
      row_i <- cbind(summary_cols, xiyp_summary, xixp_summary, yiyp_summary, yiyp_summary_adj, xiyi_summary, xpyp_summary)

      summarized_result <- rbind(summarized_result, row_i)
    }


    summarized_result <- as_tibble(summarized_result) %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)

    #expsoure_sex_temp <- summarized_result$exposure_sex
    #outcome_sex_temp <- ifelse(expsoure_sex_temp=="male", "female", "male")

    ## outcome sex previously was the same as exposure sex coming from the standard MR (where exposure and outcome sex would be the same)
    #summarized_result$outcome_sex <- outcome_sex_temp

    # Calc products and their SEs

    summarized_result <- summarized_result %>% #dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, contains(var), contains("N_outcome_GWAS"), contains("N_snps")) %>% rename_all(~stringr::str_replace(., paste0("_", var, "_"),"_")) %>%
      mutate(xixp_xpyp_beta = xixp_IVW_beta*xpyp_IVW_beta) %>%
      mutate(xixp_xpyp_se = sqrt(variance_of_product(xixp_IVW_beta, xixp_IVW_se, xpyp_IVW_beta, xpyp_IVW_se))) %>%

      mutate(xiyi_yiyp_beta = xiyi_IVW_beta*yiyp_IVW_beta_adj) %>%
      mutate(xiyi_yiyp_se = sqrt(variance_of_product(xiyi_IVW_beta, xiyi_IVW_se, yiyp_IVW_beta_adj, yiyp_IVW_se_adj))) %>%

      mutate(gam_beta = xixp_xpyp_beta) %>% mutate(gam_se = xixp_xpyp_se) %>%
      mutate(rho_beta = xiyi_yiyp_beta) %>% mutate(rho_se = xiyi_yiyp_se) %>%
      mutate(omega_beta = xiyp_IVW_beta) %>% mutate(omega_se = xiyp_IVW_se) %>%
      mutate(gam_rho_beta = gam_beta + rho_beta) %>%
      mutate(gam_rho_se = sqrt(variance_of_sum(gam_se, rho_se)))

    summarized_result_prev <- summarized_result

    # Calc difference between estimates (gamma, rho and omega)
    prod_diff_list <- list()

    prod_diff_result <- summarized_result %>%
      mutate(omega_vs_gam_pval = z.test_p(gam_beta, gam_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_rho_pval = z.test_p(rho_beta, rho_se, omega_beta, omega_se)) %>%
      mutate(omega_vs_gam_rho_pval = z.test_p(gam_rho_beta, gam_rho_se, omega_beta, omega_se)) %>%
      mutate(gam_vs_rho_pval = z.test_p(rho_beta, rho_se, gam_beta, gam_se)) %>%

      dplyr::select(exposure_ID, outcome_ID, exposure_description, outcome_description, exposure_sex, outcome_sex, starts_with("gam"), starts_with("rho"), starts_with("omega"))


    summarized_result_sub <- summarized_result %>% dplyr::select(-starts_with("gam")) %>% dplyr::select(-starts_with("omega")) %>% dplyr::select(-starts_with("rho")) %>% dplyr::select(-starts_with("xiyi_yiyp")) %>% dplyr::select(-starts_with("xixp_xpyp"))
    output <- list(MR_results = summarized_result_sub, proxy_MR_comparison = prod_diff_result)

  } else output <- NULL

  return(output)


}

summarize_proxyMR_paths <- function(proxyMR_comparison){

  MR_paths_result <- bind_rows(lapply(proxyMR_comparison, function(x) {x[[1]]}))
  return(MR_paths_result)

}

summarize_proxyMR_comparison <- function(proxyMR_comparison, traits_corr2_filled){

  comparison_result <- bind_rows(lapply(proxyMR_comparison, function(x) {x[[2]]}))

  num_result <- dim(comparison_result)[1]/2


  comparison_result <- comparison_result %>%
    mutate(omega_vs_gam_BF_sig = case_when(TRUE ~ omega_vs_gam_pval < 0.05/num_result,
                                                        TRUE ~ TRUE)) %>%
    mutate(omega_vs_rho_BF_sig = case_when(TRUE ~ omega_vs_rho_pval < 0.05/num_result,
                                                        TRUE ~ TRUE)) %>%
    mutate(gam_vs_rho_BF_sig = case_when(TRUE ~ gam_vs_rho_pval < 0.05/num_result,
                                                      TRUE ~ TRUE)) %>%
    mutate(omega_vs_gam_rho_BF_sig = case_when(TRUE ~ omega_vs_gam_rho_pval < 0.05/num_result,
                                              TRUE ~ TRUE))


  comparison_result_plus_cat <- left_join(comparison_result, traits_corr2_filled %>% dplyr::select(Neale_pheno_ID, category) %>% rename(exposure_category = category) %>% unique(), by = c("exposure_ID" = "Neale_pheno_ID")) %>%
    left_join(traits_corr2_filled %>% dplyr::select(Neale_pheno_ID, category) %>% rename(outcome_category = category) %>% unique(), by = c("outcome_ID" = "Neale_pheno_ID"))

  return(comparison_result_plus_cat)

}

summarize_proxyMR_comparison_SNPmeta <- function(proxyMR_comparison, traits_corr2_filled){

  comparison_result <- bind_rows(lapply(proxyMR_comparison, function(x) {x[[2]]}))

  num_result <- dim(comparison_result)[1]


  comparison_result <- comparison_result %>%
    mutate(omega_vs_gam_BF_sig = case_when(TRUE ~ omega_vs_gam_pval < 0.05/num_result,
                                                        TRUE ~ TRUE)) %>%
    mutate(omega_vs_rho_BF_sig = case_when(TRUE ~ omega_vs_rho_pval < 0.05/num_result,
                                                        TRUE ~ TRUE)) %>%
    mutate(gam_vs_rho_BF_sig = case_when(TRUE ~ gam_vs_rho_pval < 0.05/num_result,
                                                      TRUE ~ TRUE)) %>%
    mutate(omega_vs_gam_rho_BF_sig = case_when(TRUE ~ omega_vs_gam_rho_pval < 0.05/num_result,
                                                    TRUE ~ TRUE))

  comparison_result_plus_cat <- left_join(comparison_result, traits_corr2_filled %>% dplyr::select(Neale_pheno_ID, category) %>% rename(exposure_category = category) %>% unique(), by = c("exposure_ID" = "Neale_pheno_ID")) %>%
    left_join(traits_corr2_filled %>% dplyr::select(Neale_pheno_ID, category) %>% rename(outcome_category = category) %>% unique(), by = c("outcome_ID" = "Neale_pheno_ID"))

  data <- comparison_result_plus_cat

  model <- lm(rho_beta ~ gam_beta, data = data)

  data$rho_resid <- resid(model)

  a <- summary(model)$coefficients[2,1]
  b <- summary(model)$coefficients[1,1]
  # below is the same thing as `data$rho_resid`
  # rho_resid <- data$rho_beta - data$gam_beta*a - b
  data$gam_rho_resid <- data$rho_resid + data$gam_beta

  cor_rho_gamma <- cor(data$rho_beta,data$gam_beta)

  data$rho_var <- data$rho_se^2
  data$gam_var <- data$gam_se^2

  data$gam_rho_resid <- data$rho_resid + data$gam_beta
  data$gam_rho_resid_var <- data$rho_var + data$gam_var *(1-a)^2 + 2*data$rho_se*(1-a)*data$gam_se*cor_rho_gamma
  data$gam_rho_resid_se <- sqrt(data$gam_rho_resid_var)

  output <- data %>% mutate(omega_vs_gam_rho_resid_pval = z.test_p(gam_rho_resid, gam_rho_resid_se, omega_beta, omega_se)) %>%
    mutate(omega_vs_gam_rho_resid_BF_sig = case_when(TRUE ~ omega_vs_gam_rho_resid_pval < 0.05/num_result,
                                           TRUE ~ TRUE))

  return(output)

}

prune_proxyMR_comparison_SNPmeta <- function(proxyMR_comparison_summary_yiyp_adj_SNPmeta, household_MR_summary_BF_sig, corr_mat_traits, corr_trait_threshold){

  data <- proxyMR_comparison_summary_yiyp_adj_SNPmeta
  data <- data %>%
    mutate(omega_pval = pnorm(-abs(omega_beta/omega_se))*2) %>%
    mutate(gam_pval = pnorm(-abs(gam_beta/gam_se))*2) %>%
    mutate(rho_pval = pnorm(-abs(rho_beta/rho_se))*2) %>%
    mutate(gam_rho_pval = pnorm(-abs(gam_rho_beta/gam_rho_se))*2) %>%
    left_join(household_MR_summary_BF_sig %>% dplyr::select(exposure_ID, outcome_ID, corr_traits),  by = c("exposure_ID" = "exposure_ID", "outcome_ID" = "outcome_ID")) %>%
    arrange(omega_pval)

  ## prune as follows:
  ## if there is a pair A-B and C-D and if max(corr(A,C)*corr(B,D),corr(A,D)*corr(B,C)) > 0.8

  output_pruned <- data
  counter <- 1

  corr_traits_long <- corr_mat_traits %>% cor_gather()

  while (counter < dim(output_pruned)[1]){

    index_trait_A <- output_pruned$exposure_ID[counter]
    index_trait_B <- output_pruned$outcome_ID[counter]

    index_trait_A_phes_ID <- gsub("_irnt", "", index_trait_A)
    index_trait_B_phes_ID <- gsub("_irnt", "", index_trait_B)

    corr_trait_A <- corr_traits_long %>% filter(var1 == index_trait_A_phes_ID) %>% dplyr::select(var2, cor)
    corr_trait_B <- corr_traits_long %>% filter(var1 == index_trait_B_phes_ID) %>% dplyr::select(var2, cor)

    temp <- output_pruned
    temp <- temp %>%
      mutate(exposure_ID_phes = gsub("_irnt", "", exposure_ID)) %>%
      mutate(outcome_ID_phes = gsub("_irnt", "", outcome_ID)) %>%
      left_join(corr_trait_A, by = c("exposure_ID_phes" = "var2")) %>% rename(cor_A_C = cor) %>%
      left_join(corr_trait_A, by = c("outcome_ID_phes" = "var2")) %>% rename(cor_A_D = cor) %>%
      left_join(corr_trait_B, by = c("exposure_ID_phes" = "var2")) %>% rename(cor_B_C = cor) %>%
      left_join(corr_trait_B, by = c("outcome_ID_phes" = "var2")) %>% rename(cor_B_D = cor) %>%
      mutate(cor_A_C_by_cor_B_D = cor_A_C * cor_B_D) %>%
      mutate(cor_A_D_by_cor_B_C = cor_A_D * cor_B_C) %>%
      mutate(max_cor = pmax(cor_A_C_by_cor_B_D, cor_A_D_by_cor_B_C))

    remove_rows <- which(temp[["max_cor"]] >= corr_trait_threshold)[-1]

    if(length(remove_rows)!=0){
      output_pruned <- output_pruned[-remove_rows,]
    }

    counter <- 1 + counter

  }

  return(output_pruned)


}



create_proxy_prod_comparison_fig_ind <- function(data, exposure_sex, x, y, overlay_var, count){

  custom_col <- c("#a6cee3",
                  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a", "#ffff99", "#b15928")


  # exposure_sex <- "male"
  # data <- comparison_result
  # overlay_var <- "omega_vs_gam_BF_sig_sex_specific"
  # x <- "gam_beta"
  # y <- "omega_beta"

  if(!is.na(exposure_sex)){
    if(exposure_sex!="meta"){
      fig_data <- data[which(data[["exposure_sex"]]==exposure_sex),]
    } else fig_data <- data[which(data[["exposure_sex"]]=="male"),]
  } else fig_data <- data


  fig_data <- fig_data %>% rename(x := !!x) %>% rename(y := !!y) %>% rename(overlay_var := !!overlay_var) %>%
    mutate(x_plot = case_when(x < 0 ~ abs(x), TRUE ~ x)) %>%
    mutate(y_plot = case_when(x < 0 ~ -1*y, TRUE ~ y))



  overlay_data <- fig_data[which(fig_data[["overlay_var"]]==TRUE),]

  if(x == "gam_beta" | x == "gam_meta_beta"){x_label <- "\u03B3"}
  if(x == "rho_beta" | x == "rho_meta_beta"){x_label <- "\u03C1"}
  if(x == "gam_rho_beta" | x == "gam_rho_meta_beta"){x_label <- "\u03B3 + \u03C1"}
  if(x == "gam_rho_resid_beta"){x_label <- expression( paste(rho ["resid"], " + ", gamma))}

  if(y == "omega_beta" | y == "omega_meta_beta"){y_label <- "\u03C9"}
  if(y == "rho_beta" | y == "rho_meta_beta"){y_label <- "\u03C1"}

  xmax <- 1.3
  #if(x == "rho_beta" | x == "rho_meta_beta"){xmax = 0.9}


  ylim <- -1.05
  #if(x == "rho_beta" | x == "rho_meta_beta"){ylim = -0.55}

  ymax <- 1.5
  #if((x == "gam_beta" | x == "gam_meta_beta") & (y == "omega_beta" | y == "omega_meta_beta")){ymax = 1.5}
  #if((x == "gam_beta" | x == "gam_meta_beta") & (y == "rho_beta" | y == "rho_meta_beta")){ymax = 1.05}


  fig <- ggplot(fig_data, aes(x=x_plot, y=y_plot, label=exposure_description, label2=outcome_description, label3=exposure_sex, color=overlay_var)) +
    geom_point(alpha = 3/4) +
    #geom_smooth(mapping = aes(x = x_plot, y = y_plot),inherit.aes = FALSE, method=lm, formula=y~0+x, se=FALSE, color = custom_col[4], fullrange=TRUE) +
    #geom_smooth(data = overlay_data, mapping = aes(x = x_plot, y = y_plot),inherit.aes = FALSE, method=lm, formula=y~0+x, se=FALSE, color = custom_col[2], fullrange=TRUE) +

    theme_half_open(12) +
    scale_fill_manual(values = custom_col) +
    scale_colour_manual(values = custom_col, labels=c("Non-BF significant\ndifference","BF significant difference")) +

    geom_point(data = overlay_data, aes(x=x_plot, y=y_plot), color = custom_col[2]) +
    theme(legend.position="top") + geom_abline(slope=1, intercept=0) + theme(legend.title = element_blank()) +
    #xlab(paste0(x_label, " ", exposure_sex)) + ylab(paste0(y_label, " ", exposure_sex)) +
    xlab(paste0(x_label)) + ylab(paste0(y_label)) +

    scale_x_continuous(breaks = seq(0, xmax, by = 0.2), limits = c(0, xmax)) +
    scale_y_continuous(breaks = seq(-2, ymax, by = 0.5)[seq(-2, ymax, by = 0.5) > ylim], limits = c(ylim, ymax))


  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    fig + theme(legend.box.margin = margin(0, 0, 0, 12))
  )

  fig_no_legend <- fig + theme(legend.position="none")
  #if(count==9){fig_no_legend <- fig_no_legend + theme(legend.position = c(0.5, 0.15))}


  return(list(fig_no_legend = fig_no_legend, legend = legend))

}

create_proxy_sex_comparison_fig_ind <- function(data, var, count){



  custom_col <- c("#a6cee3",
                  "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a", "#ffff99", "#b15928")


  # data <- comparison_result
  # var <- "rho

  num_result <- dim(data)[1]/2


  vars_to_select <- paste0(var, c("_beta", "_sex_het_pval"))
  sex_het_p_var <- paste0(var, "_sex_het_pval")
  beta_var <- paste0(var, "_beta")

  fig_data <- data %>% dplyr::select("exposure_ID", "outcome_ID", "exposure_sex", vars_to_select) %>%
    rename(sex_het_p := !!sex_het_p_var) %>% rename(beta := !!beta_var) %>%
    mutate(sex_het_p_BF = case_when(TRUE ~ sex_het_p < 0.05/num_result,
                                    TRUE ~ TRUE)) %>%
    pivot_wider(names_from = exposure_sex, values_from = beta)

  overlay_data <- fig_data[which(fig_data[["sex_het_p_BF"]]==TRUE),]

  if(var == "gam"){symbol <- "\u03B3"}
  if(var == "rho"){symbol <- "\u03C1"}
  if(var == "omega"){symbol <- "\u03C9"}

  ymin <- -1.15
  ymax <- 1.25

  xmin <- -1.1
  xmax <- 1.25

  fig <- ggplot(fig_data, aes(x = male, y = female, color = sex_het_p_BF)) +
    geom_point(alpha = 3/4) +
    geom_smooth(mapping = aes(x = male, y = female),inherit.aes = FALSE, method=lm, se=FALSE, color = custom_col[4], fullrange=TRUE) +
    geom_smooth(data = overlay_data, mapping = aes(x = male, y = female),
                inherit.aes = FALSE, method=lm, se=FALSE, color = custom_col[2], fullrange=TRUE) +

    theme_half_open(12) +
    scale_fill_manual(values = custom_col) +
    scale_colour_manual(values = custom_col, labels=c("Non-BF significant\ndifference","BF significant\ndifference")) +

    geom_point(data = overlay_data, aes(x = male, y = female), color = custom_col[2]) +
    theme(legend.position="top") + geom_abline(slope=1, intercept=0) + theme(legend.title = element_blank()) +
    xlab(paste0(symbol, " ", "male")) + ylab(paste0(symbol, " ", "female")) +
    scale_y_continuous(breaks = seq(-2, ymax, by = 0.5)[seq(-2, ymax, by = 0.5) > ymin], limits = c(ymin, ymax)) +
    scale_x_continuous(breaks = seq(-2, xmax, by = 0.5)[seq(-2, xmax, by = 0.5) > xmin], limits = c(xmin, xmax))


  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    fig + theme(legend.box.margin = margin(0, 0, 0, 12))
  )

  fig_no_legend <- fig + theme(legend.position="none")
  #if(count==3){fig_no_legend <- fig_no_legend + theme(legend.position = c(0.5, 0.08))}


  return(list(fig_no_legend = fig_no_legend, legend = legend))

}

create_proxy_prod_comparison_fig_all <- function(proxyMR_figure_data){

  figures <- list()

  count <- 1
  for(panel in 1:3){

    if(panel==1){
      x_start <- "gam"
      y_start <- "omega"
      overlay_var_start <- "omega_vs_gam_BF_sig"
    }
    if(panel==2){
      x_start <- "rho"
      y_start <- "omega"
      overlay_var_start <- "omega_vs_rho_BF_sig"
    }
    if(panel==3){
      x_start <- "gam"
      y_start <- "rho"
      overlay_var_start <- "gam_vs_rho_BF_sig"
    }

    for(sex in c("male", "female", "meta")){

      if(sex == "male" | sex=="female"){
        fig_temp <- create_proxy_prod_comparison_fig_ind(proxyMR_figure_data, exposure_sex=sex, x = paste0(x_start, "_beta"), y = paste0(y_start, "_beta"), overlay_var = paste0(overlay_var_start, "_sex_specific"), count)
      } else fig_temp <- create_proxy_prod_comparison_fig_ind(proxyMR_figure_data, exposure_sex=sex, x = paste0(x_start, "_meta_beta"), y = paste0(y_start, "_meta_beta"), overlay_var = paste0(overlay_var_start, "_meta"), count)

      figures[[count]] <- fig_temp[["fig_no_legend"]]
      if(count==1) {legend <- fig_temp[["legend"]]}
      count <- count + 1

    }


  }

  figures_grid <- plot_grid(plotlist = figures, nrow=3, ncol = 3, byrow = F, labels="AUTO", rel_heights = c(1, 1, 1))

  figures_grid_plus_legend <- plot_grid(figures_grid, legend, ncol = 1, nrow =2, rel_heights = c(1, 0.05))


  return(figures_grid_plus_legend)
}


create_proxy_prod_comparison_fig <- function(proxyMR_figure_data){

  figures <- list()

  proxyMR_figure_data %>% rename(gam_rho_resid_beta = gam_rho_resid)
  count <- 1
  for(panel in 1:4){

    if(panel==1){
      x_start <- "gam"
      y_start <- "omega"
      overlay_var_start <- "omega_vs_gam_BF_sig"
    }
    if(panel==2){
      x_start <- "rho"
      y_start <- "omega"
      overlay_var_start <- "omega_vs_rho_BF_sig"
    }
    if(panel==3){
      x_start <- "gam"
      y_start <- "rho"
      overlay_var_start <- "gam_vs_rho_BF_sig"
    }

    if(panel==4){
      x_start <- "gam_rho_resid"
      y_start <- "omega"
      overlay_var_start <- "omega_vs_gam_rho_resid_BF_sig"
    }

    #for(sex in c("male", "female", "meta")){

    fig_temp <- create_proxy_prod_comparison_fig_ind(proxyMR_figure_data, exposure_sex=NA, x = paste0(x_start, "_beta"), y = paste0(y_start, "_beta"), overlay_var = paste0(overlay_var_start), count)

    figures[[count]] <- fig_temp[["fig_no_legend"]]
    if(count==1) {legend <- fig_temp[["legend"]]}
    count <- count + 1

    #}


  }

  figures_grid <- plot_grid(plotlist = figures, nrow=2, ncol = 2, byrow = F, labels="AUTO", rel_heights = c(1, 1, 1))

  figures_grid_plus_legend <- plot_grid(figures_grid, legend, ncol = 1, nrow =2, rel_heights = c(1, 0.05))


  return(figures_grid_plus_legend)
}

create_proxy_sex_comparison_fig <- function(proxyMR_figure_data){

  panel <- 1
  sex_het_figures <- list()
  for(panel in 1:3){

    if(panel==1){
      var <- "gam"

    }
    if(panel==2){
      var <- "rho"

    }
    if(panel==3){
      var <- "omega"

    }

    fig_temp <- create_proxy_sex_comparison_fig_ind(proxyMR_figure_data, var = var, panel)
    sex_het_figures[[panel]] <- fig_temp[["fig_no_legend"]]


    if(panel==1) {legend <- fig_temp[["legend"]]}
    panel <- panel + 1

  }

  sex_het_figures_grid <- plot_grid(plotlist = sex_het_figures, nrow=1, ncol = 3, byrow = F, labels="AUTO", rel_heights = c(1, 1, 1))

  figures_grid_plus_legend <- plot_grid(sex_het_figures_grid, legend, ncol = 1, nrow =2, rel_heights = c(1, 0.03))


  return(figures_grid_plus_legend)
}

find_MV_z <- function(household_MR_summary_BF_sig, standard_MR_summary){

  #exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- household_MR_summary_BF_sig %>% filter(exposure_ID!=outcome_ID)

  MR_sub$MV_z <- NA
  exposure_ID_prev <- ""
  outcome_ID_prev <- ""
  for(i in 1:dim(MR_sub)[1]){


    exposure_ID <- MR_sub$exposure_ID[i]
    outcome_ID <- MR_sub$outcome_ID[i]
    exposure_sex <- MR_sub$exposure_sex[i]

    cat(paste0("Identifying Z's MV MR for  `", exposure_ID, "` as exposure and `", outcome_ID, "` as outcome in ", exposure_sex, "s as exposure sex.\n\n"))

    if(exposure_ID!=exposure_ID_prev | outcome_ID!=outcome_ID_prev){
      z_vs_x <- standard_MR_summary %>% filter(exposure_ID==!!exposure_ID) %>% filter(outcome_ID!=!!outcome_ID)
      y_vs_z <- standard_MR_summary %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID!=!!exposure_ID)

      num_tests <- (dim(z_vs_x)[1] + dim(y_vs_z)[1])/2

      z_vs_x_BF <- z_vs_x %>% filter(IVW_meta_pval < 0.05/num_tests) %>% pull(outcome_ID)
      y_vs_z_BF <- y_vs_z %>% filter(IVW_meta_pval < 0.05/num_tests) %>% pull(exposure_ID)


      union_z <- intersect(z_vs_x_BF, y_vs_z_BF)

      y_vs_z_order <- standard_MR_summary %>% filter(outcome_ID==!!outcome_ID) %>% filter(exposure_ID!=!!exposure_ID) %>% arrange(IVW_meta_pval) %>% pull(exposure_ID) %>% unique()

      prioritze_z <- y_vs_z_order[which(y_vs_z_order %in% union_z)]

    }

    exposure_ID_prev <- exposure_ID
    outcome_ID_prev <- outcome_ID

    MR_sub$MV_z[i] <- list(prioritze_z)

  }

  return(MR_sub)
}

corr_filter_MV_z <- function(MV_z_data, corr_mat_traits, z_prune_threshold){

  MV_z_data$MV_z_prune <- NA
  exposure_ID_prev <- ""
  outcome_ID_prev <- ""

  for(i in 1:dim(MV_z_data)[1]){


    exposure_ID <- MV_z_data$exposure_ID[i]
    outcome_ID <- MV_z_data$outcome_ID[i]
    exposure_sex <- MV_z_data$exposure_sex[i]
    cat(paste0("Pruning Z's MV MR for  `", exposure_ID, "` as exposure and `", outcome_ID, "` as outcome in ", exposure_sex, "s as exposure sex.\n\n"))


    # each exposure_ID/outcome_ID is in data twice (one for each sex), since we are prioritzing by meta-anlyzed results with y, we don't need to find z's for each sex - they will be the same.
    if(exposure_ID!=exposure_ID_prev | outcome_ID!=outcome_ID_prev){

      z_list <- MV_z_data$MV_z[i][[1]]
      all_phenos <- c(z_list, exposure_ID, outcome_ID)
      all_phenos_fix <- gsub("_irnt", "", all_phenos)
      which_fixed <- which(str_detect(all_phenos, "_irnt"))
      corr_mat_z <- corr_mat_traits[all_phenos_fix, all_phenos_fix]

      # first remove any z's that are in high correaltion with either x or y

      corr_mat_z_prune <- corr_mat_z

      row_x <- which(colnames(corr_mat_z)==gsub("_irnt", "", exposure_ID))
      row_y <- which(colnames(corr_mat_z)==gsub("_irnt", "", outcome_ID))


      remove_rows <- c(which(corr_mat_z[row_x,] > z_prune_threshold), which(corr_mat_z[row_y,] > z_prune_threshold))
      remove_rows <- remove_rows[-which(remove_rows %in% c(row_x, row_y))]

      if(length(remove_rows)!=0){
        corr_mat_z_prune <- corr_mat_z_prune[-remove_rows, -remove_rows]
      }

      j <- 1
      current_col <- colnames(corr_mat_z)[j]
      while(current_col != gsub("_irnt", "", exposure_ID)){

        remove_rows <- c(which(corr_mat_z_prune[j,] > z_prune_threshold))
        remove_rows <- remove_rows[-which(remove_rows %in% c(j))]

        if(length(remove_rows)!=0){
          corr_mat_z_prune <- corr_mat_z_prune[-remove_rows, -remove_rows]
        }
        j <- j+1
        current_col <- colnames(corr_mat_z_prune)[j]

      }

      pruned_pheno <- colnames(corr_mat_z_prune)
      prune_pheno_to_fix <- which(which(all_phenos_fix %in% pruned_pheno) %in% which_fixed)
      pruned_pheno[prune_pheno_to_fix] <- paste0(pruned_pheno[prune_pheno_to_fix], "_irnt")

      pruned_pheno_z <- pruned_pheno[-which(pruned_pheno %in% c(exposure_ID, outcome_ID))]

    }

    exposure_ID_prev <- exposure_ID
    outcome_ID_prev <- outcome_ID



    MV_z_data$MV_z_prune[i] <- list(pruned_pheno_z)

  }

  return(MV_z_data)

}


pull_z_summ_stats <- function(MV_z_data){

  output_list <- list()
  for(i in 1:dim(MV_z_data)[1]){

    exposure_ID <- MV_z_data$exposure_ID[i]
    outcome_ID <- MV_z_data$outcome_ID[i]
    exposure_sex <- MV_z_data$exposure_sex[i]
    z_list <- MV_z_data$MV_z_prune[i][[1]]
    all_phenos <- c(z_list, exposure_ID, outcome_ID)
    data_i <- tibble()

    cat(paste0("Extracting summary statistics for MV MR for  `", exposure_ID, "` as exposure and `", outcome_ID, "` as outcome in ", exposure_sex, "s as exposure sex.\n"))

    for(k in all_phenos){

      if(k != outcome_ID){
        data_k <- tibble(IVs_from = k)
        for(j in all_phenos){

          GWAS_file <- paste0("analysis/traitMR/standard_GWAS/", j, "/", j, "_vs_", k, "_GWAS.csv")
          ## this will pull the GWAS results for phenotype `j` for IVs from `k`
          GWAS_results_j_vs_k <- fread(GWAS_file, data.table = F)
          GWAS_results_j_vs_k_sex <- GWAS_results_j_vs_k[which(GWAS_results_j_vs_k$sex==exposure_sex),]
          ## get the results for both sexes for pruning
          GWAS_results_j_vs_k_both <- GWAS_results_j_vs_k[which(GWAS_results_j_vs_k$sex=="both_sexes"),c("SNP", "pval")]
          colnames(GWAS_results_j_vs_k_both) <- c("SNP", "pval_both_sexes")
          GWAS_results_j_vs_k_sex <- merge(GWAS_results_j_vs_k_sex, GWAS_results_j_vs_k_both) %>% as_tibble %>% type_convert() %>% mutate_at(c("exposure_ID", "outcome_ID"), as.character)
          data_k[[paste0("GWAS_", j, "_results")]] <- list(GWAS_results_j_vs_k_sex)

        }
        data_i <- bind_rows(data_i, data_k)

      }

    }
    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_", exposure_sex)]] <- (data_i)
    cat(paste0("Finished extraction summary statistics for ", i, " of ", dim(MV_z_data)[1], " traits.\n\n"))


  }
  return(output_list)

}

prune_z_summ_stats <- function(MV_z, z_summ_stats, prune_threshold){

  output_list <- list()
  for(i in 1:length(z_summ_stats)){

    exposure_ID <- MV_z$exposure_ID[i]
    outcome_ID <- MV_z$outcome_ID[i]
    exposure_sex <- MV_z$exposure_sex[i]

    z_summ_stats_i <- z_summ_stats[[i]]
    cat(paste0("Pruning summary statistics for MV MR for  `", exposure_ID, "` as exposure and `", outcome_ID, "` as outcome in ", exposure_sex, "s as exposure sex.\n"))

    # extract all SNPs
    col <- which(colnames(z_summ_stats_i)==paste0("GWAS_", exposure_ID, "_results"))
    SNPs_i <- numeric()
    chr_i <- numeric()
    pvals_i <- numeric()

    for(k in 1:dim(z_summ_stats_i)[1]){
      SNPs_k <- z_summ_stats_i[k, col][[1]][[1]]$SNP
      chr_k <- z_summ_stats_i[k, col][[1]][[1]]$chr
      pvals_k <- z_summ_stats_i[k, col][[1]][[1]]$pval_both_sexes

      SNPs_i <- c(SNPs_i, SNPs_k)
      chr_i <- c(chr_i, chr_k)
      pvals_i <- c(pvals_i, pvals_k)

    }

    ## This should be based on meta-analyzed pvals
    to_prune_mat <- tibble(SNP = SNPs_i,
                        chr = chr_i,
                        pval = pvals_i) %>% unique()

    pruned_mat <- IV_clump(to_prune_mat, prune_threshold)

    full_dat <- tidyr::unnest(z_summ_stats_i, starts_with("GWAS_"), names_sep="_")

    ## remove duplicate columns, like SNP, chr, allele, etc.
    full_dat_shrink <- full_dat[!duplicated(as.list(full_dat))]

    SNP_col <- which(grepl("_results_SNP", colnames(full_dat_shrink)))

    ## restrict to only pruned SNPs
    full_dat_sub <- full_dat_shrink[which(full_dat_shrink[[SNP_col]] %in% pruned_mat$SNP),-1]

    full_dat_pruned <- full_dat_sub %>% unique()


    bx_dat <- list()
    bxse_dat <- list()
    mv_dat_cols <- colnames(z_summ_stats_i)[-1]
    by_col <- which(mv_dat_cols==paste0("GWAS_", outcome_ID, "_results"))
    by_col_name <- mv_dat_cols[by_col]
    bx_cols <- mv_dat_cols[-by_col]

    count <- 1
    for(mv_name in bx_cols){

      beta_col <- which(colnames(full_dat_pruned)==paste0(mv_name, "_beta"))
      se_col <- which(colnames(full_dat_pruned)==paste0(mv_name, "_se"))
      bx_dat[[count]] <- full_dat_pruned[,beta_col]
      bxse_dat[[count]] <- full_dat_pruned[,se_col]
      count <- count + 1

    }

    by_beta_col <- which(colnames(full_dat_pruned)==paste0("GWAS_", outcome_ID, "_results", "_beta"))
    by_se_col <- which(colnames(full_dat_pruned)==paste0("GWAS_", outcome_ID, "_results", "_se"))


    mv_data_format <- mr_mvinput(bx = as.matrix(data.frame(bx_dat)), bxse = as.matrix(data.frame(bxse_dat)),
                                 by = unlist(full_dat_pruned[,by_beta_col]), byse = unlist(full_dat_pruned[,by_se_col]))


    #mv_result <- mr_mvivw(mv_data_format) # see `str(mv_result)` to know how to access the results

    output_list[[paste0(outcome_ID, "_vs_", exposure_ID, "_", exposure_sex)]] <- (mv_data_format)
    cat(paste0("Finished pruning summary statistics for ", i, " of ", length(z_summ_stats), " traits.\n\n"))



  }

  return(output_list)

}


IV_clump <- function(data_IV, prune_threshold){
  data_prune <- numeric()
  for(chr in 1:22)
  {
    data_IV_temp <- data_IV[which(data_IV[["chr"]]==chr),]
    fn <- tempfile(tmpdir = tempdir())
    snps <- data_IV_temp$SNP
    if(length(snps)==0) next
    pvals <- data_IV_temp$pval
    write.table(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F)
    refdat=paste0("/data/sgg2/jenny/data/1000G/chr",chr,"/1000G_EUR_chr",chr,"_filt")

    snp_clump <- plink_clump(bfile = refdat, filename = fn, prune_threshold = prune_threshold)

    data_prune_temp <- data_IV_temp[which(data_IV_temp$SNP %in% snp_clump),]
    data_prune <- rbind(data_prune, data_prune_temp)
  }
  return(data_prune)
}

plink_clump <- function(bfile, filename, clump_p1 = 1, clump_p2 = 1, prune_threshold = 0.001, clump_kb = 10000, threads = 1){

  fun2 <- paste0(
    "plink",
    " --bfile ", bfile,
    " --clump ", filename,
    " --clump-p1 ", clump_p1,
    " --clump-p2 ", clump_p2,
    " --clump-r2 ", prune_threshold,
    " --clump-kb ", clump_kb,
    " --threads ", threads,
    " --out ", filename
  )
  system(fun2)
  a <- read.table(paste(filename, ".clumped", sep=""), he=T)
  unlink(paste(filename, "*", sep=""))
  a_out <- a[,"SNP"]
  return(a_out)
}

find_proxyMR_IV_overlap <- function(exposure_info, proxyMR_paths_summary, LD_threshold = 0.9){

  exposure_ID <- exposure_info %>% filter(Value=="trait_ID") %>% pull(Info)
  ## only run this for those where omega is significant and where exposure and outcome ID are different.
  MR_sub <- proxyMR_paths_summary %>% filter(exposure_ID==!!exposure_ID)
  MR_sub$YiXi_IV_sig_overlap <- NA
  MR_sub$YiXi_IV_exact_overlap <- NA

  summarized_result <- as_tibble(numeric())

  if(dim(MR_sub)[1]!=0){

    outcome_ID_prev <- ""
    for(i in 1:dim(MR_sub)[1]){

      outcome_ID <-  MR_sub$outcome_ID[i]
      exposure_sex <- "male" #MR_sub$exposure_sex[i]

      if(outcome_ID!=outcome_ID_prev){

        ## same IVs are used for both sexes so it doens't matter which sex we use.

        cat(paste0("Assessing IV overlap between `", exposure_ID, "` as exposure and `", outcome_ID, "` as outcome.\n\n"))
        IV_file_exposure <- paste0("analysis/traitMR/IVs/Neale/", exposure_ID, "/", exposure_sex, "_IVs.txt")
        IV_data_exposure <- fread(IV_file_exposure, data.table = F)

        IV_file_outcome <- paste0("analysis/traitMR/IVs/Neale/", outcome_ID, "/", exposure_sex, "_IVs.txt")
        IV_data_outcome <- fread(IV_file_outcome, data.table = F)

        exposure_snps <- IV_data_exposure$rsid
        outcome_snps <- IV_data_outcome$rsid

        exact_overlap <- length(which(outcome_snps %in% exposure_snps))

        signal_overlap <- 0

        for(chr in 1:22)
        {

          exposure_snps_chr <- IV_data_exposure[which(IV_data_exposure$chr==chr),"rsid"]
          outcome_snps_chr <- IV_data_outcome[which(IV_data_outcome$chr==chr),"rsid"]
          exact_overlap_chr <- length(which(outcome_snps_chr %in% exposure_snps_chr))

          fn <- tempfile(tmpdir = tempdir())
          snps <- unique(c(exposure_snps_chr, outcome_snps_chr))

          if(length(outcome_snps_chr)<=1 | length(exposure_snps_chr)==0) next

          write.table(data.frame(SNP=snps), file=fn, row=F, col=F, qu=F)

          refdat=paste0("/data/sgg2/jenny/data/1000G/chr",chr,"/1000G_EUR_chr",chr,"_filt")

          r2_dat <- plink_r2_dat(bfile = refdat, filename = fn)


          if(length(which(outcome_snps_chr %in% c(r2_dat$SNP_A, r2_dat$SNP_B)))==0){
            signal_overlap_chr <- exact_overlap_chr
          } else{

            signal_overlap_chr <- exact_overlap_chr
            for(snp_outcome in outcome_snps_chr){

              rows_y <- which(r2_dat$SNP_A %in% snp_outcome | r2_dat$SNP_B %in% snp_outcome)
              if(length(rows_y)!=0){

                r2_dat_sub <- r2_dat[rows_y,]
                rows_x <- which(r2_dat_sub$SNP_A %in% exposure_snps_chr | r2_dat_sub$SNP_B %in% exposure_snps_chr)
                r2_dat_sub2 <- r2_dat_sub[rows_x,]
                signal_overlap_chr <- exact_overlap_chr + length(which(r2_dat_sub2$R2 > LD_threshold))

              }

            }

          }

          signal_overlap <- signal_overlap + signal_overlap_chr

        }

        percent_overlap <- signal_overlap / length(outcome_snps)

      }

      MR_sub$YiXi_IV_sig_overlap[i] <- percent_overlap
      MR_sub$YiXi_IV_exact_overlap[i] <- exact_overlap / length(outcome_snps)
      outcome_ID_prev <- outcome_ID

    }
    output <- MR_sub

  } else output <- NULL


  return(output)

}


plink_r2_dat <- function(bfile, filename, threads = 1){

  # If you want the r2 matrix

  # fun2 <- paste0(
  #   "plink",
  #   " --recodeAD ",
  #   " --extract ", filename, # create a list of snps to calculated r2 matrix between 'mysnps.txt', see: https://zzz.bwh.harvard.edu/plink/ld.shtml
  #   " --bfile ", bfile,
  #   " --threads ", threads,
  #   " --out ", filename
  # )
  #

  # d <- read.table(paste0(filename, ".raw"),header=T)
  # remove_cols <- c(1:6, which(endsWith(colnames(d), "_HET")))
  # d_sub <- d[,-remove_cols]
  #
  # r2_mat <- cor(d_sub, use = "complete.obs")^2

  fun2 <- paste0(
    "plink",
    " --r2 ",
    " --extract ", filename, # create a list of snps to calculate r2 matrix between 'mysnps.txt', see: https://zzz.bwh.harvard.edu/plink/ld.shtml
    " --bfile ", bfile,
    " --threads ", threads,
    " --out ", filename,
    " --ld-window-kb 10000",
    " --ld-window 99999"

  )

  system(fun2)

  t <- read.table(paste0(filename, ".ld"), he=T)

  return(t)
}


#saveRDS(out, "code/shiny/mr_summary.rds")


household_MR_plot <-  function (dat, original_MR){

  exposure_ID <- dat$exposure[1]
  outcome_ID <- dat$outcome[1]
  exposure_sex <- dat$exposure_sex[1]
  outcome_sex <- dat$outcome_sex[1]


  exposure_description <- dat$exposure_description[1]
  outcome_description <- dat$outcome_description[1]

  mr_title <- bquote(atop(.(paste0("Estimate of the assortative mating effect of ")),
                          italic(.(outcome_description)) ~ 'on' ~ italic(.(exposure_description)) ~ .(paste0(' for ', exposure_sex, "s on ", outcome_sex,"s"))))



  xlab <- paste0("SNP effect of ", exposure_ID, " in ", exposure_sex, "s (general population)")
  ylab <- paste("SNP effect of ", outcome_ID, " in ",
                outcome_sex, "partner")

  if(is.null(exposure_sex)){
    mr_title <- bquote(atop(.(paste0("Estimate of the assortative mating effect of ")),
                            italic(.(outcome_description)) ~ 'on' ~ italic(.(exposure_description))))

    xlab <- paste0("SNP effect of ", exposure_ID, " in (general population)")
    ylab <- paste("SNP effect of ", outcome_ID, " in partner")

  }

  mr_results <- original_MR

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
                                                                                                                     x = xlab, y = ylab) + ggplot2::theme(legend.position = "bottom",
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


ivt <- function(x){
  out <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(out)
}


order_bgen <- function(bgen_file, data){

  colnames(bgen_file) <- c("ID_1", "ID_2", "missing")
  ord = match(bgen_file$ID_1[-1], data$eid)

  bgen_merge <- left_join(bgen_file, data, by=c("ID_1" = "userId"))

  add_na <- bgen_merge %>% mutate_at(vars(-ID_1, -ID_2, -missing), ~replace_na(., -999))

  output <- add_na %>% dplyr::select(-ID_1, -ID_2, -missing)
  return(output)
}

prep_PC_GWAS <- function(data_id_age, data_id_sex, sqc_munge, bgen_file = data_UKBB_sample){

  merge_sex_age <- merge(data_id_age, data_id_sex)
  out_list <- list()
  for(PC in 1:40){

    PC_i_data <- sqc_munge[,c("ID", paste0("PC_", PC))]
    colnames(PC_i_data) <- c("userId", "PCi")
    resid_data <- merge(PC_i_data, merge_sex_age)
    resid_data <- resid_data %>% mutate(PCi_ivt = ivt(PCi))

    for(sex_var in c(0,1)){

      resid_data_sub <- resid_data[which(resid_data$sex==sex_var),]
      sex <- ifelse(sex_var==0, "female", "male")
      PCi_resid <- resid(lm(PCi_ivt ~ ., data = subset(resid_data_sub, select=c( -userId, -PCi, -sex) ), na.action = na.exclude))
      out <- as.data.frame(cbind(resid_data_sub[["userId"]], PCi_resid))
      colnames(out) <- c("userId", paste0("PC_", PC, "_", sex, "_resid"))
      out_list[[paste0("PC_", PC, "_", sex, "_resid")]] <- out

    }


    sex <- "both_sexes"
    PCi_resid <- resid(lm(PCi_ivt ~ ., data = subset(resid_data, select=c( -userId, -PCi) ), na.action = na.exclude))
    out <- as.data.frame(cbind(resid_data[["userId"]], PCi_resid))
    colnames(out) <- c("userId", paste0("PC_", PC, "_", sex, "_resid"))
    out_list[[paste0("PC_", PC, "_", sex, "_resid")]] <- out


  }

  full_output <- reduce(out_list, full_join)


  ordered_output <- order_bgen(bgen_file = bgen_file, full_output)

  return(ordered_output)


}


write_PC_gwas_input <- function(PC_gwas_input){
  write.table(PC_gwas_input, "data/processed/PC_UKBB_GWAS_input", sep=" ", quote=F, row.names=F, col.names = T)
  return("data/processed/PC_UKBB_GWAS_input")
}


create_UKBB_v2_snp_file_list <- function(UKBB_processed){
  chr_num = tibble(chr = 1:22)
  output <- chr_num %>% mutate(v2_snp_list_files = paste0(UKBB_processed, "/v2_snp_list/", "snp_list_chr", chr, ".txt"))

  return(output$v2_snp_list_files)
}


make_ukbb_chunks <- function(v2_snp_list_file, chunk_size=1e6){

  v2_snp_list <- fread(v2_snp_list_file, data.table=F)
  min_positon <- min(v2_snp_list$position)
  max_position <- max(v2_snp_list$position)
  chr <- as.character(v2_snp_list$chr[1])
  if(nchar(as.character(v2_snp_list$chr[1]))==1){
    chr_char <- paste0("0", as.character(v2_snp_list$chr[1])) } else chr_char <- as.character(v2_snp_list$chr[1])

  out <- numeric()

  basepositions <- sample(v2_snp_list$position)

  num_chunks <- ceiling(length(basepositions)[1]/chunk_size)

  chunks <- split(basepositions, ceiling(rank(basepositions)/chunk_size))

  for(i in 1:length(chunks)){
    chunk_num <- i
    start <- min(chunks[[i]])
    end <- max(chunks[[i]])

    out_i <- cbind(chunk_num, chr, chr_char, start, end)
    out <- rbind(out, out_i)
  }


  out <- as_tibble(out) %>% mutate_all(as.character)

}

create_bgenie_GWAS_dir <- function(){

  dir.create("analysis/bgenie_GWAS", showWarnings = FALSE)
  return("analysis/bgenie_GWAS")
}

launch_bgenie <- function(chr, phenofile, UKBB_dir, chr_char, start_pos, end_pos, chunk_num, output_dir, output_prefix){


  output_file <- paste0(output_dir, "/", output_prefix, "_chr", chr, "_chunk", chunk_num, ".out")
  cat(paste0("Running chr ", chr, ", chunk ", chunk_num, ".\n"))
  system(paste0("/data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 ",
                "--bgen ", UKBB_dir, "/imp/_001_ukb_imp_chr", chr, "_v2.bgen ",
                "--pheno ", phenofile, " ",
                "--range ", chr_char, " ", start_pos, " ", end_pos, " ",
                "--pvals --out ", output_file))
  if(!file.exists(paste0(output_file, ".gz"))){
    file.create(paste0(output_file, ".gz"))
  }

  return(paste0(output_file, ".gz"))

  # Eleonora's command line
  # /data/sgg3/jonathan/bgenie_v1.3/bgenie_v1.3_static1 --bgen
  # /data/sgg3/eleonora/projects/UKBB_GWAS/UK10K_SNPrs/CHR11/chr11.bgen ## this script would not include SNPs specific to HRC_list
  # --pheno ../phenofile --pvals --out chr11.out
}


unzip_bgenie <- function(bgenie_file){
  # system(paste0("gzip -dk analysis/GWAS/UKBB/chr", chr, ".out.gz")) ## keep flag doesn't exist on this system

  unzipped_file <- gsub(".gz", "", bgenie_file)
  #cancel_if(file.exists(file))

  system(paste0("gunzip < ",bgenie_file, " > ", unzipped_file))
  return(unzipped_file)
}


process_bgenie <- function(directory, extension=".out", HRC_panel){

  current_dir <- getwd()
  setwd(directory)
  system(paste0("tail -q -n +2 *", extension, " > temp"))

  header_file <- fread(list.files(pattern=paste0("\\", extension, "$"))[1], nrows=2, data.table=F)
  bgenie_columns <- colnames(header_file)
  full_data <- fread("temp", data.table=F, header=F)
  setwd(current_dir)
  colnames(full_data) <- bgenie_columns

  HRC_list <- fread(HRC_panel, data.table=F)
  HRC_filtered_data <- full_data %>% filter(rsid %in%  HRC_list$ID)
  return(HRC_filtered_data)

}

mr_scatter_plot_custom <-  function (mr_results, dat, mr_title, exposure_sex, outcome_sex)
{
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


mr_plot_binned_sex_specific <- function(harmonised_data, MR_binned, exposure_sex, group){

  if(exposure_sex=="male"){outcome_sex="female"}
  if(exposure_sex=="female"){outcome_sex="male"}

  harmonised_data.sex_spec <- harmonised_data[[paste0("exp_", exposure_sex, "_harmonised_data_filter")]]
  bin_summary <- MR_binned %>% dplyr::filter(grouping_var==!!group) %>% dplyr::filter(exposure_sex==!!exposure_sex) %>% dplyr::filter(bin != "all") %>% mutate_at('bin', as.factor)

  harmonise_dat_plot <- harmonised_data.sex_spec %>% dplyr::filter(grouping_var == !!group) %>% dplyr::filter(bin != "all") %>% mutate_at('bin', as.factor)

  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])

  outcome_levels <- levels(harmonise_dat_plot$bin)
  outcome_levels_num <- str_first_number(outcome_levels)

  harmonise_dat_plot$bin <- factor(harmonise_dat_plot$bin, levels = outcome_levels[order(outcome_levels_num)])

  bin_summary$bin <- factor(bin_summary$bin, levels = outcome_levels[order(outcome_levels_num)])



  legend_title <- ifelse(group=="age_even_bins", "Median age \n of couples (years)", "Time together \n in same household (years)")

  index <- harmonise_dat_plot$beta.exposure < 0
  harmonise_dat_plot$beta.exposure[index] <- harmonise_dat_plot$beta.exposure[index] *-1
  harmonise_dat_plot$beta.outcome[index] <- harmonise_dat_plot$beta.outcome[index] *-1


  plot <- ggplot2::ggplot(data = harmonise_dat_plot, ggplot2::aes(x = beta.exposure,
                                                                  y = beta.outcome)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                           colour = "grey", width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(bin))) +
    ggplot2::geom_abline(data = bin_summary, ggplot2::aes(intercept = 0,
                                                          slope = IVW_beta, colour = factor(bin)), show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = brewer.pal(n = 9, name = "Blues")[-c(1:2)]) +
    theme_minimal() +
    ggplot2::labs(colour = legend_title, x = paste0("SNP effect on ", tolower(trait_description), "\n(estimated in ", exposure_sex, "s)"),
                  y = paste0("SNP effect on ", outcome_sex, " partner")) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))


  return(plot)

}

mr_plot_sex_specific <- function(harmonised_data, MR_binned, MR_binned_meta, custom_col){

  harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])


  ## which grouping_var we choose doesn't matter
  bin_summary <- MR_binned %>% dplyr::filter(grouping_var=="age_even_bins") %>% dplyr::filter(bin == "all") %>% mutate_at('exposure_sex', as.factor)
  bin_summary_meta <- MR_binned_meta %>% dplyr::filter(grouping_var=="age_even_bins") %>% dplyr::filter(bin == "all")


  outcome_levels <- levels(bin_summary$exposure_sex)

  legend_title <- "Exposure sex"

  index <- harmonise_dat_plot$beta.exposure < 0
  harmonise_dat_plot$beta.exposure[index] <- harmonise_dat_plot$beta.exposure[index] *-1
  harmonise_dat_plot$beta.outcome[index] <- harmonise_dat_plot$beta.outcome[index] *-1



  plot <- ggplot2::ggplot(data = harmonise_dat_plot, ggplot2::aes(x = beta.exposure,
                                                                  y = beta.outcome)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                           colour = "grey", width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(exposure_sex))) +
    ggplot2::geom_abline(data = bin_summary, ggplot2::aes(intercept = 0,
                                                          slope = IVW_beta, colour = factor(exposure_sex)), show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = custom_col[c(2,4)]) +

    theme_minimal() +
    ggplot2::labs(colour = legend_title, x = paste0("SNP effect on ", tolower(trait_description), "\n(sex-specific)"),
                  y = paste0("SNP effect on ", "partner\n(sex-specific)")) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::geom_abline(intercept = 0, slope = bin_summary$IVW_meta_beta, color="black", size=1)   ### FIX THIS LINE.


  return(plot)
}

forest_plot_binned_sex_specific <- function(harmonised_data, MR_binned, group, custom_col){

  harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])


  bin_summary <- MR_binned %>% dplyr::filter(grouping_var==!!group) %>% dplyr::filter(bin != "all")

  outcome_levels <- levels(bin_summary$exposure_sex)

  xlab <- paste0("AM MR estimate for ", tolower(trait_description), "\n(sex-specific)")

  xmin <- min(bin_summary$IVW_beta - bin_summary$IVW_se)
  xmax <- max(bin_summary$IVW_beta + bin_summary$IVW_se)

  xmin_lab <- floor(xmin*10)/10
  xmax_lab <- ceiling(xmax*10)/10


  fp_data <-  bin_summary %>% mutate(upper = IVW_beta + IVW_se) %>% mutate(lower = IVW_beta - IVW_se)

  bins <- c(bin_summary %>% filter(outcome_sex == "male") %>% pull(bin))
  first_label <- ifelse(group=="age_even_bins", "Median age \nof couples (years)", "Time together \nin same household (years)")
  bin_labels <- c(first_label, "")

  for(i in 1:length(bins)){
    bin <- bins[i]
    bin_labels <- c(bin_labels, bin, "   Female", "   Male")

  }


  mean <- c(NA, NA)
  upper <- c(NA, NA)
  lower <- c(NA, NA)

  for(bin in bins){

    mean <- unlist(c(mean, NA, fp_data[which(fp_data$bin==bin),"IVW_beta"]))
    upper <- unlist(c(upper, NA, fp_data[which(fp_data$bin==bin),"IVW_beta"] + fp_data[which(fp_data$bin==bin),"IVW_se"]))
    lower <- unlist(c(lower, NA, fp_data[which(fp_data$bin==bin),"IVW_beta"] - fp_data[which(fp_data$bin==bin),"IVW_se"]))

  }

  MR_text <- c("MR estimate (95% CI)", "")
  for(bin in bins){

    means <- unlist(fp_data[which(fp_data$bin==bin),"IVW_beta"])
    ses <- unlist(fp_data[which(fp_data$bin==bin),"IVW_se"])

    MR_text <- c(MR_text, "")
    for(i in 1:length(means)){
      MR_text <- c(MR_text, make_beta_95ci(means[i], ses[i]))

    }

    MR_text <- c(MR_text)

  }

  upper <- unlist(c(upper, fp_data[which(fp_data$bin==bin),"IVW_beta"] + fp_data[which(fp_data$bin==bin),"IVW_se"], NA))
  lower <- unlist(c(lower, fp_data[which(fp_data$bin==bin),"IVW_beta"] - fp_data[which(fp_data$bin==bin),"IVW_se"], NA))


  summary_rows <- rep(FALSE, length(bin_labels))
  index <- which(bin_labels!="")
  summary_rows[index] <- TRUE



  fn <- local({
    i = 0
    no_lines <- sum(!is.na(mean))
    b_clrs=rep(c(custom_col[2], custom_col[4]), length(which(!is.na(mean)))/2)
    l_clrs = rep("grey", length(which(!is.na(mean))))
    function(..., clr.line, clr.marker){
      i <<- i + 1
      fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    }
  })

  fp_plot <- forestplot(
    labeltext = cbind(bin_labels, MR_text), fn.ci_norm = fn,
    mean = mean,
    upper = upper,
    lower = lower,
    is.summary = summary_rows,
    boxsize = .3, # We set the box size to better visualize the type

    col = fpColors(box = custom_col[c(2,4)],  lines=c("gray50", "gray50")), ci.vertices=TRUE,
    hrzl_lines=list("3" = gpar(lwd=1, col="#99999922"),
                    "6"  = gpar(lwd=1, col="#99999922"),
                    "9"  = gpar(lwd=1, col="#99999922"),
                    "12"  = gpar(lwd=1, col="#99999922"),
                    "15"  = gpar(lwd=1, col="#99999922")),

    #col=fpColors(box="black", lines="black", zero = "gray50"),
    xlab = paste0("\n",xlab),


    txt_gp = fpTxtGp(xlab = gpar(cex = 0.9),
                     ticks = gpar(cex = 0.8),
                     label = gpar(cex = 0.9)),

    legend_args = fpLegend(pos = "right"), xticks = seq(xmin_lab, xmax_lab, by=0.1))

  return(fp_plot)

}

forest_plot_binned_meta <- function(harmonised_data, MR_binned, group, custom_col){

  harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])


  ## which exposure_sex we choose doesn't matter
  bin_summary <- MR_binned %>% dplyr::filter(grouping_var==!!group) %>% dplyr::filter(bin != "all") %>% dplyr::filter(exposure_sex=="male")

  outcome_levels <- levels(bin_summary$exposure_sex)

  xlab <- paste0("AM MR estimate for ", tolower(trait_description))

  xmin <- min(bin_summary$IVW_meta_beta - bin_summary$IVW_meta_se)
  xmax <- max(bin_summary$IVW_meta_beta + bin_summary$IVW_meta_se)

  xmin_lab <- floor(xmin*10)/10
  xmax_lab <- ceiling(xmax*10)/10


  fp_data <-  bin_summary %>% mutate(upper = IVW_meta_beta + IVW_meta_se) %>% mutate(lower = IVW_meta_beta - IVW_meta_se)

  bins <- c(bin_summary %>% pull(bin))
  first_label <- ifelse(group=="age_even_bins", "Median age\n of couples (years)", "Time together\nin same household (years)")
  bin_labels <- c(first_label)

  for(i in 1:length(bins)){
    bin <- bins[i]
    bin_labels <- c(bin_labels, bin)

  }


  mean <- c(NA)
  upper <- c(NA)
  lower <- c(NA)

  for(bin in bins){

    mean <- unlist(c(mean, fp_data[which(fp_data$bin==bin),"IVW_meta_beta"]))
    upper <- unlist(c(upper, fp_data[which(fp_data$bin==bin),"IVW_meta_beta"] + fp_data[which(fp_data$bin==bin),"IVW_meta_se"]))
    lower <- unlist(c(lower, fp_data[which(fp_data$bin==bin),"IVW_meta_beta"] - fp_data[which(fp_data$bin==bin),"IVW_meta_se"]))

  }

  MR_text <- c("MR estimate (95% CI)")
  for(bin in bins){

    means <- unlist(fp_data[which(fp_data$bin==bin),"IVW_meta_beta"])
    ses <- unlist(fp_data[which(fp_data$bin==bin),"IVW_meta_se"])

    MR_text <- c(MR_text)
    for(i in 1:length(means)){
      MR_text <- c(MR_text, make_beta_95ci(means[i], ses[i]))

    }

    MR_text <- c(MR_text)

  }

  upper <- unlist(c(upper, fp_data[which(fp_data$bin==bin),"IVW_meta_beta"] + fp_data[which(fp_data$bin==bin),"IVW_meta_se"], NA))
  lower <- unlist(c(lower, fp_data[which(fp_data$bin==bin),"IVW_meta_beta"] - fp_data[which(fp_data$bin==bin),"IVW_meta_se"], NA))


  summary_rows <- rep(FALSE, length(bin_labels))
  index <- 1
  summary_rows[index] <- TRUE



  fp_plot <- forestplot(
    labeltext = cbind(bin_labels, MR_text), fn.ci_norm = fpDrawNormalCI,
    mean = mean,
    upper = upper,
    lower = lower,
    is.summary = summary_rows,
    boxsize = .1, # We set the box size to better visualize the type

    col = fpColors(box = custom_col[c(2,4)],  lines=c("gray50", "gray50")), ci.vertices=TRUE,

    #col=fpColors(box="black", lines="black", zero = "gray50"),
    xlab = paste0("\n",xlab),


    txt_gp = fpTxtGp(xlab = gpar(cex = 0.9),
                     ticks = gpar(cex = 0.8),
                     label = gpar(cex = 0.9)),

    legend_args = fpLegend(pos = "right"), xticks = seq(xmin_lab, xmax_lab, by=0.1))

  return(fp_plot)



}

xy_plot_binned_sex_specifc <- function(harmonised_data, MR_binned, group, custom_col){

  harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])

  fig_data <- MR_binned %>% filter(bin!="all") %>% separate(bin, c("bin_start_temp", "bin_stop_temp"), ",", remove = F) %>% filter(grouping_var==!!group) %>% mutate(bin_start = substring(bin_start_temp, 2)) %>%
    mutate(bin_stop = str_sub(bin_stop_temp,1,nchar(bin_stop_temp)-1)) %>% rowwise() %>%
    mutate(bin_median = median(c(as.numeric(bin_start), as.numeric(bin_stop)))) %>%
    mutate(bin_median_plot = case_when(exposure_sex=="female" ~ bin_median + 0.4,
                                       exposure_sex=="male" ~ bin_median - 0.4))

  x_ticks <- unique(fig_data %>% pull(bin_median))
  x_labels <- unique(fig_data %>% pull(bin))

  xlab <- ifelse(group == "time_together_even_bins", "Time together in same household (years)", "Median age of couples (years)")

  legend_title <- "Exposure sex"

  plot <- ggplot2::ggplot(data = fig_data, ggplot2::aes(x = bin_median_plot,
                                                        y = IVW_beta, colour = factor(exposure_sex))) +
    geom_smooth(method="lm",formula=y~x, se = F, fullrange = T, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IVW_beta - IVW_se, ymax = IVW_beta + IVW_se),
                           colour = "grey", width = 0) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(exposure_sex))) +

    ggplot2::scale_colour_manual(values = custom_col[c(2,4)]) +

    scale_x_continuous(breaks=x_ticks, limits = c(min(x_ticks)-2, max(x_ticks)+2),
                       labels = x_labels) +

    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line=element_blank(),
                       axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
                       axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"))) +

    ggplot2::labs(colour = legend_title, x =xlab,
                  y = paste0("AM MR estimate for\n", tolower(trait_description) ,"\n(sex-specific)")) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))

  return(plot)

}

xy_plot_binned_meta <- function(harmonised_data, MR_binned, group, custom_col){

  #harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  #harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  #harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonised_data[1, "exposure_description"])

  fig_data <- MR_binned %>% filter(bin!="all") %>% separate(bin, c("bin_start_temp", "bin_stop_temp"), ",", remove = F) %>% filter(grouping_var==!!group) %>% mutate(bin_start = substring(bin_start_temp, 2)) %>%
    mutate(bin_stop = str_sub(bin_stop_temp,1,nchar(bin_stop_temp)-1)) %>% rowwise() %>%
    mutate(bin_median = median(c(as.numeric(bin_start), as.numeric(bin_stop))))


  # %>%
  #   mutate(bin_median_plot = case_when(exposure_sex=="female" ~ bin_median + 0.4,
  #                                      exposure_sex=="male" ~ bin_median - 0.4))

  x_ticks <- unique(fig_data %>% pull(bin_median))
  x_labels <- unique(fig_data %>% pull(bin))

  xlab <- ifelse(group == "time_together_even_bins", "Time together in same household (years)", "Median age of couples (years)")

  legend_title <- "Exposure sex"

  plot <- ggplot2::ggplot(data = fig_data, ggplot2::aes(x = bin_median,
                                                              y = IVW_beta)) +
    geom_smooth(method="lm",formula=y~x, se = F, fullrange = T, color = custom_col[2], linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IVW_beta - IVW_se, ymax = IVW_beta + IVW_se),
                           colour = "grey", width = 0) +
    ggplot2::geom_point(color = custom_col[2]) +
    scale_x_continuous(breaks=x_ticks, limits = c(min(x_ticks)-2, max(x_ticks)+2),
                       labels = x_labels) +

    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line=element_blank(),
                       axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
                       axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"))) +

    ggplot2::labs(colour = legend_title, x =xlab,
                  y = paste0("AM MR estimate for\n", tolower(trait_description) ,"\n(meta-analyzed across sexes)"))

  return(plot)

}

xy_plot_binned_single_sex <- function(harmonised_data, MR_binned, exposure_sex, group, custom_col){

  harmonised_data.sex_male <- harmonised_data[[paste0("exp_", "male", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)
  harmonised_data.sex_female <- harmonised_data[[paste0("exp_", "female", "_harmonised_data_filter")]] %>% dplyr::filter(grouping_var == "age_even_bins") %>% dplyr::filter(bin== "all") %>% mutate_at('bin', as.factor)

  fig_data <- MR_binned %>% filter(bin!="all") %>% separate(bin, c("bin_start_temp", "bin_stop_temp"), ",", remove = F) %>% filter(grouping_var==!!group) %>% mutate(bin_start = substring(bin_start_temp, 2)) %>%
    mutate(bin_stop = str_sub(bin_stop_temp,1,nchar(bin_stop_temp)-1)) %>% rowwise() %>%
    mutate(bin_median = median(c(as.numeric(bin_start), as.numeric(bin_stop)))) %>%
    mutate(bin_median_plot = case_when(exposure_sex=="female" ~ bin_median + 0.4,
                                       exposure_sex=="male" ~ bin_median - 0.4))

  harmonise_dat_plot <- rbind(harmonised_data.sex_male, harmonised_data.sex_female) %>% mutate_at('exposure_sex', as.factor)
  trait_description <- as.character(harmonise_dat_plot[1, "exposure_description"])

  x_ticks <- unique(fig_data %>% pull(bin_median))
  x_labels <- unique(fig_data %>% pull(bin))

  xlab <- ifelse(group == "time_together_even_bins", "Time together in\nsame household (years)", "Median age of\ncouples (years)")

  legend_title <- "Exposure sex"

  if(exposure_sex=="male"){outcome_sex="female"}
  if(exposure_sex=="female"){outcome_sex="male"}

  fig_data_sex_spec <- fig_data_sex_spec <- fig_data %>% dplyr::filter(exposure_sex==!!exposure_sex)

  plot <- ggplot2::ggplot(data = fig_data_sex_spec, ggplot2::aes(x = bin_median,
                                                                 y = IVW_beta, colour = factor(exposure_sex))) +
    geom_smooth(method="lm",formula=y~x, se = F, fullrange = T, linetype = "dashed", color = custom_col[2]) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IVW_beta - IVW_se, ymax = IVW_beta + IVW_se),
                           colour = "grey", width = 0) +
    ggplot2::geom_point(color = custom_col[2]) +


    scale_x_continuous(breaks=x_ticks, limits = c(min(x_ticks)-2, max(x_ticks)+2),
                       labels = x_labels) +

    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line=element_blank(),
                       axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
                       axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"))) +

    ggplot2::labs(colour = legend_title, x =xlab,
                  y = paste0("AM MR estimate for\n", tolower(trait_description) ,"\n(", exposure_sex, "s to ", outcome_sex, "s)")) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))


  return(plot)


}

create_MR_binned_AM_figs <- function(household_harmonised_data_meta_reverse_filter, household_harmonised_data_reverse_filter,
                                     household_MR_binned_SNPmeta, household_MR_binned_MRmeta, custom_col){

  names_data <- names(household_harmonised_data_meta_reverse_filter)

  for(j in 1:length(names_data)){
    outcome_ID <- unlist(str_split(names_data[j], "_vs_"))[1]
    exposure_ID <- unlist(str_split(unlist(str_split(names_data[j], "_vs_"))[2], "_harmonised_data_meta_filter"))[1]
    if(outcome_ID==exposure_ID){
      index_same_trait <- j
    }
  }


  household_harmonised_data_meta.AM  <- household_harmonised_data_meta_reverse_filter[[index_same_trait]]
  household_harmonised_data.AM  <- household_harmonised_data_reverse_filter[[index_same_trait]]

  household_MR_binned_SNPmeta.AM <- household_MR_binned_SNPmeta[[index_same_trait]]
  household_MR_binned_MRmeta.AM <- household_MR_binned_MRmeta[[index_same_trait]]

  for(group in c("age_even_bins", "time_together_even_bins")){

    p_male <- mr_plot_binned_sex_specific(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, "male", group)
    p_female <- mr_plot_binned_sex_specific(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, "female", group)

    prow <- plot_grid(
      p_male + theme(legend.position="none"),
      p_female + theme(legend.position="none"),
      hjust = -1,
      nrow = 1
    )

    legend_b <- get_legend(
      p_male +
        guides(color = guide_legend(nrow = 1, title.position = "left")) +

        theme(legend.position = "bottom", legend.title=element_text(size=10))
    )

    # add the legend underneath the row we made earlier. Give it 10%
    # of the height of one plot (via rel_heights).


    p_combined <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
    assign(paste0("MR_sex_specific_", group, "_fig"), p_combined)
  }

  for(group in c("age_even_bins", "time_together_even_bins")){

    #fp_sex_spec <- forest_plot_binned_sex_specific(household_harmonised_data.AM, household_MR_binned_meta.AM, group, custom_col)
    #fp_meta <- forest_plot_binned_meta(household_harmonised_data.AM, household_MR_binned_meta.AM, group, custom_col)

    xy_sex_spec <- xy_plot_binned_sex_specifc(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, group, custom_col)
    xy_meta <- xy_plot_binned_meta(household_harmonised_data_meta.AM, household_MR_binned_SNPmeta.AM, group, custom_col)


    #assign(paste0("FP_sex_specific_", group, "_fig"), fp_sex_spec)
    #assign(paste0("FP_", group, "_fig"), fp_meta)
    assign(paste0("XY_sex_specific_", group, "_fig"), xy_sex_spec)

    assign(paste0("XY_", group, "_fig"), xy_meta)

  }

  for(group in c("age_even_bins", "time_together_even_bins")){

    xy_binM <- xy_plot_binned_single_sex(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, "male", group, custom_col)
    xy_binF <- xy_plot_binned_single_sex(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, "female", group, custom_col)

    xy_bin_side_by_side <- plot_grid(xy_binM, xy_binF)
    assign(paste0("XY_sex_sidebyside_", group, "_fig"), xy_bin_side_by_side)
  }

  overall_MR_fig <- mr_plot_sex_specific(household_harmonised_data.AM, household_MR_binned_MRmeta.AM, household_MR_binned_SNPmeta.AM, custom_col)

  return(list(overall_MR_fig = overall_MR_fig,
              MR_sex_specific_age_bins_fig = MR_sex_specific_age_even_bins_fig,
              MR_sex_specific_time_together_bins_fig = MR_sex_specific_time_together_even_bins_fig,

              XY_sex_specific_age_bins_fig = XY_sex_specific_age_even_bins_fig,
              XY_sex_specific_time_together_bins_fig = XY_sex_specific_time_together_even_bins_fig,
              XY_sex_sidebyside_age_bins_fig = XY_sex_sidebyside_age_even_bins_fig,
              XY_sex_sidebyside_time_together_bins_fig = XY_sex_sidebyside_time_together_even_bins_fig,

              XY_age_bins_fig = XY_age_even_bins_fig,
              XY_time_together_bins_fig = XY_time_together_even_bins_fig

              #FP_sex_specific_age_bins_fig = FP_sex_specific_age_even_bins_fig,
              #FP_sex_specific_time_together_bins_fig = FP_sex_specific_time_together_even_bins_fig,
              #FP_age_bins_fig = FP_age_even_bins_fig,
              #FP_time_together_bins_fig = FP_time_together_even_bins_fig
              ))

}

create_household_MR_AM_FvsM_fig <- function(household_MR_binned_het, custom_col){

  AM_result <- household_MR_binned_het %>% filter(same_trait)

  plot_data <- AM_result %>% dplyr::select(exposure_ID, exposure_sex, IVW_beta, IVW_se, IVW_pval, sex_het_pval) %>% unique() %>%
    pivot_wider(names_from = "exposure_sex", values_from = c("IVW_beta", "IVW_se", "IVW_pval")) %>%
    mutate(sex_het_sig = case_when(sex_het_pval < 0.05 ~ TRUE, TRUE ~ FALSE))


  legend_title <- "Sex heterogeneity"
  plot <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = IVW_beta_male,
                                                         y = IVW_beta_female, color = factor(sex_het_sig), label = exposure_ID)) +
    geom_point(alpha = 3/4) +

    geom_smooth(mapping = aes(x = IVW_beta_male, y = IVW_beta_female), method = "lm", se=FALSE, formula = y~x+0, fullrange=TRUE, color = custom_col[4]) +
    geom_abline(slope=1, intercept=0, color = "black") +
    scale_color_manual(values = custom_col[c(1,2)], labels=c("p >= 0.01","p < 0.01")) +
    theme_minimal() +
    ggplot2::labs(colour = legend_title, x = paste0("AM MR estimate (male to female)"),
                  y = paste("AM MR estimate (female to male)"))

  plot2 <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x = IVW_beta_female,
                                                         y = IVW_beta_male, color = factor(sex_het_sig), label = exposure_ID)) +
    geom_point(alpha = 3/4) +

    geom_smooth(mapping = aes(x = IVW_beta_female, y = IVW_beta_male), method = "lm", se=FALSE, formula = y~x+0, fullrange=TRUE, color = custom_col[4]) +
    geom_abline(slope=1, intercept=0, color = "black") +
    scale_color_manual(values = custom_col[c(1,2)], labels=c("p >= 0.01","p < 0.01")) +
    theme_minimal() +
    ggplot2::labs(colour = legend_title, x = paste0("AM MR estimate (female to male)"),
                  y = paste("AM MR estimate (male to female)"))

  return(list(plot_male_to_female = plot, plot_female_to_male = plot2))

}

