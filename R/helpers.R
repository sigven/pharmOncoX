#' Utility function to get on-label/off-label drugs (not complete)
#'
#' @param cache_dir Local cache directory
#' @keywords internal
#'
#'

get_targeted_drugs <- function(cache_dir = NA) {


  onco_drugs <- 
    pharmOncoX::get_drugs(
      treatment_category = c("targeted_therapy_classified",
                             "targeted_therapy_unclassified",
                             "hormone_therapy_classified"),
              cache_dir = cache_dir,
              drug_classified_cancer = TRUE,
              drug_cancer_indication = TRUE,
              drug_action_inhibition = TRUE,
              drug_source_opentargets = TRUE)
  
  metadata <- onco_drugs$metadata
    
    
  targeted_onco_inhibitors <- as.data.frame(
    onco_drugs$records |>
      dplyr::filter(!is.na(.data$drug_max_ct_phase)) |>
      dplyr::filter(.data$drug_max_ct_phase >= 1) |>
      dplyr::filter(!is.na(.data$molecule_chembl_id)) |>
      dplyr::mutate(
        drug_link =
          paste0("<a href='https://www.targetvalidation.org/summary?drug=",
                 .data$molecule_chembl_id, "' target='_blank'>",
                 .data$drug_name, "</a>")) |>
      dplyr::select(c("target_symbol",
                    "target_genename",
                    "molecule_chembl_id",
                    "drug_max_phase_indication",
                    "drug_approved_indication",
                    "drug_cancer_relevance",
                    "drug_frac_cancer_indications",
                    "drug_year_first_approval",
                    "drug_clinical_id",
                    "disease_efo_id",
                    "disease_efo_label",
                    "primary_site",
                    "drug_type",
                    "drug_action_type",
                    "nci_concept_definition",
                    "drug_name",
                    "drug_link",
                    "atc_level3",
                    "atc_treatment_category")) |>
      dplyr::mutate(drug_max_phase_indication = dplyr::if_else(
        is.na(.data$drug_max_phase_indication) |
          .data$drug_max_phase_indication == "",
        as.integer(0),
        as.integer(.data$drug_max_phase_indication)
      )) |>
      dplyr::rename(symbol = target_symbol,
                    genename = target_genename) |>
      tidyr::separate_rows(
        drug_clinical_id, sep = ",") |>
      dplyr::group_by(.data$symbol,
                      .data$genename,
                      .data$molecule_chembl_id,
                      .data$drug_name,
                      .data$drug_action_type,
                      .data$drug_type,
                      .data$primary_site,
                      .data$drug_link,
                      .data$atc_level3,
                      .data$atc_treatment_category,
                      .data$drug_cancer_relevance,
                      .data$drug_frac_cancer_indications,
                      .data$drug_year_first_approval,
                      .data$nci_concept_definition) |>
      dplyr::summarise(
        drug_clinical_id =
          paste(sort(unique(.data$drug_clinical_id)),
                collapse = "|"),
        max_phase = max(.data$drug_max_phase_indication, na.rm = T),
        max_all_phase = paste(sort(unique(.data$drug_max_phase_indication)),
                              collapse = "|"),
        approved_indication = 
          paste(sort(unique(.data$drug_approved_indication)),
                collapse = "|"),
        disease_efo_label = paste(sort(unique(.data$disease_efo_label)),
                                  collapse = "|"),
        disease_efo_id = paste(sort(unique(.data$disease_efo_id)),
                               collapse = "|"),
        .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(
        num_clinical_id =
          stringr::str_count(.data$drug_clinical_id, pattern = "|")) |>
      dplyr::rename(drug_indication_label =
                      disease_efo_label,
                    drug_indication_id =
                      disease_efo_id) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        !is.na(.data$primary_site) &
          .data$primary_site == "Other/Unknown"
        & stringr::str_detect(
          .data$drug_indication_label,"carcinosarcoma"),
        "Soft Tissue",
        as.character(.data$primary_site))) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        !is.na(.data$primary_site) &
          .data$primary_site == "Other/Unknown"
        & stringr::str_detect(
          .data$drug_indication_label,"breast"),
        "Breast",
        as.character(.data$primary_site))) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        is.na(.data$primary_site),
        "Any",
        as.character(.data$primary_site))) |>
      dplyr::filter(
        !stringr::str_detect(
          tolower(.data$drug_name),"^yttrium")) |>
      dplyr::mutate(approved_indication = dplyr::if_else(
        !is.na(approved_indication) &
          stringr::str_detect(
            approved_indication,"\\|"),
        as.logical(TRUE),
        as.logical(approved_indication)
      )) |>     
      dplyr::select(-"genename") |>
      dplyr::distinct() |>
      dplyr::select(c("primary_site",
                    "molecule_chembl_id",
                    "symbol",
                    "drug_name",
                    "drug_type",
                    "drug_action_type",
                    "drug_link",
                    "atc_level3",
                    "atc_treatment_category",
                    "max_phase",
                    "max_all_phase",
                    "drug_frac_cancer_indications",
                    "drug_cancer_relevance",
                    "drug_year_first_approval",
                    "approved_indication",
                    "drug_indication_label")) |>
      dplyr::rename(drug_primary_site = primary_site,
                    drug_max_phase_indication = max_phase,
                    drug_max_phase_indication_all = max_all_phase)
  )

  targeted_drugs_per_site <- list()
  targeted_drugs_non_site_specific <- targeted_onco_inhibitors |>
    dplyr::filter(.data$drug_primary_site == "Any") |>
    dplyr::distinct()


  unique_sites <- unique(targeted_onco_inhibitors$drug_primary_site)

  all_tt_records <- data.frame()

  for (t in unique(targeted_onco_inhibitors$drug_primary_site)) {
    #if (!is.na(t)) {
    targeted_drugs_per_site[[t]] <- list()
    targeted_drugs_per_site[[t]][['off_label']] <- data.frame()
    targeted_drugs_per_site[[t]][['other_any_phase']] <- data.frame()
    targeted_drugs_per_site[[t]][['on_label']] <- list()
    targeted_drugs_per_site[[t]][['on_label']][['early_phase']] <- data.frame()
    targeted_drugs_per_site[[t]][['on_label']][['late_phase']] <- data.frame()
    
    targeted_drugs_per_site[[t]][['on_label']][['late_phase']] <- 
      targeted_onco_inhibitors |>
      dplyr::filter(.data$drug_primary_site == t &
                      .data$drug_max_phase_indication > 2) |>
      dplyr::distinct() |>
      dplyr::arrange(dplyr::desc(.data$drug_max_phase_indication),
                     dplyr::desc(.data$drug_year_first_approval),
                     .data$atc_treatment_category,
                     dplyr::desc(.data$drug_frac_cancer_indications),
                     .data$symbol) |>
      dplyr::filter(!(.data$drug_cancer_relevance == "by_cancer_target_nci" &
                        .data$atc_treatment_category == "cancer_unclassified"))
    
    targeted_drugs_per_site[[t]][['on_label']][['early_phase']] <- 
      targeted_onco_inhibitors |>
      dplyr::filter(.data$drug_primary_site == t &
                      .data$drug_max_phase_indication <= 2) |>
      dplyr::distinct() |>
      dplyr::arrange(dplyr::desc(.data$drug_max_phase_indication),
                     dplyr::desc(.data$drug_year_first_approval),
                     .data$atc_treatment_category,
                     dplyr::desc(.data$drug_frac_cancer_indications),
                     .data$symbol) |>
      dplyr::filter(!(.data$drug_cancer_relevance == "by_cancer_target_nci" &
                        .data$atc_treatment_category == "cancer_unclassified"))
    
    targeted_drugs_per_site[[t]][['off_label']]  <- 
      targeted_onco_inhibitors |>
      dplyr::filter(
        .data$drug_primary_site != t &
          .data$drug_primary_site != "Any" &
          .data$drug_max_phase_indication > 2) |>
      dplyr::anti_join(
        targeted_drugs_per_site[[t]][['on_label']][['late_phase']],
        by = c("molecule_chembl_id","drug_name","symbol")) |>
      dplyr::anti_join(
        targeted_drugs_per_site[[t]][['on_label']][['early_phase']],
        by = c("molecule_chembl_id","drug_name","symbol")) |>
      dplyr::distinct() |>
      dplyr::arrange(dplyr::desc(.data$drug_max_phase_indication),
                     dplyr::desc(.data$drug_year_first_approval),
                     .data$atc_treatment_category,
                     dplyr::desc(.data$drug_frac_cancer_indications),
                     .data$symbol) |>
      dplyr::filter(!(.data$drug_cancer_relevance == "by_cancer_target_nci" &
                        .data$atc_treatment_category == "cancer_unclassified"))
    
    other_any_phase <- 
      targeted_drugs_non_site_specific |>
      dplyr::anti_join(
        targeted_drugs_per_site[[t]][['off_label']],
        by = c("molecule_chembl_id","drug_name","symbol")) |>
      dplyr::anti_join(
        targeted_drugs_per_site[[t]][['on_label']][['late_phase']],
        by = c("molecule_chembl_id","drug_name","symbol")) |>
      dplyr::anti_join(
        targeted_drugs_per_site[[t]][['on_label']][['early_phase']],
        by = c("molecule_chembl_id","drug_name","symbol"))
    
    targeted_drugs_per_site[[t]][['other_any_phase']] <- other_any_phase |>
      dplyr::arrange(dplyr::desc(.data$drug_max_phase_indication),
                     dplyr::desc(.data$drug_year_first_approval),
                     .data$atc_treatment_category,
                     dplyr::desc(.data$drug_frac_cancer_indications),
                     .data$symbol) |>
      dplyr::filter(!(.data$drug_cancer_relevance == "by_cancer_target_nci" &
                        .data$atc_treatment_category == "cancer_unclassified"))
    
    
    lgr::lgr$info(
      paste0(
        "Found n = ",
        length(unique(
          targeted_drugs_per_site[[t]][['on_label']][['late_phase']]$drug_name)),
        " on-label targeted drugs for primary site: ",
        t))
    
    tmp <- list()
    ## On label - late phase
    tmp[['on_label_lp']] <- 
      targeted_drugs_per_site[[t]][['on_label']][['late_phase']] |>
      dplyr::mutate(
        off_label = FALSE,
        drug_label = "ON_LABEL",
        query_site = t,
        drug_clinical_phase = "late") |>
      dplyr::select("query_site", dplyr::everything())
    
    ## On label early phase
    tmp[['on_label_ep']] <- 
      targeted_drugs_per_site[[t]][['on_label']][['early_phase']] |>
      dplyr::mutate(
        off_label = FALSE,
        drug_label = "ON_LABEL",
        query_site = t,
        drug_clinical_phase = "early") |>
      dplyr::select("query_site", dplyr::everything())
    
    ## Off label (late phase)
    tmp[['off_label']] <- 
      targeted_drugs_per_site[[t]][['off_label']] |>
      dplyr::mutate(
        off_label = TRUE,
        drug_label = "OFF_LABEL",
        query_site = t,
        drug_clinical_phase = "late") |>
      dplyr::select("query_site", dplyr::everything())
    
    ## Other drugs
    tmp[['other']] <- 
      other_any_phase |>
      dplyr::mutate(
        drug_label = "OTHER",
        query_site = t,
        drug_clinical_phase = dplyr::if_else(
          .data$drug_max_phase_indication <= 2,
          "early",
          "late"
        )) |>
      dplyr::select("query_site", dplyr::everything())
    
    all_tt_records <- all_tt_records |>
      dplyr::bind_rows(
        tmp[['on_label_lp']],
        tmp[['on_label_ep']],
        tmp[['off_label']],
        tmp[['other']]
      )

  }
  
  all_tt_records <- all_tt_records |>
    dplyr::filter(!(.data$atc_treatment_category == "cancer_unclassified" & 
                    .data$drug_cancer_relevance == "by_cancer_condition_otp" & 
                    .data$drug_frac_cancer_indications < 0.6 & 
                    .data$approved_indication == F)) |>
    dplyr::filter(!stringr::str_detect(.data$symbol,"^RRM"))
  
  return(list('records' = all_tt_records, 
              'metadata' = metadata))

}



#' Function that retrieves pharmOncoX data from Google Drive
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should force downloaded
#' (i.e. set to TRUE to re-download even if data exists in cache)
#'
#' @keywords internal
#'
#'
get_drug_records <- function(cache_dir = NA,
                             force_download = F
                             ) {


  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  if (is.na(cache_dir)) {
    lgr::lgr$fatal(paste0("Argument cache_dir = '",
                          cache_dir, "' is not defined"))
    stop()
  }

  if (!dir.exists(cache_dir)) {
    lgr::lgr$fatal(paste0("Argument cache_dir = '",
                          cache_dir, "' does not exist"))
    stop()
  }


  drug_datasets <- list()
  file_maps <- c('drug_map_name',
                 'drug_map_alias',
                 'drug_map_basic',
                 'drug_map_target',
                 'drug_map_indication')

  for (elem in file_maps) {

    fname_local <- file.path(
      cache_dir,
      paste0(elem,"_v",
             db_id_ref[db_id_ref$name == elem,]$pVersion,
             '.rds')
    )

    fname_gd <- googledrive::as_id(
      db_id_ref[db_id_ref$name == elem,]$gid)

    md5checksum_package <-
      db_id_ref[db_id_ref$name == elem,]$md5Checksum

    if (length(fname_local) > 1) {
      fname_local <- fname_local[-1]
      lgr::lgr$warn(
        paste0(
          "Multiple data repos entries for ", elem,
          " - using latest entry")
        
      )
    }
    
    #dat <- NULL
    if (file.exists(fname_local) && force_download == F) {
      drug_datasets[[elem]] <- readRDS(fname_local)
      drug_datasets[[elem]][['fpath']] <- fname_local
      if (!is.null(drug_datasets[[elem]][['records']])) {
        lgr::lgr$info(paste0(
          "Reading from cache_dir = '", 
          cache_dir, "', argument force_download = F"))
        lgr::lgr$info(
          paste0("Object '", 
                 elem,"' sucessfully loaded"))
        lgr::lgr$info(paste0(
          "Retrieved n = ", 
          nrow(drug_datasets[[elem]][['records']]), " records"))
      }
    }else{
      
      googledrive::drive_deauth()

      lgr::lgr$info(
        "Downloading remote dataset from Google Drive to cache_dir")
      dl <- googledrive::with_drive_quiet(
        googledrive::drive_download(
          fname_gd,
          path = fname_local,
          overwrite = TRUE)
      )

      md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
      md5checksum_local <- tools::md5sum(fname_local)
      names(md5checksum_local) <- NULL

      if (md5checksum_remote == md5checksum_local) {
        drug_datasets[[elem]] <- readRDS(fname_local)
        drug_datasets[[elem]]$fpath <- fname_local
        if (!is.null(drug_datasets[[elem]][['records']]) &
           !is.null(drug_datasets[[elem]][['metadata']])) {

          lgr::lgr$info(paste0(
            "Reading from cache_dir = ' (", 
            cache_dir, "'), argument force_download = F"))
          lgr::lgr$info(paste0("Object '", elem, "' sucessfully loaded"))
          lgr::lgr$info(paste0("md5 checksum is valid: ", md5checksum_remote))
          lgr::lgr$info(paste0(
            "Retrieved ", nrow(drug_datasets[[elem]][['records']]), " records"))

        }
      }else{
        lgr::lgr$error(
          paste0("md5 checksum of local file (", md5checksum_local,
                 ") is inconsistent with remote file (",
                 md5checksum_remote,")"))
        stop()
      }
      
    }
  }

  drug_data <- list()

  drug_aliases_collapsed <- as.data.frame(
    drug_datasets[['drug_map_alias']][['records']] |>
      dplyr::group_by(.data$drug_id) |>
      dplyr::summarise(
        drug_alias = paste(sort(unique(.data$alias)), collapse = "|"),
        .groups = "drop")
  )


  drug_data[['records']] <- drug_datasets[['drug_map_name']][['records']] |>
    dplyr::left_join(drug_datasets[['drug_map_basic']][['records']],
                     by = "drug_id")  |>
    dplyr::left_join(drug_aliases_collapsed, by = "drug_id") |>
    dplyr::left_join(drug_datasets[['drug_map_target']][['records']],
                     by = "drug_id", multiple = "all",
                     relationship = "many-to-many") |>
    dplyr::left_join(drug_datasets[['drug_map_indication']][['records']],
                     by = "drug_id", multiple = "all",
                     relationship = "many-to-many")

  drug_data[['metadata']] <- drug_datasets[['drug_map_name']][['metadata']]

  return(drug_data)

}



