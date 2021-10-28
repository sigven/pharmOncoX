#' Retrieve antineoplastic drugs and drug regimens
#'
#' @param drug_is_targeted logical indicating if resulting drug records should contain molecularly targeted drugs only
#' @param drug_is_approved logical indicating if resulting drug records should contain approved drugs only
#' @param drug_target character vector with drug targets (gene symbols) for drug records included in results
#' @param drug_action_type character vector with drug action types to include in drug record list -
#' possible values "INHIBITOR","AGONIST","MODULATOR","ANTAGONIST","BLOCKER","ACTIVATOR","BINDING AGENT","OPENER",
#' "STABILISER","CROSS-LINKING AGENT",DISRUPTING AGENT","OTHER"
#' @param drug_indication_main character vector with main tumor types for which drug(s) are indicated
#' possible values "Adrenal Gland","Biliary Tract","Bladder/Urinary Tract","Bone","Breast","Cervix","CNS/Brain",
#' "Colon/Rectum","Esophagus/Stomach","Eye","Head and Neck","Kidney","Liver","Lung","Lymphoid","Myeloid","Ovary/Fallopian Tube",
#' "Pancreas","Penis","Peripheral Nervous System","Peritoneum","Pleura","Prostate","Skin","Soft Tissue","Testis","Thymus",
#' "Thyroid","Uterus","Vulva/Vagina"
#' @param source_opentargets_only logical indicating if resulting drug records should contain drug records from Open Targets Platform only
#' @param drug_has_blackbox_warning logical indicating if resulting drug records should contain drugs with black box warnings only
#' @param drug_approved_later_than only include records for drugs approved later than this date (year)
#' @param drug_minimum_phase_any_indication only include drug records that are in a clinical phase (any indication) greater or equal than this phase
#' @param output_style comprehensiveness of drug records ('extensive' - all annotations, 'narrow' - main annotations)
#' @param output_resolution output record resolution ('drug','drug2target','drug2target2indication')
#' @param list_per_drug_synonym list resulting drug records per drug synonym
#' @param is_alkylating_agent logical indicating if only this drug class is wanted
#' @param is_angiogenesis_inhibitor logical indicating if only this drug class is wanted
#' @param is_anthracycline logical indicating if only this drug class is wanted
#' @param is_antimetabolite logical indicating if only this drug class is wanted
#' @param is_ar_antagonist logical indicating if only this drug class is wanted
#' @param is_bet_inhibitor logical indicating if only this drug class is wanted
#' @param is_hedgehog_antagonist logical indicating if only this drug class is wanted
#' @param is_hdac_inhibitor logical indicating if only this drug class is wanted
#' @param is_hormone_therapy logical indicating if only this drug class is wanted
#' @param is_immune_checkpoint_inhibitor logical indicating if only this drug class is wanted
#' @param is_kinase_inhibitor logical indicating if only this drug class is wanted
#' @param is_monoclonal_antibody logical indicating if only this drug class is wanted
#' @param is_parp_inhibitor logical indicating if only this drug class is wanted
#' @param is_platinum_compound logical indicating if only this drug class is wanted
#' @param is_proteasome_inhibitor logical indicating if only this drug class is wanted
#' @param is_topoisomerase_inhibitor logical indicating if only this drug class is wanted
#' @param is_tubulin_inhibitor logical indicating if only this drug class is wanted
#'
#' @export
get_onco_drugs <- function(drug_is_targeted = F,
                     drug_is_approved = F,
                     drug_target = NULL,
                     drug_action_type = NULL,
                     drug_indication_main = NULL,
                     source_opentargets_only = F,
                     drug_has_blackbox_warning = F,
                     drug_approved_later_than = 1939,
                     drug_minimum_phase_any_indication = 0,
                     list_per_drug_synonym = F,
                     output_resolution = "drug2target2indication",
                     output_style = "extensive",
                     is_alkylating_agent = F,
                     is_angiogenesis_inhibitor = F,
                     is_anthracycline = F,
                     is_antimetabolite = F,
                     is_ar_antagonist = F,
                     is_bet_inhibitor = F,
                     is_hedgehog_antagonist = F,
                     is_hdac_inhibitor = F,
                     is_hormone_therapy = F,
                     is_immune_checkpoint_inhibitor = F,
                     is_kinase_inhibitor = F,
                     is_monoclonal_antibody = F,
                     is_parp_inhibitor = F,
                     is_platinum_compound = F,
                     is_proteasome_inhibitor = F,
                     is_topoisomerase_inhibitor = F,
                     is_tubulin_inhibitor = F
                     ){

  drug_records <- oncoPharmaDB::oncopharmadb

  valid_output_styles <- c("extensive","narrow")
  valid_output_resolutions <- c("drug","drug2target","drug2target2indication")
  valid_drug_action_types <- c("INHIBITOR","AGONIST","MODULATOR","ANTAGONIST",
                               "BLOCKER","ACTIVATOR","BINDING AGENT","OPENER",
                               "STABILISER","CROSS-LINKING AGENT",
                               "DISRUPTING AGENT","OTHER")
  valid_indications_main <-
    c("Adrenal Gland","Biliary Tract","Bladder/Urinary Tract",
      "Bone","Breast","Cervix","CNS/Brain",
      "Colon/Rectum","Esophagus/Stomach","Eye","Head and Neck",
      "Kidney","Liver","Lung","Lymphoid",
      "Myeloid","Ovary/Fallopian Tube",
      "Pancreas","Penis","Peripheral Nervous System","Peritoneum","Pleura",
      "Prostate","Skin","Soft Tissue","Testis","Thymus",
      "Thyroid","Uterus","Vulva/Vagina")

  # valid_tumor_types <- sort(
  #   unique(drug_records$primary_site[!is.na(drug_records$primary_site)]))

  arg_validation_messages <- list()

  arg_validation_messages[[1]] <-
     assertthat::validate_that(
       is.numeric(drug_approved_later_than),
       msg = "ERROR: Argument 'drug_approved_later_than' must be of type 'numeric'")
  arg_validation_messages[[2]] <-
    assertthat::validate_that(
      drug_approved_later_than >= 1939 & drug_approved_later_than <= 2021,
      msg = "ERROR: Argument 'drug_approved_later_than' must be larger than 1939 and less than 2021")
  arg_validation_messages[[3]] <-
    assertthat::validate_that(
      is.logical(drug_is_approved),
      msg = "ERROR: Argument 'drug_is_approved' must be of type 'logical'")
  arg_validation_messages[[4]] <-
    assertthat::validate_that(
      is.logical(drug_is_targeted),
      msg = "ERROR: Argument 'drug_is_targeted' must be of type 'logical'")
  arg_validation_messages[[5]] <-
    assertthat::validate_that(
      is.logical(drug_has_blackbox_warning),
      msg = "ERROR: Argument 'drug_has_blackbox_warning' must be of type 'logical'")
  arg_validation_messages[[6]] <-
    assertthat::validate_that(
      is.logical(list_per_drug_synonym),
      msg = "ERROR: Argument 'list_per_drug_synonym' must be of type 'logical'")
  arg_validation_messages[[7]] <-
    assertthat::validate_that(
      output_resolution %in% valid_output_resolutions,
      msg = "ERROR: Argument 'output_resolution' must be either 'drug','drug2target', or 'drug2target2indication'")
  arg_validation_messages[[8]] <-
    assertthat::validate_that(
      output_style %in% valid_output_styles,
      msg = "ERROR: Argument 'output_style' must be either 'extensive' or 'narrow'")

  arg_counter <- 9
  if(!is.null(drug_action_type)){
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        drug_is_targeted == T,
        msg = "ERROR: Argument 'drug_is_targeted' must be set to TRUE when 'drug_action_type' is non-NULL")

    arg_counter <- arg_counter + 1
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_action_type),
        msg = "ERROR: Argument 'drug_action_type' must be a vector of type 'character'")

    arg_counter <- arg_counter + 1

    if(length(unique(drug_action_type %in% valid_drug_action_types)) > 1){
      arg_validation_messages[[arg_counter]] <-
        assertthat::validate_that(
          F,
          msg = paste0("ERROR: Argument 'drug_action_type' must be a character vector",
                       " with any of the following action types: ",
                       paste(valid_drug_action_types, collapse=", ")))
      arg_counter <- arg_counter + 1

    }

  }

  if(!is.null(drug_indication_main)){

    #arg_counter <- arg_counter + 1
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_indication_main),
        msg = "ERROR: Argument 'drug_indication_main' must be a vector of type 'character'")

    arg_counter <- arg_counter + 1

    if(length(unique(drug_indication_main %in% valid_indications_main)) > 1){
      arg_validation_messages[[arg_counter]] <-
        assertthat::validate_that(
          F,
          msg = paste0("ERROR: Argument 'drug_indication_main' must be a character vector",
                       " with any of the following action types: ",
                       paste(valid_indications_main, collapse=", ")))
      arg_counter <- arg_counter + 1

    }

  }

  if(!is.null(drug_target)){
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        drug_is_targeted == T,
        msg = "ERROR: Argument 'drug_is_targeted' must be set to TRUE when 'drug_target' is non-NULL")
    arg_counter <- arg_counter + 1

    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_target),
        msg = "ERROR: Argument 'drug_target' must be a vector of type 'character'")
    arg_counter <- arg_counter + 1

  }

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_alkylating_agent),
      msg = "ERROR: Argument 'is_alkylating_agent' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_angiogenesis_inhibitor),
      msg = "ERROR: Argument 'is_angiogenesis_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_anthracycline),
      msg = "ERROR: Argument 'is_anthracycline' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_antimetabolite),
      msg = "ERROR: Argument 'is_antimetabolite' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_ar_antagonist),
      msg = "ERROR: Argument 'is_ar_antagonist' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_bet_inhibitor),
      msg = "ERROR: Argument 'is_bet_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_hdac_inhibitor),
      msg = "ERROR: Argument 'is_hdac_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_hedgehog_antagonist),
      msg = "ERROR: Argument 'is_hedgehog_antagonist' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_hormone_therapy),
      msg = "ERROR: Argument 'is_hormone_therapy' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_immune_checkpoint_inhibitor),
      msg = "ERROR: Argument 'is_immune_checkpoint_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_kinase_inhibitor),
      msg = "ERROR: Argument 'is_kinase_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_monoclonal_antibody),
      msg = "ERROR: Argument 'is_monoclonal_antibody' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_parp_inhibitor),
      msg = "ERROR: Argument 'is_parp_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_platinum_compound),
      msg = "ERROR: Argument 'is_platinum_compound' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_proteasome_inhibitor),
      msg = "ERROR: Argument 'is_proteasome_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_topoisomerase_inhibitor),
      msg = "ERROR: Argument 'is_topoisomerase_inhibitor' must be logical (TRUE/FALSE)")
  arg_counter <- arg_counter + 1

  arg_validation_messages[[arg_counter]] <-
    assertthat::validate_that(
      is.logical(is_tubulin_inhibitor),
      msg = "ERROR: Argument 'is_tubulin_inhibitor' must be logical (TRUE/FALSE)")


  i <- 1
  error_messages <- c()
  while(i <= length(arg_validation_messages)){
    if(!is.logical(arg_validation_messages[[i]])){
      error_messages <- c(error_messages, arg_validation_messages[[i]])
    }
    i <- i + 1
  }
  if(length(error_messages) > 0){
    cat('\n')
    cat(error_messages, sep="\n")
    return()
  }

  if(drug_is_targeted == T){
    drug_records <- drug_records %>%
      dplyr::filter(!is.na(target_symbol))

  }
  if(source_opentargets_only == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(opentargets_version))
    }
  }

  if(drug_has_blackbox_warning == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(is.na(drug_blackbox_warning) |
                        drug_blackbox_warning == T)

      if(nrow(drug_records) == 0){
        cat(paste0("WARNING: For the conditions listed below, NO drugs were found with a blackbox warning\n"))
        cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        cat('\n')
      }
    }
  }

  if(nrow(drug_records) > 0){
    drug_records <- drug_records %>%
      dplyr::filter(is.na(drug_year_first_approval) |
                      drug_year_first_approval >= drug_approved_later_than)

    if(nrow(drug_records) == 0){
      cat(paste0("WARNING: For the conditions listed below, NO drugs were found with an approval date greater than or equal to: ",
                 drug_approved_later_than),"\n")
      cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
      cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
      cat(paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))

      cat('\n')
    }
  }

  if(drug_minimum_phase_any_indication > 0){
    if(nrow(drug_records) > 0){

      ## include indications customly retrieved in DailyMed (assuming they are all max phase)
      drug_records <- drug_records %>%
        dplyr::filter(
          (!is.na(drug_max_ct_phase) &
             drug_max_ct_phase >= drug_minimum_phase_any_indication) |
            (is.na(drug_max_ct_phase) &
               drug_clinical_source == "DailyMedParseSN"))

      if(nrow(drug_records) == 0){
        cat(paste0("WARNING: For the conditions listed below, NO drugs were found with a clinical phase greater or equal than: ",
                   drug_minimum_phase_any_indication),"\n")
        cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        cat(paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        cat(paste0("Condition 4: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
        cat('\n')
      }

    }
  }

  if(drug_is_approved == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(drug_approved_indication) &
                        drug_approved_indication == T)
    }
    if(nrow(drug_records) == 0){
      cat(paste0("WARNING: For the conditions listed below, NO approved drugs were found\n"))
      cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
      cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
      cat(paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
      cat(paste0("Condition 4: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
      cat(paste0("Condition 5: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))

    }
  }

  all_drug_targets <- data.frame()

  if(!is.null(drug_target)){
    for(target in drug_target){
      df <- data.frame('target_symbol' = target,
                       stringsAsFactors = F)
      all_drug_targets <- dplyr::bind_rows(
        all_drug_targets, df
      )
    }
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::inner_join(all_drug_targets, by = "target_symbol")

      if(nrow(drug_records) == 0){
        cat(paste0("WARNING: For the conditions listed below, NO molecularly targeted drugs were found for the target proteins provided in the 'drug_target' argument: ",
                    paste(all_drug_targets$target_symbol, collapse=", "),"\n"))
        cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        cat(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        cat(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        cat(paste0("Condition 5: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
        cat(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        cat('\n')
      }else{
        drug_records <- drug_records %>%
          dplyr::arrange(desc(drug_max_ct_phase), desc(drug_year_first_approval))
      }
    }

  }

  if(!is.null(drug_action_type)){

    if(nrow(drug_records) > 0){
      drug_records_action_type <- data.frame()
      drug_records$atype <- drug_records$drug_action_type

      for(actype in drug_action_type){

        recs <- drug_records %>%
          dplyr::filter(
            !is.na(atype) &
              stringr::str_detect(atype, actype)
          )
        if(nrow(recs) > 0){
          drug_records_action_type <- drug_records_action_type %>%
            dplyr::bind_rows(recs)
        }
      }


      if(nrow(drug_records_action_type) == 0){
        cat(paste0("WARNING: For the conditions listed below, no drug records were found for the action types provided in the argument 'drug_action_type': "
                   , paste(drug_action_type, collapse=", ")),"\n")
        cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        cat(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        cat(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        cat(paste0("Condition 5: Drugs with approval later than: ", drug_approved_later_than, "\n"))
        cat(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        if(nrow(all_drug_targets) > 0){
          cat(paste0("Condition 7: Drugs targeted towards the following targets: ", paste(all_drug_targets$target_symbol, collapse=", "), "\n"))
        }
        cat('\n')
      }else{
        drug_records <- drug_records_action_type %>%
          dplyr::select(-atype)
      }
    }

  }


  if(!is.null(drug_indication_main)){

    if(nrow(drug_records) > 0){
      drug_records_indication <- data.frame()

      for(ind in drug_indication_main){

        recs <- drug_records %>%
          dplyr::filter(
            !is.na(primary_site) &
              stringr::str_detect(primary_site, ind)
          )
        if(nrow(recs) > 0){
          drug_records_indication <- drug_records_indication %>%
            dplyr::bind_rows(recs)
        }
      }


      if(nrow(drug_records_indication) == 0){
        cat(paste0("WARNING: For the conditions listed below, no drug records were found for the indications (main tumor types) provided in the argument 'drug_indication_main': "
                   , paste(drug_indication_main, collapse=", ")),"\n")
        cat(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        cat(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        cat(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        cat(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        cat(paste0("Condition 5: Drugs with approval later than: ", drug_approved_later_than, "\n"))
        cat(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        if(nrow(all_drug_targets) > 0){
          cat(paste0("Condition 7: Drugs targeted towards the following targets only: ", paste(all_drug_targets$target_symbol, collapse=", "), "\n"))
        }
        if(!is.null(drug_action_type)){
          cat(paste0("Condition 8: Drugs with the following action types only: ", paste(drug_action_type, collapse=", "), "\n"))
        }
        cat('\n')
      }else{
        drug_records <- drug_records_indication
      }
    }

  }


  if(list_per_drug_synonym == F){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::select(-nci_concept_synonym) %>%
        dplyr::distinct()
    }
  }

  if(output_resolution == "drug"){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::select(-c(target_symbol,
                         target_ensembl_gene_id,
                         target_entrezgene,
                         target_type,
                         target_genename,
                         target_uniprot_id,
                         disease_efo_id,
                         cui,
                         cui_name,
                         drug_clinical_id,
                         comb_regimen_indication,
                         drug_approved_indication,
                         drug_clinical_source)) %>%
        dplyr::distinct()

      drug_records <- as.data.frame(drug_records %>%
        dplyr::group_by(dplyr::across(-c(disease_efo_label, primary_site, drug_max_phase_indication))) %>%
        dplyr::summarise(disease_indication = paste(unique(sort(disease_efo_label)), collapse="|"),
                         disease_indicaton_max_phase = paste(unique(sort(drug_max_phase_indication)), collapse="|"),
                         disease_main_group = paste(unique(sort(primary_site)), collapse="|"),
                         .groups = "drop")
      )


    }
  }

  if(output_resolution == "drug2target"){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::select(-c(disease_efo_id,
                         cui,
                         cui_name,
                         drug_clinical_id,
                         comb_regimen_indication,
                         drug_approved_indication,
                         drug_clinical_source)) %>%
        dplyr::distinct()

      drug_records <- as.data.frame(
        drug_records %>%
          dplyr::group_by(dplyr::across(-c(disease_efo_label, primary_site, drug_max_phase_indication))) %>%
          dplyr::summarise(disease_indication = paste(unique(sort(disease_efo_label)), collapse="|"),
                           disease_indicaton_max_phase = paste(unique(sort(drug_max_phase_indication)), collapse="|"),
                           disease_main_group = paste(unique(sort(primary_site)), collapse="|"),
                           .groups = "drop")
      )
    }
  }

  if(output_resolution == "drug2target2indication"){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::select(-c(drug_clinical_id)) %>%
        dplyr::distinct()
    }
  }


  if(is_alkylating_agent == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(alkylating_agent) & alkylating_agent == T)
    }
  }

  if(is_angiogenesis_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(angiogenesis_inhibitor) & angiogenesis_inhibitor == T)
    }
  }

  if(is_anthracycline == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(anthracycline) & anthracycline == T)
    }
  }

  if(is_antimetabolite == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(antimetabolite) & antimetabolite == T)
    }
  }

  if(is_ar_antagonist == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(ar_antagonist) & ar_antagonist == T)
    }
  }
  if(is_bet_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(bet_inhibitor) & bet_inhibitor == T)
    }
  }

  if(is_hedgehog_antagonist == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(hedgehog_antagonist) & hedgehog_antagonist == T)
    }
  }
  if(is_hdac_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(hdac_inhibitor) & hdac_inhibitor == T)
    }
  }

  if(is_hormone_therapy == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(hormone_therapy) & hormone_therapy == T)
    }
  }

  if(is_immune_checkpoint_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(immune_checkpoint_inhibitor) & immune_checkpoint_inhibitor == T)
    }
  }

  if(is_kinase_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(kinase_inhibitor) & kinase_inhibitor == T)
    }
  }

  if(is_monoclonal_antibody == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(monoclonal_antibody) & monoclonal_antibody == T)
    }
  }

  if(is_parp_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(parp_inhibitor) & parp_inhibitor == T)
    }
  }

  if(is_platinum_compound == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(platinum_compound) & platinum_compound == T)
    }
  }

  if(is_proteasome_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(proteasome_inhibitor) & proteasome_inhibitor == T)
    }
  }

  if(is_topoisomerase_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(topoisomerase_inhibitor) & topoisomerase_inhibitor == T)
    }
  }

  if(is_tubulin_inhibitor == T){
    if(nrow(drug_records) > 0){
      drug_records <- drug_records %>%
        dplyr::filter(!is.na(tubulin_inhibitor) & tubulin_inhibitor == T)
    }
  }


  if(output_style == "narrow"){
    cols_for_stripping_classes <-
      c('angiogenesis_inhibitor','hormone_therapy','tubulin_inhibitor',
        'bet_inhibitor','proteasome_inhibitor','kinase_inhibitor',
        'parp_inhibitor','hdac_inhibitor','alkylating_agent',
        'antimetabolite','ar_antagonist','immune_checkpoint_inhibitor',
        'topoisomerase_inhibitor','monoclonal_antibody',
        'anthracycline','hedgehog_antagonist','platinum_compound')

    cols_for_stripping_target <-
      c('target_genename','target_ensembl_gene_id','target_type',
        'target_entrezgene','target_uniprot_id')

    cols_for_stripping_drug <-
      c('drug_name',
        'drug_clinical_source','drug_clinical_id',
        'nci_t','nci_concept_synonym','nci_concept_synonym_all',
        'drug_max_ct_phase, comb_regimen_indication')

    cols_for_stripping_disease <-
      c('cui','cui_name')

    for(c in cols_for_stripping_classes){
      if(c %in% colnames(drug_records)){
        drug_records[,c] <- NULL
      }
    }
    for(c in cols_for_stripping_target){
      if(c %in% colnames(drug_records)){
        drug_records[,c] <- NULL
      }
    }
    for(c in cols_for_stripping_drug){
      if(c %in% colnames(drug_records)){
        drug_records[,c] <- NULL
      }
    }
    for(c in cols_for_stripping_disease){
      if(c %in% colnames(drug_records)){
        drug_records[,c] <- NULL
      }
    }
  }

  if(nrow(drug_records) > 0){
    drug_records <- drug_records %>%
      dplyr::select(nci_concept_display_name,
                    drug_type,
                    molecule_chembl_id,
                    drug_action_type,
                    nci_concept_definition,
                    dplyr::everything()
                    ) %>%
      dplyr::arrange(desc(drug_year_first_approval),
                     desc(drug_max_ct_phase),
                     desc(opentargets_version),
                     nchar(nci_concept_display_name))

    drug_records_content <- drug_records %>%
      dplyr::select(-c(nci_version, chembl_version, opentargets_version))
    drug_records_version <- drug_records %>%
      dplyr::select(nci_version, chembl_version, opentargets_version)

    drug_records <- dplyr::bind_cols(drug_records_content,
                                     drug_records_version)
  }


  return(drug_records)


}
