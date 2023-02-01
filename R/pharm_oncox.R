#' Get antineoplastic drugs and drug regimens
#'
#' @description
#' Downloads preprocessed datasets to a local cache directory and returns a selected set of drugs
#' based on various criteria set by the user.
#'
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding drug resources used
#' * `records` - a data frame with drug records
#'
#' @param cache_dir local cache directory for data retrieval
#' @param force_download force download data from remote repository even if data exists in cache
#' @param exclude_salt_forms exclude salt forms of drugs
#' @param exclude_adc exclude antibody-drug conjugates (ADCs)
#' @param drug_is_targeted logical indicating if resulting drug records should contain drugs with known drug targets only
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
#' @param output_resolution dictate output record resolution ('drug','drug2target','drug2target2indication')
#' @param inhibitor_only logical indicating to return only drugs with inhibitory mechanism-of-action 
#' @param is_alkylating_agent logical indicating if only this drug class is to be considered
#' @param is_angiogenesis_inhibitor logical indicating if only this drug class is to be considered
#' @param is_anthracycline logical indicating if only this drug class is to be considered
#' @param is_antimetabolite logical indicating if only this drug class is to be considered
#' @param is_ar_antagonist logical indicating if only this drug class is to be considered
#' @param is_bet_inhibitor logical indicating if only this drug class is to be considered
#' @param is_hedgehog_antagonist logical indicating if only this drug class is to be considered
#' @param is_hdac_inhibitor logical indicating if only this drug class is to be considered
#' @param is_hormone_therapy logical indicating if only this drug class is to be considered
#' @param is_iap_inhibitor logical indicating if only this drug class is to be considered
#' @param is_immune_checkpoint_inhibitor logical indicating if only this drug class is to be considered
#' @param is_kinase_inhibitor logical indicating if only this drug class is to be considered
#' @param is_monoclonal_antibody logical indicating if only this drug class is to be considered
#' @param is_parp_inhibitor logical indicating if only this drug class is to be considered
#' @param is_platinum_compound logical indicating if only this drug class is to be considered
#' @param is_proteasome_inhibitor logical indicating if only this drug class is to be considered
#' @param is_topoisomerase_inhibitor logical indicating if only this drug class is to be considered
#' @param is_tubulin_inhibitor logical indicating if only this drug class is to be considered
#'
#' @return
#'
#' The `records` data frame contains the following columns (only selected columns
#' will be shown based on the value of `output_resolution`)
#'
#' \itemize{
#'   \item \emph{drug_id} - drug identifier (pharmaOncoX)
#'   \item \emph{drug_name} - primary drug name (upper case, NCI Thesaurus)
#'   \item \emph{drug_type} - type of drug molecule (Antibody, small molecule etc)
#'   \item \emph{molecule_chembl_id} - ChEMBL compound identifier
#'   \item \emph{drug_action_type} - main action elicited by drug (antagonist, inhibitor, stabiliser etc)
#'   \item \emph{drug_alias} - collection of unambiguous drug aliases (separated by '|')
#'   \item \emph{nci_concept_definition} - detailed description of drug mechanism-of-action (NCI Thesaurus)
#'   \item \emph{opentargets} - logical - drug is found in the Open Targets Platform resource
#'   \item \emph{is_salt} - logical - drug record represents a salt form (excluded by default)
#'   \item \emph{is_adc} - logical - drug record represents an antibody-drug conjugate (ADC - excluded by default)
#'   \item \emph{drug_blacbox_warning} - logical indicating if drug has blackbox warning
#'   \item \emph{nci_t} - NCI thesaurus identifier
#'   \item \emph{immune_checkpoint_inhibitor} - logical indicating if drug is an immune checkpoint inhibitor
#'   \item \emph{topoisomerase_inhibitor} - logical indicating if drug is a topoisomerase inhibitor
#'   \item \emph{iap_inhibitor} - logical indicating if drug is an IAP inhibitor
#'   \item \emph{tubulin_inhibitor} - logical indicating if drug is a tubulin inhibitor
#'   \item \emph{kinase_inhibitor} - logical indicating if drug is a kinase inhibitor
#'   \item \emph{hdac_inhibitor} - logical indicating if drug is a HDAC inhibitor
#'   \item \emph{parp_inhibitor} - logical indicating if drug is a PARP inhibitor
#'   \item \emph{bet_inhibitor} - logical indicating if drug is a BET inhibitor
#'   \item \emph{ar_antagonist} - logical indicating if drug is an androgen receptor antagonist
#'   \item \emph{monoclonal_antibody} - logical indicating if drug is an antibody
#'   \item \emph{antimetabolite} - logical indicating if drug is an antimetabolite
#'   \item \emph{angiogenesis_inhibitor} - logical indicating if drug is an angiogenesis inhibitor
#'   \item \emph{alkylating_agent} - logical indicating if drug is an alkylating agent
#'   \item \emph{platinum_compound} - logical indicating if drug is a platinum compound
#'   \item \emph{anthracycline} - logical indicating if drug is an anthracycline
#'   \item \emph{proteasome_inhibitor} - logical indicating if drug is a proteasome inhibitor
#'   \item \emph{hormone_therapy} - logical indicating if drug is a hormone therapy
#'   \item \emph{hedgehog_antagonist} - logical indicating if drug is a hedgehog pathway antagonist
#'   \item \emph{fda_epc_category} - FDA Established Pharmacologic Class (EPC) text phrase
#'   \item \emph{target_symbol} - gene symbol of drug target
#'   \item \emph{target_entrezgene} - Entrez gene identifier of drug target
#'   \item \emph{target_genename} - gene name/description of drug target
#'   \item \emph{target_ensembl_gene_id} - Ensembl gene identifier of drug target
#'   \item \emph{target_type} - type of drug target (single protein, protein family etc.)
#'   \item \emph{drug_max_phase_indication} - maximum clinical phase for drug (given indication)
#'   \item \emph{drug_approved_indication} - logical indicating if drug has an approved indication
#'   \item \emph{drug_year_first_approval} - year drug was first approved
#'   \item \emph{drug_max_ct_phase} - maximum clinical phase for drug (any indication)
#'   \item \emph{disease_efo_id} - EFO (Experimental Factor Ontology) identifier for drug indication
#'   \item \emph{disease_efo_label} - EFO (Experimental Factor Ontology) label for drug indication
#'   \item \emph{primary_site} - primary tumor site/type (obtained through https://github.com/sigven/oncoPhenoMap)
#'   \item \emph{drug_clinical_source} - underlying source for drug entry (DailyMed, clinicaltrials.gov, FDA etc.)
#'   \item \emph{drug_clinical_id} - drug clinical identifier (clinicaltrials.gov, DailyMed, FDA etc.)
#' }
#'
#'
#'
#' @export
#'

get_drugs <- function(
    cache_dir = NA,
    force_download = F,
    exclude_salt_forms = T,
    exclude_adc = T,
    drug_is_targeted = F,
    drug_is_approved = F,
    drug_target = NULL,
    drug_action_type = NULL,
    drug_indication_main = NULL,
    source_opentargets_only = F,
    drug_has_blackbox_warning = F,
    drug_approved_later_than = 1939,
    drug_minimum_phase_any_indication = 0,
    output_resolution = "drug2target2indication",
    inhibitor_only = F,
    is_alkylating_agent = F,
    is_angiogenesis_inhibitor = F,
    is_anthracycline = F,
    is_antimetabolite = F,
    is_ar_antagonist = F,
    is_bet_inhibitor = F,
    is_hedgehog_antagonist = F,
    is_hdac_inhibitor = F,
    is_hormone_therapy = F,
    is_iap_inhibitor = F,
    is_immune_checkpoint_inhibitor = F,
    is_kinase_inhibitor = F,
    is_monoclonal_antibody = F,
    is_parp_inhibitor = F,
    is_platinum_compound = F,
    is_proteasome_inhibitor = F,
    is_topoisomerase_inhibitor = F,
    is_tubulin_inhibitor = F) {

  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  valid_output_resolutions <- c("drug","drug2target","drug2target2indication")
  valid_drug_action_types <-
    c("INHIBITOR",
      "AGONIST",
      "MODULATOR",
      "ANTAGONIST",
      "BLOCKER",
      "ACTIVATOR",
      "BINDING AGENT",
      "OPENER",
      "STABILISER",
      "CROSS-LINKING AGENT",
      "DISRUPTING AGENT",
      "OTHER")

  valid_indications_main <-
    c("Adrenal Gland",
      "Biliary Tract",
      "Bladder/Urinary Tract",
      "Bone",
      "Breast",
      "Cervix",
      "CNS/Brain",
      "Colon/Rectum",
      "Esophagus/Stomach",
      "Eye",
      "Head and Neck",
      "Kidney",
      "Liver",
      "Lung",
      "Lymphoid",
      "Myeloid",
      "Ovary/Fallopian Tube",
      "Pancreas",
      "Penis",
      "Peripheral Nervous System",
      "Peritoneum",
      "Pleura",
      "Prostate",
      "Skin",
      "Soft Tissue",
      "Testis",
      "Thymus",
      "Thyroid",
      "Uterus",
      "Vulva/Vagina")

  arg_validation_messages <- list()

  arg_validation_messages[[1]] <-
     assertthat::validate_that(
       is.numeric(drug_approved_later_than),
       msg = "ERROR: Argument 'drug_approved_later_than' must be of type 'numeric'")
  arg_validation_messages[[2]] <-
    assertthat::validate_that(
      drug_approved_later_than >= 1939 & drug_approved_later_than <= 2022,
      msg = "ERROR: Argument 'drug_approved_later_than' must be larger than 1939 and less than or equal to 2022")
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
      output_resolution %in% valid_output_resolutions,
      msg = "ERROR: Argument 'output_resolution' must be either 'drug','drug2target', or 'drug2target2indication'")

  arg_counter <- 7
  if (!is.null(drug_action_type)) {
    # arg_validation_messages[[arg_counter]] <-
    #   assertthat::validate_that(
    #     drug_is_targeted == T,
    #     msg = "ERROR: Argument 'drug_is_targeted' must be set to TRUE when 'drug_action_type' is non-NULL")

    arg_counter <- arg_counter + 1
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_action_type),
        msg = "ERROR: Argument 'drug_action_type' must be a vector of type 'character'")

    arg_counter <- arg_counter + 1

    if (length(unique(drug_action_type %in% valid_drug_action_types)) > 1) {
      arg_validation_messages[[arg_counter]] <-
        assertthat::validate_that(
          F,
          msg = paste0(
            "ERROR: Argument 'drug_action_type' must be a character vector",
            " with any of the following action types: ",
            paste(valid_drug_action_types, collapse = ", ")))
      arg_counter <- arg_counter + 1

    }

  }

  if (!is.null(drug_indication_main)) {

    #arg_counter <- arg_counter + 1
    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_indication_main),
        msg = "ERROR: Argument 'drug_indication_main' must be a vector of type 'character'")

    arg_counter <- arg_counter + 1

    if (length(unique(drug_indication_main %in% valid_indications_main)) > 1) {
      arg_validation_messages[[arg_counter]] <-
        assertthat::validate_that(
          F,
          msg = paste0(
            "ERROR: Argument 'drug_indication_main' must be a character vector",
            " with any of the following action types: ",
            paste(valid_indications_main, collapse = ", ")))
      arg_counter <- arg_counter + 1

    }

  }

  if (!is.null(drug_target)) {
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
      is.logical(is_iap_inhibitor),
      msg = "ERROR: Argument 'is_inhibitor_therapy' must be logical (TRUE/FALSE)")
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
  while (i <= length(arg_validation_messages)) {
    if (!is.logical(arg_validation_messages[[i]])) {
      error_messages <- c(error_messages, arg_validation_messages[[i]])
    }
    i <- i + 1
  }
  if (length(error_messages) > 0) {
    lgr::lgr$info('\n')
    lgr::lgr$info(error_messages, sep = "\n")
    return()
  }


  ## Get full list of drug records
  all_drug_recs <- get_drug_records(cache_dir, force_download)
  drug_records <- as.data.frame(all_drug_recs[['records']])
  metadata <- as.data.frame(all_drug_recs[['metadata']])

  if (exclude_salt_forms == T) {
    drug_records <- drug_records |>
      dplyr::filter(.data$is_salt == F)
  }

  if (exclude_adc == T) {
    drug_records <- drug_records |>
      dplyr::filter(.data$is_adc == F)
  }


  if (drug_is_targeted == T) {
    drug_records <- drug_records |>
      dplyr::filter(!is.na(.data$target_symbol))

  }
  if (source_opentargets_only == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(.data$opentargets == T)
    }
  }
  
  if (inhibitor_only == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(.data$inhibition_moa == T)
    }
  }

  if (drug_has_blackbox_warning == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(is.na(.data$drug_blackbox_warning) |
                        .data$drug_blackbox_warning == T)

      if (nrow(drug_records) == 0) {
        lgr::lgr$info(
          paste0("WARNING: For the conditions listed below, NO drugs were found with a blackbox warning\n"))
        lgr::lgr$info(
          paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        lgr::lgr$info(
          paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info('\n')
      }
    }
  }

  if (nrow(drug_records) > 0) {
    drug_records <- drug_records |>
      dplyr::filter(is.na(.data$drug_year_first_approval) |
                      .data$drug_year_first_approval >= drug_approved_later_than)

    if (nrow(drug_records) == 0) {
      lgr::lgr$info(
        paste0("WARNING: For the conditions listed below, NO drugs were found with an approval date greater than or equal to: ",
                 drug_approved_later_than),"\n")
      lgr::lgr$info(
        paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
      lgr::lgr$info(
        paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
      lgr::lgr$info(
        paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))

      lgr::lgr$info('\n')
    }
  }

  if (drug_minimum_phase_any_indication > 0) {
    if (nrow(drug_records) > 0) {

      ## include indications customly retrieved in DailyMed (assuming they are all max phase)
      drug_records <- drug_records |>
        dplyr::filter(
          (!is.na(.data$drug_max_ct_phase) &
             .data$drug_max_ct_phase >= drug_minimum_phase_any_indication) |
            (is.na(.data$drug_max_ct_phase) &
               .data$drug_clinical_source == "DailyMedParseSN"))

      if (nrow(drug_records) == 0) {
        lgr::lgr$info(paste0("WARNING: For the conditions listed below, NO drugs were found with a clinical phase greater or equal than: ",
                   drug_minimum_phase_any_indication),"\n")
        lgr::lgr$info(
          paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        lgr::lgr$info(
          paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info(
          paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(
          paste0("Condition 4: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
        lgr::lgr$info('\n')
      }

    }
  }

  if (drug_is_approved == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$drug_approved_indication) &
                        .data$drug_approved_indication == T)
    }
    if (nrow(drug_records) == 0) {
      lgr::lgr$info(
        paste0("WARNING: For the conditions listed below, NO approved drugs were found\n"))
      lgr::lgr$info(
        paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
      lgr::lgr$info(
        paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
      lgr::lgr$info(
        paste0("Condition 3: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
      lgr::lgr$info(
        paste0("Condition 4: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
      lgr::lgr$info(
        paste0("Condition 5: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))

    }
  }

  all_drug_targets <- data.frame()

  if (!is.null(drug_target)) {
    for (target in drug_target) {
      df <- data.frame('target_symbol' = target,
                       stringsAsFactors = F)
      all_drug_targets <- dplyr::bind_rows(
        all_drug_targets, df
      )
    }
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::inner_join(
          all_drug_targets, by = "target_symbol")

      if (nrow(drug_records) == 0) {
        lgr::lgr$info(paste0("WARNING: For the conditions listed below, NO molecularly targeted drugs were found for the target proteins provided in the 'drug_target' argument: ",
                    paste(all_drug_targets$target_symbol, collapse = ", "),"\n"))
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than only: ", drug_approved_later_than, "\n"))
        lgr::lgr$info(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        lgr::lgr$info('\n')
      }else{
        drug_records <- drug_records |>
          dplyr::arrange(
            dplyr::desc(.data$drug_max_ct_phase),
            dplyr::desc(.data$drug_year_first_approval))
      }
    }

  }

  if (!is.null(drug_action_type)) {

    if (nrow(drug_records) > 0) {
      drug_records_action_type <- data.frame()
      drug_records$atype <- drug_records$drug_action_type

      for (actype in drug_action_type) {

        recs <- drug_records |>
          dplyr::filter(
            !is.na(.data$atype) &
              stringr::str_detect(.data$atype, actype)
          )
        if (nrow(recs) > 0) {
          drug_records_action_type <- drug_records_action_type |>
            dplyr::bind_rows(recs)
        }
      }


      if (nrow(drug_records_action_type) == 0) {
        lgr::lgr$info(paste0("WARNING: For the conditions listed below, no drug records were found for the action types provided in the argument 'drug_action_type': "
                   , paste(drug_action_type, collapse = ", ")),"\n")
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than: ", drug_approved_later_than, "\n"))
        lgr::lgr$info(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        if (nrow(all_drug_targets) > 0) {
          lgr::lgr$info(paste0("Condition 7: Drugs targeted towards the following targets: ", paste(all_drug_targets$target_symbol, collapse = ", "), "\n"))
        }
        lgr::lgr$info('\n')
      }else{
        drug_records <- drug_records_action_type |>
          dplyr::select(-.data$atype)
      }
    }

  }


  if (!is.null(drug_indication_main)) {

    if (nrow(drug_records) > 0) {
      drug_records_indication <- data.frame()

      for (ind in drug_indication_main) {

        recs <- drug_records |>
          dplyr::filter(
            !is.na(.data$primary_site) &
              stringr::str_detect(.data$primary_site, ind)
          )
        if (nrow(recs) > 0) {
          drug_records_indication <- drug_records_indication |>
            dplyr::bind_rows(recs)
        }
      }


      if (nrow(drug_records_indication) == 0) {
        lgr::lgr$info(paste0("WARNING: For the conditions listed below, no drug records were found for the indications (main tumor types) provided in the argument 'drug_indication_main': "
                   , paste(drug_indication_main, collapse = ", ")),"\n")
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_is_targeted, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than: ", drug_approved_later_than, "\n"))
        lgr::lgr$info(paste0("Condition 6: Drugs with minimum clinical phase (any indication) only: ", drug_minimum_phase_any_indication, "\n"))
        if (nrow(all_drug_targets) > 0) {
          lgr::lgr$info(paste0("Condition 7: Drugs targeted towards the following targets only: ", paste(all_drug_targets$target_symbol, collapse = ", "), "\n"))
        }
        if (!is.null(drug_action_type)) {
          lgr::lgr$info(paste0(
            "Condition 8: Drugs with the following action types only: ",
            paste(drug_action_type, collapse = ", "), "\n"))
        }
        lgr::lgr$info('\n')
      }else{
        drug_records <- drug_records_indication
      }
    }

  }

  if (is_alkylating_agent == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$alkylating_agent) &
                        .data$alkylating_agent == T)
    }
  }

  if (is_angiogenesis_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$angiogenesis_inhibitor) &
                        .data$angiogenesis_inhibitor == T)
    }
  }

  if (is_anthracycline == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$anthracycline) &
                        .data$anthracycline == T)
    }
  }

  if (is_antimetabolite == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$antimetabolite) &
                        .data$antimetabolite == T)
    }
  }

  if (is_ar_antagonist == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$ar_antagonist) &
                        .data$ar_antagonist == T)
    }
  }
  if (is_bet_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$bet_inhibitor) &
                        .data$bet_inhibitor == T)
    }
  }

  if (is_hedgehog_antagonist == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$hedgehog_antagonist) &
                        .data$hedgehog_antagonist == T)
    }
  }
  if (is_hdac_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$hdac_inhibitor) &
                        .data$hdac_inhibitor == T)
    }
  }

  if (is_hormone_therapy == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$hormone_therapy) &
                        .data$hormone_therapy == T)
    }
  }
  
  if (is_iap_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$iap_inhibitor) &
                        .data$iap_inhibitor == T)
    }
  }

  if (is_immune_checkpoint_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$immune_checkpoint_inhibitor) &
                        .data$immune_checkpoint_inhibitor == T)
    }
  }

  if (is_kinase_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$kinase_inhibitor) &
                        .data$kinase_inhibitor == T)
    }
  }

  if (is_monoclonal_antibody == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$monoclonal_antibody) &
                        .data$monoclonal_antibody == T)
    }
  }

  if (is_parp_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$parp_inhibitor) &
                        .data$parp_inhibitor == T)
    }
  }

  if (is_platinum_compound == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$platinum_compound) &
                        .data$platinum_compound == T)
    }
  }

  if (is_proteasome_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$proteasome_inhibitor) &
                        .data$proteasome_inhibitor == T)
    }
  }

  if (is_topoisomerase_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$topoisomerase_inhibitor) &
                        .data$topoisomerase_inhibitor == T)
    }
  }

  if (is_tubulin_inhibitor == T) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$tubulin_inhibitor) &
                        .data$tubulin_inhibitor == T)
    }
  }

  lgr::lgr$info(
    paste0("Record set satisfying user-defined criteria: n = ",
           nrow(drug_records)))

  if (output_resolution == "drug") {
    lgr::lgr$info(
      "Collapsing record set - providing output on a 'per_drug' resolution")
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::select(-c(.data$target_symbol,
                         .data$target_ensembl_gene_id,
                         .data$target_entrezgene,
                         .data$target_type,
                         .data$target_genename,
                         .data$disease_efo_id,
                         .data$cui,
                         .data$cui_name,
                         #.data$comb_regimen_indication,
                         .data$drug_approved_indication,
                         .data$drug_clinical_source)) |>
        dplyr::distinct()

      drug_records <- as.data.frame(
        drug_records |>
          dplyr::group_by(
            dplyr::across(-c(.data$disease_efo_label,
                             .data$primary_site,
                             .data$drug_clinical_id,
                             .data$drug_max_phase_indication))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            disease_indication = paste(
              unique(sort(.data$disease_efo_label)), collapse = "|"),
            disease_indicaton_max_phase = paste(
              unique(sort(.data$drug_max_phase_indication)), collapse = "|"),
            disease_main_group = paste(
              unique(sort(.data$primary_site)), collapse = "|"),
            .groups = "drop")
      )


    }
    lgr::lgr$info(
      paste0("Final record set: n = ", nrow(drug_records), " records"))

  }

  if (output_resolution == "drug2target") {
    lgr::lgr$info("Collapsing record set - providing output on a 'per_drug_target' resolution")
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::select(-c(.data$disease_efo_id,
                         .data$cui,
                         .data$cui_name,
                         #.data$comb_regimen_indication,
                         .data$drug_approved_indication,
                         .data$drug_clinical_source)) |>
        dplyr::distinct()

      drug_records <- as.data.frame(
        drug_records |>
          dplyr::group_by(
            dplyr::across(
              -c(.data$disease_efo_label,
                 .data$primary_site,
                 .data$drug_clinical_id,
                 .data$drug_max_phase_indication))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            disease_indication = paste(
              unique(sort(.data$disease_efo_label)), collapse = "|"),
            disease_indicaton_max_phase = paste(
              unique(sort(.data$drug_max_phase_indication)), collapse = "|"),
            disease_main_group = paste(
              unique(sort(.data$primary_site)), collapse = "|"),
            .groups = "drop")
      )
    }
    lgr::lgr$info(
      paste0("Final record set: n = ", nrow(drug_records), " records"))
  }

  if (output_resolution == "drug2target2indication") {
    lgr::lgr$info(
      "Collapsing record set - providing output on a 'per_drug_target_indication' resolution")
    if (nrow(drug_records) > 0) {
      drug_records <- as.data.frame(
        drug_records |> 
          dplyr::select(-c(.data$cui, .data$cui_name)) |>
          dplyr::distinct() |>
          dplyr::group_by(
            dplyr::across(
              -c(.data$drug_clinical_id,
                 .data$drug_clinical_source))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            drug_clinical_source = paste(
              unique(sort(.data$drug_clinical_source)), collapse = "|"),
            .groups = "drop"
          )
      )

    }
    lgr::lgr$info(
      paste0("Final record set: n = ", nrow(drug_records), " records"))
  }

  if (nrow(drug_records) > 0) {

    drug_records <- drug_records |>
      dplyr::select(
        .data$drug_id,
        .data$drug_name,
        .data$drug_type,
        .data$molecule_chembl_id,
        .data$drug_action_type,
        .data$drug_alias,
        .data$nci_concept_definition,
        .data$opentargets,
        dplyr::everything()
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$drug_year_first_approval),
        dplyr::desc(.data$drug_max_ct_phase),
        dplyr::desc(.data$opentargets),
        nchar(.data$drug_name))

  }

  oncodrugs <- list()
  oncodrugs[['records']] <- drug_records
  oncodrugs[['metadata']] <- metadata

  return(oncodrugs)

}

get_drug_biomarkers <- function(
    drug_query = "Vemurafenib",
    force_download = FALSE,
    cache_dir = NA){
  
  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  
  all_drug_recs <- get_drug_records(cache_dir, force_download)
  
  drugname_primary <- all_drug_recs$records |>
    dplyr::select("drug_name", 
                  "drug_id",
                  "molecule_chembl_id", 
                  "drug_alias") |>
    dplyr::distinct() |>
    dplyr::mutate(alias = tolower(.data$drug_name))
  
  drugname_alias <- drugname_primary |>
    dplyr::select(c("drug_id","drug_alias")) |>
    tidyr::separate_rows(.data$drug_alias, sep = "\\|") |>
    dplyr::mutate(alias = tolower(.data$drug_alias)) |>
    dplyr::select(c("drug_id", "alias")) |>
    dplyr::distinct()
  
  drugname_chembl <- drugname_primary |>
    dplyr::select("drug_id", "molecule_chembl_id") |>
    dplyr::filter(!is.na(.data$molecule_chembl_id)) |>
    dplyr::rename(alias = .data$molecule_chembl_id) |>
    dplyr::select(c("drug_id","alias")) |>
    dplyr::distinct()
  
  drugname_primary_lc <- drugname_primary |>
    dplyr::select(c("drug_id", "alias")) |>
    dplyr::distinct()
  
  all_aliases <- drugname_primary_lc |>
    dplyr::bind_rows(drugname_alias) |>
    dplyr::bind_rows(drugname_chembl) |>
    dplyr::mutate(alias = stringr::str_trim(.data$alias)) |>
    dplyr::filter(nchar(.data$alias) > 2) |>
    dplyr::distinct()
  
  drug_dictionary <- drugname_primary |>
    dplyr::select(c("drug_id", "drug_name")) |>
    dplyr::left_join(all_aliases, by = "drug_id") |>
    dplyr::filter(
      !stringr::str_detect(tolower(.data$drug_name),"regimen")
    )
  
  query_hit <- drug_dictionary |>
    dplyr::filter(.data$alias == tolower(drug_query))
  
  
  #biomarkers <- get_biomarkers(cache_dir = cache_dir,
                               #force_download = force_download)
  
  # therapy_biomarkers <- biomarkers$records |>
  #   dplyr::filter(!is.na(therapeutic_context)) |>
  #   dplyr::filter(nchar(therapeutic_context) > 0) |>
  #   dplyr::select(evidence_id, symbol, variant, 
  #                 variant_id, therapeutic_context,
  #                 alteration_type, cancer_type, 
  #                 molecule_chembl_id, evidence_level,
  #                 evidence_direction, evidence_description,
  #                 clinical_significance) |>
  #   dplyr::distinct() |>
  #   dplyr::filter(alteration_type != "BIA" & 
  #                   alteration_type != "METHYL" &
  #                   alteration_type != "CODON") |>
  #   dplyr::filter(!stringr::str_detect(
  #     therapeutic_context, "(hibitors|gens)$"
  #   )) |>
  #   dplyr::filter(!stringr::str_detect(
  #     therapeutic_context, "alone or in combination"
  #   ))
  #   
  # 
  # aa_markers <- therapy_biomarkers |>
  #   dplyr::filter(stringr::str_detect(variant, "p\\.[A-Z][0-9]{1,}"))
  # 
  # 
  
  ## check if hits
  
}


#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

utils::globalVariables(c("."))

