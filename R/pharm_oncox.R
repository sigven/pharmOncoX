#' Get antineoplastic drugs and drug regimens
#'
#' @description
#' Downloads preprocessed datasets to a local cache directory and returns a 
#' selected set of drugs based on various criteria set by the user.
#'
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding drug resources used
#' * `records` - a data frame with drug records
#'
#' @param cache_dir local cache directory for data retrieval
#' @param force_download force download data from remote repository even if 
#' data exists in cache
#' @param exclude_salt_forms exclude salt forms of drugs
#' @param exclude_adc exclude antibody-drug conjugates (ADCs)
#' @param drug_targeted_agent logical indicating if resulting drug records 
#' should only contain drugs that are considered targeted agents 
#' (i.e. not chemotherapeutic drugs/immunosuppressants)
#' @param drug_is_approved logical indicating if resulting drug records 
#' should contain approved drugs only
#' @param drug_target character vector with drug targets (gene symbols) 
#' for drug records included in results
#' @param drug_action_type character vector with drug action types to 
#' include in drug record list -
#' possible values "INHIBITOR","AGONIST","MODULATOR","ANTAGONIST",
#' "BLOCKER","ACTIVATOR","BINDING AGENT","OPENER",
#' "STABILISER","CROSS-LINKING AGENT",DISRUPTING AGENT","OTHER"
#' @param drug_action_inhibition logical indicating to only return drug records 
#' with inhibitory mechanism-of-action 
#' @param drug_indication_main character vector with main tumor types for 
#' which drug(s) are indicated. Possible values: "Adrenal Gland","Biliary Tract",
#' "Bladder/Urinary Tract","Bone","Breast","Cervix","CNS/Brain",
#' "Colon/Rectum","Esophagus/Stomach","Eye","Head and Neck","Kidney",
#' "Liver","Lung","Lymphoid","Myeloid","Ovary/Fallopian Tube",
#' "Pancreas","Penis","Peripheral Nervous System","Peritoneum",
#' "Pleura","Prostate","Skin","Soft Tissue","Testis","Thymus",
#' "Thyroid","Uterus","Vulva/Vagina"
#' @param drug_source_opentargets logical indicating if resulting drug records 
#' should contain drug records from Open Targets Platform/ChEMBL only
#' @param drug_cancer_indication logical indicating if resulting drug records 
#' should be for those indicated for cancer conditions only (from approved 
#' conditions, found in clinical trials etc.) 
#' @param drug_classified_cancer logical indicating if resulting drug 
#' records should be for those classified only in the "L" class of ATC (
#' "ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS") only
#' @param drug_has_blackbox_warning logical indicating if resulting drug 
#' records should contain drugs with black box warnings only
#' @param drug_approval_year only include records for drugs approved 
#' later than this date (year)
#' @param drug_minimum_phase_any_indication only include drug records that 
#' are in a clinical phase (any indication) greater or equal than this phase
#' @param output_resolution dictate output record resolution 
#' ('drug','drug2target','drug2target2indication')
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
#'   \item \emph{drug_action_type} - main action elicited by drug 
#'   (antagonist, inhibitor, stabiliser etc)
#'   \item \emph{drug_alias} - collection of unambiguous drug aliases 
#'   (separated by '|')
#'   \item \emph{nci_concept_definition} - detailed description of drug 
#'   mechanism-of-action (NCI Thesaurus)
#'   \item \emph{opentargets} - logical - drug is found in the Open Targets 
#'   Platform resource
#'   \item \emph{is_salt} - logical - drug record represents a salt form 
#'   (excluded by default)
#'   \item \emph{is_adc} - logical - drug record represents an 
#'   antibody-drug conjugate (ADC - excluded by default)
#'   \item \emph{drug_blacbox_warning} - logical indicating if drug has 
#'   blackbox warning
#'   \item \emph{nci_t} - NCI thesaurus identifier
#'   \item \emph{target_symbol} - gene symbol of drug target
#'   \item \emph{target_entrezgene} - Entrez gene identifier of drug target
#'   \item \emph{target_genename} - gene name/description of drug target
#'   \item \emph{target_ensembl_gene_id} - Ensembl gene identifier of drug target
#'   \item \emph{target_type} - type of drug target (single protein, protein 
#'   family etc.)
#'   \item \emph{drug_max_phase_indication} - maximum clinical phase for drug 
#'   (given indication)
#'   \item \emph{drug_approved_indication} - logical indicating if drug has 
#'   an approved indication
#'   \item \emph{drug_frac_cancer_indications} - fraction of drug indications that
#'   are for cancers
#'   \item \emph{drug_approved_noncancer} - logical indicating if drug
#'   is approved for a non-cancer disease
#'   \item \emph{drug_n_indications} - number of indications for the
#'   given drug (from approved indications, clinical trials etc)
#'   \item \emph{drug_year_first_approval} - year drug was first approved
#'   \item \emph{drug_max_ct_phase} - maximum clinical phase for drug 
#'   (any indication)
#'   \item \emph{disease_efo_id} - EFO (Experimental Factor Ontology) 
#'   identifier for drug indication
#'   \item \emph{disease_efo_label} - EFO (Experimental Factor Ontology) 
#'   label for drug indication
#'   \item \emph{primary_site} - primary tumor site/type (obtained through 
#'   https://github.com/sigven/oncoPhenoMap)
#'   \item \emph{drug_clinical_id} - drug clinical identifier 
#'   (clinicaltrials.gov, DailyMed, FDA etc.)
#'   \item \emph{drug_clinical_source} - underlying source for drug entry 
#'   (DailyMed, clinicaltrials.gov, FDA etc.)
#'   \item \emph{atc_code_level1} - drug identifier ATC (level 1)
#'   \item \emph{atc_level1} - drug label ATC (level 1)
#'   \item \emph{atc_code_level1} - drug identifier ATC (level 2)
#'   \item \emph{atc_level2} - drug label ATC (level 2)
#'   \item \emph{atc_code_level3} - drug identifier ATC (level 3)
#'   \item \emph{atc_level3} - drug label ATC (level 3)
#' }
#'
#'
#'
#' @export
#'

get_drugs <- function(
    cache_dir = NA,
    force_download = FALSE,
    exclude_salt_forms = TRUE,
    exclude_adc = TRUE,
    drug_targeted_agent = FALSE,
    drug_is_approved = FALSE,
    drug_target = NULL,
    drug_action_type = NULL,
    drug_indication_main = NULL,
    drug_source_opentargets = FALSE,
    drug_cancer_indication = TRUE, 
    drug_classified_cancer = TRUE,
    drug_has_blackbox_warning = FALSE,
    drug_approval_year = 1939,
    drug_minimum_phase_any_indication = 0,
    output_resolution = "drug2target2indication",
    drug_action_inhibition = F){

  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  valid_output_resolutions <- 
    c("drug","drug2target","drug2target2indication")
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
       is.numeric(drug_approval_year),
       msg = "ERROR: Argument 'drug_approval_year' must be of type 'numeric'")
  arg_validation_messages[[2]] <-
    assertthat::validate_that(
      drug_approval_year >= 1939 & drug_approval_year <= 2023,
      msg = "ERROR: Argument 'drug_approval_year' must be larger than 1939 and less than or equal to 2023")
  arg_validation_messages[[3]] <-
    assertthat::validate_that(
      is.logical(drug_is_approved),
      msg = "ERROR: Argument 'drug_is_approved' must be of type 'logical'")
  arg_validation_messages[[4]] <-
    assertthat::validate_that(
      is.logical(drug_source_opentargets),
      msg = "ERROR: Argument 'drug_source_opentargets' must be of type 'logical'")
  arg_validation_messages[[5]] <-
    assertthat::validate_that(
      is.logical(drug_classified_cancer),
      msg = "ERROR: Argument 'drug_classified_cancer' must be of type 'logical'")
  arg_validation_messages[[6]] <-
    assertthat::validate_that(
      is.logical(drug_cancer_indication),
      msg = "ERROR: Argument 'drug_cancer_indication' must be of type 'logical'")
  arg_validation_messages[[7]] <-
    assertthat::validate_that(
      is.logical(drug_targeted_agent),
      msg = "ERROR: Argument 'drug_targeted_agent' must be of type 'logical'")
  arg_validation_messages[[8]] <-
    assertthat::validate_that(
      is.logical(drug_has_blackbox_warning),
      msg = "ERROR: Argument 'drug_has_blackbox_warning' must be of type 'logical'")
  arg_validation_messages[[9]] <-
    assertthat::validate_that(
      output_resolution %in% valid_output_resolutions,
      msg = "ERROR: Argument 'output_resolution' must be either 'drug','drug2target', or 'drug2target2indication'")

  arg_counter <- 10
  if (!is.null(drug_action_type)) {
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
        drug_targeted_agent == TRUE,
        msg = "ERROR: Argument 'drug_targeted_agent' must be set to TRUE when 'drug_target' is non-NULL")
    arg_counter <- arg_counter + 1

    arg_validation_messages[[arg_counter]] <-
      assertthat::validate_that(
        is.character(drug_target),
        msg = "ERROR: Argument 'drug_target' must be a vector of type 'character'")
    arg_counter <- arg_counter + 1

  }

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

  if (exclude_salt_forms == TRUE) {
    drug_records <- drug_records |>
      dplyr::filter(.data$is_salt == FALSE)
  }

  if (exclude_adc == TRUE) {
    drug_records <- drug_records |>
      dplyr::filter(.data$is_adc == FALSE)
  }


  if (drug_targeted_agent == TRUE) {
    drug_records <- drug_records |>
      dplyr::filter(!is.na(.data$target_symbol)) |>
      
      ## exclude chemotherapeutic agents/immunosuppressants
      dplyr::filter(
        (is.na(.data$atc_level2) | 
           (!is.na(.data$atc_level2) & !stringr::str_detect(
             .data$atc_level2, "PLANT ALKALOIDS"
           ))) &
          (is.na(.data$atc_level3) |
             (!is.na(.data$atc_level3) &
                .data$atc_level3 != "Taxanes" &
                .data$atc_level3 != "Other immunosuppressants" &
                .data$atc_level3 != "Other cytotoxic antibiotics" &
                .data$atc_level3 != "Anthracyclines and related substances" &
                .data$atc_level3 != "Colony stimulating factors" &
                !stringr::str_detect(
                  tolower(.data$atc_level3), "analogues") &
                .data$atc_level3 != "Interferons" &
                .data$atc_level3 != "Interleukins" &
                .data$atc_level3 != "Podophyllotoxin derivatives" &
                .data$atc_level3 != "Nitrosoureas")) &
          !stringr::str_detect(
            .data$target_symbol,"TUBB|TUBA")
      )
  }
  if (drug_source_opentargets == TRUE) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(.data$opentargets == TRUE)
    }
  }
  
  if (drug_cancer_indication == TRUE){
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(stringr::str_detect(
          drug_cancer_relevance, "by_cancer_condition_otp"
        ))
    }
  }
  
  if (drug_classified_cancer == TRUE){
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(
          !is.na(atc_code_level1) &
          stringr::str_detect(
          atc_code_level1, "^L"
        ))
    }
  }
  
  if (drug_action_inhibition == TRUE) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(.data$inhibition_moa == TRUE)
    }
  }

  if (drug_has_blackbox_warning == TRUE) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(is.na(.data$drug_blackbox_warning) |
                        .data$drug_blackbox_warning == TRUE)

      if (nrow(drug_records) == 0) {
        lgr::lgr$info(
          paste0("WARNING: For the conditions listed below, NO drugs were found with a blackbox warning\n"))
        lgr::lgr$info(
          paste0("Condition 1: Molecularly targeted drugs only: ", drug_targeted_agent, "\n"))
        lgr::lgr$info(
          paste0("Condition 2: Open Targets Platform only: ", source_opentargets_only, "\n"))
        lgr::lgr$info('\n')
      }
    }
  }

  if (nrow(drug_records) > 0) {
    drug_records <- drug_records |>
      dplyr::filter(is.na(.data$drug_year_first_approval) |
                      .data$drug_year_first_approval >= drug_approval_year)

    if (nrow(drug_records) == 0) {
      lgr::lgr$info(
        paste0("WARNING: For the conditions listed below, NO drugs were found with an approval date greater than or equal to: ",
                 drug_approval_year),"\n")
      lgr::lgr$info(
        paste0("Condition 1: Molecularly targeted drugs only: ", drug_targeted_agent, "\n"))
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
             .data$drug_max_ct_phase >= drug_minimum_phase_any_indication))

      if (nrow(drug_records) == 0) {
        lgr::lgr$info(paste0(
          "WARNING: For the conditions listed below, NO drugs were found",
          "with a clinical phase greater or equal than: ",
          drug_minimum_phase_any_indication),"\n")
        lgr::lgr$info(
          paste0("Condition 1: Molecularly targeted drugs only: ", 
                 drug_targeted_agent, "\n"))
        lgr::lgr$info(
          paste0("Condition 2: Open Targets Platform only: ", 
                 drug_source_opentargets, "\n"))
        lgr::lgr$info(
          paste0("Condition 3: Drugs with blackbox warnings only: ", 
                 drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(
          paste0("Condition 4: Drugs with approval later than only: ", 
                 drug_approval_year, "\n"))
        lgr::lgr$info('\n')
      }

    }
  }

  if (drug_is_approved == TRUE) {
    if (nrow(drug_records) > 0) {
      drug_records <- drug_records |>
        dplyr::filter(!is.na(.data$drug_approved_indication) &
                        .data$drug_approved_indication == TRUE)
    }
    if (nrow(drug_records) == 0) {
      lgr::lgr$info(
        paste0("WARNING: For the conditions listed below, NO approved drugs were found\n"))
      lgr::lgr$info(
        paste0("Condition 1: Molecularly targeted drugs only: ", 
               drug_targeted_agent, "\n"))
      lgr::lgr$info(
        paste0("Condition 2: Open Targets Platform only: ", 
               drug_source_opentargets, "\n"))
      lgr::lgr$info(
        paste0("Condition 3: Drugs with blackbox warnings only: ", 
               drug_has_blackbox_warning, "\n"))
      lgr::lgr$info(
        paste0("Condition 4: Drugs with approval later than only: ", 
               drug_approval_year, "\n"))
      lgr::lgr$info(
        paste0("Condition 5: Drugs with minimum clinical phase (any indication) only: ", 
               drug_minimum_phase_any_indication, "\n"))

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
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_targeted_agent, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", drug_source_opentargets, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than only: ", drug_approval_year, "\n"))
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
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_targeted_agent, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", drug_source_opentargets, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than: ", drug_approval_year, "\n"))
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
        lgr::lgr$info(paste0("Condition 1: Molecularly targeted drugs only: ", drug_targeted_agent, "\n"))
        lgr::lgr$info(paste0("Condition 2: Open Targets Platform only: ", drug_source_opentargets, "\n"))
        lgr::lgr$info(paste0("Condition 3: Drugs with approved indications only: ", drug_is_approved, "\n"))
        lgr::lgr$info(paste0("Condition 4: Drugs with blackbox warnings only: ", drug_has_blackbox_warning, "\n"))
        lgr::lgr$info(paste0("Condition 5: Drugs with approval later than: ", drug_approval_year, "\n"))
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
                         .data$drug_approved_indication,
                         .data$drug_clinical_source)) |>
        dplyr::distinct()

      drug_records <- as.data.frame(
        drug_records |>
          dplyr::group_by(
            dplyr::across(-c(
              .data$atc_code_level1,
              .data$atc_level1,
              .data$atc_code_level2,
              .data$atc_level2,
              .data$atc_code_level3,
              .data$atc_level3,
              .data$disease_efo_label,
              .data$primary_site,
              .data$drug_clinical_id,
              .data$drug_max_phase_indication))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            disease_indication = paste(
              unique(sort(.data$disease_efo_label)), collapse = "|"),
            disease_indication_max_phase = paste(
              sort(unique(.data$drug_max_phase_indication),
                     decreasing = T), collapse = "|"),
            disease_main_group = paste(
              unique(sort(.data$primary_site)), collapse = "|"),
            atc_code_level1 = paste(
              unique(sort(atc_code_level1)),
              collapse = "|"
            ),
            atc_level1 = paste(
              unique(sort(atc_level1)),
              collapse = "|"
            ),
            atc_code_level2 = paste(
              unique(sort(atc_code_level2)),
              collapse = "|"
            ),
            atc_level2 = paste(
              unique(sort(atc_level2)),
              collapse = "|"
            ),
            atc_code_level3 = paste(
              unique(sort(atc_code_level3)),
              collapse = "|"
            ),
            atc_level3 = paste(
              unique(sort(atc_level3)),
              collapse = "|"
            ), .groups = "drop")
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
                         .data$drug_approved_indication,
                         .data$drug_clinical_source)) |>
        dplyr::distinct()

      drug_records <- as.data.frame(
        drug_records |>
          dplyr::group_by(
            dplyr::across(
              -c(.data$disease_efo_label,
                 .data$atc_code_level1,
                 .data$atc_level1,
                 .data$atc_code_level2,
                 .data$atc_level2,
                 .data$atc_code_level3,
                 .data$atc_level3,
                 .data$primary_site,
                 .data$drug_clinical_id,
                 .data$drug_max_phase_indication))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            disease_indication = paste(
              unique(sort(.data$disease_efo_label)), collapse = "|"),
            disease_indication_max_phase = paste(
              sort(unique(.data$drug_max_phase_indication),
                          decreasing = F), collapse = "|"),
            disease_main_group = paste(
              unique(sort(.data$primary_site)), collapse = "|"),
            atc_code_level1 = paste(
              unique(sort(.data$atc_code_level1)),
              collapse = "|"
            ),
            atc_level1 = paste(
              unique(sort(.data$atc_level1)),
              collapse = "|"
            ),
            atc_code_level2 = paste(
              unique(sort(.data$atc_code_level2)),
              collapse = "|"
            ),
            atc_level2 = paste(
              unique(sort(.data$atc_level2)),
              collapse = "|"
            ),
            atc_code_level3 = paste(
              unique(sort(.data$atc_code_level3)),
              collapse = "|"
            ),
            atc_level3 = paste(
              unique(sort(.data$atc_level3)),
              collapse = "|"
            ), .groups = "drop")
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
                 .data$atc_code_level1,
                 .data$atc_level1,
                 .data$atc_code_level2,
                 .data$atc_level2,
                 .data$atc_code_level3,
                 .data$atc_level3,
                 .data$drug_clinical_source))) |>
          dplyr::summarise(
            drug_clinical_id = paste(
              unique(sort(.data$drug_clinical_id)), collapse = "|"),
            drug_clinical_source = paste(
              unique(sort(.data$drug_clinical_source)), collapse = "|"),
            atc_code_level1 = paste(
              unique(sort(.data$atc_code_level1)),
              collapse = "|"
            ),
            atc_level1 = paste(
              unique(sort(.data$atc_level1)),
              collapse = "|"
            ),
            atc_code_level2 = paste(
              unique(sort(.data$atc_code_level2)),
              collapse = "|"
            ),
            atc_level2 = paste(
              unique(sort(.data$atc_level2)),
              collapse = "|"
            ),
            atc_code_level3 = paste(
              unique(sort(.data$atc_code_level3)),
              collapse = "|"
            ),
            atc_level3 = paste(
              unique(sort(.data$atc_level3)),
              collapse = "|"
            ),
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
        dplyr::desc(.data$drug_max_ct_phase),
        dplyr::desc(.data$drug_year_first_approval),
        dplyr::desc(.data$opentargets),
        nchar(.data$drug_name))
    
    if(!is.null(drug_indication_main)){
      
      if(output_resolution != "drug2target2indication"){
        drug_records <- drug_records |>
          dplyr::arrange(
            dplyr::desc(.data$disease_indication_max_phase),
            dplyr::desc(.data$drug_year_first_approval),
            dplyr::desc(.data$opentargets),
            nchar(.data$drug_name))
      }
    }

  }

  oncodrugs <- list()
  oncodrugs[['records']] <- drug_records
  oncodrugs[['metadata']] <- metadata

  return(oncodrugs)

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
