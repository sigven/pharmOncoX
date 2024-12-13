

## Read URL
#'
#' Function that retrieves PubChem compound properties using PubChem public user gateway (PUG)
#'
#' @param q_url Query URL
#' @return a data frame with the following columns of chemical compound properties:
#'
#'
readUrl <- function(q_url) {
  out <- tryCatch(
    {
      readLines(con=q_url, warn=FALSE)
    },
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      #message(paste("Processed URL:", q_url))
      #message("Some other message at the end")
    }
  )
  return(out)
}


### NCI DRUG DISPLAY LABELS

process_nci_labels <- function(path_data_raw, overwrite = F) {

  nci_thesaurus_labels_filepath = file.path(path_data_raw, "nci_thesaurus", "ThesaurusLabels.tsv")
  nci_thesaurus_owl_filepath <- file.path(path_data_raw, "nci_thesaurus", "ThesaurusInferred.owl")
  if(!is.null(nci_thesaurus_owl_filepath) & 
     file.exists(nci_thesaurus_owl_filepath) & 
     (!file.exists(nci_thesaurus_labels_filepath) | overwrite == T)){
    system(paste0('/Users/sigven/miniconda3/bin/python data-raw/parse_owl_labels.py ',nci_thesaurus_owl_filepath,' > ',file.path(path_data_raw,'nci_thesaurus','ThesaurusLabels.tsv')))
  }
  nci_labels <- read.table(file = nci_thesaurus_labels_filepath, sep="\t",comment.char="",
                           stringsAsFactors = F,quote="",header=F) |>
    magrittr::set_colnames(c('nci_t','nci_cd_name','cui')) |>
    dplyr::filter(!is.na(nci_t))

  return(nci_labels)
}

### CHEMBL-PUBCHEM COMPOUND CROSS-REFERENCE

get_chembl_pubchem_xref <- function(datestamp = '20241024',
                                             chembl_release = "v34",
                                             path_data_raw = NULL,
                                             update = F){
  chembl_pubchem_xref_fname <- file.path(
    path_data_raw, "chembl",
    paste0("chembl_pubchem_mapping_", datestamp,".txt.gz"))
  if(!file.exists(chembl_pubchem_xref_fname) | update == T){
    download.file(
      "ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz",
      destfile = chembl_pubchem_xref_fname)
  }
  chembl_pubchem_compound_xref <- read.table(
    gzfile(chembl_pubchem_xref_fname), header = F, skip = 1,
    stringsAsFactors = F,sep="\t",quote="") |>
    magrittr::set_colnames(c('molecule_chembl_id','pubchem_cid')) |>
    dplyr::mutate(chembl_db_version = chembl_release)

  return(chembl_pubchem_compound_xref)
}



### TARGETED ANTICANCER COMPOUNDS FROM OPEN TARGETS
get_otp_cancer_drugs <-
  function(path_data_raw = NULL,
           ot_version = "2024.09"){

    cancer_terms <- list()
    cancer_terms[['all']] <- phenOncoX::get_terms(
      cache_dir = path_data_raw
    )
    
    cancer_terms[['efo']] <- cancer_terms[['all']]$records |>
      dplyr::select(cui, cui_name, efo_id, primary_site) |> 
      dplyr::filter(!is.na(efo_id)) |>
      dplyr::distinct() |>
      dplyr::rename(disease_efo_id = efo_id)

  fname <- paste0(path_data_raw,
                  paste0("/opentargets/opentargets_drugs_",
                         ot_version,".rds"))
  ot_compounds <- as.data.frame(
    readRDS(file = fname) |>
      dplyr::select(target_genename,
                    target_symbol,
                    target_type,
                    target_ensembl_gene_id,
                    target_entrezgene,
                    disease_efo_id,
                    disease_efo_label,
                    drug_name,
                    drug_type,
                    drug_max_ct_phase,
                    drug_max_phase_indication,
                    molecule_chembl_id,
                    drug_action_type,
                    drug_moa,
                    drug_tradenames,
                    drug_synonyms,
                    drug_clinical_source,
                    drug_clinical_id,
                    drug_description,
                    drug_withdrawn,
                    drug_approved_indication,
                    drug_blackbox_warning,
                    drug_year_first_approval)) |>
    dplyr::mutate(
      disease_efo_id = stringr::str_replace_all(
        disease_efo_id, "_", ":")) |>
    dplyr::distinct()  |>
    dplyr::left_join(
      cancer_terms[['efo']], 
      by = "disease_efo_id", 
      multiple = "all", 
      relationship = "many-to-many") |>
    ## general cancer ontology terms (neoplasm, cancer, carcinoma, squamous cell carcinoma)
    dplyr::mutate(
      drug_cancer_relevance =
        dplyr::if_else(
          !is.na(primary_site) |
            disease_efo_id == 'EFO:0000616' |
            disease_efo_id == 'EFO:0000311' |
            disease_efo_id == 'EFO:0000313' |
            disease_efo_id == 'MONDO:0004992' |
            disease_efo_id == 'MONDO:0002120' |
            disease_efo_id == 'EFO:0003826' |
            disease_efo_id == 'EFO:0000707' |
            disease_efo_id == 'EFO:0000228',
          "by_cancer_condition_otp", "by_other_condition_otp")
    ) |>
    dplyr::mutate(
      drug_cancer_relevance = dplyr::if_else(
        is.na(disease_efo_id) &
          (drug_type == "Antibody" | 
             drug_type == "Unknown" |
             drug_type == "Small molecule") &
          !is.na(target_genename) &
          stringr::str_detect(
            tolower(target_genename),
            paste0(
              "(kinase|epidermal growth factor|tubulin|",
              "oncogene|fibroblast growth factor|",
              "dna topoisomerase|tnf receptor superfamily|",
              "protein tyrosine phosphatase|",
              "receptor tyrosine kinase|",
              "platelet derived growth factor)"
              )) &
          (drug_action_type == "INHIBITOR" | 
             drug_action_type == "BINDING AGENT" |
             drug_action_type == "ANTAGONIST"),
        "by_cancer_target_otp",
        as.character(drug_cancer_relevance)
      )
    ) |>
    dplyr::mutate(
      drug_type = dplyr::if_else(
        drug_type == "unknown",
        "Unknown",
        as.character(drug_type)
      )
    ) |>
    # dplyr::filter(
    #   !is.na(drug_action_type)
    # ) |>
    # dplyr::filter(
    #   !stringr::str_detect(
    #     tolower(drug_action_type),
    #     "allosteric|modulator|substrate|releasing|opener"
    #   )
    # ) |>
    dplyr::filter(
      drug_type != "Cell" &
        drug_type != "Enzyme" &
        drug_type != "Gene" &
        drug_type != "Oligonucleotide" &
        drug_type != "Oligosaccharide" 
    )
  
  cancer_indication_stats <- ot_compounds |> 
    dplyr::filter(
      !is.na(drug_cancer_relevance) & 
        drug_cancer_relevance == "by_cancer_condition_otp") |> 
    dplyr::select(drug_name, disease_efo_label) |> 
    dplyr::filter(!is.na(disease_efo_label)) |>
    dplyr::distinct() |> 
    dplyr::group_by(drug_name) |> 
    dplyr::summarise(n_cancer_indications = dplyr::n())
  
  non_cancer_indication_stats <- ot_compounds |>
    dplyr::filter(
      !is.na(drug_cancer_relevance) & 
        drug_cancer_relevance == "by_other_condition_otp") |>
    #dplyr::filter(drug_cancer_relevance == F) |> 
    dplyr::filter(!is.na(disease_efo_label)) |>
    dplyr::select(drug_name, disease_efo_label,
                  drug_approved_indication) |> 
    dplyr::distinct() |> 
    dplyr::group_by(drug_name) |> 
    dplyr::summarise(n_noncancer_indications = dplyr::n(),
                     approved_non_cancer = paste(
                       unique(drug_approved_indication),
                       collapse = ";"
                     ), .groups = "drop") |>
    dplyr::mutate(approved_non_cancer = dplyr::if_else(
      stringr::str_detect(approved_non_cancer,"TRUE"),
      as.logical(TRUE),
      as.logical(FALSE)
    ))
  
  indication_stats <- cancer_indication_stats |>
    dplyr::full_join(non_cancer_indication_stats, by = "drug_name") |>
    dplyr::mutate(n_cancer_indications = dplyr::if_else(
      is.na(n_cancer_indications),
      as.numeric(0),
      as.numeric(n_cancer_indications)
    )) |>
    dplyr::mutate(n_noncancer_indications = dplyr::if_else(
      is.na(n_noncancer_indications),
      as.numeric(0),
      as.numeric(n_noncancer_indications)
    )) |>
    dplyr::mutate(drug_n_indications = n_cancer_indications +
                    n_noncancer_indications) |>
    dplyr::mutate(drug_frac_cancer_indications = dplyr::if_else(
      n_cancer_indications == 0,
      as.numeric(0),
      round(as.numeric(n_cancer_indications / 
                   (n_cancer_indications + n_noncancer_indications)),
            digits = 3)
    )) |>
    dplyr::mutate(drug_approved_noncancer = dplyr::if_else(
      is.na(approved_non_cancer),
      FALSE,
      as.logical(approved_non_cancer)
    )) |>
    dplyr::select(drug_name,
                  drug_n_indications,
                  drug_approved_noncancer,
                  drug_frac_cancer_indications)
  
  ot_compounds <- ot_compounds |>
    dplyr::left_join(indication_stats, by = "drug_name") |>
    dplyr::mutate(
      drug_n_indications = dplyr::if_else(
        is.na(drug_n_indications) | 
          drug_n_indications == "NA",
        0,
        as.integer(drug_n_indications)
      )
    )

  targeted_compounds <- as.data.frame(
    ot_compounds |>
      dplyr::distinct() |>
      dplyr::group_by_at(
        dplyr::vars(-c(drug_clinical_id))) |>
      dplyr::summarise(
        drug_clinical_id = paste(
          unique(drug_clinical_id),
          collapse = ","),
        drug_clinical_source = paste(
          unique(sort(drug_clinical_source)),
          collapse = ","),
        .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(
        drug_name_lc = tolower(drug_name))
  )
  
  
  ## adjust max ct phase
  drugs_with_max_phase_adj <- as.data.frame(
    targeted_compounds |>
    #dplyr::filter(drug_cancer_relevance == T) |>
    dplyr::filter(!is.na(drug_max_ct_phase) & 
                    !is.na(drug_max_phase_indication)) |>
    dplyr::group_by(drug_name, molecule_chembl_id) |>
    dplyr::summarise(drug_max_ct_phase = max(drug_max_phase_indication),
                     .groups = "drop")
  )

  targeted_compounds <- as.data.frame(targeted_compounds |>
    dplyr::select(-drug_max_ct_phase) |>
    dplyr::left_join(
      drugs_with_max_phase_adj, 
      by = c("drug_name","molecule_chembl_id"),
      multiple = "all", relationship = "many-to-many")
  )
  
  targeted_cancer_compounds <- targeted_compounds |>
    dplyr::filter(
      !is.na(drug_cancer_relevance) &
        (drug_cancer_relevance == "by_cancer_condition_otp" |
           drug_cancer_relevance == "by_cancer_target_otp")) |>
    dplyr::mutate(drug_n_indications = dplyr::if_else(
      is.na(drug_n_indications),as.numeric(0),
      as.numeric(drug_n_indications)
    )) |>
    dplyr::mutate(drug_frac_cancer_indications = dplyr::if_else(
      is.na(drug_frac_cancer_indications),as.numeric(0),
      as.numeric(drug_frac_cancer_indications)
    )) |>
    dplyr::mutate(drug_approved_noncancer = dplyr::if_else(
      is.na(drug_approved_noncancer),as.logical(FALSE),
      as.logical(drug_approved_noncancer)
    )) |>
    dplyr::filter(drug_withdrawn == F)
  
  targeted_noncancer_compounds <- targeted_compounds |>
    # dplyr::filter(
    #   !is.na(drug_cancer_relevance) & 
    #     drug_cancer_relevance == "by_other_condition_otp") |> 
    dplyr::anti_join(
      dplyr::select(
        targeted_cancer_compounds, molecule_chembl_id
      ), by = "molecule_chembl_id"
    ) |>
    dplyr::filter(drug_withdrawn == F)
  
  targeted_noncancer_compounds$drug_name_lc <- NULL
  targeted_cancer_compounds$drug_name_lc <- NULL
  
  ot_recs <- dplyr::bind_rows(
    targeted_cancer_compounds,
    targeted_noncancer_compounds) |> 
    dplyr::distinct()
  
  return(ot_recs)

}

get_atc_drug_classification <- function(
    path_data_raw = NULL){
  
  atc_codes <- as.data.frame(
    readr::read_csv(
      file = file.path(
      path_data_raw, "atc", 
      "ATC.csv.gz"), 
      show_col_types = F)) |>
    janitor::clean_names() |>
    dplyr::mutate(class_id = stringr::str_replace(
      class_id, "http://purl.bioontology.org/ontology/ATC/",""
    )) |>
    dplyr::mutate(parents = stringr::str_replace(
      parents, "http://purl.bioontology.org/ontology/ATC/",""
    )) |>
    dplyr::mutate(
      is_drug_class = dplyr::if_else(
        is.na(is_drug_class),
        FALSE,
        TRUE
      )
    ) |>
    dplyr::select(class_id, preferred_label, synonyms, obsolete,
                  parents, atc_level, is_drug_class) |>
    dplyr::arrange(atc_level, class_id, parents)
    #dplyr::filter(stringr::str_detect(class_id,"^L"))
  
  atc_level1 <- atc_codes |>
    dplyr::filter(atc_level == 2) |>
    dplyr::rename(atc_level1 = preferred_label,
                  atc_code_level1 = class_id) |>
    dplyr::select(atc_code_level1, atc_level1)
  
  atc_level2 <- atc_codes |>
    dplyr::filter(atc_level == 3) |>
    dplyr::rename(atc_level2 = preferred_label,
                  atc_code_level2 = class_id, parent_code = parents) |>
    dplyr::select(atc_code_level2, atc_level2, parent_code)
  
  
  atc_level3 <- atc_codes |>
    dplyr::filter(atc_level == 4) |>
    dplyr::rename(atc_level3 = preferred_label,
                  atc_code_level3 = class_id, 
                  parent_code = parents) |>
    dplyr::select(atc_code_level3, 
                  atc_level3, parent_code)
  
  atc_level4 <- atc_codes |>
    dplyr::filter(atc_level == 5) |>
    dplyr::rename(atc_level4 = preferred_label,
                  atc_code_level4 = class_id, 
                  parent_code = parents) |>
    dplyr::select(atc_code_level4, 
                  atc_level4, 
                  parent_code)
  
  atc_drug_classification <- atc_level1 |> 
    dplyr::left_join(
      atc_level2, by = c("atc_code_level1" = "parent_code")) |>
    dplyr::left_join(
      atc_level3, by = c("atc_code_level2" = "parent_code")
    )  |>
    dplyr::left_join(
      atc_level4, by = c("atc_code_level3" = "parent_code")
    )  |>
    dplyr::rename(atc_drug_entry = atc_level4) |>
    dplyr::mutate(atc_level2 = stringr::str_replace(
      atc_level2, 
      "in ATC|, ANTINEOPLASTIC DRUGS|, ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS",""
    )) |>
    dplyr::mutate(atc_level3 = stringr::str_replace(
      atc_level3, 
      "in ATC",""
    )) |>
    dplyr::mutate(atc_level2 = stringr::str_replace(
      atc_level2, "AND ANTIBODY DRUG CONJUGATES","AND ADCs"
    )) |>
    dplyr::mutate(atc_level3 = dplyr::case_when(
      atc_code_level3 == "L03AX" ~ "Other immunostimulants",
      atc_code_level3 == "L01CD" ~ "Taxanes",
      atc_code_level3 == "L04AX" ~ "Other immunosuppressants",
      atc_code_level3 == "L04AB" ~ "TNF-alpha inhibitors",
      atc_code_level3 == "L02BX" ~ "Other hormone antagonists",
      atc_code_level3 == "L02AE" ~ "GnRH analogs",
      atc_code_level3 == "L01XK" ~ "PARP inhibitors",
      atc_code_level3 == "L01XX" ~ "Other antineoplastic agents",
      atc_code_level3 == "L02AX" ~ "Other hormones for endocrine therapy",
      atc_code_level3 == "L02BB" ~ "Anti-androgens",
      atc_code_level3 == "L01XH" ~ "HDAC inhibitors",
      atc_code_level3 == "L01XG" ~ "Proteasome inhibitors",
      atc_code_level3 == "L01XA" ~ "Platinum compounds",
      atc_code_level3 == "L01FX" ~ "Other MaBs and ADCs",
      atc_code_level3 == "L01FF" ~ "PD-1/PDL-1 inhibitors",
      atc_code_level3 == "L01FG" ~ "VEGFR inhibitors",
      atc_code_level3 == "L01FE" ~ "EGFR inhibitors",
      atc_code_level3 == "L01FD" ~ "HER2 inhibitors",
      atc_code_level3 == "L01FC" ~ "CD38 inhibitors",
      atc_code_level3 == "L01FB" ~ "CD22 inhibitors",
      atc_code_level3 == "L01FA" ~ "CD20 inhibitors",
      atc_code_level3 == "L01EX" ~ "Other protein kinase inhibitors",
      atc_code_level3 == "L01EN" ~ "FGFR inhibitors",
      atc_code_level3 == "L01EM" ~ "PI3K inhibitors",
      atc_code_level3 == "L01EL" ~ "BTK inhibitors",
      atc_code_level3 == "L01EK" ~ "VEGFR inhibitors",
      atc_code_level3 == "L01EJ" ~ "JAK inhibitors",
      atc_code_level3 == "L01EH" ~ "HER2 inhibitors",
      atc_code_level3 == "L01EG" ~ "mTOR inhibitors",
      atc_code_level3 == "L01EF" ~ "CDK inhibitors",
      atc_code_level3 == "L01EE" ~ "MEK inhibitors",
      atc_code_level3 == "L01ED" ~ "ALK inhibitors",
      atc_code_level3 == "L01EC" ~ "BRAF inhibitors",
      atc_code_level3 == "L01EB" ~ "EGFR inhibitors",
      atc_code_level3 == "L01EA" ~ "BCR-ABL inhibitors",
      atc_code_level3 == "L01DC" ~ "Other cytotoxic antibiotics",
      atc_code_level3 == "L01CX" ~ "Other plant alkaloids and natural products",
      atc_code_level3 == "L01CE" ~ "TOP1 inhibitors",
      atc_code_level3 == "L01CC" ~ "Colchicine derivatives",
      atc_code_level3 == "L01CB" ~ "Podophyllotoxin derivatives",
      atc_code_level3 == "L01CA" ~ "Vinca alkaloids and analogues",
      atc_code_level3 == "L01BC" ~ "Pyrimidine analogues",
      atc_code_level3 == "L01BB" ~ "Purine analogues",
      atc_code_level3 == "L01BA" ~ "Folic acid analogues",
      atc_code_level3 == "L01AX" ~ "Other alkylating agents",
      atc_code_level3 == "L01XB" ~ "Methylhydrazines",
      atc_code_level3 == "L01AG" ~ "Epoxides",
      atc_code_level3 == "L01AD" ~ "Nitrosoureas",
      atc_code_level3 == "L01AC" ~ "Ethylene imines",
      TRUE ~ as.character(atc_level3)
    )) |>
    dplyr::select(-atc_code_level4) |>
    dplyr::distinct() 
    #dplyr::mutate(atc_drug_entry = dplyr::if_else(
      #atc_code_level3 == "L01EX" |
        #atc_code_level3 == "L01FX" |
        #atc_code_level3 == "L01XX",
      #as.character(NA),
      #as.character(atc_drug_entry)
    #))
  
  atc_custom <- as.data.frame(
    readr::read_csv(
      file.path(
        path_data_raw,
        "atc",
        "custom_classification.csv"
      ), col_names = T, show_col_types = F, na = c(".")
    )
  )
  
  atc_drug_classification <- atc_drug_classification |>
    dplyr::bind_rows(atc_custom) |>
    dplyr::arrange(atc_code_level3) |>
    dplyr::distinct() |>
    dplyr::mutate(atc_code_level3 = dplyr::case_when(
      atc_drug_entry == "ipilimumab" ~ "L01FXA",
      atc_drug_entry == "tremelimumab" ~ "L01FXA",
      TRUE ~ as.character(atc_code_level3)
    )) |>
    dplyr::mutate(atc_level3 = dplyr::case_when(
      atc_drug_entry == "ipilimumab" ~ "Other ICIs - CTLA4 inhibitors",
      atc_drug_entry == "tremelimumab" ~ "Other ICIs - CTLA4 inhibitors",
      TRUE ~ as.character(atc_level3)
    ))
  
  
  return(atc_drug_classification)
  
  
}



## NCI THESAURUS CANCER DRUGS/TREATMENTS
get_nci_drugs <- function(nci_db_release = nci_db_release,
                          overwrite = F,
                          path_data_raw = NULL,
                          path_data_processed = NULL){
  nci_antineo_thesaurus <- NULL
  nci_drugs <- NULL
  
  nci_ftp_base <- paste0(
    "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/",
    nci_db_release,
    "_Release/")
  
  nci_thesaurus_files <- list()
  nci_thesaurus_files[['flat']] <- 
    paste0("Thesaurus_", nci_db_release,".FLAT.zip")
  nci_thesaurus_files[['owl']] <- 
    paste0("Thesaurus_", nci_db_release,".OWL.zip")
  nci_thesaurus_files[['inf_owl']] <- 
    paste0("ThesaurusInf_", nci_db_release,".OWL.zip")
  
  
  options(timeout = 50000)
  for (elem in c('flat','owl','inf_owl')) {
    remote_file <- paste0(
      nci_ftp_base, nci_thesaurus_files[[elem]])
    local_file <- file.path(
      path_data_raw,"nci_thesaurus", 
      nci_thesaurus_files[[elem]])
    if (!file.exists(local_file)) {
      download.file(
        url = remote_file, 
        destfile = local_file, 
        quiet = T)
      system(paste0('unzip -d ',
                    file.path(path_data_raw, "nci_thesaurus"), 
                    ' -o -u ',local_file))
    }
  }
  
  antineo_agents_url <-
    'https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Drug_or_Substance/Antineoplastic_Agent.txt'
  antineo_agents_local <-
    file.path(
      path_data_raw,
      "nci_thesaurus",
      "Antineoplastic_Agent.txt")
  if (!file.exists(antineo_agents_local)) {
    download.file(
      url = antineo_agents_url, 
      destfile = antineo_agents_local, quiet = T)
  }
  
  
  if(!file.exists(
    file.path(
      path_data_processed,
      "nci_thesaurus",
      "nci_treatment_thesaurus_antineo_v2.rds")) | overwrite == T){
    
    sorafenib_definition <-
      'A synthetic compound targeting growth signaling and angiogenesis. Sorafenib blocks the enzyme RAF kinase, a critical component of the RAF/MEK/ERK signaling pathway that controls cell division and proliferation; in addition, sorafenib inhibits the VEGFR-2/PDGFR-beta signaling cascade, thereby blocking tumor angiogenesis.'
    
    nci_display_labels <- process_nci_labels(
      path_data_raw = path_data_raw,
      overwrite = overwrite)
    
    ## Agents/compounds marked as antineplastic according to NCI
    nci_antineo_agents <-
      read.table(file = file.path(path_data_raw,"nci_thesaurus","Antineoplastic_Agent.txt"),
                 header = T, stringsAsFactors = F, comment.char = "",
                 quote = "", sep = "\t") |>
      janitor::clean_names() |>
      dplyr::rename(nci_t = code) |>
      dplyr::select(nci_t) |>
      dplyr::mutate(antineoplastic_agent = TRUE) |>
      dplyr::distinct()
    
    
    ## parse all entries in nci thesaurus where the 
    ## semantic concept type is treatment-related
    nci_antineo_thesaurus_raw <- as.data.frame(
      read.table(file = file.path(path_data_raw, "nci_thesaurus", "Thesaurus.txt"),
                 header = F, stringsAsFactors = F, sep="\t",
                 comment.char="", quote = "") |>
        dplyr::rename(nci_t = V1, nci_concept_name = V2,
                      nci_t_parent = V3, concept_synonym = V4,
                      nci_concept_definition = V5,
                      nci_cd_name = V6, nci_concept_status = V7,
                      nci_concept_semantic_type = V8) |>
        
        dplyr::filter(
          stringr::str_detect(
            nci_concept_semantic_type,
            "Chemical|Substance|Therapeutic|Drug|Immunologic")) |>
        dplyr::filter(
          nci_concept_semantic_type != "Drug Delivery Device" &
          nci_concept_semantic_type != "Chemical Viewed Functionally"  
        ) |>
        dplyr::filter(
          !stringr::str_detect(
            tolower(nci_concept_definition), "surgical"
          )
        ) |>
        dplyr::filter(
          nci_concept_semantic_type != "Therapeutic or Preventive Procedure" |
            (nci_concept_semantic_type == "Therapeutic or Preventive Procedure" &
               stringr::str_detect(
                 tolower(concept_synonym), "/"
               ) &
               !stringr::str_detect(
                 tolower(nci_concept_definition), "regimen"
               )
            )
        ) |>
        dplyr::left_join(
          nci_antineo_agents, by = c("nci_t"),
          multiple = "all", relationship = "many-to-many") |>
        dplyr::select(-c(nci_concept_name, nci_cd_name)) |>
        dplyr::mutate(nci_concept_synonym_all = concept_synonym) |>
        tidyr::separate_rows(concept_synonym,sep="\\|") |>
        dplyr::mutate(nci_concept_synonym = tolower(concept_synonym)) |>
        dplyr::select(-concept_synonym) |>
        dplyr::filter(
          nci_t != 'C147908' & nci_t != 'C71622') |> ## Hormone Therapy Agent
        dplyr::mutate(
          nci_concept_semantic_type =
            dplyr::if_else(
              nci_t == 'C61948','Pharmacologic Substance',
              nci_concept_semantic_type)) |> ## redundant Sorafenib entries
        dplyr::mutate(
          nci_concept_definition =
            dplyr::if_else(
              nci_t == 'C61948', sorafenib_definition,
              nci_concept_definition)) |> ## redundant Sorafenib entries
        dplyr::distinct() |>
        dplyr::left_join(
          nci_display_labels, by = c("nci_t"), 
          multiple = "all", relationship = "many-to-many") |>
        dplyr::filter(
          !(nci_t == "C1806" & nci_concept_synonym == "gemtuzumab")) |>
        dplyr::filter(
          !(nci_t == "C405" & nci_concept_synonym == "ctx")) |>
        dplyr::mutate(nci_db_version = nci_db_release) |>
        dplyr::filter(!is.na(nci_cd_name)) |>
        dplyr::filter(!stringr::str_detect(
          tolower(nci_concept_definition), "coronavirus")) |>
        dplyr::filter(!stringr::str_detect(
          nci_concept_synonym_all, 
          "SARS-CoV-2|COVID-19|CoV-19|Coronary|Corona|Covid-19|covid-19")) |>
        #))
        dplyr::mutate(antineoplastic_agent = dplyr::if_else(
          is.na(antineoplastic_agent),
          as.logical(TRUE),
          as.logical(antineoplastic_agent)
        )) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name," (Gel|Oil|Cream|Seed|Block|Field|Supplement|Factor)$")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"(Vaccination|Lotion|Therapeutic Heat|Procedure|Rehabilitation|Prevention|Rinse)$")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"(Epitope|Exract|Influenza|Ginseng|Ointment|Management|Injection|Tool)$")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Vitamin A Compound|Inactivated Poliovirus|Antineoplastic Immune Cell|Topical)")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Sheng-Yu|Ginseng|Dry Cleaning|Boost|Tobacco|Microwave)")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Blood|Interruption of|Gum Arabic|Vaginal Cylinder|Laser Ablation|Wheatgrass)")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Chemical Challenge|Prevention of|Magic Mouthwash|Wood Dust|Soot|Cocaine)")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Antibody|Antigen|Antioxidant|Vaccination|Acetate|Antiserum|Asbestos|Aspirate|Autoantigen|Cytokine)$")) |>
        dplyr::filter(!stringr::str_detect(nci_cd_name," Spray| Extract| Antidiabetic| Implant|(Green Tea|Living Healthy|Pollutant|Probe|Protective Agent|Supportive Care|Caffe)"))
        #dplyr::filter(
        #  !stringr::str_detect(
        #    tolower(nci_concept_definition),
        #    "chinese |antidiabet|diabetes|antidepress|analgesic|pulmonary edema|nutritional|human carcinogen|anesthetic|nonsedating|sedative|antihyper|antiinflamma|antiarrythm|antiangin|antihist|muscle|neurotransmitter"))
    )
    
    
    nci_antineo_thesaurus <- nci_antineo_thesaurus_raw
    
    ## add apatinib as an alias (is missing in NCI thesaurus)
    apatinib_alias_entry <- nci_antineo_thesaurus |>
      dplyr::filter(nci_concept_synonym == "rivoceranib") |>
      dplyr::mutate(nci_concept_synonym = "apatinib") |>
      dplyr::mutate(nci_concept_synonym_all =
                      paste(nci_concept_synonym_all,"Apatinib",sep="|"))
    
    
    ## add nab-paclitaxel as an alias (is missing in NCI thesaurus)
    paclitaxel_entry <- nci_antineo_thesaurus |>
      dplyr::filter(nci_cd_name == "Paclitaxel") |>
      dplyr::select(nci_t, nci_t_parent, nci_concept_definition,
                    nci_concept_status, nci_concept_semantic_type,
                    nci_concept_synonym_all, cui) |>
      dplyr::distinct()
    
    nab_paclitaxel_entry <- nci_antineo_thesaurus |>
      dplyr::filter(nci_cd_name == "Nab-paclitaxel")
    
    nab_paclitaxel_entry$nci_t <- paclitaxel_entry$nci_t
    nab_paclitaxel_entry$nci_t_parent <- paclitaxel_entry$nci_t_parent
    nab_paclitaxel_entry$cui <- paclitaxel_entry$cui
    nab_paclitaxel_entry$nci_concept_definition <- paclitaxel_entry$nci_concept_definition
    nab_paclitaxel_entry$nci_concept_status <- paclitaxel_entry$nci_concept_status
    nab_paclitaxel_entry$nci_concept_semantic_type <- paclitaxel_entry$nci_concept_semantic_type
    nab_paclitaxel_entry$nci_cd_name <- "Paclitaxel"
    
    aliases_paclitaxel <- stringr::str_split(
      paclitaxel_entry$nci_concept_synonym_all, "\\|")[[1]]
    aliases_nab_paclitaxel <- stringr::str_split(
      nab_paclitaxel_entry$nci_concept_synonym_all[1], "\\|")[[1]]
    all_paclitaxel_aliases <- 
      paste(sort(unique(c(aliases_paclitaxel, aliases_nab_paclitaxel))),
                 collapse="|")
    
    nab_paclitaxel_entry$nci_concept_synonym_all <- 
      all_paclitaxel_aliases
    
    paclitaxel_entries <- nci_antineo_thesaurus |>
      dplyr::filter(nci_cd_name == "Paclitaxel")
    paclitaxel_entries$nci_concept_synonym_all <- 
      all_paclitaxel_aliases
    
    nci_antineo_thesaurus <- nci_antineo_thesaurus |>
      dplyr::filter(nci_cd_name != "Paclitaxel" &
                      nci_cd_name != "Nab-paclitaxel") |>
      dplyr::bind_rows(nab_paclitaxel_entry) |>
      dplyr::bind_rows(paclitaxel_entries) |>
      dplyr::bind_rows(apatinib_alias_entry) |>
      dplyr::arrange(nci_cd_name) |>
      dplyr::rename(nci_drug_name = nci_concept_synonym) |>
      dplyr::filter(!(nci_cd_name == "Sorafenib Tosylate" &
                        nci_drug_name == "sorafenib")) |>
      dplyr::select(-c(nci_t_parent,cui)) |>
      dplyr::distinct()
     
    nci_drugs <- list()
    
    ## NCI anticancer drugs (non-targeted) - lacking compound identifier (CHEMBL)
    nci_drugs[['no_chembl_id']] <- nci_antineo_thesaurus |>
      dplyr::select(nci_t,
                    nci_concept_definition,
                    nci_cd_name,
                    nci_drug_name,
                    nci_concept_synonym_all) |>
      dplyr::filter(!stringr::str_detect(
        nci_drug_name,
        "(material|medication|medicine|medicinal|medical|^was |^agent |no therapy|-line therapy)"
      )) |>
      dplyr::filter(!stringr::str_detect(
        tolower(nci_cd_name),
        "(material|medication|medicine|medicinal|medical|^was |^agent |no therapy|-line therapy)"
      )) |>
      
      ## remove general drug terms ("XX Agent" etc)
      dplyr::mutate(num_words = stringr::str_count(
        nci_cd_name, " ") + 1) |>
      dplyr::filter(
        !(num_words <= 4 & 
            stringr::str_detect(tolower(nci_cd_name), " agent") &
            !stringr::str_detect(nci_cd_name, "[0-9]{1,}")
        )) |>
      dplyr::select(-num_words) |>
      dplyr::distinct()
    
    saveRDS(nci_drugs, file = file.path(
      path_data_processed,  "nci_thesaurus",
      "nci_treatment_thesaurus_antineo_v2.rds"))
    
  }else{
    nci_drugs <- readRDS(
      file = file.path(
        path_data_processed, "nci_thesaurus",
        "nci_treatment_thesaurus_antineo_v2.rds"))
  }
  return(nci_drugs)
}


merge_nci_opentargets <- function(
    drug_sets = NULL,
    path_data_raw = NULL){
  
  
  ## X-ref Open Targets and NCI by drug name
  otp_drugs_all <- drug_sets[['otp']] |>
    dplyr::mutate(drug_name_lc = tolower(drug_name)) |>
    dplyr::left_join(
        dplyr::filter(
          drug_sets[['nci']][['no_chembl_id']],
          !is.na(nci_drug_name)),
      by = c("drug_name_lc" = "nci_drug_name"),
      multiple = "all", relationship = "many-to-many") |>
    dplyr::select(-c("drug_name_lc")) |>
    dplyr::distinct() |>
    dplyr::mutate(opentargets = TRUE)
  
  ## list all drug aliases in targeted drugs from OTP
  all_drug_aliases <- otp_drugs_all |>
    tidyr::separate_rows(drug_synonyms, sep = "\\|") |>
    dplyr::select(drug_synonyms) |>
    dplyr::mutate(nci_drug_name = tolower(drug_synonyms)) |>
    dplyr::distinct()
  
  ## Identify drugs in NCI that are in Open Targets
  found_drugs <- drug_sets[['nci']][['no_chembl_id']] |>
    dplyr::semi_join(all_drug_aliases,
                     by = "nci_drug_name") |>
    dplyr::select(nci_t, nci_drug_name) |>
    dplyr::distinct()
  
  
  
  nci_missing <- drug_sets[['nci']][['no_chembl_id']] |>
    dplyr::anti_join(found_drugs,
                     by = "nci_t") |>
    dplyr::select(-nci_drug_name) |>
    dplyr::mutate(opentargets = F) |>
    dplyr::distinct()
  
  ## ignore some drugs
  custom_name_ignore <- as.data.frame(readr::read_tsv(
    file = "data-raw/drug_names_ignore.tsv",
    col_names = F, show_col_types = F))
  colnames(custom_name_ignore) <- c("nci_cd_name")
  
  ## ignore some molecule identifiers (mislabeled)
  custom_chembl_ignore <- as.data.frame(readr::read_tsv(
    file = "data-raw/molecule_chembl_ids_ignore.tsv",
    col_names = F, show_col_types = F))
  colnames(custom_chembl_ignore) <- c("molecule_chembl_id")
  
  all_cancer_drugs <- otp_drugs_all |>
    dplyr::bind_rows(nci_missing) |>
    dplyr::anti_join(
      custom_chembl_ignore,
      by = "molecule_chembl_id") |>
    dplyr::distinct() |>
    dplyr::mutate(
      nci_cd_name =
        dplyr::if_else(
          is.na(nci_cd_name) &
            !stringr::str_detect(drug_name,"[0-9]"),
          Hmisc::capitalize(tolower(drug_name)),
          nci_cd_name)) |>
    dplyr::mutate(
      nci_cd_name = dplyr::if_else(
        is.na(nci_cd_name) &
          stringr::str_detect(drug_name,"[0-9]"),
        drug_name,
        nci_cd_name)) |>
    dplyr::mutate(
      nci_cd_name = dplyr::if_else(
        !is.na(nci_cd_name) &
          !is.na(molecule_chembl_id) &
          molecule_chembl_id == "CHEMBL535",
      "Sunitinib",
      nci_cd_name)) |>
    dplyr::mutate(
      nci_cd_name = dplyr::if_else(
        !is.na(drug_name) &
        !is.na(nci_cd_name) &
          !is.na(molecule_chembl_id) & 
          (molecule_chembl_id == "CHEMBL3545427" |
          molecule_chembl_id == "CHEMBL5315068" |
          molecule_chembl_id == "CHEMBL2109325" |
          molecule_chembl_id == "CHEMBL5316170" |
          molecule_chembl_id == "CHEMBL31965" |
          molecule_chembl_id == "CHEMBL5314999" |
          molecule_chembl_id == "CHEMBL2109402" |
          molecule_chembl_id == "CHEMBL1790041" |
          molecule_chembl_id == "CHEMBL4650278" |
          molecule_chembl_id == "CHEMBL3545401"),
        stringr::str_to_title(drug_name),
        nci_cd_name)) |>
    dplyr::mutate(drug_name = dplyr::if_else(
      is.na(drug_name),
      nci_cd_name,
      drug_name)) |>
    dplyr::anti_join(
      custom_name_ignore,
      by = "nci_cd_name")
  
  rownames(all_cancer_drugs) <- NULL
  
  salt_patterns <-
    readr::read_tsv(
      file.path(path_data_raw, "salts.tsv"),
      show_col_types = F, col_names = F)

  salt_patterns_regex <- paste0(
    "( (",
    paste(salt_patterns$X1, collapse="|"),
    "))$")
  
  salt_forms <- all_cancer_drugs |>
    dplyr::filter(
      stringr::str_detect(
        tolower(nci_cd_name), 
        salt_patterns_regex)) |>
    dplyr::filter(!is.na(opentargets_version)) |>
    dplyr::mutate(tradename = stringr::str_replace(
      nci_cd_name,
      salt_patterns_regex,
      "")) |>
    dplyr::select(tradename, nci_cd_name) |>
    dplyr::distinct() |>
    dplyr::mutate(is_salt = T) |>
    dplyr::inner_join(
      dplyr::select(all_cancer_drugs, nci_cd_name),
      by = c("tradename" = "nci_cd_name"),
      multiple = "all", relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::select(-tradename)
  
  
  all_cancer_drugs <- all_cancer_drugs |>
    dplyr::left_join(
      salt_forms, 
      by = "nci_cd_name",
      multiple = "all", 
      relationship = "many-to-many") |>
    dplyr::mutate(is_salt = dplyr::if_else(
      is.na(is_salt),
      as.logical(FALSE),
      as.logical(is_salt)
    )) |>
    dplyr::distinct()
  
  
  ## antibody drug conjugates
  adc_candidates <- all_cancer_drugs |>
    dplyr::filter(
      (!is.na(nci_cd_name) &
         stringr::str_detect(tolower(nci_cd_name), "mab ")) |
        stringr::str_detect(
          nci_concept_definition, "ADC|antibody(-| )drug conjugate")) |>
    dplyr::mutate(is_adc = TRUE) |>
    dplyr::mutate(is_adc = dplyr::if_else(
      stringr::str_detect(nci_concept_definition, "ADCC") &
        !(stringr::str_detect(
          tolower(nci_concept_definition),
          "(antibody(-| )drug conjugate)") &
            (is.na(drug_type) |
               drug_type == "Small molecule" |
               drug_type == "Antibody")),
      as.logical(FALSE),
      as.logical(is_adc)
    )) |>
    dplyr::mutate(is_adc = dplyr::if_else(
      is_adc == F &
        !is.na(nci_cd_name) &
        stringr::str_detect(
          tolower(nci_cd_name), "mab ") &
        stringr::str_detect(
          tolower(nci_cd_name), "(tin|ine|tan|can|tox)$"
        ),
      TRUE,
      as.logical(is_adc)
    )) |>
    dplyr::filter(!is.na(nci_cd_name)) |>
    dplyr::select(nci_cd_name, is_adc) |>
    dplyr::distinct()
  
  all_cancer_drugs <- all_cancer_drugs |>
    dplyr::left_join(
      adc_candidates, 
      by = "nci_cd_name",
      multiple = "all", 
      relationship = "many-to-many") |>
    dplyr::mutate(is_adc = dplyr::if_else(
      is.na(is_adc),
      as.logical(FALSE),
      as.logical(is_adc)
    )) |>
    dplyr::filter(drug_withdrawn == F | is.na(drug_withdrawn)) |>
    dplyr::mutate(nci_concept_definition = dplyr::if_else(
      !is.na(nci_concept_definition) & 
        nchar(nci_concept_definition) == 0,
      as.character(NA),
      as.character(nci_concept_definition)
    )) |>
    dplyr::mutate(drug_n_indications = dplyr::if_else(
      is.na(drug_n_indications),as.numeric(0),
      as.numeric(drug_n_indications)
    )) |>
    dplyr::mutate(drug_frac_cancer_indications = dplyr::if_else(
      is.na(drug_frac_cancer_indications),as.numeric(0),
      as.numeric(drug_frac_cancer_indications)
    )) |>
    dplyr::mutate(drug_approved_noncancer = dplyr::if_else(
      is.na(drug_approved_noncancer),as.logical(FALSE),
      as.logical(drug_approved_noncancer)
    )) |>
    dplyr::mutate(drug_cancer_relevance = dplyr::if_else(
      is.na(drug_cancer_relevance) &
        opentargets == F,
      "by_nci",
      as.character(drug_cancer_relevance)
    )) |>
    dplyr::distinct()
  
  return(all_cancer_drugs)
  
}



map_curated_targets <- function(gene_info = NULL,
                                   path_data_raw = NULL,
                                   drug_df = NULL){

  drug_df$target_entrezgene <- as.numeric(
    drug_df$target_entrezgene
  )
  
  drug_target_patterns <-
    read.table(file = file.path(
      path_data_raw,
      "custom_drug_target_regex_nci.tsv"),
      sep = "\t", header = T, 
      stringsAsFactors = F, quote = "") |>
    dplyr::inner_join(
      gene_info, by = "symbol", multiple = "all", 
      relationship = "many-to-many") |>
    dplyr::distinct()

  drugname_suffix <- 
    paste0("(ab|al|an|at|cl|co|da|de|dy|ea|ed|el|en|er|es|ex|fa|ib|ic|",
           "id|il|im|in|ir|is|le|ls|lt|me|na|ne|ns|nt|od|ol|on|or|pt|",
           "py|ra|rd|re|rm|rt|se|ta|te|ts|ue|um|us|yl)$")
  drug_inhibitor_regex <- 
    paste0(drugname_suffix,"|^(anti-|inhibitor of)|inhibitor")
  
  all_inhibitors_no_target <- drug_df |>
    dplyr::filter(is.na(target_symbol)) |>
    dplyr::filter(
      stringr::str_detect(
        tolower(nci_cd_name),
        "inhibitor|antagonist|antibody|blocker|sepantronium| mimetic") |
        stringr::str_detect(
          tolower(nci_cd_name),
          drug_inhibitor_regex) |
        (stringr::str_detect(nci_concept_definition,"KRAS") &
           stringr::str_detect(nci_concept_definition,"inhibitor"))) |>
    dplyr::filter(!stringr::str_detect(
      nci_cd_name,
      " CAR T|(T|t)herapy|SARS-CoV-2|( (R|r)egimen|(A|a)cid|Hydrochloride)$")) |>
    dplyr::filter(!stringr::str_detect(
      nci_concept_definition,
      "SARS-CoV-2")) |>
    dplyr::filter(nchar(nci_concept_definition) > 0)

  custom_nci_targeted_drugs <- data.frame()
  for(i in 1:nrow(drug_target_patterns)){
    pattern <- drug_target_patterns[i, "pattern"]
    target_symbol <- drug_target_patterns[i, "symbol"]
    target_genename <- drug_target_patterns[i, "genename"]
    target_entrezgene <- drug_target_patterns[i, "target_entrezgene"]
    target_type <- drug_target_patterns[i, "target_type"]
    target_ensembl_gene_id <- drug_target_patterns[i, "target_ensembl_gene_id"]
    #target_uniprot_id <- drug_target_patterns[i, "target_uniprot_id"]

    hits <- all_inhibitors_no_target |>
      dplyr::filter(stringr::str_detect(
        nci_cd_name,
        pattern = pattern) |
          (stringr::str_detect(
            tolower(nci_cd_name),
            drug_inhibitor_regex) &
             stringr::str_detect(
               nci_concept_definition, pattern))
      )

    if(nrow(hits) > 0){

      for(n in 1:nrow(hits)){
        hit <- hits[n,]

        if(stringr::str_detect(tolower(hit$nci_cd_name),
                               "(mab|art)$|monoclonal antibody|^anti-")){
          hit$drug_type <- "Antibody"
        }else{
          hit$drug_type <- "Small molecule"
        }

        hit$drug_action_type <- "INHIBITOR"
        if(stringr::str_detect(
          tolower(hit$nci_cd_name),
          "antagonist")){
          hit$drug_action_type <- "ANTAGONIST"
        }
        if(stringr::str_detect(
          tolower(hit$nci_cd_name),
          "blocker")){
          hit$drug_action_type <- "BLOCKER"
        }
        hit$target_symbol <- target_symbol
        hit$target_genename <- target_genename
        hit$target_type <- target_type
        hit$target_entrezgene <- target_entrezgene
        hit$target_ensembl_gene_id <- target_ensembl_gene_id
        hit$drug_clinical_source <- "nci_thesaurus_custom"
        hit$drug_cancer_relevance <- "by_cancer_target_nci"
        hit$drug_n_indications <- 0
        hit$drug_frac_cancer_indications <- 0
        hit$drug_approved_noncancer <- FALSE

        ## set general indications for unknown cases
        if(is.na(hit$disease_efo_id) & 
           is.na(hit$disease_efo_label) &
           is.na(hit$cui) & 
           is.na(hit$cui_name)){
          hit$disease_efo_id = "EFO:0000311"
          hit$disease_efo_label = "cancer"
          hit$cui = "C0006826"
          hit$cui_name = "Malignant neoplastic disease"
        }

        custom_nci_targeted_drugs <- custom_nci_targeted_drugs |>
          dplyr::bind_rows(hit)

      }
    }
  }


  ### CHECK HOW MANY TARGET-LACKING INHIBITORS ARE MISSING
  ### FROM THE CUSTOM NCI MATCHING ROUTINE

  inhibitors_no_target_nonmapped <- all_inhibitors_no_target |>
    dplyr::anti_join(custom_nci_targeted_drugs, by = "nci_cd_name") |>
    dplyr::filter(!stringr::str_detect(
      nci_concept_definition, 
      "(A|a)ntibody(-| )drug conjugate \\(ADC\\)"
    )) |>
    dplyr::filter(!stringr::str_detect(
      tolower(nci_cd_name),"^(allogeneic|regimen |copper |fluorine f |indium |iodine |carbon c|autologous |recombinant |lutetium |yttrium |y 90)|vaccine$"
    )) |>
    dplyr::filter(
      stringr::str_detect(
        nci_concept_definition, "antineoplastic|tumor|cancer"
      )
    ) |>
    dplyr::select(nci_cd_name,
                  nci_concept_definition) |>
    dplyr::distinct()


  drug_df$target_entrezgene <- 
    as.integer(drug_df$target_entrezgene)
  
  custom_nci_targeted_drugs$target_entrezgene <-
    as.integer(custom_nci_targeted_drugs$target_entrezgene)
  
  ot_nci_drugs_curated <-
    dplyr::anti_join(drug_df, custom_nci_targeted_drugs,
                     by = "nci_cd_name") |>
    dplyr::bind_rows(custom_nci_targeted_drugs) |>
    dplyr::arrange(target_symbol, nci_cd_name) |>
    dplyr::mutate(drug_action_type = dplyr::if_else(
      (stringr::str_detect(tolower(nci_cd_name),"inhibitor") &
         is.na(drug_action_type)) |
        (!is.na(nci_cd_name) &
           stringr::str_detect(nci_cd_name,"(mab|art)$") &
           is.na(drug_action_type)),
      "INHIBITOR",
      as.character(drug_action_type))) |>
    dplyr::mutate(drug_cancer_relevance = dplyr::if_else(
      is.na(drug_cancer_relevance) &
        (stringr::str_detect(
          tolower(nci_concept_definition),
          "anti-tumor|chemotherapy|cancer vaccine|immunothera|monoclonal antibody|antineoplastic|treatment of cancer|treatment of metastat") |
           stringr::str_detect(tolower(nci_cd_name)," regimen|recombinant|carcinoma|immune checkpoint|anti-programmed cell death ")),
      as.character("by_cancer_definition_nci"),
      as.character(drug_cancer_relevance)
    )) |>
    dplyr::mutate(drug_action_type = dplyr::if_else(
      is.na(drug_action_type) &
        stringr::str_detect(
          drug_action_type,
          "^(SUBSTRATE|HYDROLYTIC ENZYME|RELEASING AGENT)"),
      paste0(drug_action_type,"_OTHER"),
      as.character(drug_action_type)
    ))

  return(list('curated' = ot_nci_drugs_curated,
              'nonmapped' = inhibitors_no_target_nonmapped))
}

assign_drug_category <- function(drug_df = NULL,
                                  path_data_raw = NULL){

  atc_classification <- 
    get_atc_drug_classification(path_data_raw = path_data_raw)

  atc_classification_with_drugs <- atc_classification |>
    dplyr::mutate(atc_drug_entry = dplyr::case_when(
      !is.na(atc_drug_entry) &
      (tolower(atc_drug_entry) == "sotorasib" |
         tolower(atc_drug_entry) == "ivosidenib" |
         tolower(atc_drug_entry) == "enasidenib" |
         tolower(atc_drug_entry) == "venetoclax" |
        tolower(atc_drug_entry) == "adagrasib") ~ as.character(NA),
      TRUE ~ as.character(atc_drug_entry)
    )) |>
    dplyr::filter(!is.na(atc_drug_entry))
  
  atc_classification_clean <- atc_classification |>
    dplyr::filter(is.na(atc_drug_entry)) |>
    dplyr::bind_rows(atc_classification_with_drugs) |>
    dplyr::select(-atc_drug_entry) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(atc_code_level3))
  
  classified_drugs <- list()

  classified_drugs[['pre_classified_atc']] <- drug_df |>
    dplyr::mutate(drug_entry = tolower(nci_cd_name)) |>
    dplyr::left_join(
      dplyr::select(
        atc_classification_with_drugs, 
        atc_drug_entry, 
        atc_code_level3),
      by = c("drug_entry" = "atc_drug_entry"),
      relationship = "many-to-many") |>
    dplyr::filter(!is.na(atc_code_level3)) |>
    dplyr::distinct()
  
  drugs_non_classified <- drug_df |>
    dplyr::mutate(drug_entry = tolower(nci_cd_name)) |>
    dplyr::anti_join(
      classified_drugs[['pre_classified_atc']], by = "drug_entry") |>
    dplyr::distinct()

  custom_target_classifications <- drugs_non_classified |>
    dplyr::filter(!is.na(target_symbol)) |>
    dplyr::group_by(drug_name, nci_cd_name, 
                    drug_entry,
                    nci_concept_definition,
                    molecule_chembl_id, drug_action_type) |>
    dplyr::reframe(
      target_symbol = paste(
        sort(unique(target_symbol)), collapse="|")) |>
    dplyr::mutate(atc_code_level3 = dplyr::case_when(
      target_symbol == "ABL1" | target_symbol == "ABL1|BCR" ~ "L01EA",
      target_symbol == "BRAF" ~ "L01EC",
      stringr::str_detect(target_symbol,"TOP1(\\|TOP2)?") ~ "L01CE",
      stringr::str_detect(
        target_symbol, "^((PSM(A|B|C|D)[0-9]{1,})|ADRM1|SEM1)") ~ "L01XG",
      stringr::str_detect(
        target_symbol, "^((BRAF\\|(KDR|RAF))|(ARAF\\|BRAF))") ~ "L01EXJ",
      target_symbol == "PLK1" | target_symbol == "PLK4" ~ "L01EXK",
      target_symbol == "EGFR" &
        !is.na(drug_name) &
          stringr::str_detect(drug_name, "MAB") ~ "L01FE",
      target_symbol == "ERBB2" &
        !is.na(drug_name) &
        stringr::str_detect(drug_name, "MAB") ~ "L01FD",
      target_symbol == "ERBB2" &
        (is.na(drug_name) |
           (!is.na(drug_name) &
              !stringr::str_detect(drug_name, "MAB"))) ~ "L01EH",
      !is.na(drug_name) &
        drug_name == "SIROLIMUS" ~ "L01EG",
      target_symbol == "EGFR" &
        (is.na(drug_name) |
        (!is.na(drug_name) &
        !stringr::str_detect(drug_name, "MAB"))) ~ "L01EB",
      (!is.na(drug_name) &
         stringr::str_detect(drug_name, "ABIVERTINIB")) ~ "L01EB",
      target_symbol == "ALK" ~ "L01ED",
      (drug_action_type == "INHIBITOR" |
         drug_action_type == "ANTAGONIST") & 
        target_symbol == "CTLA4" ~ "L01FXA",
      (drug_action_type == "INHIBITOR" |
         drug_action_type == "ANTAGONIST") & 
        target_symbol == "LAG3"  ~ "L01FXB",
      (drug_action_type == "INHIBITOR" |
         drug_action_type == "ANTAGONIST") & 
           target_symbol == "TIGIT" ~ "L01FXC",
      (drug_action_type == "INHIBITOR" |
         drug_action_type == "ANTAGONIST") & 
           (target_symbol == "CD274|CTLA4" |
           target_symbol == "CTLA4|PDCD1") ~ "L01FXD",
      (!is.na(drug_name) & drug_name == "BRENTUXIMAB VEDOTIN") |
      (stringr::str_detect(
        target_symbol,"^(TUBA|TUBB)") &
        (drug_action_type == "INHIBITOR" |
           drug_action_type == "DISRUPTING AGENT")) ~ "L01XXK",
      stringr::str_detect(
        target_symbol, "^(PARP[1-9]{1}(\\|)?){1,}$") ~ "L01XK",
      stringr::str_detect(
        target_symbol, "^(HDAC[0-9]{1,}(\\|)?){1,}$") ~ "L01XH",
      stringr::str_detect(
        target_symbol, "^MAP2K[0-9]") ~ "L01EE",
      stringr::str_detect(
        target_symbol, "^CDK[0-9]{1,}") ~ "L01EF",
      stringr::str_detect(
        target_symbol, "^(BIRC|XIAP)") ~ "L01XXF",
      stringr::str_detect(
        target_symbol, "^JAK[1-2]") ~ "L01EJ",
      stringr::str_detect(
        target_symbol, "^AURK(A|B|C)") ~ "L01EXB",
      stringr::str_detect(
        target_symbol, "EGFR") &
        stringr::str_detect(
          target_symbol,"ERBB") ~ "L01EXH",
      stringr::str_detect(
        target_symbol, "^FGFR[1-4]{1}") ~ "L01EN",
      stringr::str_detect(
        target_symbol, "^IDH[1-2]{1}") ~ "L01XXD",
      stringr::str_detect(
        target_symbol, "^(K|N|H)RAS") ~ "L01XXC",
      stringr::str_detect(
        target_symbol, "^MET$") ~ "L01EXA",
      stringr::str_detect(
        target_symbol, "^(AKT[0-9](\\|)?){1,}$") |
        (!is.na(drug_name) &
           stringr::str_detect(
             drug_name, "GSK-690693")) ~ "L01EXC",
      stringr::str_detect(
        target_symbol, "^ATR$") ~ "L01XXI",
      stringr::str_detect(
        target_symbol, "^PDGFR(A|B)|^KIT") ~ "L01EXE",
      stringr::str_detect(
        target_symbol, "^(NTRK[0-9]\\|?){1,}$") |
        (!is.na(drug_name) & drug_name == "TALETRECTINIB") ~ "L01EXF",
      stringr::str_detect(
        target_symbol, "^(CHEK(1|2)\\|?){1,2}$") ~ "L01XXJ",
      stringr::str_detect(
        target_symbol, "^BCL2") ~ "L01XXM",
      stringr::str_detect(
        target_symbol, "^GNRH") ~ "L02AE",
      stringr::str_detect(
        target_symbol, "^BRD(T|[1-9]{1})") ~ "L01XXA",
      stringr::str_detect(
        target_symbol, "^MTOR\\|PIK3") ~ "L01XXG",
      stringr::str_detect(
        target_symbol, "^PIK3") |
        (!is.na(drug_name) & 
        (stringr::str_detect(drug_name,"UMBRALISIB") |
        (drug_name == "OMIPALISIB" |
           drug_name == "PF-04691502"))) ~ "L01EM",
      stringr::str_detect(
        target_symbol, "^BTK$") ~ "L01EL",
      (!is.na(drug_name) & 
         (drug_name == "PRALSETINIB" |
            drug_name == "SELPERCATINIB")) ~ "L01EXL",
      stringr::str_detect(
        target_symbol, "^(KDR|FLT1|FLT3|FLT4)") ~ "L01EK",
      stringr::str_detect(
        target_symbol, "^(MS4A1)") ~ "L01FA",
      stringr::str_detect(
        target_symbol, "^(CD38)") ~ "L01FC",
      stringr::str_detect(
        target_symbol, "^(CD22)") ~ "L01FB",
      target_symbol == "MAPK1|MAPK3" ~ "L01XXE",
      stringr::str_detect(target_symbol, "^AR$") & 
        (!is.na(nci_concept_definition) &
           stringr::str_detect(
             nci_concept_definition, "androgen receptor")) ~ "L02BB",
      (target_symbol == "PDCD1" |
         target_symbol == "CD274") ~ "L01FF",
      TRUE ~ as.character(NA)
    )) |>
    dplyr::mutate(atc_code_level3 = dplyr::if_else(
      (!is.na(drug_name) &
         stringr::str_detect(
           drug_name,
           paste0(
             "DOVITINIB|BRIVANIB|UCN-01|AT-9283|",
             "SURUFATINIB|ORANTINIB|LUCITANIB|",
             "ALTIRATINIB|CEP-11981|CRENOLANIB|",
             "FORETINIB|FAMITINIB|MOTESANIB|",
             "LESTAURTINIB|OSI-930|PACRITINIB|KW-2449|",
             "CABOZANTINIB|BMS-817378|BMS-794833|",
             "GOLVATINIB|TAK-593|XL-820|TANDUTINIB|IMATINIB|",
             "CERDULATINIB|CEP-2563|FEDRATINIB|IBCASERTIB|",
             "VANDETANIB|TESEVATINIB|TARLOXOTINIB|",
             "TAK-285|SKLB1028|PUQUITINIB|KBP5209|REPOTRECTINIB|",
             "VATALANIB|QUIZARTINIB|BMS-690514|IMATINIB|DASATINIB|",
             "CANERTINIB|CEP-32496|REGORAFENIB|GUSACITINIB|",
             "SU-014813|X-82|XL-999|LINIFANIB|NINGETINIB|",
             "PEXIDARTINIB|RG-1530|SITRAVATINIB|SORAFENIB|SUNITINIB|",
             "ENTRECTINIB"))) |
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           tolower(nci_concept_definition), 
        paste0(
          "multikinase|multi(-)?targeted|multi-kinase inhibitor|",
          "targets multiple|inhibitor of multiple|multiple-receptor")) &
        stringr::str_count(target_symbol, pattern = "\\|") >= 2),
      "L01EXI",
      as.character(atc_code_level3)
    )) |>
    dplyr::left_join(
      dplyr::select(
        atc_classification_clean,
        atc_code_level3,
      ),
      by = "atc_code_level3",
      relationship = "many-to-many"
    ) |>
    #dplyr::filter(!is.na(atc_code_level3)) |>
    dplyr::select(
      #molecule_chembl_id,
      drug_entry,
      atc_code_level3
    ) |>
    dplyr::group_by(drug_entry) |>
    dplyr::summarise(
      atc_code_level3 = paste(
        unique(atc_code_level3), collapse="|")) |>
    dplyr::distinct() |>
    dplyr::mutate(
      atc_code_level3 = stringr::str_replace(
        atc_code_level3, "^(NA\\|)|(\\|(NA|L01XXK))$",""
      )
    )
  
  classified_drugs[['classified_targeted_custom']] <- 
    drugs_non_classified |>
    dplyr::filter(!is.na(target_symbol)) |>
    dplyr::left_join(
      custom_target_classifications, by = "drug_entry",
      relationship = "many-to-many") |>
    dplyr::distinct()
  
  classified_drugs_all<- 
    classified_drugs[['pre_classified_atc']] |>
    dplyr::bind_rows(
      classified_drugs[['classified_targeted_custom']]) |>
    dplyr::bind_rows(
      dplyr::filter(drugs_non_classified, is.na(target_symbol))) |>
    dplyr::mutate(atc_code_level3 = dplyr::case_when(
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        (!is.na(nci_concept_definition) &
         stringr::str_detect(
           tolower(nci_concept_definition),
           "anthracycline|anthracenedione")) ~ "L01DB",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        stringr::str_detect(tolower(drug_entry),"hydroxyurea|leucovorin|eniluracil") ~ "L01BA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        stringr::str_detect(tolower(drug_entry),"arsenic trioxide") ~ "L01BA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        stringr::str_detect(tolower(drug_entry),
                            "thioguanine|aspacytarabine|aspacytarabine|troxacitabine") ~ "L01BB",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        !stringr::str_detect(drug_entry,"/| ") &
        stringr::str_detect(drug_entry, "xel$") ~ "L01CD",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        !stringr::str_detect(drug_entry,"/| ") &
        stringr::str_detect(drug_entry, "platin$") ~ "L01XA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        !is.na(nci_concept_definition) &
        stringr::str_detect(tolower(nci_concept_definition), "anti-estrogen") ~ "L02BA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        !is.na(nci_concept_definition) &
        stringr::str_detect(tolower(nci_concept_definition), "aromatase inhibitor") ~ "L02BG",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        !is.na(nci_concept_definition) &
        stringr::str_detect(tolower(nci_concept_definition), "nitrogen mustard") ~ "L01AA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        stringr::str_detect(
        tolower(nci_concept_definition), "purine( nucleoside)? analog") ~ "L01BB",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
        stringr::str_detect(
        tolower(nci_concept_definition), "pyrimidine( nucleoside)? analog") ~ "L01BC",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") &
        !is.na(drug_entry) &
        !stringr::str_detect(drug_entry, "/") &
        stringr::str_detect(
          tolower(nci_concept_definition), "vinca alkaloid") ~ "L01CA",
      (is.na(atc_code_level3) | atc_code_level3 == "NA") & 
      ((!is.na(nci_concept_definition) &
          stringr::str_detect(
            tolower(nci_concept_definition),
            "antineoplastic activit|anti-tumor activit"
          )) |
         (!is.na(drug_max_ct_phase) &
            stringr::str_detect(
              drug_entry,"(in|ib|ide|ine|ax|il|an|ate| alfa)$") &
            drug_max_ct_phase >= 2 &
            (!is.na(drug_n_indications) &
               drug_n_indications > 2) &
            (!is.na(drug_frac_cancer_indications) &
               drug_frac_cancer_indications > 0.4))) ~ "L01XX",
      
      
      TRUE ~ as.character(atc_code_level3)
    )) |>
    dplyr::distinct() |>
    dplyr::group_by(dplyr::across(-c("atc_code_level3"))) |>
    dplyr::summarise(atc_code_level3 = paste(unique(atc_code_level3), collapse="|"),
                     .groups = "drop")
  
  
  atc_classified_drugs <- classified_drugs_all |>
    dplyr::filter(!is.na(drug_entry)) |>
    dplyr::select(atc_code_level3, drug_entry) |>
    dplyr::filter(!is.na(atc_code_level3)) |>
    tidyr::separate_rows(atc_code_level3) |>
    dplyr::filter(!is.na(atc_code_level3)) |>
    dplyr::left_join(
      dplyr::select(
        atc_classification_clean,
        atc_code_level1,
        atc_level1,
        atc_code_level2,
        atc_level2,
        atc_code_level3,
        atc_level3
      )
    ) |>
    dplyr::group_by(drug_entry) |>
    dplyr::summarise(
      atc_code_level1 = paste(unique(atc_code_level1), collapse="|"),
      atc_level1 = paste(unique(stringr::str_trim(atc_level1)), collapse="|"),
      atc_code_level2 = paste(unique(atc_code_level2), collapse="|"),
      atc_level2 = paste(unique(stringr::str_trim(atc_level2)), collapse="|"),
      atc_code_level3 = paste(unique(atc_code_level3), collapse="|"),
      atc_level3 = paste(unique(stringr::str_trim(atc_level3)), collapse="|")) |>
    dplyr::mutate(atc_code_level1 = dplyr::if_else(
      !is.na(atc_code_level1),
      stringr::str_replace(
        atc_code_level1,"(^NA\\|)|\\|NA$",""),
      as.character(atc_code_level1)
    )) |>
    dplyr::mutate(atc_level1 = dplyr::if_else(
      !is.na(atc_level1),
      stringr::str_replace(
        atc_level1,"(^NA\\|)|\\|NA$",""),
      as.character(atc_level1)
    )) |>
    dplyr::mutate(atc_code_level2 = dplyr::if_else(
      !is.na(atc_code_level2),
      stringr::str_replace(
        atc_code_level2,"(^NA\\|)|\\|NA$",""),
      as.character(atc_code_level2)
    )) |>
    dplyr::mutate(atc_level2 = dplyr::if_else(
      !is.na(atc_level2),
      stringr::str_replace(
        atc_level2,"(^NA\\|)|\\|NA$",""),
      as.character(atc_level2)
    )) |>
    dplyr::mutate(atc_code_level3 = dplyr::if_else(
      !is.na(atc_code_level3),
      stringr::str_replace(
        atc_code_level3,"(^NA\\|)|\\|NA$",""),
      as.character(atc_code_level3)
    )) |>
    dplyr::mutate(atc_level3 = dplyr::if_else(
      !is.na(atc_level3),
      stringr::str_replace(
        atc_level3,"(^NA\\|)|\\|NA$",""),
      as.character(atc_level3)
    )) |>
    dplyr::mutate(atc_code_level3 = dplyr::if_else(
      !is.na(atc_code_level3),
      stringr::str_replace(
        atc_code_level3,"L01XX\\|L01XXK","L01XXK"),
      as.character(atc_code_level3)
    )) |>
    dplyr::mutate(atc_level3 = dplyr::if_else(
      !is.na(atc_level3),
      stringr::str_replace(
        atc_level3,
        "Other antineoplastic agents\\|Tubulin inhibitors",
        "Tubulin inhibitors"),
      as.character(atc_level3)
    )) |>
    dplyr::mutate(atc_level3 = dplyr::if_else(
      !is.na(atc_level3),
      stringr::str_replace(
        atc_level3,
        "Folic acid analogues|Other immunosuppressants",
        "Folic acid analogues"),
      as.character(atc_level3)
    ))
  

  classified_drugs_all_final <- classified_drugs_all |>
    dplyr::select(-atc_code_level3) |>
    dplyr::left_join(atc_classified_drugs, by = "drug_entry",
                     relationship = "many-to-many") |>
    dplyr::mutate(atc_treatment_category = "unknown") |>
    dplyr::mutate(atc_treatment_category = dplyr::case_when(
      !is.na(atc_code_level2) &
        stringr::str_detect(
          atc_code_level2, "(\\|L02(A|B))|(^L02(A|B))") ~ 
        "cancer_hormone_therapy",
      !is.na(target_symbol) & 
      !is.na(atc_code_level2) &
        stringr::str_detect(
          atc_code_level2, "(\\|L01(E|F))|(^L01(E|F))"
        ) ~ "cancer_targeted_therapy",
      !is.na(target_symbol) &
      !is.na(atc_code_level3) &
        stringr::str_detect(
          atc_code_level3, "(\\|L01XX[A-J])|(^L01XX[A-J])"
        ) ~ "cancer_targeted_therapy",
      !is.na(target_symbol) &
        !is.na(atc_code_level3) &
        stringr::str_detect(
          atc_code_level3, "(\\|L01XXM)|(^L01XXM)"
        ) ~ "cancer_targeted_therapy",
      !is.na(target_symbol) &
        !is.na(atc_code_level3) &
        stringr::str_detect(
          atc_code_level3, "(\\|L04)|(^L04)"
        ) ~ "cancer_immuno_suppressants",
      !is.na(target_symbol) &
        stringr::str_detect(drug_cancer_relevance,"^by_cancer") & 
        (is.na(atc_code_level3) | 
           atc_code_level3 == "NA" | 
           atc_code_level3 == "L01XX")
         ~ "cancer_unclassified",
      !is.na(target_symbol) &
      !is.na(atc_code_level3) &
        stringr::str_detect(
          atc_code_level3, "(\\|L01X[A-L])|(^L01X[A-L])"
        ) ~ "cancer_targeted_therapy",
      !is.na(atc_code_level2) &
        !stringr::str_detect(
          atc_code_level2, "(\\|L01(E|F))|(^L01(E|F))") &
        
        ((!is.na(atc_code_level3) &
            drug_type != "Antibody" &
        stringr::str_detect(atc_code_level3, "L01XXK")) |
          (!is.na(atc_code_level3) &
          stringr::str_detect(atc_code_level3,"L01X(A|B)")) |
        (stringr::str_detect(
          atc_code_level2, "(L01(A|B|C|D))"
          ))) ~ "cancer_chemo_therapy",
      (!is.na(atc_code_level3) &
          drug_type == "Antibody" &
          stringr::str_detect(
            atc_code_level3, "L01XXK")
       ) ~ "cancer_adc_targeted_therapy",
      !is.na(target_symbol) &
      !is.na(atc_code_level2) &
        atc_code_level2 != "NA" &
        !stringr::str_detect(
          atc_code_level2, "\\|L0|^L0"
          ) ~ "other_targeted_therapy",
      TRUE ~ as.character(atc_treatment_category)
    ))
  
  
  cancer_targets_classified <- classified_drugs_all_final |>
    dplyr::filter(
      (drug_cancer_relevance == "by_cancer_condition_otp" |
         drug_cancer_relevance == "by_cancer_target_otp" |
        drug_cancer_relevance == "by_cancer_target_nci") & 
        atc_treatment_category == "cancer_targeted_therapy") |>
    dplyr::select(target_symbol) |>
    dplyr::filter(!is.na(target_symbol)) |>
    dplyr::filter(!stringr::str_detect(target_symbol,"^(RP(S|L)|TUB)")) |>
    dplyr::distinct() |>
    dplyr::mutate(cancer_target_classified = TRUE)
  
  classified_drugs_all_final <- classified_drugs_all_final |>
    dplyr::left_join(cancer_targets_classified) |>
    dplyr::mutate(cancer_target_classified = dplyr::if_else(
      is.na(target_symbol) | is.na(cancer_target_classified),
      FALSE,
      as.logical(cancer_target_classified)
    ))
    
  drug_df <- remove_duplicate_chembl_ids(
    drug_df = classified_drugs_all_final)
  
  drug_df$drug_entry <- NULL
  
  return(drug_df)
  
}


remove_duplicate_chembl_ids <- function(drug_df = NULL){
  
  black_list <- 
    readr::read_tsv(
      file = "data-raw/drug_name_black_list.txt", 
      col_names = F, show_col_types = F)
  
  black_list$nci_cd_name <- black_list$X1
  
  drug_df <- drug_df |> dplyr::anti_join(
    black_list, by = "nci_cd_name"
  )
  
  return(drug_df)
  
}

clean_final_drug_list <- function(drug_df = NULL){

  pharmaoncox <- drug_df |>
    dplyr::filter(
      is.na(nci_cd_name) |
      !stringr::str_detect(
        nci_cd_name,
        "Sustained-release| Bead(s)?| Compound|Vector|Pegylated")
    ) |>
    dplyr::filter(
      is.na(nci_concept_definition) |
      !stringr::str_detect(
        tolower(nci_concept_definition),
        "a( (synthetic|diagnostic|targeted))?( radio(immuno)conjugate)")
    ) |>
    dplyr::filter(
      is.na(drug_action_type) |
      !stringr::str_detect(tolower(drug_action_type),"vaccine")
    ) |>
    dplyr::filter(
      is.na(nci_concept_definition) |
        !stringr::str_detect(
          tolower(nci_concept_definition),
          "(^(a|any)|cancer|tumor|dna|autologous|cell-based|synthetic|cell|peptide|valent)( cell)? vaccine")
    ) |>
    dplyr::filter(
      is.na(nci_cd_name) |
      !stringr::str_detect(
        nci_cd_name,
        "^(Carbon C|Fluorine F|Gallium Ga|Indium In|Iodine I|Lutetium Lu|Technetium Tc|Yttrium Y)"
      )
    ) |>
    dplyr::filter(is.na(nci_cd_name) | 
                    nci_cd_name != "Abivertinib Maleate") |>
    dplyr::mutate(nci_cd_name = dplyr::if_else(
      is.na(nci_cd_name),
      stringr::str_to_title(drug_name),
      as.character(nci_cd_name)
    )) |>
    dplyr::mutate(is_salt = dplyr::if_else(
      nci_cd_name == "Avitinib Maleate",
      as.logical(NA),
      as.logical(is_salt)
    )) |>
    dplyr::mutate(drug_name = dplyr::if_else(
      nci_cd_name == "Doxycycline" &
        !is.na(molecule_chembl_id) &
        molecule_chembl_id == "CHEMBL1433",
      "DOXYCYCLINE",
      as.character(drug_name)))


  drug_action_types <- as.data.frame(
    pharmaoncox |>
      dplyr::select(nci_cd_name, drug_action_type) |>
      dplyr::distinct() |>
      dplyr::group_by(nci_cd_name) |>
      dplyr::summarise(drug_action_type = paste(
        drug_action_type, collapse = "/"
      ))
  )

  pharmaoncox$drug_action_type <- NULL
  pharmaoncox <- pharmaoncox |>
    dplyr::left_join(
      drug_action_types,
      by = "nci_cd_name", 
      relationship = "many-to-many") |>
    dplyr::select(drug_name,
                  nci_cd_name,
                  drug_type,
                  drug_action_type,
                  molecule_chembl_id,
                  drug_moa,
                  drug_max_phase_indication,
                  dplyr::everything())


  drug_max_ct_phase <- as.data.frame(
    pharmaoncox |>
      dplyr::select(nci_cd_name, drug_max_ct_phase) |>
      dplyr::group_by(nci_cd_name) |>
      dplyr::summarise(drug_max_ct_phase = max(drug_max_ct_phase))
  )

  pharmaoncox$drug_max_ct_phase <- NULL

  pharmaoncox <- pharmaoncox |>
    dplyr::left_join(drug_max_ct_phase,
                     by = "nci_cd_name",
                     relationship = "many-to-many") |>
    dplyr::select(-c(drug_moa)) |>
    dplyr::mutate(nci_concept_synonym2 = dplyr::if_else(
      is.na(nci_concept_synonym_all) & !is.na(drug_synonyms),
      as.character(tolower(drug_synonyms)),
      as.character(tolower(nci_concept_synonym_all))
    )) |>
    dplyr::mutate(nci_concept_synonym_all2 = nci_concept_synonym_all) |>
    tidyr::separate_rows(nci_concept_synonym2,
                         sep="\\|") |>
    dplyr::rename(nci_concept_synonym = nci_concept_synonym2) |>
    dplyr::select(-c(nci_concept_synonym_all2,
                     drug_synonyms,
                     drug_tradenames,
                     drug_description)) |>
    dplyr::distinct() |>
    dplyr::select(drug_name, nci_cd_name, drug_type,
                  drug_action_type, drug_cancer_relevance,
                  molecule_chembl_id,
                  drug_max_phase_indication, drug_max_ct_phase,
                  target_genename, target_symbol,
                  target_type, target_ensembl_gene_id,
                  target_entrezgene,
                  disease_efo_id, disease_efo_label,
                  cui, cui_name, primary_site,
                  nci_concept_synonym,
                  nci_concept_synonym_all,
                  dplyr::everything()) |>
    dplyr::mutate(nci_concept_definition =
                    stringi::stri_enc_toascii(nci_concept_definition)) |>
    dplyr::mutate(nci_concept_synonym_all =
                    stringi::stri_enc_toascii(nci_concept_synonym_all)) |>
    dplyr::mutate(nci_concept_synonym =
                    stringi::stri_enc_toascii(nci_concept_synonym)) |>
    dplyr::mutate(drug_name =
                    stringi::stri_enc_toascii(drug_name)) |>
    dplyr::mutate(
      nci_cd_name =
        stringi::stri_enc_toascii(nci_cd_name)) |>
    dplyr::mutate(drug_action_type = stringr::str_replace_all(
      drug_action_type, "/NA|NA/",""
    )) |>
    dplyr::mutate(idx = dplyr::row_number())
    


  ## Simplify records with only "cancer" indications, 
  ## mapping them to a unique
  ## EFO/CUI cross-ref, avoiding similar records 
  ##with "neoplasm", "carcinoma" etc.

  pharmaoncox_cancer_no_indication <- pharmaoncox |>
    dplyr::filter(is.na(disease_efo_id))
  
  pharmaoncox_non_cancer <- pharmaoncox |>
    dplyr::filter(is.na(primary_site) & !is.na(disease_efo_id) &
                    drug_cancer_relevance == "by_other_condition_otp")

  pharmaoncox_cancer_NOS <- as.data.frame(
    pharmaoncox |>
      dplyr::filter(is.na(primary_site) & !is.na(disease_efo_id) &
                      drug_cancer_relevance != "by_other_condition_otp") |>
      dplyr::mutate(disease_efo_id = "EFO:0000311",
                    disease_efo_label = "cancer",
                    cui = "C0006826",
                    cui_name = "Malignant neoplastic disease") |>
      dplyr::group_by(
        dplyr::across(-dplyr::ends_with(c("drug_clinical_id")))) |>
      dplyr::summarise(
        drug_clinical_id = paste(unique(drug_clinical_id), collapse=","),
        .groups = "drop"
      ) |>
      dplyr::distinct()
  )

  pharmaoncox_cancer_specific <- pharmaoncox |>
    dplyr::filter(!is.na(primary_site))


  pharmaoncox <- pharmaoncox_cancer_no_indication |>
    dplyr::bind_rows(pharmaoncox_non_cancer) |>
    dplyr::bind_rows(pharmaoncox_cancer_specific) |>
    dplyr::bind_rows(pharmaoncox_cancer_NOS) |>
    dplyr::arrange(nci_cd_name) |>
    dplyr::mutate(nci_cd_name = dplyr::if_else(
      stringr::str_detect(nci_cd_name,";-"),
      drug_name,
      as.character(nci_cd_name)
    )) |>
    dplyr::filter(!is.na(nci_cd_name)) |>
    dplyr::mutate(drug_name_final = dplyr::case_when(
      !is.na(drug_name) & 
        !is.na(nci_cd_name) & 
        tolower(nci_cd_name) != tolower(drug_name) &
        stringr::str_detect(
          tolower(nci_cd_name),"(ib|mab|in)$") &
        stringr::str_detect(
          drug_name, "[0-9]{1,}") ~ nci_cd_name,
      !is.na(drug_name) & 
        !is.na(nci_cd_name) & 
        tolower(nci_cd_name) != tolower(drug_name) &
        !stringr::str_detect(
          tolower(nci_cd_name),"(ib|mab|in)$") &
        stringr::str_detect(
          drug_name, "[0-9]{1,}") ~ drug_name,
      !is.na(drug_name) & 
        !is.na(nci_cd_name) & 
        tolower(nci_cd_name) == tolower(drug_name) ~ 
        stringr::str_to_title(tolower(nci_cd_name)),
      is.na(drug_name) & !is.na(nci_cd_name) ~ nci_cd_name,
      TRUE ~ as.character(
        stringr::str_to_title(tolower(drug_name))
      )
    )) |>
    dplyr::filter(
      !(drug_name_final == "CFI-400945" &
          molecule_chembl_id == "CHEMBL3408947") &
        !(drug_name_final == "Relebactam" &
            molecule_chembl_id == "CHEMBL3112741") &
        !(drug_name_final == "SR16234" &
            molecule_chembl_id == "CHEMBL3545210")
    )


  pharmaoncox$drug_name <- NULL
  pharmaoncox$drug_name <- pharmaoncox$drug_name_final
  pharmaoncox$nci_cd_name <- NULL
  pharmaoncox$drug_name_final <- NULL
  pharmaoncox$nci_concept_synonym <- NULL

  blackbox_warnings <- as.data.frame(pharmaoncox |>
    dplyr::select(drug_name, drug_blackbox_warning) |>
    dplyr::distinct() |>
      dplyr::group_by(drug_name) |>
    dplyr::summarise(drug_blackbox_warning = paste(
      unique(drug_blackbox_warning), collapse=","))
  )

  nci_t_map <- as.data.frame(pharmaoncox |>
    dplyr::select(drug_name, nci_t, nci_concept_definition) |>
    dplyr::distinct() |>
    dplyr::group_by(drug_name) |>
    dplyr::summarise(
      nci_t = paste(
        unique(nci_t), collapse=","),
      nci_concept_definition = paste(
        unique(nci_concept_definition),
        collapse = ",")
    ), .groups = "drop"
  )

  pharmaoncox$nci_t <- NULL
  pharmaoncox$nci_concept_definition <- NULL
  pharmaoncox$drug_blackbox_warning <- NULL
  pharmaoncox2 <- pharmaoncox |>
    dplyr::left_join(
      nci_t_map, by = "drug_name", 
      relationship = "many-to-many") |>
    dplyr::left_join(
      blackbox_warnings, by = "drug_name", 
      relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::mutate(drug_action_type = dplyr::if_else(
      drug_action_type == "NA" &
        stringr::str_detect(tolower(drug_name),
                            "ib |ib$") &
        stringr::str_detect(tolower(nci_concept_definition),
                            "inhibitor"),
      "INHIBITOR",
      as.character(drug_action_type)
    )) |>
    dplyr::mutate(drug_type = dplyr::if_else(
      drug_action_type == "NA" &
        stringr::str_detect(tolower(drug_name),
                            "ib |ib$") &
        stringr::str_detect(tolower(nci_concept_definition),
                            "inhibitor"),
      "Small molecule",
      as.character(drug_type)
    )) |>
    
    dplyr::mutate(inhibition_moa = dplyr::if_else(
      !stringr::str_detect(
        drug_action_type,
        paste0(
          "^(POSITIVE|NEGATIVE|AGONIST|HYDROLYTIC|",
          "CHELATING|PARTIAL|NA|INVERSE|OPENER|",
          "RELEASING|ANTISENSE|STABILISER|MODULATOR|",
          "SEQUESTERING|OXIDATIVE|OTHER|ACTIVATOR|RNA|",
          "VACCINE|DEGRADER|PROTEOLYTIC|SUBSTRATE)")) &
      (is.na(drug_type) | 
         (drug_type != "Protein" & 
            drug_type != "Gene" & 
            drug_type != "Cell")) &
      (is.na(target_type) | 
         target_type != "protein_protein_interaction"),
      as.logical(TRUE),
      as.logical(FALSE),
      as.logical(FALSE)
    )) |>
    dplyr::mutate(inhibition_moa = dplyr::if_else(
      inhibition_moa == F &
        drug_action_type == "OTHER" &
        (stringr::str_detect(
          tolower(drug_name),"mab$") |
           drug_type == "Antibody"),
      as.logical(TRUE),
      as.logical(inhibition_moa)
    )) 
  
  inhibition_moa_df <- as.data.frame(pharmaoncox2 |>
    dplyr::select(drug_name, inhibition_moa) |>
    dplyr::distinct() |>
    dplyr::group_by(drug_name) |>
    dplyr::summarise(inhibition_moa = paste(
      unique(inhibition_moa), collapse="|"
    ), .groups = "drop")
  ) |>
    dplyr::mutate(inhibition_moa = dplyr::if_else(
      stringr::str_detect(inhibition_moa, "\\|"),
      as.logical(TRUE),
      as.logical(inhibition_moa)
    ))
  
  pharmaoncox2$inhibition_moa <- NULL
  pharmaoncox2 <- pharmaoncox2 |> 
    dplyr::left_join(
      inhibition_moa_df, 
      by = "drug_name", 
      relationship = "many-to-many") |>
    dplyr::filter(!((!is.na(molecule_chembl_id) &
                       molecule_chembl_id == "CHEMBL1742994") &
                      drug_name == "Brentuximab vedotin")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL1908394") &
                      drug_name == "GSK-461364")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) &
                       molecule_chembl_id == "CHEMBL2158685") &
                    drug_name == "ABC-294640")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL2331680") &
                      drug_name == "RG-7603")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL253969") &
                      drug_name == "OSI-632")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL283120") &
                    drug_name == "AXL-1717")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL52885") &
                    drug_name == "ENMD-981693")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4650827") &
                      drug_name == "SNDX-5613 FREE BASE")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL1645462") &
                      drug_name == "AC-480")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4303241") &
                    drug_name == "BAY-1161909")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4298098") &
                    drug_name == "SAR-408701")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4091801") &
                    drug_name == "APG115")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL3623290") &
                    drug_name == "AZD-3759")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL103") &
                      drug_name == "Mycophenolic Acid")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL3545007") &
                    drug_name == "RG-7602")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                      molecule_chembl_id == "CHEMBL3545003") &
                    drug_name == "MSC-2363318A")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4297310") &
                    drug_name == "ONC-201")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4297930") &
                    drug_name == "CYC-065")) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL3899477") &
                    (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL3126004") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL376408") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL4297489") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL17157") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL1201231") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!((!is.na(molecule_chembl_id) & 
                       molecule_chembl_id == "CHEMBL1433") &
                      (nci_t == "NA" | is.na(nci_t)))) |>
    dplyr::filter(!(drug_name == "Risedronate Sodium" | 
                      drug_name == "Reminertant" |
                      drug_name == 'Ibandronate Sodium'))

  
  

  drug_maps <- list()
  drug_maps[['id2name']] <- pharmaoncox2 |>
    dplyr::select(drug_name) |>
    dplyr::distinct() |>
    dplyr::mutate(
      drug_id = dplyr::row_number())

  if(nrow(drug_maps[['id2name']]) != length(
    unique(pharmaoncox2$drug_name))){
    lgr::lgr$fatal("Ambiguous drug names")
    return(0)
  }

  pharmaoncox2 <- pharmaoncox2 |>
    dplyr::left_join(
      drug_maps[['id2name']], 
      by = "drug_name", 
      multiple = "all")

  drug_maps[['id2target']] <- pharmaoncox2 |>
    dplyr::select(drug_id,
                  target_symbol,
                  target_entrezgene,
                  target_genename,
                  target_ensembl_gene_id,
                  target_type) |>
    dplyr::distinct()

  drug_maps[['id2indication']] <- pharmaoncox2 |>
    dplyr::select(drug_id,
                  drug_max_phase_indication,
                  drug_approved_indication,
                  drug_frac_cancer_indications,
                  drug_approved_noncancer,
                  drug_n_indications,
                  drug_clinical_source,
                  drug_year_first_approval,
                  drug_clinical_id,
                  drug_max_ct_phase,
                  disease_efo_id,
                  disease_efo_label,
                  cui,
                  cui_name,
                  primary_site) |>
    dplyr::distinct()

  drug_maps[['id2synonym']] <- pharmaoncox2 |>
    dplyr::select(drug_id,
                  nci_concept_synonym_all) |>
    dplyr::distinct()

  drug_maps[['id2basic']] <- as.data.frame(
    pharmaoncox2 |>
    dplyr::select(
      drug_id,
      molecule_chembl_id,
      drug_type,
      drug_action_type,
      drug_cancer_relevance,
      inhibition_moa,
      is_salt,
      is_adc,
      drug_blackbox_warning,
      nci_concept_definition,
      nci_t,
      opentargets,
      atc_code_level1,
      atc_level1,
      atc_code_level2,
      atc_level2,
      atc_code_level3,
      atc_level3,
      atc_treatment_category) |>
      dplyr::distinct())

  return(drug_maps)

}


expand_drug_aliases <- function(drug_index_map = NULL,
                                path_data_raw = NULL,
                                chembl_pubchem_datestamp = chembl_pubchem_datestamp){

  chembl_pubchem_xref <-
    get_chembl_pubchem_xref(
      datestamp = chembl_pubchem_datestamp,
      path_data_raw = path_data_raw)
  
  drugAliasPrimary <- 
    drug_index_map[['id2name']] |>
    dplyr::select(drug_id, drug_name) |>
    dplyr::mutate(alias = drug_name) |>
    dplyr::select(-drug_name) |>
    dplyr::mutate(alias_source = "primaryName")
  
  
  drugAliasNCI <- as.data.frame(
    drug_index_map[['id2synonym']] |>
      dplyr::select(drug_id, nci_concept_synonym_all) |>
      tidyr::separate_rows(nci_concept_synonym_all, sep="\\|") |>
      dplyr::distinct() |>
      dplyr::rename(alias = nci_concept_synonym_all) |>
      dplyr::filter(nchar(alias) > 3) |>
      dplyr::mutate(alias_source = "nci")
  )

  ## Extend aliases with those found in PubChem

  ## get drug set that contains PubChem cross-references
  unique_chembl_pubchem <- drug_index_map[['id2basic']] |>
    dplyr::select(drug_id, molecule_chembl_id) |>
    dplyr::filter(!is.na(molecule_chembl_id)) |>
    dplyr::distinct() |>
    dplyr::left_join(
      chembl_pubchem_xref, 
      by = "molecule_chembl_id",
      multiple = "all", relationship = "many-to-many") |>
    dplyr::filter(!is.na(pubchem_cid)) |>
    dplyr::select(-c(chembl_db_version))


  ## Retrieve aliases for drugs with PubChem x-refs
  lgr::lgr$info("Retrieving additional drug aliases/synonyms for PubChem-mapped compounds")
  pubchem_synonym_files <-
    sort(list.files(path = file.path(path_data_raw, "pubchem"),
                    pattern = "CID-Synonym-filtered_",
                    full.names = T))

  drugAliasPubchem <- data.frame()
  for(f in pubchem_synonym_files){
    synonym_data <- as.data.frame(readr::read_tsv(
      f, col_names = c('pubchem_cid','alias'),
      col_types = "dc",
      progress = F
    ))

    pubchem_alias_df <- synonym_data |>
      dplyr::inner_join(
        unique_chembl_pubchem,
        by = "pubchem_cid", 
        multiple = "all", relationship = "many-to-many")
    
    if(nrow(pubchem_alias_df) > 0){
      pubchem_alias_df <- pubchem_alias_df |>
        dplyr::select(-c(pubchem_cid, molecule_chembl_id)) |>
        dplyr::mutate(alias_source = "pubchem")
      drugAliasPubchem <-
        drugAliasPubchem |>
        dplyr::bind_rows(pubchem_alias_df)
    }
    rm(synonym_data)
  }
  
  ## ignore some drugs
  custom_name_ignore <- as.data.frame(readr::read_tsv(
    file = "data-raw/drug_names_ignore.tsv",
    col_names = F, show_col_types = F))
  colnames(custom_name_ignore) <- c("alias_lc")
  custom_name_ignore$alias_lc <- tolower(
    custom_name_ignore$alias_lc)
  
  
  drugAliasAll <- as.data.frame(drugAliasPrimary |>
    dplyr::bind_rows(drugAliasNCI) |>
    dplyr::bind_rows(drugAliasPubchem) |>
    dplyr::filter(nchar(alias) > 3 & 
                    nchar(alias) < 250) |>
    #tidyr::separate_rows(alias, sep = "\\|") |>
    dplyr::distinct() |>
    dplyr::left_join(drug_index_map[['id2name']],
                     by = "drug_id") |>
      
    ## avoid drug aliases for a particular drug that are identical 
    ## to the primary drug name of another drug
    dplyr::filter(
      !(tolower(alias) != tolower(drug_name) &
          tolower(alias) %in% tolower(drugAliasPrimary$alias))) |>
      
    ## consider unambiguous drug aliases only
    dplyr::group_by(alias) |>
    dplyr::summarise(
      alias_source = paste(unique(alias_source), collapse="|"),
      drug_id = paste(unique(drug_id), collapse=","),
      .groups = "drop") |>
      dplyr::filter(!stringr::str_detect(drug_id,",")) |>
      dplyr::distinct() |>
      dplyr::mutate(
        alias = stringi::stri_enc_toascii(alias)
      ) |>
      dplyr::mutate(
        alias_lc = tolower(alias)) |> 
      
      ## remove aliases that co-incide with ordinary
      ## english words (pubchem)
      dplyr::left_join(
        words::words, 
        by = c("alias_lc" = "word"),
        multiple = "all", relationship = "many-to-many") |> 
      dplyr::filter(
        is.na(word_length) |
          (!is.na(word_length) &
             word_length > 7 &
             (alias_source != "pubchem" ))) |>
      dplyr::anti_join(
        custom_name_ignore, by = "alias_lc") |>
      dplyr::select(-c(word_length,alias_lc)) |>
      dplyr::distinct() |>
      dplyr::select(drug_id, alias, alias_source) |>
      dplyr::mutate(
        drug_id = as.integer(drug_id)
      )
      
  )

  return(drugAliasAll)


}


get_mesh_drug_categories <- function(path_data_raw = NULL){
  
  mesh_drug_data <- whatamesh::read_mesh_file(
    file.path(path_data_raw,"data-raw",
              "mesh","d2022.bin")
  )
  i <- 1
  
  all_drugs <- data.frame()
  while(i <= nrow(mesh_drug_data)){
    drugname <- mesh_drug_data[i,]$MH
    if(!is.null(mesh_drug_data[i,]$PA)){
      drug_category <- paste(unlist(mesh_drug_data[i,]$PA),
                             collapse=";")
      
      df <- data.frame(
        'drugname' = drugname,
        'category' = drug_category,
        stringsAsFactors = F
      )
      
      all_drugs <- all_drugs |>
        dplyr::bind_rows(df)
    }
    i <- i + 1
  }
  
}
