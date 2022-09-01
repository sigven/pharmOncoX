

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

### PUBCHEM/CHEMBL - WEB SERVICE FUNCTIONS FOR COMPOUND INFORMATION

## Retrieval of compound properties from PubChem
#'
#' Function that retrieves PubChem compound properties using PubChem public user gateway (PUG)
#'
#' @param pubchem_cid PubChem compound identifier
#' @return a data frame with the following columns of chemical compound properties:
#' \itemize{
#'   \item pubchem_isomeric_SMILES
#'   \item pubchem_TPSA
#'   \item pubchem_IUPAC_name
#'   \item pubchem_complexity
#'   \item pubchem_hbond_donor_count
#'   \item pubchem_hbond_acceptor_count
#'   \item pubchem_molecular_weight
#'   \item pubchem_cid
#'   \item pubchem_name
#'   \item molecule_chembl_id
#' }
#' @examples
#' ## Retrieve compound properties for Azacitidine (Pubchem compound ID = 9444)
#' \dontrun{
#' pubchem_compound_props(pubchem_cid = 9444)
#' }
#'
pubchem_compound_properties <- function(PUBCHEM_PUG_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/', pubchem_cid = 9444){
  compound_description_list <- jsonlite::fromJSON(paste0(PUBCHEM_PUG_URL,'compound/cid/',pubchem_cid,'/description/JSON'))
  compound_description <- dplyr::filter(compound_description_list$InformationList$Information, !is.na(Title))
  if(nrow(compound_description) > 1){
    cat("More than a single description line\n",pubchem_cid)
    cat('\n')
    return(NULL)
  }
  compound_properties <- read.csv(paste0(PUBCHEM_PUG_URL,'compound/cid/',pubchem_cid,'/property/IsomericSMILES,TPSA,IUPACName,Complexity,HBondDonorCount,HBondAcceptorCount,MolecularWeight/CSV'),stringsAsFactors = F)

  compound_properties$pubchem_name <- compound_description$Title
  chembl_id <- NA
  synonym_url <- paste0(PUBCHEM_PUG_URL,'compound/cid/',pubchem_cid,'/synonyms/JSON')
  pubchem_compound_synonyms <- NULL
  if(!is.null(readUrl(synonym_url))){
    compound_synonyms <- jsonlite::fromJSON(synonym_url)
    all_synonyms <- compound_synonyms$InformationList$Information$Synonym[[1]]
    pubchem_compound_synonyms <- paste(all_synonyms,collapse="@@@@")
    j <- 1
    while(j <= length(all_synonyms)){
      if(stringr::str_detect(all_synonyms[j],"^CHEMBL[0-9]{1,}")){
        chembl_id <- all_synonyms[j]
        break
      }
      j <- j + 1
    }
  }

  compound_properties$molecule_chembl_id <- chembl_id
  compound_properties$pubchem_synonyms <- pubchem_compound_synonyms
  compound_properties <- dplyr::rename(compound_properties, pubchem_isomeric_SMILES = IsomericSMILES, pubchem_TPSA = TPSA, pubchem_IUPAC_name = IUPACName, pubchem_complexity = Complexity, pubchem_hbond_donor_count = HBondDonorCount, pubchem_hbond_acceptor_count = HBondAcceptorCount, pubchem_molecular_weight = MolecularWeight, pubchem_id = CID)

  return(compound_properties)
}

## Retrieval of molecular compound properties from ChEMBL
#'
#' Function that accepts a ChEMBL molecule identifer and uses the ChEMBL web service API to return a range of compound properties
#'
#' @param molecule_chembl_id ChEMBL molecule identifier
#' @return a data frame with the following columns of chemical compound properties:
#' \itemize{
#'   \item molecule_chembl_id
#'   \item parent_chembl_id
#'   \item chembl_acd_logp
#'   \item chembl_acd_logd
#'   \item chembl_acd_most_apka
#'   \item chembl_acd_most_bpka
#'   \item chembl_alogp
#'   \item chembl_aromatic_rings
#'   \item chembl_full_molformula
#'   \item chembl_full_mwt
#'   \item chembl_hba
#'   \item chembl_hba_lipinski
#'   \item chembl_hbd
#'   \item chembl_hbd_lipinski
#'   \item chembl_heavy_atoms
#'   \item chembl_molecular_species
#'   \item chembl_mw_freebase
#'   \item chembl_mw_monoisotopic
#'   \item chembl_num_alerts
#'   \item chembl_num_lipinski_ro5_violations
#'   \item chembl_num_ro5_violations
#'   \item chembl_psa
#'   \item chembl_qed_weighted
#'   \item chembl_ro3_pass
#'   \item chembl_rtb
#'   \item chembl_canonical_smiles
#'   \item chembl_standard_inchi
#'   \item chembl_standard_inchi_key
#'   \item chembl_pref_name
#'   \item lincs_id
#'   \item drugbank_id
#'   \item pharmgkb_id
#' }
#' @examples
#' ## Retrieve compound properties for IRINOTECAN (ChEMBL molecular compound ID = 'CHEMBL481')
#' \dontrun{
#' chembl_compound_props(molecule_chembl_id = 'CHEMBL481')
#' }
#'
chembl_compound_properties <- function(chembl_ws_base_url = 'https://www.ebi.ac.uk/chembl/api/data', molecule_chembl_id = NA){
  all_molecule_properties <- data.frame("molecule_chembl_id" = molecule_chembl_id, "parent_chembl_id" = NA, "chembl_acd_logp" = NA, "chembl_acd_logd" = NA, "chembl_acd_most_apka" = NA,
                                        "chembl_acd_most_bpka" = NA, "chembl_alogp" = NA, "chembl_aromatic_rings" = NA, "chembl_full_molformula" = NA,
                                        "chembl_full_mwt" = NA, "chembl_hba" = NA, "chembl_hba_lipinski" = NA, "chembl_hbd" = NA, "chembl_hbd_lipinski" = NA,
                                        "chembl_heavy_atoms" = NA, "chembl_molecular_species" = NA, "chembl_mw_freebase" = NA, "chembl_mw_monoisotopic" = NA,
                                        "chembl_num_alerts" = NA, "chembl_num_lipinski_ro5_violations" = NA, "chembl_num_ro5_violations" = NA, "chembl_psa" = NA,
                                        "chembl_qed_weighted" = NA, "chembl_ro3_pass" = NA, "chembl_rtb" = NA, "chembl_canonical_smiles" = NA,
                                        "chembl_standard_inchi" = NA, "chembl_standard_inchi_key" = NA, "chembl_pref_name" = NA, "lincs_id" = NA, "drugbank_id" = NA,
                                        "pharmgkb_id" = NA, stringsAsFactors = F)

  molecule_url <- paste0(chembl_ws_base_url,'/molecule?molecule_chembl_id=',molecule_chembl_id)
  molecule_hierarchy <- NULL
  molecule_properties <- NULL
  molecule_structures <- NULL
  props_all <- NULL
  lincs_id <- NA
  drugbank_id <- NA
  pharmgkb_id <- NA
  if(!is.null(readUrl(molecule_url))){
    raw_xml_string <- rawToChar(httr::GET(molecule_url)$content)
    doc <- XML::xmlParse(raw_xml_string)
    if(stringr::str_detect(raw_xml_string, '<molecule_hierarchy>')){
      molecule_hierarchy <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/molecules/molecule/molecule_hierarchy"), collectNames = T, stringsAsFactors = F)
    }
    if(stringr::str_detect(raw_xml_string, '<molecule_properties>')){
      molecule_properties <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/molecules/molecule/molecule_properties"), collectNames = T, stringsAsFactors = F)
    }
    if(!is.null(molecule_properties)){
      colnames(molecule_properties) <- paste0('chembl_',colnames(molecule_properties))
    }
    if(stringr::str_detect(raw_xml_string, '<molecule_structures>')){
      molecule_structures <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/molecules/molecule/molecule_structures"), collectNames = T, stringsAsFactors = F)
    }
    if(!is.null(molecule_structures)){
      colnames(molecule_structures) <- paste0('chembl_',colnames(molecule_structures))
    }
    if(stringr::str_detect(raw_xml_string, '<molecule>')){
      props_all <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/molecules/molecule"), collectNames = T, stringsAsFactors = F)
      props_all <- dplyr::select(props_all, pref_name) |> dplyr::rename(chembl_pref_name = pref_name)
    }

    if(!is.null(molecule_structures)){
      unichem_url <- paste0('https://www.ebi.ac.uk/unichem/rest/verbose_inchikey/',molecule_structures$chembl_standard_inchi_key)
      if(!is.null(readUrl(unichem_url))){
        unichem_xrefs <- httr::content(httr::GET(unichem_url))
        i <- 1
        while(i <= length(unichem_xrefs)){
          list_item <- unichem_xrefs[[i]]
          if(list_item$name == 'lincs'){
            lincs_id <- list_item$src_compound_id[[1]]
          }
          if(list_item$name == 'drugbank'){
            drugbank_id <- list_item$src_compound_id[[1]]
          }
          if(list_item$name == 'pharmgkb'){
            pharmgkb_id <- list_item$src_compound_id[[1]]
          }
          i <- i + 1
        }
      }
    }

    if(!is.null(molecule_properties)){
      for(c in colnames(molecule_properties)){
        all_molecule_properties[,c] <- molecule_properties[,c]
      }
    }
    if(!is.null(molecule_structures)){
      for(c in colnames(molecule_structures)){
        all_molecule_properties[,c] <- molecule_structures[,c]
      }
    }
    if(!is.null(props_all)){
      for(c in colnames(props_all)){
        all_molecule_properties[,c] <- props_all[,c]
      }
    }

    all_molecule_properties$lincs_id <- lincs_id
    all_molecule_properties$drugbank_id <- drugbank_id
    all_molecule_properties$pharmgkb_id <- pharmgkb_id
    if(!is.null(molecule_hierarchy)){
      all_molecule_properties$parent_chembl_id <- molecule_hierarchy$parent_chembl_id
    }
  }
  return(all_molecule_properties)
}


### NCI DRUG DISPLAY LABELS

process_nci_labels <- function(path_data_raw, overwrite = F) {

  nci_thesaurus_labels_filepath = file.path(path_data_raw, "nci_thesaurus", "ThesaurusLabels.tsv")
  nci_thesaurus_owl_filepath <- file.path(path_data_raw, "nci_thesaurus", "ThesaurusInferred.owl")
  if(!is.null(nci_thesaurus_owl_filepath) & file.exists(nci_thesaurus_owl_filepath) & (!file.exists(nci_thesaurus_labels_filepath) | overwrite == T)){
    system(paste0('python data-raw/parse_owl_labels.py ',nci_thesaurus_owl_filepath,' > ',file.path(path_data_raw,'nci_thesaurus','ThesaurusLabels.tsv')))
  }
  nci_labels <- read.table(file = nci_thesaurus_labels_filepath, sep="\t",comment.char="",
                           stringsAsFactors = F,quote="",header=F) |>
    magrittr::set_colnames(c('nci_t','nci_cd_name','cui')) |>
    dplyr::filter(!is.na(nci_t))

  return(nci_labels)
}

### CHEMBL-PUBCHEM COMPOUND CROSS-REFERENCE

get_chembl_pubchem_compound_xref <- function(datestamp = '20220429',
                                             chembl_release = "v31",
                                             path_data_raw = NULL,
                                             update = F){
  chembl_pubchem_xref_fname <- file.path(path_data_raw, "chembl",
                                         paste0("chembl_pubchem_mapping_", datestamp,".txt.gz"))
  if(!file.exists(chembl_pubchem_xref_fname) | update == T){
    download.file("ftp://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz",
                  destfile = chembl_pubchem_xref_fname)
  }
  chembl_pubchem_compound_xref <- read.table(gzfile(chembl_pubchem_xref_fname), header = F, skip = 1,
                                             stringsAsFactors = F,sep="\t",quote="") |>
    magrittr::set_colnames(c('molecule_chembl_id','pubchem_cid')) |>
    dplyr::mutate(chembl_db_version = chembl_release)

  return(chembl_pubchem_compound_xref)
}



### TARGETED ANTICANCER COMPOUNDS FROM OPEN TARGETS
get_opentargets_cancer_drugs <-
  function(path_data_raw = NULL,
           ot_version = "2022.06"){

  phenotype_cancer_efo <- oncoPhenoMap::auxiliary_maps$umls$concept |>
    dplyr::filter(main_term == T) |>
    dplyr::select(cui, cui_name) |>
    dplyr::distinct() |>
    dplyr::inner_join(dplyr::select(oncoPhenoMap::oncotree_expanded_full,
                                    efo_id, cui,
                                    cui_name, primary_site),
                      by = c("cui", "cui_name")) |>
    dplyr::rename(disease_efo_id = efo_id) |>
    dplyr::filter(!is.na(disease_efo_id)) |>
    dplyr::distinct()


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
                    drug_approved_indication,
                    drug_blackbox_warning,
                    drug_year_first_approval)) |>
    dplyr::mutate(
      disease_efo_id = stringr::str_replace_all(disease_efo_id, "_", ":")) |>
    dplyr::distinct()  |>
    dplyr::left_join(phenotype_cancer_efo, by = "disease_efo_id") |>
    ## general cancer ontology terms (neoplasm, cancer, carcinoma, squamous cell carcinoma)
    dplyr::mutate(
      cancer_drug =
        dplyr::if_else(
          !is.na(primary_site) |
            disease_efo_id == 'EFO:0000616' |
            disease_efo_id == 'EFO:0000311' |
            disease_efo_id == 'EFO:0000313' |
            disease_efo_id == 'EFO:0000707' |
            disease_efo_id == 'EFO:0000228',
          TRUE, FALSE)
    )

  targeted_compounds <- ot_compounds |>
    dplyr::filter(!is.na(target_symbol))

  untargeted_cancer_compounds <- as.data.frame(ot_compounds |>
    dplyr::filter(is.na(target_symbol)) |>
    dplyr::group_by_at(dplyr::vars(-c(drug_clinical_id))) |>
    dplyr::summarise(
      drug_clinical_id = paste(
        unique(drug_clinical_id),
        collapse = ","),
      drug_clinical_source = paste(
        unique(sort(drug_clinical_source)),
        collapse = ","),
      .groups = "drop") |>
    dplyr::ungroup() |>
    dplyr::filter(cancer_drug == T)
  )

  ## adjust max ct phase
  drugs_with_max_phase_adj <- as.data.frame(
    untargeted_cancer_compounds |>
      dplyr::filter(!is.na(drug_max_ct_phase) &
                      !is.na(drug_max_phase_indication)) |>
      dplyr::group_by(drug_name, molecule_chembl_id) |>
      dplyr::summarise(drug_max_ct_phase = max(drug_max_phase_indication),
                       .groups = "drop")
  )

  untargeted_cancer_compounds <- as.data.frame(
    untargeted_cancer_compounds |>
      dplyr::select(-drug_max_ct_phase) |>
      dplyr::left_join(drugs_with_max_phase_adj,
                       by = c("drug_name","molecule_chembl_id"))) |>
    dplyr::mutate(drug_name_lc = tolower(drug_name))

  targeted_compounds <- as.data.frame(
    targeted_compounds |>
      dplyr::distinct() |>
      dplyr::group_by_at(dplyr::vars(-c(drug_clinical_id))) |>
      dplyr::summarise(
        drug_clinical_id = paste(
          unique(drug_clinical_id),
          collapse = ","),
        drug_clinical_source = paste(
          unique(sort(drug_clinical_source)),
          collapse = ","),
        .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(target_symbol) &
                      !is.na(target_ensembl_gene_id)) |>
      dplyr::mutate(drug_name_lc = tolower(drug_name))
  )
  
  
  ## adjust max ct phase
  drugs_with_max_phase_adj <- as.data.frame(
    targeted_compounds |>
    dplyr::filter(!is.na(drug_max_ct_phase) & 
                    !is.na(drug_max_phase_indication)) |>
    dplyr::group_by(drug_name, molecule_chembl_id) |>
    dplyr::summarise(drug_max_ct_phase = max(drug_max_phase_indication),
                     .groups = "drop")
  )

  targeted_compounds <- as.data.frame(targeted_compounds |>
    dplyr::select(-drug_max_ct_phase) |>
    dplyr::left_join(drugs_with_max_phase_adj, 
                     by = c("drug_name","molecule_chembl_id"))
  )
  
  targeted_cancer_compounds <- targeted_compounds |>
    dplyr::filter(cancer_drug == T)
  targeted_noncancer_compounds <- targeted_compounds |>
    dplyr::filter(cancer_drug == F)
  
  return(list('targeted' = targeted_cancer_compounds,
              'targeted_noncancer' <- targeted_noncancer_compounds,
              'untargeted' = untargeted_cancer_compounds))

}

get_fda_ndc_mapping <- function(
  path_data_raw = NULL){

  fda_ndc_fname <- file.path(
    path_data_raw,"national_drug_code_fda","product.txt")

  drug2epc <- readr::read_tsv(fda_ndc_fname, show_col_types = F) |>
    janitor::clean_names() |>
    ## ignore drug regimens/combos
    dplyr::filter(!stringr::str_detect(substancename,"; ")) |>
    ## separate entries of pharm classes
    tidyr::separate_rows(pharm_classes, sep=", ") |>
    dplyr::select(substancename, pharm_classes) |>
    dplyr::filter(!is.na(pharm_classes)) |>
    dplyr::filter(stringr::str_detect(pharm_classes," \\[EPC\\]")) |>
    dplyr::rename(drug = substancename,
                  fda_epc_category = pharm_classes) |>
    dplyr::distinct()

  return(drug2epc)
}


#### NCI THESAURUS CANCER DRUGS/TREATMENTS
get_nci_drugs <- function(nci_db_release = nci_db_release,
                          overwrite = F,
                          path_data_raw = NULL,
                          path_data_processed = NULL){
  nci_antineo_thesaurus <- NULL
  nci_drugs <- NULL

  if(!file.exists(file.path(path_data_processed,"nci_thesaurus","nci_treatment_thesaurus_antineo.rds")) | overwrite == T){
    sorafenib_definition <-
      'A synthetic compound targeting growth signaling and angiogenesis. Sorafenib blocks the enzyme RAF kinase, a critical component of the RAF/MEK/ERK signaling pathway that controls cell division and proliferation; in addition, sorafenib inhibits the VEGFR-2/PDGFR-beta signaling cascade, thereby blocking tumor angiogenesis.'

    nci_display_labels <- process_nci_labels(path_data_raw = path_data_raw,
                                             overwrite = overwrite)

    drug2chembl <-
      read.table(file = file.path(path_data_raw,"dgidb","dgidb.202202.tsv"),
                 header=T, quote="", comment.char="", sep="\t",
                 stringsAsFactors = F) |>
      dplyr::filter(nchar(drug_name) > 0) |>
      dplyr::select(-drug_claim_source) |>
      dplyr::rename(molecule_chembl_id = concept_id) |>
      dplyr::mutate(molecule_chembl_id =
                      stringr::str_replace(molecule_chembl_id, "chembl:","")) |>
      dplyr::filter(!startsWith(molecule_chembl_id,"wiki")) |>
      dplyr::distinct() |>
      dplyr::mutate(drug_claim_name = tolower(drug_claim_name)) |>
      dplyr::filter(!startsWith(drug_claim_name,"chembl")) |>
      dplyr::filter(!stringr::str_detect(drug_claim_name,"^[0-9]{1,}$")) |>
      dplyr::mutate(drug_name = tolower(drug_name)) |>

      ### remove and correct wrong drug to molecule_chembl_id associations
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        drug_name == "calcifediol",
        as.character("CHEMBL1040"),
        as.character(molecule_chembl_id)
      )) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        drug_name == "leucovorin",
        as.character("CHEMBL1040"),
        as.character(molecule_chembl_id)
      )) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        drug_name == "adl-5747",
        as.character("CHEMBL561339"),
        as.character(molecule_chembl_id)
      )) |>
      dplyr::filter(molecule_chembl_id != "CHEMBL1200796") |> #cyclophosphamide
      dplyr::filter(molecule_chembl_id != "CHEMBL3989496") |> #tezacitabine
      dplyr::filter(molecule_chembl_id != "CHEMBL1588") |> #thiamine
      dplyr::filter(molecule_chembl_id != "CHEMBL1200751") |> #mercaptopurine
      dplyr::filter(molecule_chembl_id != "CHEMBL541887") |> #irinotecan
      dplyr::filter(molecule_chembl_id != "CHEMBL513000") |> #emetine hydrochloride
      dplyr::filter(molecule_chembl_id != "CHEMBL3989727") |> #vitamine e
      dplyr::filter(molecule_chembl_id != "CHEMBL1889436") |> #ouabain
      dplyr::filter(molecule_chembl_id != "CHEMBL1201113") |> #cetirizine hydrochloride
      dplyr::filter(molecule_chembl_id != "CHEMBL1200792") |> #fluphenazine hydrochloride
      dplyr::filter(molecule_chembl_id != "CHEMBL935") |> #pentamidine isethionate
      dplyr::filter(molecule_chembl_id != "CHEMBL482811") |> #u-50488 methane sulfonate
      dplyr::distinct()

    drug2chembl_all <-
      dplyr::select(drug2chembl, molecule_chembl_id, drug_claim_name) |>
      dplyr::rename(drug_name = drug_claim_name) |>
      dplyr::bind_rows(dplyr::select(drug2chembl, drug_name, molecule_chembl_id)) |>
      dplyr::filter(!stringr::str_detect(drug_name,"(^[0-9]{1,}$)|^chembl[0-9]{1,}")) |>
      dplyr::distinct() |>
      dplyr::bind_rows(data.frame('drug_name' = 'gemtuzumab',
                                  'molecule_chembl_id' = 'CHEMBL2108342',
                                  stringsAsFactors = F)) |>
      dplyr::arrange(drug_name)


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


    ## parse all entries in nci thesaurus where the semantic concept type is treatment-related
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
          nci_concept_semantic_type,"Chemical|Substance|Therapeutic|Drug|Immunologic")) |>
      dplyr::left_join(nci_antineo_agents, by = c("nci_t")) |>
      dplyr::select(-c(nci_concept_name, nci_cd_name)) |>
      # dplyr::filter(
      #   nci_concept_status != 'Obsolete_Concept' &
      #     nci_concept_status != "Retired_Concept") |>
      #dplyr::select(-nci_concept_status) |>
      dplyr::mutate(nci_concept_synonym_all = concept_synonym) |>
      tidyr::separate_rows(concept_synonym,sep="\\|") |>
      dplyr::mutate(nci_concept_synonym = tolower(concept_synonym)) |>
      dplyr::select(-concept_synonym) |>
      dplyr::filter(nci_t != 'C147908' & nci_t != 'C71622') |> ## Hormone Therapy Agent
      #dplyr::mutate(nci_t = dplyr::if_else(nci_t == 'C2194','C61948',nci_t)) |> ## redundant Sorafenib entries
      #dplyr::filter(!(nnci))
      dplyr::mutate(
        nci_concept_semantic_type =
          dplyr::if_else(nci_t == 'C61948','Pharmacologic Substance',
                         nci_concept_semantic_type)) |> ## redundant Sorafenib entries
      dplyr::mutate(
        nci_concept_definition =
          dplyr::if_else(nci_t == 'C61948', sorafenib_definition,
                         nci_concept_definition)) |> ## redundant Sorafenib entries
      dplyr::distinct() |>
      dplyr::left_join(nci_display_labels, by = c("nci_t")) |>
      dplyr::filter(!(nci_t == "C1806" & nci_concept_synonym == "gemtuzumab")) |>
      dplyr::filter(!(nci_t == "C405" & nci_concept_synonym == "ctx")) |>

      dplyr::left_join(drug2chembl_all, by=c("nci_concept_synonym" = "drug_name")) |>
      dplyr::mutate(nci_db_version = nci_db_release) |>
      dplyr::filter(!is.na(nci_cd_name)) |>
      dplyr::filter(!stringr::str_detect(
          tolower(nci_concept_definition), "coronavirus")) |>
      dplyr::filter(!stringr::str_detect(
        nci_concept_synonym_all, "SARS-CoV-2|COVID-19|CoV-19|Coronary|Corona|Covid-19|covid-19")) |>
      #))
      dplyr::mutate(antineoplastic_agent = dplyr::if_else(
        is.na(antineoplastic_agent),
        as.logical(TRUE),
        as.logical(antineoplastic_agent)
      )) |>
      # dplyr::filter(antineoplastic_agent == T |
      #                 stringr::str_detect(nci_concept_synonym,"chemother|immunother|radiation|regimen") |
      #                 stringr::str_detect(tolower(nci_cd_name),"chemother|immunother|radiation|regimen") |
      #                 stringr::str_detect(tolower(nci_cd_name),"(ib|mab)$") |
      #                 stringr::str_detect(tolower(nci_concept_definition),"immuno|b-cell|t-cell|melanoma|myelomoa|neuroblastoma|leukemia|lymphoma|thymoma|sarcoma|paraganglioma|regimen|proliferation|cancer|tumor|carcino|radiation|block|radio|chemotherapy|immune|antigen|antineo")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name," (Gel|Oil|Cream|Seed|Block|Field|Supplement|Factor)$")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"(Vaccination|Lotion|Therapeutic Heat|Procedure|Rehabilitation|Prevention|Rinse)$")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"(Epitope|Exract|Influenza|Ginseng|Ointment|Management|Injection|Tool)$")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Vitamin A Compound|Inactivated Poliovirus|Antineoplastic Immune Cell|Topical)")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Sheng-Yu|Ginseng|Dry Cleaning|Boost|Tobacco|Microwave)")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Blood|Interruption of|Gum Arabic|Vaginal Cylinder|Laser Ablation|Wheatgrass)")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Chemical Challenge|Prevention of|Magic Mouthwash|Wood Dust|Soot|Cocaine)")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name,"^(Antibody|Antigen|Antioxidant|Vaccination|Acetate|Antiserum|Asbestos|Aspirate|Autoantigen|Cytokine)$")) |>
      dplyr::filter(!stringr::str_detect(nci_cd_name," Spray| Extract| Antidiabetic| Implant|(Green Tea|Living Healthy|Pollutant|Probe|Protective Agent|Supportive Care|Caffe)")) |>
      dplyr::filter(
        !stringr::str_detect(
          tolower(nci_concept_definition),
          "chinese |antidiabet|diabetes|antidepress|analgesic|pulmonary edema|nutritional|human carcinogen|anesthetic|nonsedating|sedative|antihyper|antiinflamma|antiarrythm|antiangin|antihist|muscle|neurotransmitter"))
    )

    #nci_antineo_thesaurus
    nci2chembl <- as.data.frame(nci_antineo_thesaurus_raw |>
      dplyr::select(molecule_chembl_id, nci_cd_name) |>
      dplyr::filter(!is.na(molecule_chembl_id)) |>
      dplyr::distinct()
    )

    j <- 1
    nci2chembl_dict <- list()
    while(j <= nrow(nci2chembl)){
      nci2chembl_dict[[nci2chembl[j,"nci_cd_name"]]] <-
        nci2chembl[j,"molecule_chembl_id"]
      j <- j + 1
    }

    i <- 1
    while(i <= nrow(nci_antineo_thesaurus_raw)){
      nci_cd_name <-
        nci_antineo_thesaurus_raw[i,"nci_cd_name"]
      if(nci_cd_name %in% names(nci2chembl_dict)){
        nci_antineo_thesaurus_raw[i,"molecule_chembl_id"] <-
          nci2chembl_dict[[nci_cd_name]]
      }
      i <- i + 1
    }

    nci_antineo_thesaurus <- nci_antineo_thesaurus_raw |>
      dplyr::filter(!(molecule_chembl_id == "CHEMBL1569487" & nci_t == "C405")) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Sorafenib Tosylate", "CHEMBL1200485", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Leucovorin", "CHEMBL1679", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Masoprocol", "CHEMBL313972", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "7-Hydroxystaurosporine", "CHEMBL1236539", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Dromostanolone", "CHEMBL1201048", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Patritumab", "CHEMBL2109406", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Rivoceranib", "CHEMBL3186534", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Rivoceranib Mesylate", "CHEMBL3545414", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Tricribine Phosphate", "CHEMBL462018", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Umbrasilib" |
          nci_cd_name == "Umbrasilib Tosylate", "CHEMBL3948730",
        as.character(molecule_chembl_id))) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Brigatinib") & molecule_chembl_id != "CHEMBL3545311")) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Brentuximab") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Afatinib") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Dacomitinib") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Calcipotriene") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Balixafortide") & molecule_chembl_id != "CHEMBL2109665")) |>
      dplyr::filter(!(stringr::str_detect(nci_cd_name,"^Antineoplast(on|ic) Agent|^Support| Gel | Ointment|Caffeine|^Acetate$")))
      #CHEMBL3545055

    ## add apatinib as an alias (is missing in NCI thesaurus)
    apatinib_alias_entry <- nci_antineo_thesaurus |>
      dplyr::filter(nci_concept_synonym == "rivoceranib") |>
      dplyr::mutate(nci_concept_synonym = "apatinib") |>
      dplyr::mutate(nci_concept_synonym_all =
                      paste(nci_concept_synonym_all,"Apatinib",sep="|"))

    nci_antineo_thesaurus <- nci_antineo_thesaurus |>
      dplyr::bind_rows(apatinib_alias_entry) |>
      dplyr::arrange(nci_cd_name)

    nci_with_chembl <- nci_antineo_thesaurus |>
      dplyr::filter(!is.na(molecule_chembl_id)) |>
      dplyr::select(nci_t, molecule_chembl_id) |>
      dplyr::distinct()

    nci_2 <- nci_antineo_thesaurus |>
      dplyr::filter(is.na(molecule_chembl_id)) |>
      dplyr::select(-molecule_chembl_id) |>
      dplyr::left_join(nci_with_chembl, by = c("nci_t"))

    nci_antineo_thesaurus <- nci_antineo_thesaurus |>
      dplyr::filter(!is.na(molecule_chembl_id)) |>
      dplyr::bind_rows(nci_2) |>
      dplyr::rename(nci_drug_name = nci_concept_synonym) |>
      dplyr::filter(!(nci_cd_name == "Sorafenib Tosylate" &
                        nci_drug_name == "sorafenib")) |>
      dplyr::select(-c(nci_t_parent,cui)) |>
      dplyr::distinct() |>
      ## Removing duplicate/erroneous NCI/CHEMBL cross-ref identifiers
      dplyr::filter(!((nci_t == "C137804" | nci_t == "C137803") &
                        molecule_chembl_id == "CHEMBL1201567")) |> #FILGRASTIM
      dplyr::filter(!(nci_t == "C88325" & molecule_chembl_id == "CHEMBL2109653")) |> #BERMEKIMAB
      dplyr::filter(!(nci_t == "C2602" & molecule_chembl_id == "CHEMBL513")) |> #CARMUSTINE
      dplyr::filter(!(nci_t == "C80867" & molecule_chembl_id == "CHEMBL491473")) |> #CEDIRANIB
      dplyr::filter(!(nci_t == "C225" & molecule_chembl_id == "CHEMBL1201577")) |> #CETUXIMAB
      dplyr::filter(!(nci_t == "C2213" & molecule_chembl_id == "CHEMBL178")) |> #DAUNORUBICIN
      dplyr::filter(!(nci_t == "C128039" & molecule_chembl_id == "CHEMBL3137331")) |> #DEFACTINIB
      dplyr::filter(!(nci_t == "C62435" & molecule_chembl_id == "CHEMBL522892")) |> #DOVITINIB
      dplyr::filter(!(nci_t == "C2693" & molecule_chembl_id == "CHEMBL553")) |> #ERLOTINIB
      dplyr::filter(!(nci_t == "C478" & molecule_chembl_id == "CHEMBL135")) |> #ESTRADIOL
      dplyr::filter(!(nci_t == "C1687" & molecule_chembl_id == "CHEMBL941")) |> #IMATINIB
      dplyr::filter(!(nci_t == "C29165" & molecule_chembl_id == "CHEMBL191")) |> #LOSARTAN
      dplyr::filter(!(nci_t == "C1561" & molecule_chembl_id == "CHEMBL2109447")) |> #LYM-1
      dplyr::filter(!(nci_t == "C1155" & molecule_chembl_id == "CHEMBL717")) |> #MEDROXYPROGESTERONE ACETATE
      dplyr::filter(!(nci_t == "C9678" & molecule_chembl_id == "CHEMBL1456")) |> #MYCOPHENOLATE MOFETIL
      dplyr::filter(!(nci_t == "C20513" & molecule_chembl_id == "CHEMBL1201573")) |> #OPRELVEKIN
      dplyr::filter(!(nci_t == "C95230" & molecule_chembl_id == "CHEMBL1201421")) |> #PEGAPTANIB SODIUM
      dplyr::filter(!(nci_t == "C176878" & molecule_chembl_id == "CHEMBL3545154")) |> #POZIOTINIB
      dplyr::filter(!(nci_t == "C2297" & molecule_chembl_id == "CHEMBL103")) |> #PROGESTERONE
      dplyr::filter(!(nci_t == "C148170" & molecule_chembl_id == "CHEMBL225071")) |> #RALTITREXED
      dplyr::filter(!(nci_t == "C66506" & molecule_chembl_id == "CHEMBL1790041")) |> #RANITIDINE
      dplyr::filter(!(nci_t == "C82693" & molecule_chembl_id == "CHEMBL1738757")) |> #REBASTINIB
      dplyr::filter(!(nci_t == "C1492" & molecule_chembl_id == "CHEMBL1201670")) |> #SARGRAMOSTIM
      dplyr::filter(!(nci_t == "C88337" & molecule_chembl_id == "CHEMBL2105737")) |> #SONIDEGIB
      dplyr::filter(!(nci_t == "C80631" & molecule_chembl_id == "CHEMBL565612")) |> #SOTRASTAURIN
      dplyr::filter(!(nci_t == "C104057" & molecule_chembl_id == "CHEMBL2105694")) |> #TELAPRISTONE ACETATE
      dplyr::filter(!(nci_t == "C29523" & molecule_chembl_id == "CHEMBL1201334")) |> #TRIPTORELIN
      dplyr::filter(!(nci_t == "C80049" & molecule_chembl_id == "CHEMBL3545218")) |> #VORUCICLIB
      dplyr::filter(!(nci_t == "C133021" & molecule_chembl_id == "CHEMBL3188386")) #WNT-974



    i <- 1
    nci_compounds_no_chembl <- as.data.frame(
      nci_antineo_thesaurus |>
      dplyr::select(nci_cd_name, molecule_chembl_id) |>
      dplyr::distinct() |>
      dplyr::filter(!stringr::str_detect(tolower(nci_cd_name),
                                         "( vaccine)|^[0-9]")) |>
      dplyr::filter(is.na(molecule_chembl_id)) |>
      dplyr::select(nci_cd_name) |>
      dplyr::mutate(num_spaces = stringr::str_count(nci_cd_name," ")) |>
      dplyr::filter(num_spaces <= 1 & !stringr::str_detect(tolower(nci_cd_name),"regimen|&|/|;|,")) |>
      dplyr::filter(stringr::str_detect(nci_cd_name,"(mab|cin|ide|ib|im|bine|tin|om|lin|stat|one|ate|ole|ane|ine|xel|rol)$")) |>
      dplyr::distinct()
    )



    nci_compounds_chembl_match <- data.frame()

    ## Retrieve aliases for drugs with PubChem x-refs
    pubchem_synonym_files <-
      sort(list.files(path = file.path(here::here(), "data-raw","pubchem"),
                      pattern = "CID-Synonym-filtered_",
                      full.names = T))

    lgr::lgr$info("Mapping ChEMBL identifiers for NCI compounds")
    i <- 1
    for(f in pubchem_synonym_files){
      lgr::lgr$info("Mapping iteration..", i)
      synonym_data <- as.data.frame(readr::read_tsv(
        f, col_names = c('pubchem_cid','alias'),
        col_types = "dc",
        progress = F
      ))

      chembl2pubchem <- synonym_data |>
        dplyr::filter(stringr::str_detect(
          alias, "^CHEMBL")) |>
        dplyr::rename(molecule_chembl_id = alias) |>
        dplyr::mutate(pubchem_cid = as.integer(pubchem_cid))

      chembl2alias <- synonym_data |>
        dplyr::filter(!stringr::str_detect(
          alias, "^CHEMBL")) |>
        dplyr::mutate(pubchem_cid = as.integer(pubchem_cid)) |>
        dplyr::mutate(alias = tolower(alias))

      hits <- nci_compounds_no_chembl |>
        dplyr::mutate(nci_cd_name_lc =
                        tolower(nci_cd_name)) |>
        dplyr::inner_join(
          chembl2alias, by = c("nci_cd_name_lc" = "alias"))


      rm(chembl2alias)

      if(nrow(hits) > 0){
        hits <- hits |>
          dplyr::inner_join(
            chembl2pubchem, by = "pubchem_cid") |>
          dplyr::select(nci_cd_name,
                        molecule_chembl_id)

        lgr::lgr$info(paste0("Found ", nrow(hits), " ChEMBL identifiers"))

        nci_compounds_chembl_match <- nci_compounds_chembl_match |>
          dplyr::bind_rows(hits)
      }
      rm(chembl2pubchem)

      i <- i + 1
    }

    nci_compounds_chembl_match_unique <- nci_compounds_chembl_match |>
      dplyr::group_by(nci_cd_name) |>
      dplyr::summarise(
        n_identifiers = dplyr::n(),
        molecule_chembl_id = paste(unique(molecule_chembl_id), collapse="&"),
        .groups = "drop") |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Goserelin Acetate",
        "CHEMBL1200501",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Abacavir Sulfate",
        "CHEMBL1200666",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Roniciclib",
        "CHEMBL4442620",
        as.character(molecule_chembl_id))) |>
      dplyr::filter(!stringr::str_detect(molecule_chembl_id,"&")) |>
      dplyr::select(-n_identifiers)


    # i <- 1
    # while(i <= nrow(nci_compounds_no_chembl)){
    #   name <- nci_compounds_no_chembl[i, "nci_cd_name"]
    #   chembl_hit <- get_chembl_compound_by_name(name)
    #   if(!is.null(chembl_hit)){
    #     chembl_hit <- chembl_hit |>
    #       dplyr::rename(nci_cd_name = name)
    #     nci_compounds_chembl_match <- nci_compounds_chembl_match |>
    #       dplyr::bind_rows(chembl_hit)
    #
    #   }
    #   if(i %% 10 == 0){
    #     lgr::lgr$info("Done with querying ChEMBL for ", i, " compound names")
    #   }
    #   i <- i + 1
    # }

    nci_antineo_thesaurus_chembl <- nci_antineo_thesaurus |>
      dplyr::anti_join(nci_compounds_no_chembl,
                       by = "nci_cd_name")

    nci_antineo_thesaurus_no_chembl <- nci_antineo_thesaurus |>
      dplyr::inner_join(nci_compounds_no_chembl,
                        by = "nci_cd_name") |>
      dplyr::select(-c(num_spaces, molecule_chembl_id)) |>
      dplyr::left_join(nci_compounds_chembl_match_unique,
                       by = "nci_cd_name")

    nci_antineo_thesaurus <-
      nci_antineo_thesaurus_chembl |>
      dplyr::bind_rows(nci_antineo_thesaurus_no_chembl) |>
      dplyr::arrange(nci_cd_name) |>
      dplyr::distinct() |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Vorolanib",
        "CHEMBL3545427",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Gemtuzumab Ozogamicin",
        "CHEMBL1201506",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Ublituximab",
        "CHEMBL2108354",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Etoposide",
        "CHEMBL44657",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Lorvotuzumab Mertansine",
        "CHEMBL1743037",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Epratuzumab",
        "CHEMBL2108404",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Ocaratuzumab",
        "CHEMBL2109665",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Medroxyprogesterone Acetate",
        "CHEMBL717",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Brigatinib",
        "CHEMBL3545311",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Brentuximab",
        "CHEMBL1742994",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Afatinib",
        "CHEMBL1173655",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Spebrutinib",
        "CHEMBL3301625",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Dacomitinib",
        "CHEMBL2110732",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Lucitanib",
        "CHEMBL2220486",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Tofacitinib",
        "CHEMBL221959",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Mifamurtide",
        "CHEMBL2111100",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Umbralisib",
        "CHEMBL3948730",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Onvansertib",
        "CHEMBL1738758",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Sulindac",
        "CHEMBL15770",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_cd_name == "Calcipotriene",
        "CHEMBL1200666",
        as.character(molecule_chembl_id))) |>
      dplyr::distinct()


    nci_drugs <- list()
    ## NCI anticancer drugs (targeted) - including compound identifier (CHEMBL)
    nci_drugs[['with_chembl_id']] <- nci_antineo_thesaurus |>
      dplyr::select(nci_t,
                    nci_concept_definition,
                    nci_cd_name,
                    molecule_chembl_id,
                    nci_drug_name,
                    nci_concept_synonym_all) |>
      dplyr::filter(!is.na(molecule_chembl_id)) |>
      #dplyr::mutate(nci_drug_name_lc = tolower(nci_drug_name)) |>
      dplyr::distinct()

    ## NCI anticancer drugs (non-targeted) - lacking compound identifier (CHEMBL)
    nci_drugs[['no_chembl_id']] <- nci_antineo_thesaurus |>
      dplyr::filter(is.na(molecule_chembl_id)) |>
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
      "nci_treatment_thesaurus_antineo.rds"))

  }else{
    nci_drugs <- readRDS(
      file = file.path(
        path_data_processed, "nci_thesaurus",
        "nci_treatment_thesaurus_antineo.rds"))
  }
  return(nci_drugs)
}


get_compound_properties <- function(molecule_chembl_id = NULL,
                                    pchem_cid = NA,
                                    version_chembl = "ChEMBL_29",
                                    version_dgidb = "4.0.2"){

  compound_properties <- list()

  #compound_info_dgidb <- NULL
  compound_moa_chembl <- NULL
  compound_properties_pubchem <- NULL
  compound_properties_chembl <- NULL
  compound_indications_chembl <- NULL

  if(is.null(molecule_chembl_id)){
    return(NULL)
  }

  if(!is.na(molecule_chembl_id) & stringr::str_detect(string = molecule_chembl_id, pattern = "^CHEMBL")){

    ## fetch compound indications (EFO) from ChEMBL
    #compound_indications_chembl <- chembl_compound_indications(molecule_chembl_id = molecule_chembl_id, min_maxphase = 3) |>
      #dplyr::distinct()

    ## fetch compound MOA (mechanism of action) from ChEMBL - ignore compounds with no recorded target
    #compound_moa_chembl <- chembl_compound_moa(molecule_chembl_id = molecule_chembl_id) |>
      #dplyr::filter(!is.na(target_chembl_id))

    ## fetch compound properties from ChEMBL
    compound_properties_chembl <- chembl_compound_properties(molecule_chembl_id = molecule_chembl_id)
  }
  if(!is.na(pchem_cid)){
    ## fetch compound properties from PubChem
    compound_properties_pubchem <- pubchem_compound_properties(pubchem_cid = pchem_cid)
  }

  compound_properties$molecule_chembl_id <- molecule_chembl_id
  compound_properties$chembl_db_version <- version_chembl
  compound_properties$pubchem_cid <- pchem_cid
  compound_properties$indications_chembl <- compound_indications_chembl
  compound_properties$properties_pubchem <- compound_properties_pubchem
  compound_properties$properties_chembl <- compound_properties_chembl
  compound_properties$moa_chembl <- compound_moa_chembl

  return(compound_properties)

}


# get_dailymed_drug_indications <- function(update = F,
#                                           path_data_raw = NULL){
#
#   existing_dataset_fname <-
#     paste0(path_data_raw,"/dailymed/dm_cancerdrug_indications.rds")
#   if(update == F){
#     if(file.exists(existing_dataset_fname)){
#       dm_drugs_with_indications <- readRDS(file = existing_dataset_fname)
#       return(dm_drugs_with_indications)
#     }
#   }
#
#   cancer_drug_names <-
#     unique(toupper(oncoPharmaDB::oncopharmadb$nci_cd_name))
#
#   all_drug_indications <- data.frame()
#   for(i in 16:21){
#     year = paste0("20", i)
#     if(i < 10){
#       year = paste0("200",i)
#     }
#
#     index <- 1
#     for(f in list.files(
#       path = paste0("data-raw/dailymed/prescription/", year),
#       pattern = ".xml", full.names = T)){
#
#       index = index + 1
#
#       #cat(f,'\n')
#       if(index %% 50 == 0){
#         cat(year,' -- processed: ', index, ' indication records\n')
#       }
#       drug_indications <-
#         get_drug_indications_dm(xml_fname = f,
#                                 cancer_drug_names = cancer_drug_names)
#       if(nrow(drug_indications) > 0){
#         all_drug_indications <- all_drug_indications |>
#           dplyr::bind_rows(drug_indications)
#       }
#     }
#     cat(year,' --- processed: ', index, ' indication records\n')
#   }
#
#
#   comb_regimen_indications <- all_drug_indications |>
#     dplyr::filter(cancer_indication == T) |>
#     dplyr::mutate(comb_regimen_indication = dplyr::if_else(
#       !is.na(indication_criteria) &
#       stringr::str_detect(tolower(indication_criteria), "in combination with") &
#         !stringr::str_detect(tolower(indication_criteria), "alone (or|and) in combination"),
#       TRUE, FALSE, FALSE
#     )) |>
#     dplyr::select(drugname_generic, indication, indication_criteria,
#                   comb_regimen_indication, xml_fname) |>
#     dplyr::filter(comb_regimen_indication == T) |>
#     dplyr::distinct()
#
#
#   indications_otm_input1 <- all_drug_indications |>
#     dplyr::filter(cancer_indication == T) |>
#     dplyr::select(indication,
#                   xml_fname) |>
#     dplyr::rename(value = indication, nct_id = xml_fname) |>
#     dplyr::mutate(property = "condition") |>
#     dplyr::filter(!is.na(value))
#
#   indications_otm_input2 <- all_drug_indications |>
#     dplyr::filter(cancer_indication == T) |>
#     dplyr::select(indication_criteria,
#                   xml_fname) |>
#     dplyr::rename(value = indication_criteria,
#                   nct_id = xml_fname) |>
#     dplyr::mutate(property = "condition") |>
#     dplyr::filter(!is.na(value))
#
#   indications_otm_input <- indications_otm_input1 |>
#     dplyr::bind_rows(indications_otm_input2)
#
#   indications_mapped <- oncoTrialMiner::index_phenotype(
#     ct_data_raw = indications_otm_input
#   )
#
#   tmp <- indications_mapped$per_xref |>
#     dplyr::select(nct_id, cui, primary_site, efo_id,
#                   cui_name, ct_condition_raw) |>
#     dplyr::rename(xml_fname = nct_id,
#                   indication = ct_condition_raw) |>
#     dplyr::distinct()
#
#   dm_drugs_with_mapped_indications1 <- all_drug_indications |>
#     dplyr::select(xml_fname,
#                   drugname_generic,
#                   drugname_trade,
#                   indication) |>
#     dplyr::filter(!is.na(indication)) |>
#     dplyr::inner_join(tmp, by = c("xml_fname", "indication")) |>
#     dplyr::distinct()
#
#   dm_drugs_with_mapped_indications2 <- all_drug_indications |>
#     dplyr::select(xml_fname,
#                   drugname_generic,
#                   drugname_trade,
#                   indication_criteria) |>
#     dplyr::rename(indication = indication_criteria) |>
#     dplyr::filter(!is.na(indication)) |>
#     dplyr::inner_join(tmp, by = c("xml_fname", "indication")) |>
#     #dplyr::select(-c(xml_fname, indication)) |>
#     dplyr::distinct()
#
#   dm_drugs_with_indications <- as.data.frame(
#     dm_drugs_with_mapped_indications2 |>
#       dplyr::bind_rows(dm_drugs_with_mapped_indications1) |>
#
#       dplyr::left_join(dplyr::select(comb_regimen_indications,
#                                      drugname_generic,
#                                      indication_criteria,
#                                      comb_regimen_indication,
#                                      xml_fname),
#                        by = c("xml_fname" = "xml_fname",
#                               "indication" = "indication_criteria",
#                               "drugname_generic" = "drugname_generic")) |>
#       dplyr::left_join(dplyr::select(comb_regimen_indications,
#                                      drugname_generic,
#                                      indication,
#                                      comb_regimen_indication,
#                                      xml_fname),
#                        by = c("xml_fname" = "xml_fname",
#                               "indication" = "indication",
#                               "drugname_generic" = "drugname_generic")) |>
#       dplyr::mutate(comb_regimen_indication = dplyr::if_else(
#         !is.na(comb_regimen_indication.x) & comb_regimen_indication.x == TRUE,
#         TRUE, as.logical(FALSE)
#       )) |>
#       dplyr::mutate(comb_regimen_indication = dplyr::if_else(
#         !is.na(comb_regimen_indication.y) &
#           comb_regimen_indication.y == TRUE,
#         TRUE, as.logical(comb_regimen_indication)
#       )) |>
#       dplyr::select(-c(comb_regimen_indication.x,
#                        comb_regimen_indication.y)) |>
#       dplyr::distinct() |>
#       dplyr::mutate(nci_drug_name = tolower(drugname_generic)) |>
#       dplyr::select(-drugname_generic) |>
#       dplyr::group_by(nci_drug_name, cui, efo_id,
#                       cui_name, primary_site,
#                       comb_regimen_indication) |>
#       dplyr::summarise(drugname_trade =
#                          paste(unique(drugname_trade), collapse="|"),
#                        drug_clinical_id = paste(unique(
#                          stringr::str_replace(
#                            xml_fname,"\\.xml$",""
#                          )), collapse=","),
#                        .groups = "drop") |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "renal cell carcinoma"),
#         "EFO:0000681", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "urothelial carcinoma"),
#         "EFO:0008528", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "mantle cell lymphoma"),
#         "EFO:1001469", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "non-small cell lung cancer"),
#         "EFO:0003060", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "gastrointestinal stromal"),
#         "Orphanet:44890", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "mediastinal large b-cell lymphoma"),
#         "MONDO:0004021", as.character(efo_id)
#       )) |>
#       dplyr::mutate(efo_id = dplyr::if_else(
#         (is.na(efo_id) | efo_id == "NA") &
#           stringr::str_detect(tolower(cui_name), "pancreatic neuroendocrine tumor"),
#         "EFO:1000045", as.character(efo_id)
#       )) |>
#       dplyr::filter(!is.na(cui) & efo_id != "NA")
#   )
#
#   dm_drugs_with_indications <- as.data.frame(
#     dm_drugs_with_indications |>
#       tidyr::separate_rows(efo_id,sep="\\|") |>
#       dplyr::left_join(dplyr::select(
#         oncoPhenoMap::auxiliary_maps$efo$efo2name, efo_id, efo_name
#       ), by = "efo_id") |>
#       dplyr::rename(disease_efo_id = efo_id,
#                     disease_efo_label = efo_name) |>
#       #dplyr::mutate(nci_drug_name = tolower(drugname_generic)) |>
#       dplyr::mutate(drug_max_phase_indication = NA,
#                     drug_clinical_source = "DailyMedParseSN",
#                     #drug_clinical_id = NA,
#                     drug_max_ct_phase = NA,
#                     drug_approved_indication = TRUE) |>
#       dplyr::distinct()
#   )
#
#   saveRDS(dm_drugs_with_indications, existing_dataset_fname)
#
#   return(dm_drugs_with_indications)
# }
#
#
#
#
# get_drug_indications_dm <- function(xml_fname = NULL,
#                                     cancer_drug_names = NULL){
#
#   invisible(assertthat::assert_that(
#     file.exists(xml_fname),
#     msg = paste0("File name ", xml_fname," does not exist")))
#   con <- file(description = xml_fname)
#   tmp <- xml2::read_xml(con, options = c("HUGE","NOBLANKS"))
#   xml2::xml_ns_strip(tmp)
#
#
#   drugname_trade <-  stringr::str_trim(xml2::xml_text(xml2::xml_find_all(
#     tmp, "//manufacturedProduct/name")[1]))
#   if(length(drugname_trade) == 0){
#     drugname_trade <- NA
#   }
#
#   drug_active_moiety <- xml2::xml_text(xml2::xml_find_all(
#     tmp, "//manufacturedProduct/ingredient/ingredientSubstance/activeMoiety/activeMoiety/name")[1])
#   if(length(drug_active_moiety) == 0){
#     drug_active_moiety <- NA
#   }
#
#   initial_approval_date <- NA
#   full_title <- ""
#   title_elements <- xml2::as_list(xml2::xml_find_first(tmp,".//title"))
#
#   if(is.list(title_elements)){
#     if(length(title_elements) > 0){
#       for(i in 1:length(title_elements)){
#         title_data <- stringr::str_squish(stringr::str_trim(unlist(title_elements[[i]])))
#         if(is.list(title_data)){
#           cat('BALLE\n')
#         }
#         if(is.null(title_data)){
#           next
#         }
#         if(!is.character(title_data)){
#           next
#         }
#         full_title <- paste0(full_title, title_data)
#       }
#     }
#   }
#
#   #cat(full_title, '\n')
#   if(full_title != ""){
#     t <- stringr::str_match(full_title,"Initial U.S. Approval: [0-9]{4}")
#     if(length(t) == 1){
#       initial_approval_date <- stringr::str_trim(
#         stringr::str_replace(
#           t, "Initial U.S. Approval: ",""
#         ))
#     }
#   }
#
#   drugname_generic <- xml2::xml_text(xml2::xml_find_all(
#     tmp, "//manufacturedProduct/asEntityWithGeneric/genericMedicine/name")[1])
#   if(length(drugname_generic) == 0){
#     drugname_generic <- NA
#   }
#   sections <- xml2::xml_find_all(tmp,".//component/section")
#
#   indication_section <- 0
#   section_counter <- 0
#   all_indication_entries <- data.frame()
#   all_main_indication_entries <- data.frame()
#   effective_time <- NA
#
#   is_cancer_drug <- F
#   if(!is.null(cancer_drug_names) &
#      !is.na(drugname_trade)){
#     if(toupper(drugname_trade) %in% cancer_drug_names){
#       is_cancer_drug <- T
#     }
#   }
#
#   if(!is.null(cancer_drug_names) &
#      !is.na(drugname_generic)){
#     if(toupper(drugname_generic) %in% cancer_drug_names){
#       is_cancer_drug <- T
#     }
#   }
#
#   if(!is.null(cancer_drug_names) &
#      !is.na(drug_active_moiety)){
#     if(toupper(drug_active_moiety) %in% cancer_drug_names){
#       is_cancer_drug <- T
#     }
#   }
#
#   if(is_cancer_drug == F){
#     return(all_indication_entries)
#   }
#
#   #cat(drugname_generic, '\n')
#
#   for(s in sections){
#     section_counter <- section_counter + 1
#     section_titles <- xml2::xml_text(xml2::xml_find_first(s, ".//title"))
#     for(sec_title in unique(section_titles)){
#       sec_title <- stringr::str_squish(stringr::str_trim(sec_title))
#       if(!is.na(sec_title) & nchar(sec_title) > 0){
#
#         if(stringr::str_detect(sec_title,"INDICATIONS") &
#            stringr::str_detect(sec_title,"USAGE")){
#           indication_section <- indication_section + 1
#
#           indication_info_main <- xml2::as_list(sections[section_counter][[1]])
#           effective_time <- attributes(
#             indication_info_main$effectiveTime)$value
#
#           if("excerpt" %in% names(indication_info_main)){
#             if("highlight" %in% names(indication_info_main$excerpt)){
#               if("text" %in% names(indication_info_main$excerpt$highlight)){
#
#                 all_indication_contents <-
#                   unlist(indication_info_main$excerpt$highlight$text)
#
#                 names_indication_contents <- names(all_indication_contents)
#                 names(all_indication_contents) <- NULL
#
#                 indication <- ""
#                 indication_criteria <- NA
#                 prev_cname <- ""
#                 if(length(all_indication_contents) > 0 &
#                    length(names_indication_contents) > 0 &
#                   length(all_indication_contents) == length(names_indication_contents)){
#                   for(i in 1:length(all_indication_contents)){
#                     cname <- names_indication_contents[i]
#
#                     df <- data.frame('cname' = cname,
#                                      'indication' = all_indication_contents[i],
#                                      'indication_criteria' = NA,
#                                      'drugname_generic' = drugname_generic,
#                                      'drugname_trade' = drugname_trade,
#                                      'drug_active_moiety' = drug_active_moiety,
#                                      'effective_time' = effective_time,
#                                      stringsAsFactors = F)
#
#                     all_main_indication_entries <- all_main_indication_entries |>
#                       dplyr::bind_rows(df)
#                   }
#                 }
#               }
#             }
#           }else{
#             if("text" %in% names(indication_info_main)){
#
#               all_indication_contents <-
#                 unlist(indication_info_main$text)
#
#               if(!is.null(all_indication_contents)){
#                 names_indication_contents <- names(all_indication_contents)
#                 names(all_indication_contents) <- NULL
#                 if(length(all_indication_contents) > 0 &
#                    length(names_indication_contents) > 0 &
#                    length(all_indication_contents) == length(names_indication_contents)){
#                   for(i in 1:length(all_indication_contents)){
#                     cname <- names_indication_contents[i]
#
#                     df <- data.frame('cname' = cname,
#                                      'indication' = all_indication_contents[i],
#                                      'indication_criteria' = NA,
#                                      'drugname_generic' = drugname_generic,
#                                      'drugname_trade' = drugname_trade,
#                                      'drug_active_moiety' = drug_active_moiety,
#                                      'effective_time' = effective_time,
#                                      stringsAsFactors = F)
#
#                     all_main_indication_entries <- all_main_indication_entries |>
#                       dplyr::bind_rows(df)
#                   }
#                 }
#               }
#             }
#           }
#         }
#         if(indication_section > 0 &
#            !(stringr::str_detect(
#              sec_title, "INDICATIONS")) &
#            !stringr::str_detect(
#              sec_title, "DOSAGE"
#            )){
#
#           df_sec_title <- data.frame('cname' = 'section_title',
#                                      'indication' = sec_title,
#                                      'indication_criteria' = NA,
#                                      'drugname_generic' = drugname_generic,
#                                      'drugname_trade' = drugname_trade,
#                                      'drug_active_moiety' = drug_active_moiety,
#                                      'effective_time' = effective_time,
#                                      stringsAsFactors = F)
#
#
#           section_indication_details_found <- 0
#
#           indication_details <- xml2::as_list(sections[section_counter][[1]])
#
#           if("text" %in% names(indication_details)){
#
#             if(length(indication_details$text) == 1){
#               if("paragraph" %in% names(indication_details$text)){
#
#                 ind_details <- unlist(indication_details$text$paragraph)
#                 attributes(ind_details) <- NULL
#                 ind_details <- ind_details[!nchar(ind_details) == 1]
#                 ind_details <- stringr::str_squish(paste(ind_details, collapse=""))
#
#                 df <- data.frame('cname' = 'paragraph',
#                                  'indication' = sec_title,
#                                  'indication_criteria' = ind_details,
#                                  'drugname_generic' = drugname_generic,
#                                  'drugname_trade' = drugname_trade,
#                                  'drug_active_moiety' = drug_active_moiety,
#                                  'effective_time' = effective_time,
#                                  stringsAsFactors = F)
#
#                 all_indication_entries <- all_indication_entries |>
#                   dplyr::bind_rows(df)
#
#                 section_indication_details_found <- 1
#
#               }
#             }else{
#               if(!("list" %in% names(indication_details$text))){
#                 if(length(indication_details$text) == 0){
#                   next
#                 }
#                 for(n in 1:length(indication_details$text)){
#                   if("content" %in% names(indication_details$text[[n]])){
#
#                     if(is.list(indication_details$text[[n]]$content) &
#                        length(indication_details$text[[n]]$content) == 1){
#
#                       if(is.character(indication_details$text[[n]]$content[[1]])){
#
#                         df <- data.frame(
#                           'cname' = 'section_title',
#                           'indication' = sec_title,
#                           'indication_criteria' = indication_details$text[[n]]$content[[1]],
#                           'drugname_generic' = drugname_generic,
#                           'drugname_trade' = drugname_trade,
#                           'drug_active_moiety' = drug_active_moiety,
#                           'effective_time' = effective_time,
#                           stringsAsFactors = F)
#
#                         all_indication_entries <- all_indication_entries |>
#                           dplyr::bind_rows(df)
#                       }
#
#                       section_indication_details_found <- 1
#
#                     }else{
#                       crit <- unlist(indication_details$text[[n]])
#                       if(is.character(crit) & length(crit) == 1){
#                         names(crit) <- NULL
#                         if(nchar(crit[1]) > 1 &
#                            !(stringr::str_detect(
#                              crit[1],"Clinical Studies|Dosage and Administration"))){
#                           df <- data.frame(
#                             'cname' = 'section_title',
#                             'indication' = sec_title,
#                             'indication_criteria' = crit[1],
#                             'drugname_generic' = drugname_generic,
#                             'drugname_trade' = drugname_trade,
#                             'drug_active_moiety' = drug_active_moiety,
#                             'effective_time' = effective_time,
#                             stringsAsFactors = F)
#
#                           all_indication_entries <- all_indication_entries |>
#                             dplyr::bind_rows(df)
#
#                           section_indication_details_found <- 1
#
#                         }
#                       }
#                     }
#                   }else{
#
#                     if(is.list(indication_details$text[[n]]) & length(indication_details$text[[n]]) > 0){
#                       #cat(section_counter, '--', n, '---', length(indication_details$text[[n]]), '\n')
#                       for(i in 1:length(indication_details$text[[n]])){
#                         if(is.character(indication_details$text[[n]][[i]])){
#                           df <- data.frame(
#                             'cname' = 'section_title',
#                             'indication' = sec_title,
#                             'indication_criteria' = indication_details$text[[n]][[i]],
#                             'drugname_generic' = drugname_generic,
#                             'drugname_trade' = drugname_trade,
#                             'drug_active_moiety' = drug_active_moiety,
#                             'effective_time' = effective_time,
#                             stringsAsFactors = F)
#
#                           all_indication_entries <- all_indication_entries |>
#                             dplyr::bind_rows(df)
#
#                           section_indication_details_found <- 1
#
#                         }
#                       }
#                     }
#                   }
#                 }
#               }
#             }
#
#             if("paragraph" %in% names(indication_details$text) &
#                length(indication_details$text) > 1){
#
#               indication_title <- ""
#               indication_title_elements <- c()
#               if(length(indication_details$text$paragraph) > 1){
#                 for(j in 1:length(indication_details$text$paragraph)){
#                   if(is.character(indication_details$text$paragraph[[j]])){
#                     indication_title_elements <-
#                       c(indication_title_elements,
#                         indication_details$text$paragraph[[j]])
#                   }
#                 }
#                 indication_title <- paste(indication_title_elements, collapse="")
#               }else{
#                 if(length(indication_details$text$paragraph) == 1){
#                   if(is.character(indication_details$text$paragraph[[1]])){
#                     indication_title <- as.character(indication_details$text$paragraph[[1]])
#                   }
#                 }
#               }
#               indication_criterions <- c()
#
#               if("list" %in% names(indication_details$text)){
#
#                 if("item" %in% names(indication_details$text$list)){
#
#                   if(length(indication_details$text$list) == 1){
#                     indication_items <-
#                       unlist(indication_details$text$list$item)
#
#                     if(!is.null(indication_items) &
#                        is.vector(indication_items)){
#                       ind_crit <- ""
#                       #cat('--', unlist(names(indication_items)),'\n')
#
#                       if(length(indication_items) > 0){
#                         for(m in 1:length(indication_items)){
#                           if(is.null(names(indication_items[m]))){
#                             next
#                           }
#                           if(names(indication_items[m]) == "" |
#                              stringr::str_detect(names(indication_items[m]),"content([0-9])?$")){
#                             ind_crit <- paste(ind_crit, stringr::str_trim(indication_items[m]))
#                           }
#                         }
#                       }
#                       indication_criterions <- c(indication_criterions,
#                                                  paste(indication_title,
#                                                        ind_crit, sep=" "))
#                     }
#                   }else{
#
#                     for(m in 1:length(indication_details$text$list)){
#                       #if(nchar(unlist(indication_details$text$list[[m]])) > 1){
#                       indication_criterions <-
#                         c(indication_criterions,
#                           paste(indication_title,
#                                 unlist(indication_details$text$list[[m]]),
#                                 sep = " - ")
#                         )
#                       #}
#                     }
#                   }
#                 }
#                 else{
#                   for(j in 1:length(indication_details$text$list)){
#
#                     if(is.character(unlist(indication_details$text$list[[j]][[1]]))){
#                       if(nchar(unlist(indication_details$text$list[[j]][[1]])) > 1){
#                         indication_criterions <-
#                           c(indication_criterions,
#                             paste(
#                               indication_title,
#                               unlist(indication_details$text$list[[j]][[1]]),
#                               sep = " - ")
#                           )
#                       }
#                     }
#                   }
#                 }
#               }
#
#               for(c in indication_criterions){
#
#                 df <- data.frame(
#                   'cname' = 'section_title',
#                   'indication' = sec_title,
#                   'indication_criteria' = c,
#                   'drugname_generic' = drugname_generic,
#                   'drugname_trade' = drugname_trade,
#                   'drug_active_moiety' = drug_active_moiety,
#                   'effective_time' = effective_time,
#                   stringsAsFactors = F)
#
#                 all_indication_entries <- all_indication_entries |>
#                   dplyr::bind_rows(df)
#
#                 section_indication_details_found <- 1
#
#               }
#             }
#
#             else{
#
#               if(length(indication_details$text) >= 1){
#                 type <- unique(names(indication_details$text))
#                 if(length(type) != 1){
#                   break
#                 }
#                 if(unique(names(indication_details$text)) == "list"){
#                   all_criterions <- c()
#
#                   for(j in 1:length(indication_details$text$list)){
#                     if(length(indication_details$text$list[[j]]) == 0){
#                       next
#                     }
#                     for(k in 1:length(indication_details$text$list[[j]])){
#                       indication_criterion <- ""
#                       if(is.list(indication_details$text$list[[j]][[k]])){
#
#                         if(length(indication_details$text$list[[j]][[k]]) == 1){
#                           if(is.character(unlist(indication_details$text$list[[j]][[k]])) &
#                              nchar(unlist(indication_details$text$list[[j]][[k]])) > 1){
#                             indication_criterion <- as.character(
#                               unlist(indication_details$text$list[[j]][[k]]))
#
#                             all_criterions <- c(all_criterions, indication_criterion)
#                           }
#                         }else{
#                           indication_criterion <- ""
#                           if(length(indication_details$text$list[[j]][[k]]) == 0){
#                             next
#                           }
#                           for(l in 1:length(indication_details$text$list[[j]][[k]])){
#                             if(is.character(indication_details$text$list[[j]][[k]][[l]])){
#                               indication_criterion <- paste0(
#                                 indication_criterion,
#                                 as.character(
#                                   indication_details$text$list[[j]][[k]][[l]]))
#                             }
#                           }
#                           all_criterions <- c(all_criterions, indication_criterion)
#
#
#                         }
#                       }else{
#                         indication_criterion <-
#                           unlist(indication_details$text$list[[j]][[k]])
#                         all_criterions <- c(all_criterions, indication_criterion)
#
#                       }
#                     }
#                   }
#
#                   for(c in all_criterions){
#                     df <- data.frame(
#                       'cname' = 'section_list',
#                       'indication' = sec_title,
#                       'indication_criteria' = c,
#                       'drugname_generic' = drugname_generic,
#                       'drugname_trade' = drugname_trade,
#                       'drug_active_moiety' = drug_active_moiety,
#                       'effective_time' = effective_time,
#                       stringsAsFactors = F)
#
#                     all_indication_entries <- all_indication_entries |>
#                       dplyr::bind_rows(df)
#
#                     section_indication_details_found <- 1
#
#                   }
#                 }
#               }
#             }
#           }
#           if(section_indication_details_found == 0){
#             all_indication_entries <- all_indication_entries |>
#               dplyr::bind_rows(df_sec_title)
#           }
#         }
#
#         if(stringr::str_detect(sec_title,"CONTRAINDICATIONS") |
#            stringr::str_detect(sec_title,"WARNINGS|DOSAGE")){
#           indication_section <- 0
#           break
#         }
#       }
#     }
#   }
#
#
#
#   if(nrow(all_indication_entries) == 0 &
#      nrow(all_main_indication_entries) > 0){
#     all_indication_entries <- all_indication_entries |>
#       dplyr::bind_rows(all_main_indication_entries)
#   }
#
#   if(nrow(all_indication_entries) == 0){
#     return(all_indication_entries)
#   }
#
#   #cancer_regex <- "cancer|carcinoma|wilms|esophageal|melanoma|myeloma"
#   #cancer_regex <- paste0(cancer_regex,"|leukemia|lymphom|sarcoma|thymom|thyroid cancer")
#   #cancer_regex <- paste0(cancer_regex,"|head and neck cancer|")
#
#   cancer_tissue_regex <- "endometrial|prostate|rectal|myeloid|anal|colorectal|wilms|kidney|ovarian|"
#   cancer_tissue_regex <- paste0(
#     cancer_tissue_regex,"thyroid|esophag|cervical|cervix|gastric|stomach|head and neck|ovary|")
#   cancer_tissue_regex <- paste0(
#    cancer_tissue_regex,"bone|testis|testicular|lung|liver|skin|pancreas|pancreatic|brain|bladder|")
#   cancer_tissue_regex <- paste0(
#     cancer_tissue_regex,"adrenal gland|renal cell|urothelial|lympho|hepatocellular|breast")
#
#   cancer_other_regex <- "melanoma|myeloma|glioblastoma|sarcoma|leukemia|lymphoma|"
#   cancer_other_regex <- paste0(
#     cancer_other_regex, "cholangiocarcinom|mesotheliom|glioma|neuroblastoma|tumor")
#
#   rownames(all_indication_entries) <- NULL
#   all_indication_entries <- all_indication_entries |>
#     dplyr::mutate(xml_fname = basename(xml_fname)) |>
#     dplyr::mutate(initial_approval = initial_approval_date) |>
#     dplyr::mutate(indication = stringr::str_trim(stringr::str_squish(indication))) |>
#     dplyr::mutate(indication_criteria =
#                     stringr::str_trim(stringr::str_squish(indication_criteria))) |>
#     dplyr::mutate(cancer_indication = dplyr::if_else(
#       !is.na(indication) &
#         (stringr::str_detect(indication,"NSCLC|AML|CML|BRAF|KRAS|EGFR|HER2|ROS1|BRCA|ABL") |
#       (stringr::str_detect(tolower(indication), cancer_tissue_regex) &
#         stringr::str_detect(tolower(indication), "cancer|tumor|carcinoma|neoplasm|neoplastic"))),
#       TRUE,FALSE
#     )) |>
#     dplyr::mutate(cancer_indication = dplyr::if_else(
#       !is.na(indication) &
#       stringr::str_detect(tolower(indication), cancer_other_regex),
#       TRUE,as.logical(cancer_indication)
#     )) |>
#     dplyr::mutate(cancer_indication = dplyr::if_else(
#       !is.na(indication_criteria) &
#         (stringr::str_detect(indication_criteria,"NSCLC|AML|CML|BRAF|KRAS|EGFR|HER2|ROS1|BRCA|ABL") |
#       (stringr::str_detect(tolower(indication_criteria), cancer_tissue_regex) &
#         stringr::str_detect(tolower(indication_criteria), "cancer|tumor|carcinoma|neoplasm|neoplastic"))),
#       TRUE,as.logical(cancer_indication)
#     )) |>
#     dplyr::mutate(cancer_indication = dplyr::if_else(
#       !is.na(indication_criteria) &
#       stringr::str_detect(tolower(indication_criteria), cancer_other_regex),
#       TRUE,as.logical(cancer_indication)
#     )) |>
#     dplyr::filter(!stringr::str_detect(indication,"indicated( for)?:$")) |>
#     dplyr::filter(!(stringr::str_detect(indication,"indicated for") &
#                       stringr::str_detect(indication,":$"))) |>
#     dplyr::filter(!stringr::str_detect(indication,"^Limitations of (U|u)se")) |>
#     dplyr::filter(!stringr::str_detect(indication,"Limitations of (U|u)se")) |>
#     dplyr::filter(!stringr::str_detect(indication,"Select patients for therapy based on")) |>
#     dplyr::mutate(indication = Hmisc::capitalize(
#       stringr::str_replace(
#         indication,"^( )?:( )?",""))) |>
#     dplyr::mutate(indication_criteria = Hmisc::capitalize(
#       stringr::str_replace(
#         indication_criteria,"^( )?:( )?",""))) |>
#     dplyr::filter(is.na(indication_criteria) |
#                     (!is.na(indication_criteria) & nchar(indication_criteria) >= 3)) |>
#     dplyr::distinct()
#
#   return(all_indication_entries)
#
# }


merge_nci_open_targets <- function(ot_drugs = NULL,
                                 nci_antineo_all = NULL){

  ot_nci_matched <- list()

  ## X-ref Open Targets and NCI by molecule id
  ot_nci_matched[['targeted_by_id']] <- ot_drugs$targeted |>
    dplyr::left_join(nci_antineo_all[['with_chembl_id']],
                     by = c("molecule_chembl_id")) |>
    dplyr::filter(!is.na(nci_drug_name)) |>
    dplyr::select(-c(drug_name_lc, nci_drug_name))

  ## X-ref Open Targets and NCI by drug name
  ot_nci_matched[['targeted_by_name']] <- ot_drugs$targeted |>
    dplyr::left_join(
      dplyr::select(nci_antineo_all[['with_chembl_id']],
                    -molecule_chembl_id),
                     by = c("drug_name_lc" = "nci_drug_name")) |>
    dplyr::anti_join(
      ot_nci_matched[['targeted_by_id']]) |>
    dplyr::filter(!is.na(drug_name_lc)) |>
    dplyr::select(-c(drug_name_lc))


  ## X-ref Open Targeets (no target) by molecule identifier
  ot_nci_matched[['untargeted_by_id']] <- ot_drugs$untargeted |>
    dplyr::left_join(nci_antineo_all[['with_chembl_id']],
                     by = c("molecule_chembl_id")) |>
    dplyr::filter(!is.na(nci_drug_name)) |>
    dplyr::select(-c(drug_name_lc, nci_drug_name))

  ## X-ref Open Targeets (no target) by drug name
  ot_nci_matched[['untargeted_by_name']] <- ot_drugs$untargeted |>
    dplyr::left_join(
      dplyr::select(nci_antineo_all[['with_chembl_id']],
                    -molecule_chembl_id),
      by = c("drug_name_lc" = "nci_drug_name")) |>
    dplyr::anti_join(
      ot_nci_matched[['untargeted_by_id']]) |>
    dplyr::filter(!is.na(drug_name_lc)) |>
    dplyr::select(-drug_name_lc)


  ot_nci_matched_all <- do.call(rbind, ot_nci_matched) |>
    dplyr::distinct()

  ot_targeted_remain <- ot_drugs$targeted |>
    dplyr::anti_join(ot_nci_matched_all) |>
    dplyr::left_join(
      nci_antineo_all[['no_chembl_id']],
      by = c("drug_name_lc" = "nci_drug_name")) |>
    dplyr::select(-drug_name_lc)

  ot_drugs_all <- ot_nci_matched_all |>
    dplyr::bind_rows(ot_targeted_remain) |>
    dplyr::mutate(
      nci_cd_name =
        dplyr::if_else(is.na(nci_cd_name) &
                         !stringr::str_detect(drug_name,"[0-9]"),
                       Hmisc::capitalize(tolower(drug_name)),
                       nci_cd_name)) |>
    dplyr::mutate(
      nci_cd_name = dplyr::if_else(
        is.na(nci_cd_name) &
          stringr::str_detect(drug_name,"[0-9]"),
        drug_name,
        nci_cd_name)) |>
    dplyr::mutate(opentargets = TRUE)

    #dplyr::mutate(nci_version = nci_db_release) |>
    #dplyr::mutate(chembl_version = chembl_db_release) |>
    #dplyr::mutate(opentargets_version = opentargets_version)

  nci_missing_1 <- nci_antineo_all$no_chembl_id |>
    dplyr::anti_join(ot_drugs_all,
                     by = "nci_cd_name") |>
    dplyr::select(-nci_drug_name)

  nci_missing_2 <- nci_antineo_all$with_chembl_id |>
    dplyr::anti_join(ot_drugs_all,
                     by = "molecule_chembl_id") |>
    dplyr::select(-nci_drug_name)

  nci_missing <- dplyr::bind_rows(
    nci_missing_1,
    nci_missing_2) |>
    dplyr::distinct() |>
    dplyr::mutate(opentargets = FALSE)
  
  all_cancer_drugs <- ot_drugs_all |>
    dplyr::bind_rows(nci_missing) |>
    dplyr::distinct()


  ## Figure out cases where a single drug maps to multiple
  ## molecule chembl identifiers
  rownames(all_cancer_drugs) <- NULL

  custom_chembl_map <- readr::read_tsv(
    file = "data-raw/custom_chembl_map.tsv",
    col_names = F, show_col_types = F)

  colnames(custom_chembl_map) <-
    c("nci_cd_name", "molecule_chembl_id")

  custom_name_ignore <- readr::read_tsv(
      file = "data-raw/drug_names_ignore.tsv",
    col_names = F, show_col_types = F)
  colnames(custom_name_ignore) <- c("nci_cd_name")

  name2chembl_id <- all_cancer_drugs |>
    dplyr::filter(!is.na(molecule_chembl_id)) |>
    dplyr::group_by(nci_cd_name) |>
    dplyr::summarise(m = paste(
      sort(unique(molecule_chembl_id)), collapse=";"))

  ## Drugs that map to a single identifier
  name2chembl_unique <- name2chembl_id |>
    dplyr::filter(!stringr::str_detect(m,";")) |>
    dplyr::inner_join(all_cancer_drugs) |>
    dplyr::select(-m)

  ## Drugs that map to multiple identifiers
  name2chembl_ambiguous_curated <- name2chembl_id |>
    dplyr::filter(stringr::str_detect(m, ";")) |>
    dplyr::inner_join(custom_chembl_map) |>
    dplyr::inner_join(all_cancer_drugs) |>
    dplyr::select(-m)

  all_drugs <- dplyr::bind_rows(
    dplyr::filter(all_cancer_drugs, is.na(molecule_chembl_id)),
    name2chembl_unique,
    name2chembl_ambiguous_curated
  )

  salt_forms <- all_drugs |>
    dplyr::filter(stringr::str_detect(
      tolower(nci_cd_name), " (maleate|tosylate|dimaleate|succinate|meglumine|hydrobromide|sulfate|ditosylate|fumarate|dimesylate|lactate|potassium|dihydrochloride|anhydrous|hydrochloride|chloride|acetate|disodium|hydrochloride monohydrate|sodium|phosphate|dimeglumine|diphosphate|camsylate|s-malate|citrate|mesylate|acetate|calcium)$")) |>
    dplyr::filter(!is.na(opentargets_version)) |>
    dplyr::mutate(tradename = stringr::str_replace(
      nci_cd_name,
      " ((M|m)aleate|(A|a)nhydrous|(T|t)osylate|(D|d)imesylate|(D|d)imeglumine|(S|s)uccinate|(M|m)eglumine|(H|h)ydrobromide|(S|s)ulfate|(F|f)umarate|(D|d)imaleate|(D|d)itosylate|(L|l)actate|(P|p)otassium|(H|h)ydrochloride (M|m)onohydrate|(D|d)ihydrochloride|(H|h)ydrochloride|(S|s)-malate|(D|d)isodium|(C|c)amsylate|(C|c)hloride|(A|a)cetate|(S|s)odium|(P|p)hosphate|(D|d)iphosphate|(C|c)itrate|(M|m)esylate|(A|a)cetate|(C|c)alcium)$", "")) |>
    dplyr::select(tradename, nci_cd_name) |>
    dplyr::distinct() |>
    dplyr::mutate(is_salt = T) |>
    dplyr::inner_join(
      dplyr::select(all_drugs, nci_cd_name),
      by = c("tradename" = "nci_cd_name")) |>
    dplyr::distinct() |>
    dplyr::select(-tradename)


  all_drugs_final <- all_drugs |>
    dplyr::left_join(salt_forms) |>
    dplyr::mutate(is_salt = dplyr::if_else(
      is.na(is_salt),
      as.logical(FALSE),
      as.logical(is_salt)
    )) |>
    dplyr::anti_join(custom_name_ignore)

  ## antibody drug conjugates
  adc_candidates <- all_drugs_final |>
    dplyr::filter(stringr::str_detect(
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
    dplyr::select(nci_cd_name, is_adc) |>
    dplyr::distinct()

  all_drugs_final <- all_drugs_final |>
    dplyr::left_join(adc_candidates) |>
    dplyr::mutate(is_adc = dplyr::if_else(
      is.na(is_adc),
      as.logical(FALSE),
      as.logical(is_adc)
    ))


  return(all_drugs_final)

  # ambig <- all_drugs_final |>
  #   dplyr::filter(!is.na(molecule_chembl_id)) |>
  #   #dplyr::filter(!is.na(opentargets_version)) |>
  #   #dplyr::filter(is.na(is_salt)) |>
  #   dplyr::group_by(molecule_chembl_id) |>
  #   dplyr::summarise(m = paste(unique(nci_cd_name),
  #                              collapse = "@")) |>
  #   dplyr::distinct() |>
  #   dplyr::

}

map_custom_nci_targets <- function(gene_info = NULL,
                                   path_data_raw = NULL,
                                   drug_df = NULL){

  drug_target_patterns <-
    read.table(file = file.path(
      path_data_raw,
      "custom_drug_target_regex_nci.tsv"),
      sep = "\t", header = T, stringsAsFactors = F, quote = "") |>
    dplyr::inner_join(gene_info, by = "symbol") |>
    dplyr::distinct()


  all_inhibitors_no_target <- ot_nci_drugs |>
    dplyr::filter(is.na(target_symbol)) |>
    dplyr::filter(stringr::str_detect(
      tolower(nci_cd_name),
      "inhibitor|antagonist|antibody|blocker") |
        stringr::str_detect(
          tolower(nci_cd_name),
          "ib$|mab$|mab/|^anti-") |
        (stringr::str_detect(nci_concept_definition,"KRAS") &
           stringr::str_detect(nci_concept_definition,"inhibitor"))) |>
    dplyr::filter(!stringr::str_detect(
      nci_cd_name,
      " CAR T|SARS-CoV-2| Regimen$")) |>
    dplyr::filter(!stringr::str_detect(
      nci_concept_definition,
      "SARS-CoV-2"))

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
            nci_cd_name,
            "^(Inhibitor of|Anti-)|ib$|Inhibitor|targeting|ine$|ate$|ide$|mab$|antibody|ant$|mab/") &
             stringr::str_detect(nci_concept_definition, pattern))
      )

    if(nrow(hits) > 0){

      for(n in 1:nrow(hits)){
        hit <- hits[n,]

        if(stringr::str_detect(hit$nci_cd_name,
                               "mab$|monoclonal antibody")){
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
        #hit$target_uniprot_id <- target_uniprot_id
        hit$drug_clinical_source <- "nci_thesaurus_custom"
        hit$cancer_drug <- TRUE

        ## set general indications for unknown cases
        if(is.na(hit$disease_efo_id) & is.na(hit$disease_efo_label) &
           is.na(hit$cui) & is.na(hit$cui_name)){
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
      nci_concept_definition, "(A|a)ntibody(-| )drug conjugate \\(ADC\\)"
    )) |>
    dplyr::select(nci_cd_name,
                  nci_concept_definition) |>
    dplyr::distinct()


  ot_nci_drugs_curated <-
    dplyr::anti_join(ot_nci_drugs, custom_nci_targeted_drugs,
                     by = "nci_cd_name") |>
    dplyr::bind_rows(custom_nci_targeted_drugs) |>
    dplyr::arrange(target_symbol, nci_cd_name) |>
    dplyr::mutate(drug_action_type = dplyr::if_else(
      (stringr::str_detect(tolower(nci_cd_name),"inhibitor") &
         is.na(drug_action_type)) |
        (!is.na(nci_cd_name) &
           stringr::str_detect(nci_cd_name,"mab$") &
           is.na(drug_action_type)),
      "INHIBITOR",
      as.character(drug_action_type))) |>
    dplyr::mutate(cancer_drug = dplyr::if_else(
      is.na(cancer_drug) &
        (stringr::str_detect(
          tolower(nci_concept_definition),
          "anti-tumor|chemotherapy|cancer vaccine|immunothera|monoclonal antibody|antineoplastic|treatment of cancer|treatment of metastat") |
           stringr::str_detect(tolower(nci_cd_name)," regimen|recombinant|carcinoma|immune checkpoint|anti-programmed cell death ")),
      as.logical(TRUE),
      as.logical(cancer_drug)
    )) |>
    dplyr::mutate(drug_action_type = dplyr::if_else(
      is.na(drug_action_type) &
        stringr::str_detect(drug_action_type,"^(SUBSTRATE|HYDROLYTIC ENZYME|RELEASING AGENT)"),
      paste0(drug_action_type,"_OTHER"),
      as.character(drug_action_type)
    )) |>
    dplyr::mutate(comb_regimen_indication = F)


  return(ot_nci_drugs_curated)
}

assign_drug_category <- function(drug_df = NULL,
                                 path_data_raw = NULL){

  drug_df <- drug_df |>
    dplyr::distinct() |>
    dplyr::mutate(antimetabolite = dplyr::if_else(
      !is.na(nci_concept_definition) &
        stringr::str_detect(
          tolower(nci_concept_definition),
          "antimetabol|anti-metabol|nucleoside analog"),TRUE,FALSE)
    ) |>
    dplyr::mutate(topoisomerase_inhibitor = dplyr::if_else(
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           nci_concept_definition,
           "(T|t)opoisomerase II-mediated|(T|t)opoisomerase( I|II )? \\(.*\\) inhibitor|inhibit(ion|or) of (T|t)opoisomerase|(stabilizes|interrupts|binds to|interacts with|inhibits( the activity of)?)( the)?( DNA)? (t|T)opoisomerase|(T|t)opoisomerase( (I|II))? inhibitor")) |
        (!is.na(target_genename) &
           stringr::str_detect(target_genename,"topoisomerase")),TRUE,FALSE)
    ) |>
    dplyr::mutate(hedgehog_antagonist = dplyr::if_else(
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           nci_concept_definition,
           "Hedgehog") & stringr::str_detect(
             nci_cd_name,"Smoothened Antagonist|(ate|ib)$")) |
        (!is.na(nci_cd_name) &
           stringr::str_detect(
             nci_cd_name,"Hedgehog Inhibitor|SMO Protein Inhibitor")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(hdac_inhibitor = dplyr::if_else(
      (!is.na(target_symbol) &
         stringr::str_detect(
           target_symbol,
           "^HDAC")) |
        (!is.na(nci_concept_definition) &
           stringr::str_detect(nci_concept_definition,"inhibitor of histone deacetylase")) |
        (!is.na(nci_cd_name) &
           stringr::str_detect(nci_cd_name,"HDAC Inhibitor")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(alkylating_agent = dplyr::if_else(
      is.na(drug_moa) &
        !stringr::str_detect(nci_cd_name,
                             "antiangiogenic") &
        !is.na(nci_concept_definition) &
        stringr::str_detect(
          tolower(nci_concept_definition),
          "alkylating agent|alkylating activities"),TRUE,FALSE)
    ) |>
    dplyr::mutate(parp_inhibitor = dplyr::if_else(
      !is.na(target_symbol) &
        stringr::str_detect(
          target_symbol,
          "^PARP[0-9]{1}"),TRUE,FALSE)
    ) |>
    dplyr::mutate(bet_inhibitor = dplyr::if_else(
      !is.na(target_symbol) &
        stringr::str_detect(
          target_symbol,
          "^BRD(T|[1-9]{1})") |
        (!is.na(nci_cd_name) &
           stringr::str_detect(
             nci_cd_name,"BET( Bromodomain)? Inhibitor")),TRUE,FALSE)
    ) |>
    dplyr::mutate(tubulin_inhibitor = dplyr::if_else(
      (!is.na(drug_action_type) &
         drug_action_type != "STABILISER" &
         !is.na(target_genename) &
         stringr::str_detect(
           tolower(target_genename),
           "tubulin")) |
        (!is.na(nci_concept_definition) & stringr::str_detect(
          tolower(nci_concept_definition),
          "binds to tubulin|disrupts microtubule|microtubule disrupt")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(ar_antagonist = dplyr::if_else(
      (!is.na(target_genename) &
         stringr::str_detect(
           tolower(target_genename),
           "androgen receptor")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(kinase_inhibitor = dplyr::if_else(
      (!is.na(target_symbol) & stringr::str_detect(target_symbol,"EGFR|PTPN11|ABL1|FGFR|PDGFR|CSF1R")) |
        (((!is.na(drug_action_type) &
             stringr::str_detect(tolower(drug_action_type),"blocker|inhibitor|antagonist")) |
            stringr::str_detect(tolower(nci_cd_name),"ib$")) &
           (!is.na(target_genename) &
              stringr::str_detect(tolower(target_genename),"kinase|eph receptor"))) |
        (!is.na(nci_concept_definition) &
           stringr::str_detect(nci_concept_definition,"kinase inhibit(or|ion)")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(angiogenesis_inhibitor = dplyr::if_else(
      stringr::str_detect(tolower(drug_action_type),"blocker|inhibitor|antagonist") &
        (!is.na(nci_cd_name) &
           stringr::str_detect(tolower(nci_cd_name),
                               "antiangiogenic|angiogenesis inhibitor")) |
        (!is.na(nci_concept_definition) &
           stringr::str_detect(
             tolower(nci_concept_definition),
             "antiangiogenic activities|angiogenesis inhibitor|(inhibiting|blocking)( tumor)? angiogenesis|anti(-)?angiogenic|(inhibits|((inhibition|reduction) of))( .*) angiogenesis")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(monoclonal_antibody = dplyr::if_else(
      (!is.na(drug_type) & drug_type == "Antibody") |
        (stringr::str_detect(tolower(nci_cd_name),
                             "^anti-|mab |mab$|monoclonal antibody") &
           (!is.na(nci_concept_definition) &
              stringr::str_detect(nci_concept_definition,"monoclonal antibody"))),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(proteasome_inhibitor = dplyr::if_else(
      (stringr::str_detect(tolower(nci_cd_name),
                           "^proteasome") &
         !stringr::str_detect(tolower(nci_cd_name),"vaccine")) |
        (!is.na(nci_concept_definition) &
           stringr::str_detect(
             tolower(nci_concept_definition),"proteasome inhibitor|inhibits the proteasome|inhibition of proteasome")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(hormone_therapy = dplyr::if_else(
      stringr::str_detect(tolower(nci_cd_name),
                          "aromatase inhib|estrogen receptor (inhibitor|degrader|modulator)") |
        (!is.na(nci_concept_definition) &
           stringr::str_detect(
             tolower(nci_concept_definition),"inhibitor of estrogen|estrogen receptor (modulator|inhibitor|degrader)|antiestrogen|aromatase inhibit(or|ion)") &
           !stringr::str_detect(nci_concept_definition,"antiestrogen resistance")) |
        (!is.na(target_symbol) & stringr::str_detect(target_symbol,"ESR[0-9]|GNRHR")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(anthracycline = dplyr::if_else(
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           tolower(nci_concept_definition),
           "anthracycline|anthracenedione")),
      TRUE, FALSE)
    ) |>
    dplyr::mutate(immune_checkpoint_inhibitor = dplyr::if_else(
      (!is.na(nci_concept_definition) &
         !stringr::str_detect(
           tolower(nci_cd_name), "oncolytic|pentoxifylline|vaccine") &
         stringr::str_detect(
           tolower(nci_concept_definition),
           "immune checkpoint inhib")) |
        stringr::str_detect(nci_cd_name,"(Anti-(PD|CTLA)-)") |
        (stringr::str_detect(nci_cd_name,
                             "Tremelimumab|Milatuzumab")) |
        (!is.na(target_symbol) &
           !stringr::str_detect(
             tolower(nci_cd_name), "oncolytic|pentoxifylline") &
           (target_symbol == "CD274" |
              target_symbol == "CTLA4" |
              target_symbol == "PDCD1" |
              target_symbol == "TIGIT")) |
        (!is.na(nci_concept_definition) &
           !is.na(target_symbol) &
           stringr::str_detect(nci_concept_definition,
                               "immunemodulating|immune response") &
           !stringr::str_detect(
             tolower(nci_cd_name), "oncolytic|pentoxifylline|vaccine") &
           (target_symbol == "ADORA2A" |
              target_symbol == "ADORA2B")),
      TRUE,FALSE)
    ) |>
    dplyr::mutate(immune_checkpoint_inhibitor = dplyr::if_else(
      !is.na(nci_cd_name) &
        immune_checkpoint_inhibitor == T &
        stringr::str_detect(nci_cd_name,"NLM-001|CEA-MUC-1|Oncolytic|Vaccine"),
      as.logical(FALSE),
      as.logical(immune_checkpoint_inhibitor)
    )) |>
    dplyr::mutate(platinum_compound = dplyr::if_else(
      !is.na(drug_name) &
        stringr::str_detect(tolower(drug_name),"platin$"),
      as.logical(TRUE),
      as.logical(FALSE)
    ))

  ## Make sure each drug is assigned an unambiguous value for each category
  nciCDN2Category <- list()
  for(c in c('immune_checkpoint_inhibitor',
             'topoisomerase_inhibitor',
             'tubulin_inhibitor',
             'kinase_inhibitor',
             'hdac_inhibitor',
             'parp_inhibitor',
             'bet_inhibitor',
             'ar_antagonist',
             'monoclonal_antibody',
             'antimetabolite',
             'angiogenesis_inhibitor',
             'alkylating_agent',
             'anthracycline',
             'platinum_compound',
             'proteasome_inhibitor',
             'hormone_therapy',
             'hedgehog_antagonist')){

    cat <- drug_df[,c]
    name <- drug_df$nci_cd_name

    nciCDN2Category[[c]] <- as.data.frame(
      data.frame(
        'nci_cd_name' = name,
        stringsAsFactors = F
      ) |>
        dplyr::mutate(!!c := cat) |>
        dplyr::distinct() |>
        dplyr::group_by(nci_cd_name) |>
        dplyr::summarise(!!c := paste(!!dplyr::sym(c), collapse="/")) |>
        dplyr::mutate(!!c := dplyr::if_else(
          stringr::str_detect(!!dplyr::sym(c),"/"),
          TRUE,
          as.logical(!!dplyr::sym(c))))
    )

    drug_df[,c] <- NULL
    drug_df <- drug_df |>
      dplyr::left_join(
        nciCDN2Category[[c]], by = "nci_cd_name"
      )

  }

  fda_epc_codes <-
    as.data.frame(
      get_fda_ndc_mapping(path_data_raw = path_data_raw) |>
        dplyr::group_by(drug) |>
        dplyr::summarise(fda_epc_category = paste(fda_epc_category, collapse = "; "),
                         .groups = "drop")
    )

  drug_df <- drug_df |>
    dplyr::left_join(fda_epc_codes, by = c("drug_name" = "drug"))


  return(drug_df)

}

clean_final_drug_list <- function(drug_df = NULL){

  pharmaoncox <- drug_df |>
    #dplyr::filter(nci_cd_name != "Lapatinib Ditosylate") |>
    dplyr::filter(nci_cd_name != "Abivertinib Maleate") |>
    dplyr::mutate(nci_cd_name = dplyr::if_else(
      is.na(nci_cd_name),
      stringr::str_to_title(drug_name),
      as.character(nci_cd_name)
    )) |>
    dplyr::mutate(nci_cd_name = dplyr::if_else(
      nci_cd_name == "Cediranib Maleate",
      "Cediranib",
      as.character(nci_cd_name)
    )) |>
    # dplyr::mutate(nci_cd_name = dplyr::if_else(
    #   nci_cd_name == "Abivertinib Maleate",
    #   "Avitinib Maleate",
    #   as.character(nci_cd_name)
    # )) |>
    # dplyr::mutate(drug_name = dplyr::if_else(
    #   nci_cd_name == "Avitinib Maleate",
    #   "AVITINIB MALEATE",
    #   as.character(drug_name)
    # )) |>
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
      as.character(drug_name))) |>
    dplyr::distinct()


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
    dplyr::left_join(drug_action_types,
                     by = "nci_cd_name") |>
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

  #tmp2 <- tmp |>
  pharmaoncox <- pharmaoncox |>
    dplyr::left_join(drug_max_ct_phase,
                     by = "nci_cd_name") |>
    #dplyr::filter(!is.na(cancer_drug)) |>
    dplyr::select(-c(drug_moa, cancer_drug)) |>

    #dplyr::rename(nci_concept_synonym = nci_drug_name) |>
    dplyr::mutate(nci_concept_synonym2 = dplyr::if_else(
      is.na(nci_concept_synonym_all) & !is.na(drug_synonyms),
      as.character(tolower(drug_synonyms)),
      as.character(tolower(nci_concept_synonym_all))
    )) |>
    dplyr::mutate(nci_concept_synonym_all2 = nci_concept_synonym_all) |>
    #dplyr::rename(nci_concept_synonym_old = nci_concept_synonym) |>
    tidyr::separate_rows(nci_concept_synonym2,
                         sep="\\|") |>
    dplyr::rename(nci_concept_synonym = nci_concept_synonym2) |>
    dplyr::select(-c(nci_concept_synonym_all2,
                     drug_synonyms,
                     drug_tradenames,
                     drug_description)) |>
    dplyr::distinct() |>
    dplyr::select(drug_name, nci_cd_name, drug_type,
                  drug_action_type, molecule_chembl_id,
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
    #dplyr::filter(!stringr::str_detect(nci_concept_synonym,"^([a-z]{3,4})$")) |>
    dplyr::mutate(drug_action_type = stringr::str_replace_all(
      drug_action_type, "/NA|NA/",""
    ))


  ## Simplify records with only "cancer" indications, mapping them to a unique
  ## EFO/CUI cross-ref, avoiding similar records with "neoplasm", "carcinoma" etc.

  pharmaoncox_cancer_no_indication <- pharmaoncox |>
    dplyr::filter(is.na(disease_efo_id))

  pharmaoncox_cancer_NOS <- as.data.frame(
    pharmaoncox |>
      dplyr::filter(is.na(primary_site) & !is.na(disease_efo_id)) |>
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
    dplyr::bind_rows(pharmaoncox_cancer_specific) |>
    dplyr::bind_rows(pharmaoncox_cancer_NOS) |>
    dplyr::arrange(nci_cd_name) |>
    dplyr::mutate(nci_cd_name = dplyr::if_else(
      stringr::str_detect(nci_cd_name,";-"),
      drug_name,
      as.character(nci_cd_name)
    )) |>
    dplyr::filter(!is.na(nci_cd_name))


  pharmaoncox$drug_name <- NULL
  pharmaoncox$drug_name <- pharmaoncox$nci_cd_name
  pharmaoncox$nci_cd_name <- NULL
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
    dplyr::left_join(nci_t_map) |>
    dplyr::left_join(blackbox_warnings) |>
    dplyr::distinct()

  drug_maps <- list()
  drug_maps[['id2name']] <- pharmaoncox2 |>
    dplyr::select(drug_name) |>
    dplyr::distinct() |>
    dplyr::mutate(drug_id = dplyr::row_number())

  if(nrow(drug_maps[['id2name']]) != length(unique(pharmaoncox2$drug_name))){
    lgr::lgr$fatal("Ambiguous drug names")
    return(0)
  }

  pharmaoncox2 <- pharmaoncox2 |>
    dplyr::left_join(drug_maps[['id2name']], by = "drug_name")

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

  drug_maps[['id2basic']] <- pharmaoncox2 |>
    dplyr::select(drug_id,
                  molecule_chembl_id,
                  drug_type,
                  drug_action_type,
                  is_salt,
                  is_adc,
                  drug_blackbox_warning,
                  nci_concept_definition,
                  nci_t,
                  opentargets,
                  comb_regimen_indication,
                  immune_checkpoint_inhibitor,
                  topoisomerase_inhibitor,
                  tubulin_inhibitor,
                  kinase_inhibitor,
                  hdac_inhibitor,
                  parp_inhibitor,
                  bet_inhibitor,
                  monoclonal_antibody,
                  antimetabolite,
                  angiogenesis_inhibitor,
                  alkylating_agent,
                  anthracycline,
                  platinum_compound,
                  proteasome_inhibitor,
                  hormone_therapy,
                  hedgehog_antagonist,
                  fda_epc_category) |>
        dplyr::distinct()

  # tmp <- drug_maps$id2name |>
  #   dplyr::left_join(drug_maps$id2basic) |>
  #   dplyr::left_join(drug_maps$id2target) |>
  #   dplyr::left_join(drug_maps$id2indication) |>
  #   dplyr::left_join(drug_maps$id2synonym)
  #
  #
  return(drug_maps)

}


expand_drug_aliases <- function(drug_index_map = NULL,
                                path_data_raw = NULL,
                                chembl_pubchem_datestamp = chembl_pubchem_datestamp){

  chembl_pubchem_xref <-
    get_chembl_pubchem_compound_xref(
      datestamp = chembl_pubchem_datestamp,
      path_data_raw = path_data_raw)

  non_ambiguous_synonyms <- as.data.frame(
    drug_index_map[['id2synonym']] |>
      dplyr::select(drug_id, nci_concept_synonym_all) |>
      tidyr::separate_rows(nci_concept_synonym_all, sep="\\|") |>
      dplyr::distinct() |>
      dplyr::rename(nci_concept_synonym = nci_concept_synonym_all) |>
      dplyr::group_by(nci_concept_synonym) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n == 1 & nchar(nci_concept_synonym) >= 4)
  )

  ## Drug alias to NCI concept display name (primary name)
  antineopharma_synonyms <-
    drug_index_map[['id2synonym']] |>
    dplyr::select(drug_id, nci_concept_synonym_all) |>
    tidyr::separate_rows(nci_concept_synonym_all, sep="\\|") |>
    dplyr::distinct() |>
    dplyr::rename(nci_concept_synonym = nci_concept_synonym_all) |>
    dplyr::inner_join(non_ambiguous_synonyms, by = "nci_concept_synonym") |>
    dplyr::rename(alias = nci_concept_synonym) |>
    dplyr::select(-n) |>
    dplyr::left_join(
      dplyr::select(drug_index_map[['id2basic']],
                    drug_id, molecule_chembl_id),
      by = "drug_id") |>
    dplyr::distinct()

  ## include also the primary name among aliases
  tmp <- drug_index_map[['id2name']] |>
    dplyr::select(drug_id, drug_name) |>
    dplyr::mutate(alias = drug_name) |>
    dplyr::select(-drug_name) |>
    dplyr::left_join(
      dplyr::select(drug_index_map[['id2basic']], drug_id, molecule_chembl_id),
      by = "drug_id"
    ) |>
    dplyr::distinct()

  antineopharma_synonyms <- antineopharma_synonyms |>
    dplyr::bind_rows(tmp) |>
    dplyr::arrange(drug_id) |>
    dplyr::distinct()


  ## Extend aliases with those found in PubChem

  ## get drug set that contains PubChem cross-references
  unique_chembl_pubchem <- drug_index_map[['id2basic']] |>
    dplyr::select(drug_id, molecule_chembl_id) |>
    dplyr::filter(!is.na(molecule_chembl_id)) |>
    dplyr::distinct() |>
    dplyr::left_join(chembl_pubchem_xref, by="molecule_chembl_id") |>
    dplyr::filter(!is.na(pubchem_cid)) |>
    dplyr::select(-c(chembl_db_version))


  ## Retrieve aliases for drugs with PubChem x-refs
  lgr::lgr$info("Retrieving additional drug aliases/synonyms for PubChem-mapped compounds")
  pubchem_synonym_files <-
    sort(list.files(path = file.path(path_data_raw, "pubchem"),
                    pattern = "CID-Synonym-filtered_",
                    full.names = T))

  antineopharma_synonyms_pubchem <- data.frame()
  for(f in pubchem_synonym_files){
    synonym_data <- as.data.frame(readr::read_tsv(
      f, col_names = c('pubchem_cid','alias'),
      col_types = "dc",
      progress = F
    ))

    pubchem_alias_df <- synonym_data |>
      dplyr::inner_join(unique_chembl_pubchem,
                        by = "pubchem_cid")

    if(nrow(pubchem_alias_df) > 0){
      pubchem_alias_df <- pubchem_alias_df |>
        dplyr::select(-pubchem_cid)
      antineopharma_synonyms_pubchem <-
        antineopharma_synonyms_pubchem |>
        dplyr::bind_rows(pubchem_alias_df)
    }
    rm(synonym_data)
  }

  # antineopharma_synonyms_pubchem <- antineopharma_synonyms_pubchem |>
  #   dplyr::left_join(drug_index_map[['id2name']])


  ## Only include drug aliases that are unambiguous
  unambiguous_drug_aliases <- antineopharma_synonyms |>
    dplyr::bind_rows(antineopharma_synonyms_pubchem) |>
    dplyr::left_join(drug_index_map[['id2name']], by = "drug_id") |>
    dplyr::filter(!(alias == "nab-paclitaxel" &
                      drug_name == "Paclitaxel")) |>
    dplyr::select(alias, drug_name) |>
    dplyr::distinct() |>
    dplyr::group_by(alias) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n == 1) |>
    dplyr::select(alias)

  compound_synonyms <-
    dplyr::bind_rows(antineopharma_synonyms, antineopharma_synonyms_pubchem) |>
    dplyr::left_join(drug_index_map[['id2name']], by = "drug_id") |>
    dplyr::distinct() |>
    dplyr::inner_join(unambiguous_drug_aliases, by = "alias") |>
    dplyr::distinct()  |>
    dplyr::mutate(alias = dplyr::if_else(
      alias == "nab-paclitaxel" & drug_name == "Paclitaxel",
      "paclitaxel",
      as.character(alias)
    )) |>
    dplyr::mutate(
      alias = stringi::stri_enc_toascii(alias)
    ) |>
    dplyr::mutate(
      drug_name =
        stringi::stri_enc_toascii(drug_name)
    ) |>
    dplyr::arrange(drug_id) |>
    dplyr::select(-c(molecule_chembl_id, drug_name)) |>
    dplyr::distinct() |>
    dplyr::mutate(
      alias_lc = tolower(alias)) |> 
    
    ## remove aliases that co-incide with ordinary
    ## english words
    dplyr::left_join(
      words::words, by = c("alias_lc" = "word")) |> 
    dplyr::filter(
      is.na(word_length) |
        (!is.na(word_length) &
           word_length > 6 & 
           !stringr::str_detect(
             alias_lc,"(one|i(n|l)(e|s)?|ol|id(e)?|oid|trel|opa|phen|xel|tal|rapy|lite(s)?|tomy|ase|ole|pam|fan|fen|yl|ane|ose|ium|(ph|z)ene|yde|lan|tam|fam|xal|strogen|gen|nal|xan|ene|gon|ram|glycan|prim|vir|yte(s)?|ate(s)?)$"
           )
        )
  ) |>
    dplyr::filter(alias_lc != "medium" & 
                    alias_lc != "serum" &
                    alias_lc != "bypass" &
                    alias_lc != "molecule" &
                    alias_lc != "feces" &
                    alias_lc != "fiber" &
                    alias_lc != "start" &
                    alias_lc != "other" &
                    alias_lc != "freeze" &
                    alias_lc != "transplantation" &
                    alias_lc != "hydrochloride") |>
    dplyr::select(-c(word_length,alias_lc)) |>
    dplyr::distinct()
  
  return(compound_synonyms)


}
