

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

## Retrieval of mechanism-of-action for ChEMBL molecular compounds
#'
#' Function that accepts a ChEMBL molecule identifer and uses the ChEMBL web service API to return
#' its mechanism of action (including disease efficacy, trial phase etc)
#'
#' @param molecule_chembl_id ChEMBL molecule identifier
#' @return a data frame with the following columns of chemical compound properties:
#' \itemize{
#'   \item molecule_chembl_id_salt
#'   \item molecule_chembl_id = parent_chembl_id,
#'   \item target_chembl_id
#'   \item chembl_disease_efficacy
#'   \item chembl_max_phase
#'   \item chembl_mechanism_comment
#'   \item chembl_moa
#' }
#' @examples
#' ## Retrieve mechanism of action for IRINOTECAN (ChEMBL molecular compound ID = 'CHEMBL481')
#' \dontrun{
#' chembl_compound_moa(molecule_chembl_id = 'CHEMBL481')
#' }
#'

chembl_compound_moa <- function(chembl_ws_base_url = 'https://www.ebi.ac.uk/chembl/api/data', molecule_chembl_id = NA){

  moa_all <- data.frame()

  ## get salts
  if(!is.na(molecule_chembl_id)){
    molecule_chembl_id_salt <- NA
    parent_chembl_id <- NA
    salt_df <- chembl_compound_salts(molecule_chembl_id = molecule_chembl_id)
    j <- 1
    ## loop over salts to retrieve mechanism of action (phase, targets etc.)
    while(j <= nrow(salt_df)){
      molecule_chembl_id_salt <- salt_df[j,]$molecule_chembl_id
      parent_chembl_id <- salt_df[j,]$parent_chembl_id
      chembl_moa <- NULL
      moa_url <- paste0(chembl_ws_base_url,'/mechanism?molecule_chembl_id=',molecule_chembl_id_salt)
      if(!is.null(readUrl(moa_url))){
        raw_xml_string <- rawToChar(httr::GET(moa_url)$content)
        if(stringr::str_detect(raw_xml_string,"<mechanism>")){
          doc <- XML::xmlParse(raw_xml_string)
          mechanism_df <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/mechanisms/mechanism"), collectNames = T, stringsAsFactors = F)
          if(nrow(mechanism_df) > 0){
            mechanism_df <- dplyr::select(mechanism_df, disease_efficacy, max_phase, mechanism_comment, mechanism_of_action, molecule_chembl_id, target_chembl_id)
            chembl_moa <- mechanism_df |>
              dplyr::rename(chembl_disease_efficacy = disease_efficacy, chembl_max_phase = max_phase,
                            chembl_mechanism_comment = mechanism_comment, chembl_moa = mechanism_of_action) |>
              dplyr::select(molecule_chembl_id,target_chembl_id,chembl_disease_efficacy,chembl_max_phase,chembl_mechanism_comment,chembl_moa)
            chembl_moa$molecule_chembl_id_salt <- molecule_chembl_id_salt
            chembl_moa$molecule_chembl_id <- as.character(parent_chembl_id)
            chembl_moa$chembl_disease_efficacy <- as.logical(dplyr::recode(chembl_moa$chembl_disease_efficacy, True = TRUE, False = FALSE))
          }
        }else{
          chembl_moa <- data.frame(molecule_chembl_id_salt = molecule_chembl_id_salt,
                                   molecule_chembl_id = as.character(parent_chembl_id),
                                   target_chembl_id = as.character(NA),
                                   chembl_disease_efficacy = as.logical(NA),
                                   chembl_max_phase = as.character(NA),
                                   chembl_mechanism_comment = as.character(NA),
                                   chembl_moa = as.character(NA), stringsAsFactors = F)
        }
      }else{
        chembl_moa <- data.frame(molecule_chembl_id_salt = molecule_chembl_id_salt,
                                 molecule_chembl_id = as.character(parent_chembl_id),
                                 target_chembl_id = as.character(NA),
                                 chembl_disease_efficacy = as.logical(NA),
                                 chembl_max_phase = as.character(NA),
                                 chembl_mechanism_comment = as.character(NA),
                                 chembl_moa = as.character(NA), stringsAsFactors = F)
      }
      if(!is.null(chembl_moa)){
        moa_all <- rbind(moa_all, chembl_moa)
      }
      j <- j + 1
    }
    if(nrow(salt_df) == 0){
      moa_all <- data.frame(molecule_chembl_id_salt = NA,
                            molecule_chembl_id = as.character(molecule_chembl_id),
                            target_chembl_id = as.character(NA),
                            chembl_disease_efficacy = as.logical(NA),
                            chembl_max_phase = as.character(NA),
                            chembl_mechanism_comment = as.character(NA),
                            chembl_moa = as.character(NA), stringsAsFactors = F)
    }
    else{
      moa_all <- dplyr::filter(moa_all, !is.na(target_chembl_id))
    }
  }
  return(moa_all)
}

## Retrieval of compound salts for ChEMBL molecular compounds
#'
#' Function that accepts a ChEMBL molecule identifer and uses the ChEMBL web service API to return
#' relationships between molecule parents and salts (https://www.ebi.ac.uk/chembl/ws)
#'
#' @param molecule_chembl_id ChEMBL molecule identifier
#' @return a data frame with the following columns of chemical compound properties:
#' \itemize{
#'   \item is_parent - logical indicator if a salt is considered a parent compound
#'   \item molecule_chembl_id - query compound identifier
#'   \item parent_chembl_id - parent compound identifier (salt)
#' }
#'
chembl_compound_salts <- function(chembl_ws_base_url = 'https://www.ebi.ac.uk/chembl/api/data', molecule_chembl_id = NA){
  salt_df <- data.frame()
  if(!is.na(molecule_chembl_id)){
    salt_url <- paste0(chembl_ws_base_url, '/molecule_form/', molecule_chembl_id)

    if(!is.null(readUrl(salt_url))){
      doc <- XML::xmlParse(rawToChar(httr::GET(salt_url)$content))
      salt_df <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/molecule_forms/molecule_form"), collectNames = T, stringsAsFactors = F)
    }else{
      salt_df <- data.frame('is_parent' = NA, molecule_chembl_id = molecule_chembl_id, parent_chembl_id = NA, stringsAsFactors = F)
    }
  }
  return(salt_df)
}


## Retrieval of drug target predictions for ChEMBL molecular compounds
#'
#' Function that accepts a ChEMBL molecule identifer and uses the ChEMBL web service API to return
#' targets for which there is a predictied binding of a the molecular target (https://www.ebi.ac.uk/chembl/ws)
#'
#' @param molecule_chembl_id ChEMBL molecule identifier
#' @param min_probability minimum probability of compound-target binding
#' @return a data frame with the following columns of chemical compound properties:
#' \itemize{
#'   \item target_chembl_id
#'   \item target_uniprot_acc
#'   \item molecule_chembl_id_salt
#'   \item molecule_chembl_id
#'   \item target_uniprot_id
#'   \item symbol
#'   \item entrezgene
#'   \item gene_name
#'   \item target_compound_binding_probability
#' }
#' @examples
#' ## Retrieve predicted targets for IRINOTECAN (ChEMBL molecular compound ID = 'CHEMBL481')
#' \dontrun{
#' chembl_compound_target_predictions(molecule_chembl_id = 'CHEMBL481', probability = 0.9)
#' }
#'

chembl_compound_target_predictions <- function(chembl_ws_base_url = 'https://www.ebi.ac.uk/chembl/api/data',
                                               min_probability = 0.9, molecule_chembl_id = NA){

  target_predictions_all <- data.frame()
  ## get salts
  if(!is.na(molecule_chembl_id)){
    salt_df <- chembl_compound_salts(molecule_chembl_id = molecule_chembl_id)
    j <- 1
    while(j <= nrow(salt_df)){
      molecule_chembl_id_salt <- salt_df[j,]$molecule_chembl_id
      parent_chembl_id <- salt_df[j,]$parent_chembl_id
      chembl_target_predictions <- NULL
      target_pred_url <- paste0(chembl_ws_base_url, '/target_prediction?limit=500&probability__gt=',min_probability,'&target_tax_id__exact=9606&molecule_chembl_id=', molecule_chembl_id)
      if(!is.null(readUrl(target_pred_url))){
        doc <- XML::xmlParse(rawToChar(httr::GET(target_pred_url)$content))
        target_df <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//response/target_predictions/target_prediction"), collectNames = T, stringsAsFactors = F)
        if(nrow(target_df) > 0){
          chembl_target_predictions <- dplyr::select(target_df, target_chembl_id, target_accession, probability) |>
            dplyr::rename(target_uniprot_acc = target_accession) |>
            dplyr::distinct()
          chembl_target_predictions$molecule_chembl_id_salt <- molecule_chembl_id_salt
          chembl_target_predictions$molecule_chembl_id <- parent_chembl_id
          chembl_target_predictions <- chembl_target_predictions |>
            #dplyr::left_join(dplyr::select(pharmaMine::all_protein_targets, target_uniprot_acc, target_uniprot_id, symbol, entrezgene, gene_name)) |>
            dplyr::mutate(target_compound_binding_probability = round(as.numeric(probability), digits = 6)) |>
            dplyr::arrange(desc(target_compound_binding_probability)) |>
            dplyr::select(-probability)


        }else{
          chembl_target_predictions <- data.frame(target_chembl_id = NA,
                                                  target_uniprot_acc = NA,
                                                  target_compound_binding_probability = NA,
                                                  molecule_chembl_id_salt = molecule_chembl_id_salt,
                                                  molecule_chembl_id = parent_chembl_id,
                                                  target_uniprot_id = NA,
                                                  symbol = NA,
                                                  entrezgene = NA,
                                                  gene_name = NA,
                                                  stringsAsFactors = F)
        }
      }
      if(!is.null(chembl_target_predictions)){
        target_predictions_all <- rbind(target_predictions_all, chembl_target_predictions)
      }
      j <- j + 1
    }
  }
  if(nrow(salt_df) == 0){
    target_predictions_all <- data.frame(target_chembl_id = NA,
                                         target_uniprot_acc = NA,
                                         target_compound_binding_probability = NA,
                                         molecule_chembl_id_salt = molecule_chembl_id_salt,
                                         molecule_chembl_id = parent_chembl_id,
                                         target_uniprot_id = NA,
                                         symbol = NA,
                                         entrezgene = NA,
                                         gene_name = NA,
                                         stringsAsFactors = F)
  }
  target_predictions_all <- target_predictions_all |> dplyr::distinct()

  return(target_predictions_all)
}


chembl_compound_indications <- function(chembl_ws_base_url = "https://www.ebi.ac.uk/chembl/api/data/drug_indication?molecule_chembl_id=", molecule_chembl_id = NULL, min_maxphase = 3){

  ## get drug indications (i.e. diseases, cancer types)
  indications <- NULL
  indication_url <- paste0(chembl_ws_base_url,molecule_chembl_id,"&limit=100")
  doc <- XML::xmlParse(rawToChar(httr::GET(indication_url)$content))
  efo_nodes <- XML::getNodeSet(doc, "//response/drug_indications/drug_indication/efo_id")
  efo_term_nodes <- XML::getNodeSet(doc, "//response/drug_indications/drug_indication/efo_term")
  phase_nodes <- XML::getNodeSet(doc, "//response/drug_indications/drug_indication/max_phase_for_ind")
  efo_ids <- sapply(efo_nodes, function(el) XML::xmlValue(el, "efo_id"))
  efo_terms <- sapply(efo_term_nodes, function(el) XML::xmlValue(el, "efo_term"))
  max_phase_indications <- sapply(phase_nodes, function(el) XML::xmlValue(el, "max_phase_for_id"))
  j <- 1
  while(j <= length(max_phase_indications)){
    if(max_phase_indications[j] == ""){
      max_phase_indications[j] <- 0
    }
    j <- j + 1
  }
  if(length(efo_ids) == length(efo_terms) & length(efo_terms) == length(max_phase_indications)){
    indications <- data.frame('efo_id' = efo_ids, 'efo_term' = efo_terms, 'max_phase' = max_phase_indications, stringsAsFactors = F)
    indications$max_phase <- as.integer(indications$max_phase)
    if(nrow(indications) == 0){
      indications$efo_id <- as.character(indications$efo_id)
      indications$efo_term <- as.character(indications$efo_term)
    }

    if(nrow(indications) > 0){
      if(is.null(indications[nchar(indications$efo_id) == 0,])){
        indications[nchar(indications$efo_id) == 0,] <- NA
      }
      indications <- dplyr::filter(indications, max_phase >= min_maxphase)
      indications <- indications |> dplyr::distinct()
      if(nrow(indications) > 0){
        indications$molecule_chembl_id <- molecule_chembl_id
      }
    }
  }
  if(length(indications) ==  0){
    indications <- data.frame('efo_id' = character(), 'efo_term' = character(), 'max_phase' = integer(), 'molecule_chembl_id' = character(), stringsAsFactors = F)
  }

  return(indications)
}

## Retrieval of structurally similar compounds based on SMILES
#'
#' Function that accepts a SMILEs and queries PubChem for structurally similar compounds
#'
#' @param approach 'pubchem_pug' or 'chemminer'
#' @param smiles query SMILES
#' @param manimoto_threshold minimum similarity
#' @param MaxRec maximum number of records
#' @param MaxSec maximum number of
#' @return a character vector of PubChem compound identifiers (CIDs)
#'

get_similar_compounds <- function(approach = 'pubchem_pug', smiles = NULL, manimoto_threshold = 100, MaxRec = 100, MaxSec = 10, mw = NULL){
  cids <- integer()
  smiles_set <- unlist(stringr::str_split(smiles,";"))
  filtered_cids <- integer()

  if(approach == 'pubchem_pug'){
    for(i in 1:length(smiles_set)){
      cids <- integer()
      query_url <- paste0(PUBCHEM_PUG_URL,'compound/similarity/smiles/',smiles_set[i],'/JSON?Threshold=',manimoto_threshold,'&MaxRecords=',MaxRec,'&MaxSec=',MaxSec)
      tcOut <- tryCatch(
        {
          #con1 <- query_url
          async_query <- jsonlite::fromJSON(query_url)
          Sys.sleep(10)
          async_query_results <- jsonlite::fromJSON(paste0(PUBCHEM_PUG_URL,'compound/listkey/',async_query$Waiting$ListKey,'/cids/JSON'))
          cids <- async_query_results$IdentifierList$CID
          length(cids)
        },
        error = function(cond){
          message(paste("URL does not seem to exist:", query_url))
          return(NA)
        },
        warning = function(cond){
          message(paste("URL caused a warning:", query_url))
          return(NA)
        }
      )
      if(!is.na(tcOut) && length(cids) > 0){
        cids_subset <- match_cids_by_molecularweight(cids, mw)
        if(length(cids_subset) > 0){
          filtered_cids <- append(filtered_cids, cids_subset)
        }
      }
    }
  }
  # if(approach == 'chemminer'){
  #   ## ChemmineR as an alternative option
  #   for(i in 1:length(smiles_set)){
  #     try({
  #       results <- ChemmineR::searchString(smiles_set[i])
  #       num_results <- length(results)
  #       if(num_results > MaxRec){
  #         num_results <- MaxRec
  #       }
  #       for(m in 1:num_results){
  #         cids <- append(cids, as.integer(unlist(datablock(results[m]))[paste0("CMP",m,".PUBCHEM_COMPOUND_CID")]))
  #       }
  #     })
  #   }
  #   if(length(cids) > 0){
  #     filtered_cids <- match_cids_by_molecularweight(cids, mw)
  #   }
  # }
  return(filtered_cids)
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
    magrittr::set_colnames(c('nci_t','nci_concept_display_name','cui')) |>
    dplyr::filter(!is.na(nci_t))

  return(nci_labels)
}

### CHEMBL-PUBCHEM COMPOUND CROSS-REFERENCE

get_chembl_pubchem_compound_xref <- function(datestamp = '20220429',
                                             chembl_release = NULL,
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



get_chembl_compound_by_name <- function(name = NULL){
  df <- NULL
  if(!is.null(name)){
    query <-
      paste0("https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup/search?format=json&q=",
             tolower(stringr::str_replace_all(name," ","%20")))
    result <- jsonlite::fromJSON(query)
    if(length(result$chembl_id_lookups) > 0){
      if(nrow(result$chembl_id_lookups) >= 1){
        recs <- dplyr::filter(result$chembl_id_lookups,
                              entity_type == "COMPOUND")
        if(nrow(recs) == 1){
          df <- data.frame('name' = name,
                           'molecule_chembl_id' = recs$chembl_id,
                           stringsAsFactors = F)
          rlogging::message(paste('Found',nrow(recs), 'records for',
                                  name,'-',sep=" "))
          return(df)
        }
      }
    }
  }
  return(df)
}

### CHEMBL DRUG TARGETS

get_chembl_drug_targets <- function(datestamp = pharmamine_datestamp,
                                    chembl_release = chembl_db_release,
                                    path_data_raw = NULL){

  uniprot_idmapping_url <-
    'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'
  uniprot_chembl_mapping_url <-
    "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_uniprot_mapping.txt"

  chembl_uniprot_mapping_fname <- file.path(path_data_raw, "chembl", paste0("chembl_uniprot_mapping",datestamp,".txt"))
  gene_info_fname <- file.path(path_data_raw, "Homo_sapiens.gene_info.gz")
  uniprot_id_mapping <- file.path(path_data_raw, "uniprot", "HUMAN_9606_idmapping.dat.gz")

  if(!file.exists(chembl_uniprot_mapping_fname)){
    download.file(uniprot_chembl_mapping_url,
                  destfile = chembl_uniprot_mapping_fname)
  }
  chembl_drug_targets <- read.table(file = chembl_uniprot_mapping_fname, sep="\t",
                                    skip = 1, stringsAsFactors = F,quote = "")
  colnames(chembl_drug_targets) <- c('uniprot_acc','target_chembl_id','target_chembl_name','target_chembl_type')
  chembl_drug_targets <- dplyr::select(chembl_drug_targets, -target_chembl_name)

  if(!file.exists(gene_info_fname)){
    download.file("ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                  destfile = gene_info_fname, quiet=T)
  }
  gene_info_ncbi <- read.table(gzfile(gene_info_fname), sep="\t", stringsAsFactors = F,na.strings="-", skip=1, comment.char = "#", quote = "", fill=T)
  gene_info_ncbi <- gene_info_ncbi |> dplyr::filter(V1 == 9606)
  gene_info_ncbi <- dplyr::select(gene_info_ncbi, c(V2,V3,V9,V10))
  colnames(gene_info_ncbi) <- c('entrezgene','symbol','gene_name','gene_biotype')
  gene_info_ncbi$entrezgene <- as.character(gene_info_ncbi$entrezgene)
  if(nrow(gene_info_ncbi[!is.na(gene_info_ncbi$gene_biotype) & gene_info_ncbi$gene_biotype == 'protein-coding',])>0){
    gene_info_ncbi[!is.na(gene_info_ncbi$gene_biotype) & gene_info_ncbi$gene_biotype == 'protein-coding',]$gene_biotype <- 'protein_coding'
  }

  if(!file.exists(uniprot_id_mapping)){
    download.file(uniprot_idmapping_url, destfile = uniprot_id_mapping)
  }

  idmapping <- read.table(gzfile(uniprot_id_mapping),sep="\t",
                          header = F,quote = "", stringsAsFactors = F)
  idmapping_up_kb <- dplyr::filter(idmapping, V2 == 'UniProtKB-ID' | V2 == 'GeneID')
  colnames(idmapping_up_kb) <- c('uniprot_acc','type','name')

  idmapping_up_id <- dplyr::filter(idmapping_up_kb, type == 'UniProtKB-ID')
  idmapping_up_id <- dplyr::rename(idmapping_up_id, uniprot_id = name) |> dplyr::select(uniprot_acc, uniprot_id) |> dplyr::distinct()
  idmapping_geneid <- dplyr::filter(idmapping_up_kb, type == 'GeneID')
  idmapping_geneid <- dplyr::rename(idmapping_geneid, entrezgene = name) |> dplyr::select(uniprot_acc, entrezgene) |> dplyr::distinct()

  uniprot_acc_mapping <- dplyr::inner_join(idmapping_up_id, idmapping_geneid)
  uniprot_acc_mapping <- dplyr::left_join(uniprot_acc_mapping, gene_info_ncbi)
  uniprot_acc_mapping <- uniprot_acc_mapping |> dplyr::filter(!is.na(symbol) | !is.na(uniprot_id))

  rm(idmapping_up_id)
  rm(idmapping_geneid)

  chembl_drug_targets <- dplyr::left_join(chembl_drug_targets, uniprot_acc_mapping)
  chembl_drug_targets <- dplyr::filter(chembl_drug_targets, !is.na(uniprot_id)) |>
    dplyr::distinct() |>
    dplyr::rename(target_uniprot_acc = uniprot_acc, target_uniprot_id = uniprot_id)

  rm(uniprot_acc_mapping)
  chembl_drug_targets$chembl_db_version <- chembl_release
  return(chembl_drug_targets)
}


### TARGETED ANTICANCER COMPOUNDS FROM OPEN TARGETS
get_opentargets_cancer_drugs <-
  function(path_data_raw = NULL,
           ot_version = "2022.02",
           uniprot_release = "2021_04"){

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
                       by = c("drug_name","molecule_chembl_id"))
  ) |>
    dplyr::mutate(target_uniprot_id = NA)


  up_ensembl_map <- as.data.frame(get_uniprot_ensembl_map(
    path_data_raw = path_data_raw,
    uniprot_release = uniprot_release) |>
    dplyr::rename(target_symbol = symbol,
                  target_uniprot_id = uniprot_id) |>
    dplyr::select(-c(entrezgene, target_chembl_id)) |>
    dplyr::group_by(target_symbol, target_ensembl_gene_id) |>
    dplyr::summarise(target_uniprot_id = paste(
      target_uniprot_id, collapse = "|"),
      .groups = "drop"
    ))

  targeted_cancer_compounds <- as.data.frame(
    targeted_compounds |>
      dplyr::distinct() |>
      dplyr::left_join(up_ensembl_map,
                       by = c("target_ensembl_gene_id","target_symbol")) |>
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
      dplyr::filter(cancer_drug == T &
                      !is.na(target_symbol) &
                      !is.na(target_ensembl_gene_id))
  )

  # n_entries <- as.data.frame(
  #   targeted_cancer_compounds |>
  #     dplyr::group_by(target_ensembl_gene_id,
  #                     disease_efo_id,
  #                     molecule_chembl_id) |>
  #     dplyr::summarise(n = dplyr::n(), .groups = "drop")
  # )
  #
  # targeted_cancer_compounds <- targeted_cancer_compounds |>
  #   dplyr::left_join(n_entries,
  #                    by = c("target_ensembl_gene_id",
  #                           "disease_efo_id",
  #                           "molecule_chembl_id")) |>
  #   #dplyr::filter(n == 1) |>
  #   dplyr::select(-n) |>
  #   dplyr::distinct()

  ## adjust max ct phase
  drugs_with_max_phase_adj <- as.data.frame(
    targeted_cancer_compounds |>
    dplyr::filter(!is.na(drug_max_ct_phase) & !is.na(drug_max_phase_indication)) |>
    dplyr::group_by(drug_name, molecule_chembl_id) |>
    dplyr::summarise(drug_max_ct_phase = max(drug_max_phase_indication),
                     .groups = "drop")
  )

  targeted_cancer_compounds <- as.data.frame(targeted_cancer_compounds |>
    dplyr::select(-drug_max_ct_phase) |>
    dplyr::left_join(drugs_with_max_phase_adj, by = c("drug_name","molecule_chembl_id"))
  )


  return(list('targeted' = targeted_cancer_compounds,
              'untargeted' = untargeted_cancer_compounds))

}

get_fda_ndc_mapping <- function(
  path_data_raw = NULL
){

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

  if(!file.exists(file.path(path_data_processed,"nci_thesaurus","nci_treatment_thesaurus_antineo.rds")) | overwrite == T){
    sorafenib_definition <- 'A synthetic compound targeting growth signaling and angiogenesis. Sorafenib blocks the enzyme RAF kinase, a critical component of the RAF/MEK/ERK signaling pathway that controls cell division and proliferation; in addition, sorafenib inhibits the VEGFR-2/PDGFR-beta signaling cascade, thereby blocking tumor angiogenesis.'

    nci_display_labels <- process_nci_labels(path_data_raw = path_data_raw,
                                             overwrite = overwrite)

    drug2chembl <-
      read.table(file = file.path(path_data_raw,"dgidb","dgidb.202101.tsv"),
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
                    nci_concept_display_name = V6, nci_concept_status = V7,
                    nci_concept_semantic_type = V8) |>

      dplyr::filter(
        stringr::str_detect(
          nci_concept_semantic_type,"Chemical|Substance|Therapeutic|Drug|Immunologic")) |>
      dplyr::left_join(nci_antineo_agents, by = c("nci_t")) |>
      dplyr::select(-c(nci_concept_name, nci_concept_display_name)) |>
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
      dplyr::filter(!is.na(nci_concept_display_name)) |>
      dplyr::filter(antineoplastic_agent == T |
                      stringr::str_detect(nci_concept_synonym,"chemother|immunother|radiation|regimen") |
                      stringr::str_detect(tolower(nci_concept_display_name),"chemother|immunother|radiation|regimen") |
                      stringr::str_detect(tolower(nci_concept_display_name),"(ib|mab)$") |
                      stringr::str_detect(tolower(nci_concept_definition),"immuno|b-cell|t-cell|melanoma|myelomoa|neuroblastoma|leukemia|lymphoma|thymoma|sarcoma|paraganglioma|regimen|proliferation|cancer|tumor|carcino|radiation|block|radio|chemotherapy|immune|antigen|antineo")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name," (Gel|Oil|Cream|Seed|Block|Field|Supplement|Factor)$")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"(Vaccination|Lotion|Therapeutic Heat|Procedure|Rehabilitation|Prevention|Rinse)$")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"(Epitope|Exract|Influenza|Ginseng|Ointment|Management|Injection|Tool)$")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"^(Vitamin A Compound|Inactivated Poliovirus|Antineoplastic Immune Cell|Topical)")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"^(Sheng-Yu|Ginseng|Dry Cleaning|Boost|Tobacco|Microwave)")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"^(Blood|Interruption of|Gum Arabic|Vaginal Cylinder|Laser Ablation|Wheatgrass)")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"^(Chemical Challenge|Prevention of|Magic Mouthwash|Wood Dust|Soot|Cocaine)")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name,"^(Antibody|Antigen|Antioxidant|Vaccination|Acetate|Antiserum|Asbestos|Aspirate|Autoantigen|Cytokine)$")) |>
      dplyr::filter(!stringr::str_detect(nci_concept_display_name," Spray| Extract| Antidiabetic| Implant|(Green Tea|Living Healthy|Pollutant|Probe|Protective Agent|Supportive Care|Caffe)")) |>
      dplyr::filter(
        !stringr::str_detect(
          tolower(nci_concept_definition),
          "chinese |antidiabet|diabetes|antidepress|analgesic|pulmonary edema|nutritional|human carcinogen|anesthetic|nonsedating|sedative|antihyper|antiinflamma|antiarrythm|antiangin|antihist|muscle|neurotransmitter"))
    )

    #nci_antineo_thesaurus
    nci2chembl <- as.data.frame(nci_antineo_thesaurus_raw |>
      dplyr::select(molecule_chembl_id, nci_concept_display_name) |>
      dplyr::filter(!is.na(molecule_chembl_id)) |>
      dplyr::distinct()
    )

    j <- 1
    nci2chembl_dict <- list()
    while(j <= nrow(nci2chembl)){
      nci2chembl_dict[[nci2chembl[j,"nci_concept_display_name"]]] <-
        nci2chembl[j,"molecule_chembl_id"]
      j <- j + 1
    }

    i <- 1
    while(i <= nrow(nci_antineo_thesaurus_raw)){
      nci_concept_display_name <-
        nci_antineo_thesaurus_raw[i,"nci_concept_display_name"]
      if(nci_concept_display_name %in% names(nci2chembl_dict)){
        nci_antineo_thesaurus_raw[i,"molecule_chembl_id"] <-
          nci2chembl_dict[[nci_concept_display_name]]
      }
      i <- i + 1
    }

    nci_antineo_thesaurus <- nci_antineo_thesaurus_raw |>
      dplyr::filter(!(molecule_chembl_id == "CHEMBL1569487" & nci_t == "C405")) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Sorafenib Tosylate", "CHEMBL1200485", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Leucovorin", "CHEMBL1679", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Masoprocol", "CHEMBL313972", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "7-Hydroxystaurosporine", "CHEMBL1236539", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Dromostanolone", "CHEMBL1201048", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Patritumab", "CHEMBL2109406", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Rivoceranib", "CHEMBL3186534", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Rivoceranib Mesylate", "CHEMBL3545414", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Tricribine Phosphate", "CHEMBL462018", as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Umbrasilib" |
          nci_concept_display_name == "Umbrasilib Tosylate", "CHEMBL3948730",
        as.character(molecule_chembl_id))) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Brigatinib") & molecule_chembl_id != "CHEMBL3545311")) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Brentuximab") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Afatinib") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Dacomitinib") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Calcipotriene") & molecule_chembl_id != "CHEMBL2109665")) |>
      # dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Balixafortide") & molecule_chembl_id != "CHEMBL2109665")) |>
      dplyr::filter(!(stringr::str_detect(nci_concept_display_name,"^Antineoplast(on|ic) Agent|^Support| Gel | Ointment|Caffeine|^Acetate$")))
      #CHEMBL3545055

    ## add apatinib as an alias (is missing in NCI thesaurus)
    apatinib_alias_entry <- nci_antineo_thesaurus |>
      dplyr::filter(nci_concept_synonym == "rivoceranib") |>
      dplyr::mutate(nci_concept_synonym = "apatinib") |>
      dplyr::mutate(nci_concept_synonym_all =
                      paste(nci_concept_synonym_all,"Apatinib",sep="|"))

    nci_antineo_thesaurus <- nci_antineo_thesaurus |>
      dplyr::bind_rows(apatinib_alias_entry) |>
      dplyr::arrange(nci_concept_display_name)

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
      dplyr::rename(drug_name_nci = nci_concept_synonym) |>
      dplyr::filter(!(nci_concept_display_name == "Sorafenib Tosylate" &
                        drug_name_nci == "sorafenib")) |>
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
      dplyr::select(nci_concept_display_name, molecule_chembl_id) |>
      dplyr::distinct() |>
      dplyr::filter(!stringr::str_detect(tolower(nci_concept_display_name),
                                         "( vaccine)|^[0-9]")) |>
      dplyr::filter(is.na(molecule_chembl_id)) |>
      dplyr::select(nci_concept_display_name) |>
      dplyr::mutate(num_spaces = stringr::str_count(nci_concept_display_name," ")) |>
      dplyr::filter(num_spaces <= 1 & !stringr::str_detect(tolower(nci_concept_display_name),"regimen|&|/|;|,")) |>
      dplyr::filter(stringr::str_detect(nci_concept_display_name,"(mab|cin|ide|ib|im|bine|tin|om|lin|stat|one|ate|ole|ane|ine|xel|rol)$")) |>
      dplyr::distinct()
    )



    nci_compounds_chembl_match <- data.frame()

    ## Retrieve aliases for drugs with PubChem x-refs
    pubchem_synonym_files <-
      sort(list.files(path = file.path(here::here(), "data-raw","pubchem"),
                      pattern = "CID-Synonym-filtered_",
                      full.names = T))

    rlogging::message("Mapping ChEMBL identifiers for NCI compounds")
    #antineopharma_synonyms_pubchem <- data.frame()
    i <- 1
    for(f in pubchem_synonym_files){
      rlogging::message("Mapping iteration..", i)
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
        dplyr::mutate(nci_concept_display_name_lc =
                        tolower(nci_concept_display_name)) |>
        dplyr::inner_join(
          chembl2alias, by = c("nci_concept_display_name_lc" = "alias"))


      rm(chembl2alias)

      if(nrow(hits) > 0){
        hits <- hits |>
          dplyr::inner_join(
            chembl2pubchem, by = "pubchem_cid") |>
          dplyr::select(nci_concept_display_name,
                        molecule_chembl_id)

        rlogging::message("Found ", nrow(hits), " ChEMBL identifiers")

        nci_compounds_chembl_match <- nci_compounds_chembl_match |>
          dplyr::bind_rows(hits)
      }
      rm(chembl2pubchem)

      i <- i + 1
    }

    nci_compounds_chembl_match_unique <- nci_compounds_chembl_match |>
      dplyr::group_by(nci_concept_display_name) |>
      dplyr::summarise(
        n_identifiers = dplyr::n(),
        molecule_chembl_id = paste(unique(molecule_chembl_id), collapse="&"),
        .groups = "drop") |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Goserelin Acetate",
        "CHEMBL1200501",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Abacavir Sulfate",
        "CHEMBL1200666",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Roniciclib",
        "CHEMBL4442620",
        as.character(molecule_chembl_id))) |>
      dplyr::filter(!stringr::str_detect(molecule_chembl_id,"&")) |>
      dplyr::select(-n_identifiers)


    # i <- 1
    # while(i <= nrow(nci_compounds_no_chembl)){
    #   name <- nci_compounds_no_chembl[i, "nci_concept_display_name"]
    #   chembl_hit <- get_chembl_compound_by_name(name)
    #   if(!is.null(chembl_hit)){
    #     chembl_hit <- chembl_hit |>
    #       dplyr::rename(nci_concept_display_name = name)
    #     nci_compounds_chembl_match <- nci_compounds_chembl_match |>
    #       dplyr::bind_rows(chembl_hit)
    #
    #   }
    #   if(i %% 10 == 0){
    #     rlogging::message("Done with querying ChEMBL for ", i, " compound names")
    #   }
    #   i <- i + 1
    # }

    nci_antineo_thesaurus_chembl <- nci_antineo_thesaurus |>
      dplyr::anti_join(nci_compounds_no_chembl,
                       by = "nci_concept_display_name")

    nci_antineo_thesaurus_no_chembl <- nci_antineo_thesaurus |>
      dplyr::inner_join(nci_compounds_no_chembl,
                        by = "nci_concept_display_name") |>
      dplyr::select(-c(num_spaces, molecule_chembl_id)) |>
      dplyr::left_join(nci_compounds_chembl_match_unique,
                       by = "nci_concept_display_name")

    nci_antineo_thesaurus <-
      nci_antineo_thesaurus_chembl |>
      dplyr::bind_rows(nci_antineo_thesaurus_no_chembl) |>
      dplyr::arrange(nci_concept_display_name) |>
      dplyr::distinct() |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Vorolanib",
        "CHEMBL3545427",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Gemtuzumab Ozogamicin",
        "CHEMBL1201506",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Ublituximab",
        "CHEMBL2108354",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Etoposide",
        "CHEMBL44657",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Lorvotuzumab Mertansine",
        "CHEMBL1743037",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Epratuzumab",
        "CHEMBL2108404",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Ocaratuzumab",
        "CHEMBL2109665",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Medroxyprogesterone Acetate",
        "CHEMBL717",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Brigatinib",
        "CHEMBL3545311",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Brentuximab",
        "CHEMBL1742994",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Afatinib",
        "CHEMBL1173655",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Spebrutinib",
        "CHEMBL3301625",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Dacomitinib",
        "CHEMBL2110732",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Lucitanib",
        "CHEMBL2220486",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Tofacitinib",
        "CHEMBL221959",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Mifamurtide",
        "CHEMBL2111100",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Umbralisib",
        "CHEMBL3948730",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Onvansertib",
        "CHEMBL1738758",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Sulindac",
        "CHEMBL15770",
        as.character(molecule_chembl_id))) |>
      dplyr::mutate(molecule_chembl_id = dplyr::if_else(
        nci_concept_display_name == "Calcipotriene",
        "CHEMBL1200666",
        as.character(molecule_chembl_id))) |>
      dplyr::distinct()

    saveRDS(nci_antineo_thesaurus, file = file.path(path_data_processed, "nci_thesaurus", "nci_treatment_thesaurus_antineo.rds"))
  }else{
    nci_antineo_thesaurus <- readRDS(file = file.path(path_data_processed, "nci_thesaurus", "nci_treatment_thesaurus_antineo.rds"))
  }
  return(nci_antineo_thesaurus)
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

get_uniprot_ensembl_map <-
  function(path_data_raw = NULL,
           uniprot_release = "2021_03") {

  ## get NCBI gene info (gene symbols, names, and entrez gene identifiers)
  gene_info <- get_gene_info_ncbi(path_data_raw = path_data_raw,
                                  update = F) |>
    dplyr::rename(target_ensembl_gene_id = ensembl_gene_id)

  ## read uniprot ID mapping
  idmapping <- read.table(gzfile(paste0(path_data_raw, '/uniprot/',
                                        uniprot_release, '/HUMAN_9606_idmapping.dat.gz')),
                          sep = "\t", header = F, quote = "", stringsAsFactors = F) |>
    dplyr::filter(V2 == 'ChEMBL' | V2 == 'UniProtKB-ID' |
                    V2 == 'Ensembl' | V2 == 'GeneID' |
                    V2 == 'Gene_Name') |>
    magrittr::set_colnames(c('acc', 'type', 'name')) |>
    dplyr::mutate(acc = stringr::str_replace(acc, "-[0-9]{1,}", ""))

  mappings <- list()
  for (m in c('ChEMBL', 'UniProtKB-ID', 'GeneID', 'Ensembl', 'Gene_Name')) {
    mappings[[m]] <- dplyr::filter(idmapping, type == m) |>
      dplyr::select(acc, name) |>
      dplyr::distinct()

    if(m == 'ChEMBL') {
      mappings[[m]] <- mappings[[m]] |> dplyr::rename(target_chembl_id = name)
    }
    if(m == 'Ensembl') {
      mappings[[m]] <- mappings[[m]] |> dplyr::rename(target_ensembl_gene_id = name)
    }
    if(m == 'GeneID') {
      mappings[[m]] <- mappings[[m]] |> dplyr::rename(entrezgene = name)
    }
    if(m == 'UniProtKB-ID') {
      mappings[[m]] <- mappings[[m]] |> dplyr::rename(uniprot_id = name)
    }
    if(m == 'Gene_Name') {
      mappings[[m]] <- mappings[[m]] |> dplyr::rename(symbol = name)
    }
  }

  mappings[['Ensembl']] <- mappings[['Ensembl']] |>
    dplyr::left_join(mappings[['GeneID']], by = "acc") |>
    dplyr::left_join(mappings[['UniProtKB-ID']], by = "acc") |>
    dplyr::left_join(mappings[['ChEMBL']], by = "acc") |>
    dplyr::left_join(mappings[['Gene_Name']], by = "acc") |>
    dplyr::select(-acc) |>
    dplyr::distinct() |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::inner_join(dplyr::select(gene_info, entrezgene, symbol, target_ensembl_gene_id),
                      by = c("target_ensembl_gene_id", "entrezgene", "symbol"))
  df <- mappings[['Ensembl']]
  return(df)
}


get_gene_info_ncbi <- function(path_data_raw = NULL,
                               update = T){

  invisible(assertthat::assert_that(
    update == T | update == F,
    msg = "'update' is not of type logical (TRUE/FALSE)"))
  invisible(assertthat::assert_that(
    dir.exists(path_data_raw),
    msg = paste0("Directory '",path_data_raw,"' does not exist")))

  rlogging::message("Retrieving gene_info from NCBI/Entrez")
  gene_info_fname <- file.path(path_data_raw,
                               "gene_info",
                               "Homo_sapiens.gene_info.gz")
  if(update == F){
    invisible(assertthat::assert_that(
      file.exists(gene_info_fname),
      msg = paste0("File ",gene_info_fname,
                   " does not exist")))
  }

  if(!file.exists(gene_info_fname) | update == T){
    remote_url <- "ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
    if(RCurl::url.exists(remote_url)){
      download.file(remote_url, destfile = gene_info_fname, quiet=T)
    }else{
      rlogging::message(paste0("Cannot update gene_info - download not available: ",
                               remote_url))
    }
  }
  gene_info <- read.table(
    gzfile(gene_info_fname),
    sep = "\t", stringsAsFactors = F,na.strings = "-",
    skip = 1, comment.char = "#", quote = "", fill = T) |>
    dplyr::filter(V1 == 9606) |>
    dplyr::select(c(V2, V3, V5, V9, V6, V10, V11)) |>
    dplyr::rename(entrezgene = V2, synonyms = V5,
                  symbol = V3, name = V9,
                  gene_biotype = V10, symbol_entrez = V11) |>
    dplyr::mutate(ensembl_gene_id = stringr::str_replace(
      stringr::str_match(V6,"Ensembl:ENSG[0-9]{1,}"), "Ensembl:", "")) |>
    dplyr::mutate(hgnc_id = stringr::str_replace(
      stringr::str_match(V6,"HGNC:HGNC:[0-9]{1,}"), "HGNC:HGNC:", "")) |>
    dplyr::select(-V6) |>
    dplyr::mutate(symbol_entrez = dplyr::if_else(
      is.na(symbol_entrez) | symbol_entrez == "-",
      symbol, as.character(symbol_entrez))) |>
    dplyr::mutate(gene_biotype = dplyr::if_else(
      gene_biotype ==  "protin-coding",
      "protein_coding", as.character(gene_biotype))) |>
    dplyr::filter(nchar(symbol_entrez) > 0)

  ### for genes annotated with the same ensembl gene ids, ignore this annotation (set to NA)
  ensgene_id_count <- as.data.frame(
    dplyr::filter(gene_info, !is.na(ensembl_gene_id)) |>
      dplyr::group_by(ensembl_gene_id) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
  )

  gene_info <- gene_info |>
    dplyr::left_join(ensgene_id_count,by=c("ensembl_gene_id")) |>
    dplyr::mutate(ensembl_gene_id = dplyr::if_else(
      !is.na(n) & n > 1,
      as.character(NA),
      as.character(ensembl_gene_id))) |>
    dplyr::select(-n) |>
    ## TEC
    dplyr::filter(entrezgene != 100124696) |>
    ## HBD
    dplyr::filter(entrezgene != 100187828) |>
    ## MMD2
    dplyr::filter(entrezgene != 100505381) |>
    ## MEMO1
    dplyr::filter(entrezgene != 7795)


  return(gene_info)

}


# dm_data_dir <-
#   "/Users/sigven/project_data/package__oncoPharmaDB/oncoPharmaDB/data-raw/dailymed/prescription"

# test_drugs <-
#   c(paste0(dm_data_dir,"/2021/1843c04f-c7f4-44a3-a8be-7b346e015616.xml") #olaparib
    # paste0(dm_data_dir,"/2020/9b009bb2-02c2-d844-e053-2a95a90a242c.xml"),#mitoxantrone
    # paste0(dm_data_dir,"/2021/2fbaadf6-c86c-48bb-bbba-803377841733.xml"), #ipilimumab
    # paste0(dm_data_dir,"/2021/226c861d-4f85-487e-a69d-e2a070ca65e0.xml"), #encorafenib
    # paste0(dm_data_dir,"/2021/bd6b1fd2-24b4-4204-ba3d-fd90a6fbc457.xml"), #nivolumab
    # paste0(dm_data_dir,"/2021/37fc5f22-5832-4293-8b7f-d316e0166c1e.xml"), #lenvatinib
    # paste0(dm_data_dir,"/2021/9ab5f41c-5fee-47a5-ab6d-12852c3eeff3.xml"), #sunitinib
    # paste0(dm_data_dir,"/2020/98f607d6-12a5-41e9-9003-dbfb500b3c11.xml"), #pazopanib
    # paste0(dm_data_dir,"/2021//5294834e-e905-4b80-9d4a-120005ef29f7.xml"), #crizotinib
    # paste0(dm_data_dir,"/2021/68805c84-fc83-4bc6-9311-2fec065cc3b8.xml"), #trastuzumab
    # paste0(dm_data_dir,"/2021/ab6f3cb3-34a8-4492-a5d7-fb6b055a2d6b.xml") #erlotinib
 # )


# all_drug_indications <- data.frame()
# for(f in test_drugs){
#   drug_indications <- get_drug_indications_dm(xml_fname = f,
#                                               cancer_drug_names = cancer_drug_names)
#   all_drug_indications <-
#     dplyr::bind_rows(all_drug_indications, drug_indications)
# }

get_dailymed_drug_indications <- function(update = F,
                                          path_data_raw = NULL){

  existing_dataset_fname <-
    paste0(path_data_raw,"/dailymed/dm_cancerdrug_indications.rds")
  if(update == F){
    if(file.exists(existing_dataset_fname)){
      dm_drugs_with_indications <- readRDS(file = existing_dataset_fname)
      return(dm_drugs_with_indications)
    }
  }

  cancer_drug_names <-
    unique(toupper(oncoPharmaDB::oncopharmadb$nci_concept_display_name))

  all_drug_indications <- data.frame()
  for(i in 16:21){
    year = paste0("20", i)
    if(i < 10){
      year = paste0("200",i)
    }

    index <- 1
    for(f in list.files(
      path = paste0("data-raw/dailymed/prescription/", year),
      pattern = ".xml", full.names = T)){

      index = index + 1

      #cat(f,'\n')
      if(index %% 50 == 0){
        cat(year,' -- processed: ', index, ' indication records\n')
      }
      drug_indications <-
        get_drug_indications_dm(xml_fname = f,
                                cancer_drug_names = cancer_drug_names)
      if(nrow(drug_indications) > 0){
        all_drug_indications <- all_drug_indications |>
          dplyr::bind_rows(drug_indications)
      }
    }
    cat(year,' --- processed: ', index, ' indication records\n')
  }


  comb_regimen_indications <- all_drug_indications |>
    dplyr::filter(cancer_indication == T) |>
    dplyr::mutate(comb_regimen_indication = dplyr::if_else(
      !is.na(indication_criteria) &
      stringr::str_detect(tolower(indication_criteria), "in combination with") &
        !stringr::str_detect(tolower(indication_criteria), "alone (or|and) in combination"),
      TRUE, FALSE, FALSE
    )) |>
    dplyr::select(drugname_generic, indication, indication_criteria,
                  comb_regimen_indication, xml_fname) |>
    dplyr::filter(comb_regimen_indication == T) |>
    dplyr::distinct()


  indications_otm_input1 <- all_drug_indications |>
    dplyr::filter(cancer_indication == T) |>
    dplyr::select(indication,
                  xml_fname) |>
    dplyr::rename(value = indication, nct_id = xml_fname) |>
    dplyr::mutate(property = "condition") |>
    dplyr::filter(!is.na(value))

  indications_otm_input2 <- all_drug_indications |>
    dplyr::filter(cancer_indication == T) |>
    dplyr::select(indication_criteria,
                  xml_fname) |>
    dplyr::rename(value = indication_criteria,
                  nct_id = xml_fname) |>
    dplyr::mutate(property = "condition") |>
    dplyr::filter(!is.na(value))

  indications_otm_input <- indications_otm_input1 |>
    dplyr::bind_rows(indications_otm_input2)

  indications_mapped <- oncoTrialMiner::index_phenotype(
    ct_data_raw = indications_otm_input
  )

  tmp <- indications_mapped$per_xref |>
    dplyr::select(nct_id, cui, primary_site, efo_id,
                  cui_name, ct_condition_raw) |>
    dplyr::rename(xml_fname = nct_id,
                  indication = ct_condition_raw) |>
    dplyr::distinct()

  dm_drugs_with_mapped_indications1 <- all_drug_indications |>
    dplyr::select(xml_fname,
                  drugname_generic,
                  drugname_trade,
                  indication) |>
    dplyr::filter(!is.na(indication)) |>
    dplyr::inner_join(tmp, by = c("xml_fname", "indication")) |>
    dplyr::distinct()

  dm_drugs_with_mapped_indications2 <- all_drug_indications |>
    dplyr::select(xml_fname,
                  drugname_generic,
                  drugname_trade,
                  indication_criteria) |>
    dplyr::rename(indication = indication_criteria) |>
    dplyr::filter(!is.na(indication)) |>
    dplyr::inner_join(tmp, by = c("xml_fname", "indication")) |>
    #dplyr::select(-c(xml_fname, indication)) |>
    dplyr::distinct()

  dm_drugs_with_indications <- as.data.frame(
    dm_drugs_with_mapped_indications2 |>
      dplyr::bind_rows(dm_drugs_with_mapped_indications1) |>

      dplyr::left_join(dplyr::select(comb_regimen_indications,
                                     drugname_generic,
                                     indication_criteria,
                                     comb_regimen_indication,
                                     xml_fname),
                       by = c("xml_fname" = "xml_fname",
                              "indication" = "indication_criteria",
                              "drugname_generic" = "drugname_generic")) |>
      dplyr::left_join(dplyr::select(comb_regimen_indications,
                                     drugname_generic,
                                     indication,
                                     comb_regimen_indication,
                                     xml_fname),
                       by = c("xml_fname" = "xml_fname",
                              "indication" = "indication",
                              "drugname_generic" = "drugname_generic")) |>
      dplyr::mutate(comb_regimen_indication = dplyr::if_else(
        !is.na(comb_regimen_indication.x) & comb_regimen_indication.x == TRUE,
        TRUE, as.logical(FALSE)
      )) |>
      dplyr::mutate(comb_regimen_indication = dplyr::if_else(
        !is.na(comb_regimen_indication.y) &
          comb_regimen_indication.y == TRUE,
        TRUE, as.logical(comb_regimen_indication)
      )) |>
      dplyr::select(-c(comb_regimen_indication.x,
                       comb_regimen_indication.y)) |>
      dplyr::distinct() |>
      dplyr::mutate(drug_name_nci = tolower(drugname_generic)) |>
      dplyr::select(-drugname_generic) |>
      dplyr::group_by(drug_name_nci, cui, efo_id,
                      cui_name, primary_site,
                      comb_regimen_indication) |>
      dplyr::summarise(drugname_trade =
                         paste(unique(drugname_trade), collapse="|"),
                       drug_clinical_id = paste(unique(
                         stringr::str_replace(
                           xml_fname,"\\.xml$",""
                         )), collapse=","),
                       .groups = "drop") |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "renal cell carcinoma"),
        "EFO:0000681", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "urothelial carcinoma"),
        "EFO:0008528", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "mantle cell lymphoma"),
        "EFO:1001469", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "non-small cell lung cancer"),
        "EFO:0003060", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "gastrointestinal stromal"),
        "Orphanet:44890", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "mediastinal large b-cell lymphoma"),
        "MONDO:0004021", as.character(efo_id)
      )) |>
      dplyr::mutate(efo_id = dplyr::if_else(
        (is.na(efo_id) | efo_id == "NA") &
          stringr::str_detect(tolower(cui_name), "pancreatic neuroendocrine tumor"),
        "EFO:1000045", as.character(efo_id)
      )) |>
      dplyr::filter(!is.na(cui) & efo_id != "NA")
  )

  dm_drugs_with_indications <- as.data.frame(
    dm_drugs_with_indications |>
      tidyr::separate_rows(efo_id,sep="\\|") |>
      dplyr::left_join(dplyr::select(
        oncoPhenoMap::auxiliary_maps$efo$efo2name, efo_id, efo_name
      ), by = "efo_id") |>
      dplyr::rename(disease_efo_id = efo_id,
                    disease_efo_label = efo_name) |>
      #dplyr::mutate(drug_name_nci = tolower(drugname_generic)) |>
      dplyr::mutate(drug_max_phase_indication = NA,
                    drug_clinical_source = "DailyMedParseSN",
                    #drug_clinical_id = NA,
                    drug_max_ct_phase = NA,
                    drug_approved_indication = TRUE) |>
      dplyr::distinct()
  )

  saveRDS(dm_drugs_with_indications, existing_dataset_fname)

  return(dm_drugs_with_indications)
}




get_drug_indications_dm <- function(xml_fname = NULL,
                                    cancer_drug_names = NULL){

  invisible(assertthat::assert_that(
    file.exists(xml_fname),
    msg = paste0("File name ", xml_fname," does not exist")))
  con <- file(description = xml_fname)
  tmp <- xml2::read_xml(con, options = c("HUGE","NOBLANKS"))
  xml2::xml_ns_strip(tmp)


  drugname_trade <-  stringr::str_trim(xml2::xml_text(xml2::xml_find_all(
    tmp, "//manufacturedProduct/name")[1]))
  if(length(drugname_trade) == 0){
    drugname_trade <- NA
  }

  drug_active_moiety <- xml2::xml_text(xml2::xml_find_all(
    tmp, "//manufacturedProduct/ingredient/ingredientSubstance/activeMoiety/activeMoiety/name")[1])
  if(length(drug_active_moiety) == 0){
    drug_active_moiety <- NA
  }

  initial_approval_date <- NA
  full_title <- ""
  title_elements <- xml2::as_list(xml2::xml_find_first(tmp,".//title"))

  if(is.list(title_elements)){
    if(length(title_elements) > 0){
      for(i in 1:length(title_elements)){
        title_data <- stringr::str_squish(stringr::str_trim(unlist(title_elements[[i]])))
        if(is.list(title_data)){
          cat('BALLE\n')
        }
        if(is.null(title_data)){
          next
        }
        if(!is.character(title_data)){
          next
        }
        full_title <- paste0(full_title, title_data)
      }
    }
  }

  #cat(full_title, '\n')
  if(full_title != ""){
    t <- stringr::str_match(full_title,"Initial U.S. Approval: [0-9]{4}")
    if(length(t) == 1){
      initial_approval_date <- stringr::str_trim(
        stringr::str_replace(
          t, "Initial U.S. Approval: ",""
        ))
    }
  }

  drugname_generic <- xml2::xml_text(xml2::xml_find_all(
    tmp, "//manufacturedProduct/asEntityWithGeneric/genericMedicine/name")[1])
  if(length(drugname_generic) == 0){
    drugname_generic <- NA
  }
  sections <- xml2::xml_find_all(tmp,".//component/section")

  indication_section <- 0
  section_counter <- 0
  all_indication_entries <- data.frame()
  all_main_indication_entries <- data.frame()
  effective_time <- NA

  is_cancer_drug <- F
  if(!is.null(cancer_drug_names) &
     !is.na(drugname_trade)){
    if(toupper(drugname_trade) %in% cancer_drug_names){
      is_cancer_drug <- T
    }
  }

  if(!is.null(cancer_drug_names) &
     !is.na(drugname_generic)){
    if(toupper(drugname_generic) %in% cancer_drug_names){
      is_cancer_drug <- T
    }
  }

  if(!is.null(cancer_drug_names) &
     !is.na(drug_active_moiety)){
    if(toupper(drug_active_moiety) %in% cancer_drug_names){
      is_cancer_drug <- T
    }
  }

  if(is_cancer_drug == F){
    return(all_indication_entries)
  }

  #cat(drugname_generic, '\n')

  for(s in sections){
    section_counter <- section_counter + 1
    section_titles <- xml2::xml_text(xml2::xml_find_first(s, ".//title"))
    for(sec_title in unique(section_titles)){
      sec_title <- stringr::str_squish(stringr::str_trim(sec_title))
      if(!is.na(sec_title) & nchar(sec_title) > 0){

        if(stringr::str_detect(sec_title,"INDICATIONS") &
           stringr::str_detect(sec_title,"USAGE")){
          indication_section <- indication_section + 1

          indication_info_main <- xml2::as_list(sections[section_counter][[1]])
          effective_time <- attributes(
            indication_info_main$effectiveTime)$value

          if("excerpt" %in% names(indication_info_main)){
            if("highlight" %in% names(indication_info_main$excerpt)){
              if("text" %in% names(indication_info_main$excerpt$highlight)){

                all_indication_contents <-
                  unlist(indication_info_main$excerpt$highlight$text)

                names_indication_contents <- names(all_indication_contents)
                names(all_indication_contents) <- NULL

                indication <- ""
                indication_criteria <- NA
                prev_cname <- ""
                if(length(all_indication_contents) > 0 &
                   length(names_indication_contents) > 0 &
                  length(all_indication_contents) == length(names_indication_contents)){
                  for(i in 1:length(all_indication_contents)){
                    cname <- names_indication_contents[i]

                    df <- data.frame('cname' = cname,
                                     'indication' = all_indication_contents[i],
                                     'indication_criteria' = NA,
                                     'drugname_generic' = drugname_generic,
                                     'drugname_trade' = drugname_trade,
                                     'drug_active_moiety' = drug_active_moiety,
                                     'effective_time' = effective_time,
                                     stringsAsFactors = F)

                    all_main_indication_entries <- all_main_indication_entries |>
                      dplyr::bind_rows(df)
                  }
                }
              }
            }
          }else{
            if("text" %in% names(indication_info_main)){

              all_indication_contents <-
                unlist(indication_info_main$text)

              if(!is.null(all_indication_contents)){
                names_indication_contents <- names(all_indication_contents)
                names(all_indication_contents) <- NULL
                if(length(all_indication_contents) > 0 &
                   length(names_indication_contents) > 0 &
                   length(all_indication_contents) == length(names_indication_contents)){
                  for(i in 1:length(all_indication_contents)){
                    cname <- names_indication_contents[i]

                    df <- data.frame('cname' = cname,
                                     'indication' = all_indication_contents[i],
                                     'indication_criteria' = NA,
                                     'drugname_generic' = drugname_generic,
                                     'drugname_trade' = drugname_trade,
                                     'drug_active_moiety' = drug_active_moiety,
                                     'effective_time' = effective_time,
                                     stringsAsFactors = F)

                    all_main_indication_entries <- all_main_indication_entries |>
                      dplyr::bind_rows(df)
                  }
                }
              }
            }
          }
        }
        if(indication_section > 0 &
           !(stringr::str_detect(
             sec_title, "INDICATIONS")) &
           !stringr::str_detect(
             sec_title, "DOSAGE"
           )){

          df_sec_title <- data.frame('cname' = 'section_title',
                                     'indication' = sec_title,
                                     'indication_criteria' = NA,
                                     'drugname_generic' = drugname_generic,
                                     'drugname_trade' = drugname_trade,
                                     'drug_active_moiety' = drug_active_moiety,
                                     'effective_time' = effective_time,
                                     stringsAsFactors = F)


          section_indication_details_found <- 0

          indication_details <- xml2::as_list(sections[section_counter][[1]])

          if("text" %in% names(indication_details)){

            if(length(indication_details$text) == 1){
              if("paragraph" %in% names(indication_details$text)){

                ind_details <- unlist(indication_details$text$paragraph)
                attributes(ind_details) <- NULL
                ind_details <- ind_details[!nchar(ind_details) == 1]
                ind_details <- stringr::str_squish(paste(ind_details, collapse=""))

                df <- data.frame('cname' = 'paragraph',
                                 'indication' = sec_title,
                                 'indication_criteria' = ind_details,
                                 'drugname_generic' = drugname_generic,
                                 'drugname_trade' = drugname_trade,
                                 'drug_active_moiety' = drug_active_moiety,
                                 'effective_time' = effective_time,
                                 stringsAsFactors = F)

                all_indication_entries <- all_indication_entries |>
                  dplyr::bind_rows(df)

                section_indication_details_found <- 1

              }
            }else{
              if(!("list" %in% names(indication_details$text))){
                if(length(indication_details$text) == 0){
                  next
                }
                for(n in 1:length(indication_details$text)){
                  if("content" %in% names(indication_details$text[[n]])){

                    if(is.list(indication_details$text[[n]]$content) &
                       length(indication_details$text[[n]]$content) == 1){

                      if(is.character(indication_details$text[[n]]$content[[1]])){

                        df <- data.frame(
                          'cname' = 'section_title',
                          'indication' = sec_title,
                          'indication_criteria' = indication_details$text[[n]]$content[[1]],
                          'drugname_generic' = drugname_generic,
                          'drugname_trade' = drugname_trade,
                          'drug_active_moiety' = drug_active_moiety,
                          'effective_time' = effective_time,
                          stringsAsFactors = F)

                        all_indication_entries <- all_indication_entries |>
                          dplyr::bind_rows(df)
                      }

                      section_indication_details_found <- 1

                    }else{
                      crit <- unlist(indication_details$text[[n]])
                      if(is.character(crit) & length(crit) == 1){
                        names(crit) <- NULL
                        if(nchar(crit[1]) > 1 &
                           !(stringr::str_detect(
                             crit[1],"Clinical Studies|Dosage and Administration"))){
                          df <- data.frame(
                            'cname' = 'section_title',
                            'indication' = sec_title,
                            'indication_criteria' = crit[1],
                            'drugname_generic' = drugname_generic,
                            'drugname_trade' = drugname_trade,
                            'drug_active_moiety' = drug_active_moiety,
                            'effective_time' = effective_time,
                            stringsAsFactors = F)

                          all_indication_entries <- all_indication_entries |>
                            dplyr::bind_rows(df)

                          section_indication_details_found <- 1

                        }
                      }
                    }
                  }else{

                    if(is.list(indication_details$text[[n]]) & length(indication_details$text[[n]]) > 0){
                      #cat(section_counter, '--', n, '---', length(indication_details$text[[n]]), '\n')
                      for(i in 1:length(indication_details$text[[n]])){
                        if(is.character(indication_details$text[[n]][[i]])){
                          df <- data.frame(
                            'cname' = 'section_title',
                            'indication' = sec_title,
                            'indication_criteria' = indication_details$text[[n]][[i]],
                            'drugname_generic' = drugname_generic,
                            'drugname_trade' = drugname_trade,
                            'drug_active_moiety' = drug_active_moiety,
                            'effective_time' = effective_time,
                            stringsAsFactors = F)

                          all_indication_entries <- all_indication_entries |>
                            dplyr::bind_rows(df)

                          section_indication_details_found <- 1

                        }
                      }
                    }
                  }
                }
              }
            }

            if("paragraph" %in% names(indication_details$text) &
               length(indication_details$text) > 1){

              indication_title <- ""
              indication_title_elements <- c()
              if(length(indication_details$text$paragraph) > 1){
                for(j in 1:length(indication_details$text$paragraph)){
                  if(is.character(indication_details$text$paragraph[[j]])){
                    indication_title_elements <-
                      c(indication_title_elements,
                        indication_details$text$paragraph[[j]])
                  }
                }
                indication_title <- paste(indication_title_elements, collapse="")
              }else{
                if(length(indication_details$text$paragraph) == 1){
                  if(is.character(indication_details$text$paragraph[[1]])){
                    indication_title <- as.character(indication_details$text$paragraph[[1]])
                  }
                }
              }
              indication_criterions <- c()

              if("list" %in% names(indication_details$text)){

                if("item" %in% names(indication_details$text$list)){

                  if(length(indication_details$text$list) == 1){
                    indication_items <-
                      unlist(indication_details$text$list$item)

                    if(!is.null(indication_items) &
                       is.vector(indication_items)){
                      ind_crit <- ""
                      #cat('--', unlist(names(indication_items)),'\n')

                      if(length(indication_items) > 0){
                        for(m in 1:length(indication_items)){
                          if(is.null(names(indication_items[m]))){
                            next
                          }
                          if(names(indication_items[m]) == "" |
                             stringr::str_detect(names(indication_items[m]),"content([0-9])?$")){
                            ind_crit <- paste(ind_crit, stringr::str_trim(indication_items[m]))
                          }
                        }
                      }
                      indication_criterions <- c(indication_criterions,
                                                 paste(indication_title,
                                                       ind_crit, sep=" "))
                    }
                  }else{

                    for(m in 1:length(indication_details$text$list)){
                      #if(nchar(unlist(indication_details$text$list[[m]])) > 1){
                      indication_criterions <-
                        c(indication_criterions,
                          paste(indication_title,
                                unlist(indication_details$text$list[[m]]),
                                sep = " - ")
                        )
                      #}
                    }
                  }
                }
                else{
                  for(j in 1:length(indication_details$text$list)){

                    if(is.character(unlist(indication_details$text$list[[j]][[1]]))){
                      if(nchar(unlist(indication_details$text$list[[j]][[1]])) > 1){
                        indication_criterions <-
                          c(indication_criterions,
                            paste(
                              indication_title,
                              unlist(indication_details$text$list[[j]][[1]]),
                              sep = " - ")
                          )
                      }
                    }
                  }
                }
              }

              for(c in indication_criterions){

                df <- data.frame(
                  'cname' = 'section_title',
                  'indication' = sec_title,
                  'indication_criteria' = c,
                  'drugname_generic' = drugname_generic,
                  'drugname_trade' = drugname_trade,
                  'drug_active_moiety' = drug_active_moiety,
                  'effective_time' = effective_time,
                  stringsAsFactors = F)

                all_indication_entries <- all_indication_entries |>
                  dplyr::bind_rows(df)

                section_indication_details_found <- 1

              }
            }

            else{

              if(length(indication_details$text) >= 1){
                type <- unique(names(indication_details$text))
                if(length(type) != 1){
                  break
                }
                if(unique(names(indication_details$text)) == "list"){
                  all_criterions <- c()

                  for(j in 1:length(indication_details$text$list)){
                    if(length(indication_details$text$list[[j]]) == 0){
                      next
                    }
                    for(k in 1:length(indication_details$text$list[[j]])){
                      indication_criterion <- ""
                      if(is.list(indication_details$text$list[[j]][[k]])){

                        if(length(indication_details$text$list[[j]][[k]]) == 1){
                          if(is.character(unlist(indication_details$text$list[[j]][[k]])) &
                             nchar(unlist(indication_details$text$list[[j]][[k]])) > 1){
                            indication_criterion <- as.character(
                              unlist(indication_details$text$list[[j]][[k]]))

                            all_criterions <- c(all_criterions, indication_criterion)
                          }
                        }else{
                          indication_criterion <- ""
                          if(length(indication_details$text$list[[j]][[k]]) == 0){
                            next
                          }
                          for(l in 1:length(indication_details$text$list[[j]][[k]])){
                            if(is.character(indication_details$text$list[[j]][[k]][[l]])){
                              indication_criterion <- paste0(
                                indication_criterion,
                                as.character(
                                  indication_details$text$list[[j]][[k]][[l]]))
                            }
                          }
                          all_criterions <- c(all_criterions, indication_criterion)


                        }
                      }else{
                        indication_criterion <-
                          unlist(indication_details$text$list[[j]][[k]])
                        all_criterions <- c(all_criterions, indication_criterion)

                      }
                    }
                  }

                  for(c in all_criterions){
                    df <- data.frame(
                      'cname' = 'section_list',
                      'indication' = sec_title,
                      'indication_criteria' = c,
                      'drugname_generic' = drugname_generic,
                      'drugname_trade' = drugname_trade,
                      'drug_active_moiety' = drug_active_moiety,
                      'effective_time' = effective_time,
                      stringsAsFactors = F)

                    all_indication_entries <- all_indication_entries |>
                      dplyr::bind_rows(df)

                    section_indication_details_found <- 1

                  }
                }
              }
            }
          }
          if(section_indication_details_found == 0){
            all_indication_entries <- all_indication_entries |>
              dplyr::bind_rows(df_sec_title)
          }
        }

        if(stringr::str_detect(sec_title,"CONTRAINDICATIONS") |
           stringr::str_detect(sec_title,"WARNINGS|DOSAGE")){
          indication_section <- 0
          break
        }
      }
    }
  }



  if(nrow(all_indication_entries) == 0 &
     nrow(all_main_indication_entries) > 0){
    all_indication_entries <- all_indication_entries |>
      dplyr::bind_rows(all_main_indication_entries)
  }

  if(nrow(all_indication_entries) == 0){
    return(all_indication_entries)
  }

  #cancer_regex <- "cancer|carcinoma|wilms|esophageal|melanoma|myeloma"
  #cancer_regex <- paste0(cancer_regex,"|leukemia|lymphom|sarcoma|thymom|thyroid cancer")
  #cancer_regex <- paste0(cancer_regex,"|head and neck cancer|")

  cancer_tissue_regex <- "endometrial|prostate|rectal|myeloid|anal|colorectal|wilms|kidney|ovarian|"
  cancer_tissue_regex <- paste0(
    cancer_tissue_regex,"thyroid|esophag|cervical|cervix|gastric|stomach|head and neck|ovary|")
  cancer_tissue_regex <- paste0(
   cancer_tissue_regex,"bone|testis|testicular|lung|liver|skin|pancreas|pancreatic|brain|bladder|")
  cancer_tissue_regex <- paste0(
    cancer_tissue_regex,"adrenal gland|renal cell|urothelial|lympho|hepatocellular|breast")

  cancer_other_regex <- "melanoma|myeloma|glioblastoma|sarcoma|leukemia|lymphoma|"
  cancer_other_regex <- paste0(
    cancer_other_regex, "cholangiocarcinom|mesotheliom|glioma|neuroblastoma|tumor")

  rownames(all_indication_entries) <- NULL
  all_indication_entries <- all_indication_entries |>
    dplyr::mutate(xml_fname = basename(xml_fname)) |>
    dplyr::mutate(initial_approval = initial_approval_date) |>
    dplyr::mutate(indication = stringr::str_trim(stringr::str_squish(indication))) |>
    dplyr::mutate(indication_criteria =
                    stringr::str_trim(stringr::str_squish(indication_criteria))) |>
    dplyr::mutate(cancer_indication = dplyr::if_else(
      !is.na(indication) &
        (stringr::str_detect(indication,"NSCLC|AML|CML|BRAF|KRAS|EGFR|HER2|ROS1|BRCA|ABL") |
      (stringr::str_detect(tolower(indication), cancer_tissue_regex) &
        stringr::str_detect(tolower(indication), "cancer|tumor|carcinoma|neoplasm|neoplastic"))),
      TRUE,FALSE
    )) |>
    dplyr::mutate(cancer_indication = dplyr::if_else(
      !is.na(indication) &
      stringr::str_detect(tolower(indication), cancer_other_regex),
      TRUE,as.logical(cancer_indication)
    )) |>
    dplyr::mutate(cancer_indication = dplyr::if_else(
      !is.na(indication_criteria) &
        (stringr::str_detect(indication_criteria,"NSCLC|AML|CML|BRAF|KRAS|EGFR|HER2|ROS1|BRCA|ABL") |
      (stringr::str_detect(tolower(indication_criteria), cancer_tissue_regex) &
        stringr::str_detect(tolower(indication_criteria), "cancer|tumor|carcinoma|neoplasm|neoplastic"))),
      TRUE,as.logical(cancer_indication)
    )) |>
    dplyr::mutate(cancer_indication = dplyr::if_else(
      !is.na(indication_criteria) &
      stringr::str_detect(tolower(indication_criteria), cancer_other_regex),
      TRUE,as.logical(cancer_indication)
    )) |>
    dplyr::filter(!stringr::str_detect(indication,"indicated( for)?:$")) |>
    dplyr::filter(!(stringr::str_detect(indication,"indicated for") &
                      stringr::str_detect(indication,":$"))) |>
    dplyr::filter(!stringr::str_detect(indication,"^Limitations of (U|u)se")) |>
    dplyr::filter(!stringr::str_detect(indication,"Limitations of (U|u)se")) |>
    dplyr::filter(!stringr::str_detect(indication,"Select patients for therapy based on")) |>
    dplyr::mutate(indication = Hmisc::capitalize(
      stringr::str_replace(
        indication,"^( )?:( )?",""))) |>
    dplyr::mutate(indication_criteria = Hmisc::capitalize(
      stringr::str_replace(
        indication_criteria,"^( )?:( )?",""))) |>
    dplyr::filter(is.na(indication_criteria) |
                    (!is.na(indication_criteria) & nchar(indication_criteria) >= 3)) |>
    dplyr::distinct()

  return(all_indication_entries)

}


get_uniprot_map <- function(basedir = NULL,
                            uniprot_release = "2020_03"){

  ## corum protein-complex data

  #download.file("http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip",destfile="data-raw/chorum/coreComplexes.txt.zip", quiet = T)
  # corum_complexes <-
  #   read.table(file=paste0(basedir,"/data-raw/corum/coreComplexes.txt"),
  #              sep="\t",header=T, quote="",comment.char="",stringsAsFactors = F) |>
  #   dplyr::select(c(Organism, ComplexID, subunits.UniProt.IDs.,
  #                   ComplexName,FunCat.description)) |>
  #   dplyr::filter(Organism == 'Human') |>
  #   dplyr::mutate(n_complexes = nrow(.)) |>
  #   dplyr::rename(uniprot_acc = subunits.UniProt.IDs.,
  #                 corum_id = ComplexID, corum_name = ComplexName,
  #                 corum_funcat_description = FunCat.description) |>
  #   tidyr::separate_rows(uniprot_acc, sep=";") |>
  #   dplyr::select(-Organism)

  # n_complexes <- unique(corum_complexes$n_complexes)
  # n_protein_complex_interactions <- nrow(corum_complexes)
  # rlogging::message("Retrieving protein complex annotation from CORUM")
  # rlogging::message("A total of ", n_complexes,
  #                   " protein complexes were parsed, containing ",
  #                   n_protein_complex_interactions," protein-complex interactions")
  #
  # uniprot_acc_to_corum_ids <- as.data.frame(
  #   dplyr::group_by(corum_complexes, uniprot_acc) |>
  #     dplyr::summarise(corum_id = paste(unique(corum_id), collapse = "&"), .groups = "drop"))
  rlogging::message("Retrieving UniProtKB annotation")

  ## read uniprot ID mapping
  idmapping_fname <- paste0(basedir,'/data-raw/uniprot/',uniprot_release,'/HUMAN_9606_idmapping.dat.gz')
  idmapping_up_kb <- read.table(gzfile(idmapping_fname),sep="\t",header = F,quote = "", stringsAsFactors = F) |>
    dplyr::filter(V2 == 'Ensembl_TRS' | V2 == 'UniProtKB-ID' | V2 == 'Ensembl' | V2 == 'GeneID' | V2 == 'RefSeq_NT') |>
    magrittr::set_colnames(c('acc','type','name')) |>
    dplyr::mutate(acc = stringr::str_replace(acc,"-[0-9]{1,}",""))

  reviewed_uniprot_accs <-
    readRDS(file=paste0(basedir,"/data-raw/uniprot/",uniprot_release,"/reviewed_uniprot_accessions.rds"))

  ## UniProt accession to Ensembl transcript ID
  ensembl_up_acc <- dplyr::filter(idmapping_up_kb, type == 'Ensembl_TRS') |>
    dplyr::select(acc, name) |>
    dplyr::rename(uniprot_acc = acc, ensembl_transcript_id = name) |>
    dplyr::distinct()
  ## remove transcripts mapped to more than one UniProt acc
  unique_ensembl_trans <- as.data.frame(
    dplyr::group_by(ensembl_up_acc, ensembl_transcript_id) |>
      dplyr::summarise(n_trans = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n_trans == 1))
  ensembl_up_acc <- dplyr::semi_join(
    ensembl_up_acc,
    dplyr::select(unique_ensembl_trans, ensembl_transcript_id),
    by=c("ensembl_transcript_id"))

  ## UniProt accession to UniProt ID
  up_id_acc <- as.data.frame(
    dplyr::filter(idmapping_up_kb, type == 'UniProtKB-ID') |>
      dplyr::group_by(acc) |>
      dplyr::summarise(uniprot_id = paste(unique(name), collapse="&"),
                       .groups = "drop") |>
      dplyr::rename(uniprot_acc = acc))

  ## UniProt accession to UniProt ID
  ensembl_id_acc <- as.data.frame(
    dplyr::filter(idmapping_up_kb, type == 'Ensembl') |>
      dplyr::group_by(acc) |>
      dplyr::summarise(ensembl_gene_id = paste(unique(name), collapse="&"),
                       .groups = "drop") |>
      dplyr::rename(uniprot_acc = acc))

  rlogging::message("A total of ",nrow(up_id_acc)," UniProt protein accessions were parsed")


  ## UniProt accession to RefSeq ID
  refseq_id_acc <- as.data.frame(
    dplyr::filter(idmapping_up_kb, type == 'RefSeq_NT') |>
      dplyr::mutate(name = stringr::str_replace(name,"\\.[0-9]{1,}$","")) |>
      dplyr::group_by(acc) |>
      dplyr::summarise(refseq_mrna = paste(unique(name), collapse="&"),
                       .groups = "drop") |>
      dplyr::rename(uniprot_acc = acc))


  uniprot_map <- dplyr::full_join(
    dplyr::left_join(refseq_id_acc, up_id_acc,by=c("uniprot_acc")),
    dplyr::left_join(ensembl_up_acc, up_id_acc,by=c("uniprot_acc")),
    by=c("uniprot_acc","uniprot_id")) |>
    dplyr::left_join(ensembl_id_acc)

  ## UniProt accession and ID for Ensembl transcript IDs
  uniprot_map <- uniprot_map |>
    #dplyr::left_join(uniprot_acc_to_corum_ids,by=c("uniprot_acc")) |>
    dplyr::left_join(reviewed_uniprot_accs, by=c("uniprot_acc","uniprot_id")) |>
    dplyr::rename(uniprot_reviewed = reviewed)

  return(list('uniprot_map' = uniprot_map))

}

