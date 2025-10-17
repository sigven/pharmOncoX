
clean_variant_consequence <- function(variant_df){
  
  variant_df <- variant_df |>
    dplyr::mutate(
      variant_consequence =
        dplyr::if_else(
          stringr::str_detect(variant,"Overexpression") &
            variant_consequence == "N/A",
          "transcript_amplification",
          as.character(variant_consequence))) |>
    dplyr::mutate(
      variant_consequence =
        dplyr::if_else(
          stringr::str_detect(variant,"Underexpression") &
            variant_consequence == "N/A",
          "transcript_ablation",
          as.character(variant_consequence))) |>
    dplyr::mutate(
      variant_consequence =
        dplyr::if_else(
          stringr::str_detect(
            variant_consequence,
            "^(plus_1_frameshift|minus_1_frameshift|frameshift_)"),
          "frameshift_variant",
          as.character(variant_consequence))) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      variant_consequence == "" &
        stringr::str_detect(variant,"EXON") &
        stringr::str_detect(variant,"DELETION"),
      "exon_loss_variant",
      as.character(variant_consequence)
    )) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      variant_consequence == "" &
        stringr::str_detect(variant,"Exon") &
        stringr::str_detect(variant,"Mutation"),
      "exon_variant",
      as.character(variant_consequence)
    )) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      variant_consequence == "" &
        stringr::str_detect(variant,"^(Loss|Loss-of-function)$"),
      "loss_of_function_variant",
      as.character(variant_consequence)
    )) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      stringr::str_detect(variant,"^(Activating |Gain-of)") &
        (is.na(variant_consequence) | nchar(variant_consequence) == 0),
      "gain_of_function_variant",
      as.character(variant_consequence)
    )) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      stringr::str_detect(variant,"^Wildtype") &
        (is.na(variant_consequence) | nchar(variant_consequence) == 0),
      "wild_type",
      as.character(variant_consequence)
    ))

  return(variant_df)
}

set_alteration_type <- function(variant_df){
  
  variant_df <- variant_df |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(toupper(variant),"EXPRESSION"),
          "EXP",as.character(alteration_type))) |>
    dplyr::mutate(
      alteration_type =
        dplyr::case_when(
          alteration_type == "EXP" &
            stringr::str_detect(
              toupper(variant),"OVER") ~ "EXP_OVER",
          alteration_type == "EXP" &
            stringr::str_detect(
              toupper(variant),"UNDER") ~ "EXP_UNDER",
          TRUE ~ as.character(alteration_type)
        )
    ) |>
    dplyr::mutate(alteration_type = dplyr::if_else(
      variant_consequence == "stop_gained",
      "MUT_LOF",as.character(alteration_type)
    )) |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(
            toupper(variant),"LOSS-OF-FUNCTION|LOSS|TRUNCAT|INACTIVAT") |
            stringr::str_detect(
              variant_consequence, 
              "stop_gained,loss_of_function_variant|loss_of_function_variant,stop_gained") |
            (stringr::str_detect(
              variant_consequence, "loss_of_function_variant") &
               toupper(variant) == "MUTATION") |
            stringr::str_detect(variant_consequence,"truncation"),
          "MUT_LOF",as.character(alteration_type))) |>
    
    dplyr::mutate(alteration_type = dplyr::if_else(
      toupper(variant) == "BIALLELIC INACTIVATION",
      "MUT_LOF_BIALLELIC",
      as.character(alteration_type)
    )) |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          alteration_type == 'MUT_LOF' &
            (stringr::str_detect(
              toupper(variant_consequence),"FRAMESHIFT") |
               stringr::str_detect(
                 variant_consequence,"feature_truncation")),
          "MUT_LOF_FS",as.character(alteration_type))) |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(
            toupper(variant), 
            "AMPLIFICATION|DELETION|COPY") &
            stringr::str_detect(
              variant_consequence, "_amplification|_ablation"),
          "CNA",as.character(alteration_type))) |>
    
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(
            toupper(variant_consequence),
            "TRANSCRIPT_FUSION|REGION_FUSION") |
            stringr::str_detect(toupper(variant),"FUSION") |
            (variant_consequence == "" & 
               stringr::str_detect(variant,"::")),
          "FUSION",
          as.character(alteration_type))) |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(
            toupper(variant),"REARRANGEMENT") &
            variant_consequence == "",
          "REARRANGEMENT",
          as.character(alteration_type)))
  
  return(variant_df)
  
  
}

clean_variant_naming <- function(
    variant_df, col = "molecular_profile_name"){
  
  variant_df[,"cc_clean"] <- variant_df[,col]
  variant_df <- variant_df |>
    dplyr::mutate(cc_clean = stringr::str_trim(cc_clean)) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean,"RS[0-9]{5,}"),
      stringr::str_replace_all(cc_clean,"RS","rs"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = stringr::str_replace_all(
      cc_clean, 
      paste0(
        "(ENST[0-9]{1,}|",
        "NM_[0-9]{1,})(\\.[0-9]{1,})?:c\\."),"c."
    )) |>
    dplyr::mutate(cc_clean = stringr::str_replace_all(
      cc_clean, 
      paste0(
        "(ENSP[0-9]{1,}|NP_[0-9]{1,}",
        ")(\\.[0-9]{1,})?:p\\."),"p."
    )) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "^NC_[0-9]{1,}"),
      "", as.character(cc_clean)
    )) |>
    dplyr::mutate(cc_clean = dplyr::case_when(
      cc_clean == "D835H/Y" ~ "D835H,D835Y",
      cc_clean == "G12/G13" ~  "G12,G13",
      cc_clean == "Q157P/R" ~ "Q157P,Q157R",
      cc_clean == "S310F/Y" ~ "S310F,S310Y",
      cc_clean == "S34Y/F" ~ "S34Y,S34F",
      cc_clean == "S893A/T" ~ "S893A,S893T",
      TRUE ~ as.character(cc_clean)
    )) |> 
   
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, "^(CD44s|p16|PRPS1|IKZF1|NUCLEAR) "),
      stringr::str_replace_all(
        cc_clean, "^(CD44s|p16|PRPS1|IKZF1|NUCLEAR) ", ""),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, " NUCLEAR "),
      stringr::str_replace_all(
        cc_clean, " NUCLEAR ", ""),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "IVS2\\+1G>A"),
      stringr::str_replace_all(
        cc_clean, "^IVS2\\+1G>A", "c.444+1G>A"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, " HOMOZYGOSITY"),
      stringr::str_replace_all(
        cc_clean, " HOMOZYGOSITY", ""),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, "deletion and mutation"),
      stringr::str_replace_all(
        cc_clean, "deletion and mutation", "Mutation"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = stringr::str_replace(
      cc_clean, " \\(.+\\)$","")) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "\\*$"),
      stringr::str_replace_all(cc_clean, "\\*$", "X"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "^\\*"),
      stringr::str_replace_all(cc_clean, "^\\*", "X"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "\\* "),
      stringr::str_replace_all(cc_clean, "\\* ", "X "),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "Rare "),
      stringr::str_replace_all(cc_clean, "Rare ", ""),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, " expression"),
      stringr::str_replace_all(
        cc_clean, " expression", " Expression"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "BIALLELIC"),
      stringr::str_replace_all(
        cc_clean, "BIALLELIC", "Biallelic"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "MUTATION|mutation"),
      stringr::str_replace_all(
        cc_clean, "MUTATION|mutation", "Mutation"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "DOMAIN"),
      stringr::str_replace_all(
        cc_clean, "DOMAIN", "Domain"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "CONSERVED"),
      stringr::str_replace_all(
        cc_clean, "CONSERVED", "Conserved"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "DELETION"),
      stringr::str_replace_all(
        cc_clean, "DELETION", "Deletion"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, "FRAMESHIFT|FRAME SHIFT"),
      stringr::str_replace_all(
        cc_clean, "FRAMESHIFT|FRAME SHIFT", "Frameshift"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "INACTIVATING"),
      stringr::str_replace_all(
        cc_clean, "INACTIVATING", "Inactivating"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "INACTIVATION"),
      stringr::str_replace_all(
        cc_clean, "INACTIVATION", "Inactivation"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "OVEREXPRESSION"),
      stringr::str_replace_all(
        cc_clean, "OVEREXPRESSION", "Overexpression"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "UNDEREXPRESSION"),
      stringr::str_replace_all(
        cc_clean, "UNDEREXPRESSION", "Underexpression"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "EXPRESSION"),
      stringr::str_replace_all(
        cc_clean, "EXPRESSION", "Expression"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "EXON"),
      stringr::str_replace_all(cc_clean, "EXON", "Exon"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "LOSS$"),
      stringr::str_replace_all(cc_clean, "LOSS", "Loss"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "AMPLIFICATION"),
      stringr::str_replace_all(
        cc_clean, "AMPLIFICATION", "Amplification"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "TRUNCATING"),
      stringr::str_replace_all(
        cc_clean, "TRUNCATING", "Truncating"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "TRUNCATION"),
      stringr::str_replace_all(
        cc_clean, "TRUNCATION", "Truncation"),
      as.character(cc_clean))) |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(
        cc_clean, "TRANSLOCATION|rearrangement"),
      stringr::str_replace_all(
        cc_clean, "TRANSLOCATION|rearrangement", "Rearrangement"),
      as.character(cc_clean)))
  
  if(col != "molecular_profile_name"){
    variant_df <- variant_df |>
      dplyr::mutate(cc_clean = dplyr::if_else(
        stringr::str_detect(cc_clean, "Binding Domain Mutation") &
          feature_name == "TP53",
        "aa_region:100-288",
        as.character(cc_clean))) |>
      dplyr::mutate(cc_clean = dplyr::if_else(
        stringr::str_detect(cc_clean, "SH2 Domain") &
          feature_name == "STAT3",
        "aa_region:584-651",
        as.character(cc_clean))) |>
      dplyr::mutate(cc_clean = dplyr::if_else(
        stringr::str_detect(cc_clean, "Kinase Domain") &
          feature_name == "ERBB2",
        "aa_region:722-975",
        as.character(cc_clean))) |>
      dplyr::mutate(cc_clean = dplyr::if_else(
        stringr::str_detect(cc_clean, "TKD ") &
          feature_name == "FLT3",
        "aa_region:610-942",
        as.character(cc_clean))) |>
      dplyr::mutate(cc_clean = dplyr::if_else(
        stringr::str_detect(cc_clean, "Conserved Domain") &
          feature_name == "TP53",
        paste0(
          "aa_region:117-143,aa_region:171-181,",
          "aa_region:234-258,aa_region:270-286"),
        as.character(cc_clean)))
    }
  variant_df <- variant_df |>
    dplyr::mutate(cc_clean = dplyr::if_else(
      stringr::str_detect(cc_clean, "N/A"),
      stringr::str_replace_all(cc_clean, "N/A", ""),
      as.character(cc_clean))) |>
    dplyr::distinct()
    
  variant_df[,col] <- variant_df[,"cc_clean"]
  variant_df[,"cc_clean"] <- NULL
  return(variant_df)
}

get_literature_references <- function(biomarker_items){
  
  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  
  biomarker_items[['literature']] <- as.data.frame(
    biomarker_items[['clinical']] |>
      dplyr::filter(!is.na(citation_id)) |>
      tidyr::separate_rows(citation_id, sep=";") |> 
      dplyr::mutate(
        citation_id = stringr::str_trim(citation_id)) |> 
      dplyr::filter(
        !stringr::str_detect(
          tolower(citation_id),
          "asco|aacr|abstr|suppl|confex|jci|blood|annonc|caris")) |> 
      dplyr::distinct() |>
      dplyr::mutate(citation_id = dplyr::case_when(
        citation_id == "PMC3638050" ~ "PMID:22038996",
        citation_id == "PMC3936420" ~ "PMID:24265155",
        TRUE ~ as.character(
          stringr::str_replace(citation_id,"PMID: ","PMID:"))
      )) |>
      dplyr::filter(!stringr::str_detect(citation_id, "^[0-9]{1,5}$")) |>
      dplyr::filter(citation_id != "PMID:238900088" &
                      citation_id != "PMID:24569783" &
                      citation_id != "2556735") |>
      dplyr::mutate(citation_id = dplyr::case_when(
        stringr::str_detect(citation_id, "^FDA") ~ "FDA",
        stringr::str_detect(citation_id, "NCCN") ~ "NCCN",
        TRUE ~ as.character(citation_id)
      )) |>
      dplyr::mutate(citation_id = stringr::str_replace(
        citation_id, "PMID:",""
      )) |>
      
      dplyr::select(citation_id, evidence_id) |>
      dplyr::rename(source_id = citation_id) |>
      dplyr::mutate(source_id = as.character(source_id)) |>
      dplyr::group_by(source_id) |>
      dplyr::summarise(evidence_id = paste(
        unique(sort(evidence_id)), collapse = ";"
      ), .groups = "drop") |>
      dplyr::left_join(
        literatureVault::vault,
        by = "source_id", 
        relationship = "many-to-many",
        multiple = "all"
      )
  )
  
  biomarker_items[['clinical']]$citation_id <- NULL
  
  not_found <- biomarker_items$literature |>
    dplyr::filter(is.na(link)) |>
    dplyr::select(source, source_id)
  
  if(NROW(not_found) > 0){
    lgr::lgr$info(
      "Citations missing in sigven/literatureVault::vault - update needed")
    return(not_found)
  }else{
    lgr::lgr$info(
      "Found all literature references in sigven/literatureVault::vault"
    )
  }
  
  ## ignore molecular profiles/biomarkers that lack literature support
  biomarkers_with_literature_support <- as.data.frame(
    biomarker_items[['literature']] |>
      dplyr::select(evidence_id, source_id) |>
      tidyr::separate_rows(evidence_id, sep = ";") |>
      dplyr::group_by(evidence_id) |>
      dplyr::mutate(source_id = paste(
        unique(source_id), collapse = ";"
      )) |>
      dplyr::distinct()
  )
  
  biomarker_items[['clinical']] <- biomarker_items[['clinical']] |>
    dplyr::inner_join(
      biomarkers_with_literature_support,
      by = "evidence_id", multiple = "all",
      relationship = "many-to-many")
  
  return(biomarker_items)
  
}

#' #' A function that splits an array into chunks of equal size
#' chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
#' 
#' #' A function that returns a citation with first author, journal and year for a PubMed ID
#' #'
#' #' @param pmid An array of Pubmed IDs
#' #' @param chunk_size Size of PMID chunks
#' #'
#' #' @export
#' #'
#' #' @return citation PubMed citation, with first author, journal and year
#' #'
#' get_citations_pubmed <- function(
#'     pmid,
#'     pmid_cache_dir = NA,
#'     chunk_size = 2){
#'   
#'   ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
#'   pmid_chunks <- chunk(pmid, ceiling(length(pmid)/chunk_size))
#'   j <- 0
#'   all_citations <- data.frame()
#'   cat('Retrieving PubMed citations for PMID list, total length', length(pmid))
#'   cat('\n')
#'   while(j < length(pmid_chunks)){
#'     #cat(unlist(pmid_chunks),"\n")
#'     pmid_chunk <- pmid_chunks[[as.character(j)]]
#'     #cat(unlist(pmid_chunk),"\n")
#'     #cat('Processing chunk ',j,' with ',length(pmid_chunk),'PMIDS')
#'     #cat('\n')
#'     pmid_string <- paste(pmid_chunk,collapse = " ")
#'     res <- RISmed::EUtilsGet(
#'       RISmed::EUtilsSummary(pmid_string, type="esearch", db="pubmed", retmax = 5000)
#'     )
#'     
#'     
#'     year <- RISmed::YearPubmed(res)
#'     authorlist <- RISmed::Author(res)
#'     pmid_list <- RISmed::PMID(res)
#'     i <- 1
#'     first_author <- c()
#'     while(i <= length(authorlist)){
#'       
#'       
#'       if(length(authorlist[[i]]) == 5){
#'         first_author <- c(
#'           first_author,
#'           paste(authorlist[[i]][1,]$LastName," et al.",sep=""))
#'       }else{
#'         first_author <- c(
#'           first_author, as.character("Unknown et al.")
#'         )
#'       }
#'       i <- i + 1
#'     }
#'     journal <- RISmed::ISOAbbreviation(res)
#'     citations <- data.frame('pmid' = as.integer(pmid_list),
#'                             'citation' = paste(first_author,year,journal,sep=", "),
#'                             stringsAsFactors = F)
#'     citations$link <- paste0('<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',
#'                              citations$pmid,'\' target=\'_blank\'>',
#'                              citations$citation,'</a>')
#'     all_citations <- dplyr::bind_rows(
#'       all_citations, citations)
#'     j <- j + 1
#'   }
#'   
#'   return(all_citations)
#'   
#' }
#' 
#' 
#' get_literature_links <- function(literature_df) {
#'   
#'   stopifnot(is.data.frame(literature_df))
#'   assertable::assert_colnames(
#'     literature_df, c("source","source_id","evidence_id"),
#'     quiet = T)
#'   
#'   pubmed_source_ids <- literature_df |>
#'     dplyr::filter(source == "PubMed")
#'   
#'   other_source_ids <- literature_df |>
#'     dplyr::filter(source != "PubMed")
#'   
#'   
#'   if (NROW(pubmed_source_ids) > 0) {
#'     pubmed_source_ids <- get_citations_pubmed(
#'       pmid = pubmed_source_ids$source_id, chunk_size = 100) |>
#'       dplyr::rename(
#'         name = citation,
#'       ) |>
#'       dplyr::mutate(source_id = as.character(pmid)) |>
#'       dplyr::select(-pmid) |>
#'       dplyr::left_join(
#'         pubmed_source_ids, by = "source_id")
#'   }
#'     
#'   if (NROW(other_source_ids) > 0) {
#'     other_source_ids <- other_source_ids |>
#'       dplyr::filter(source != "PubMed") |>
#'       dplyr::mutate(name = source_id) |>
#'       dplyr::mutate(link = dplyr::case_when(
#'         source == "FDA" ~
#'           "<a href='https://www.fda.gov/drugs/resources-information-approved-drugs/oncology-cancer-hematologic-malignancies-approval-notifications' target='_blank'>FDA approvals</a>",
#'         source == "NCCN" ~
#'           "<a href='https://www.nccn.org/guidelines/category_1' target='_blank'>NCCN guidelines</a>",
#'         source == "EMA" ~ "<a href='https://www.ema.europa.eu/en/medicines/field_ema_web_categories%253Aname_field/Human' target='_blank'>European Medicines Agency (EMA)</a>",
#'         source == "clinicaltrials.gov" ~ paste0(
#'           "<a href='https://clinicaltrials.gov/ct2/show/", 
#'           source_id, "' target='_blank'>", 
#'           source_id, "</a>")
#'       ))
#'   }
#'   
#'   all_literature_links <- pubmed_source_ids |>
#'     dplyr::select(source, source_id, 
#'                   evidence_id, 
#'                   dplyr::everything())
#'   
#'   if(NROW(other_source_ids) > 0){
#'     all_literature_links <- 
#'       other_source_ids |>
#'       dplyr::bind_rows(pubmed_source_ids)
#'   }
#'   
#'   return(all_literature_links)
#'   
#'   
#' }


get_amino_acid_dictionary <- function() {
  aa_df <- 
    data.frame('one_letter' = names(Biostrings::AMINO_ACID_CODE), 
               'three_letter' = Biostrings::AMINO_ACID_CODE, 
               'three_letter_uc' = toupper(Biostrings::AMINO_ACID_CODE),
               stringsAsFactors = F) |>
    dplyr::filter(one_letter != "B" &
                    one_letter != "J" &
                    one_letter != "O" &
                    one_letter != "U" &
                    one_letter != "X" &
                    one_letter != "Z")
  
  rownames(aa_df) <- NULL
  return(aa_df)
  
  
}


expand_hgvs_terms <- function(var, aa_dict, add_codon_markers = FALSE) {
  hits <- c()
  
  if (stringr::str_detect(toupper(var),"([A-Z]{3}[0-9]{1,}(([A-Z]{3}) {0,1})([0-9]{1,})?)") &
     !stringr::str_detect(toupper(var),"(DEL|INS|DUP)|=|-|\\*")) {
    aa <- stringr::str_to_title(strsplit(var, "[0-9]{1,}")[[1]])
    codon <- stringr::str_match(var, "[0-9]{1,}")[[1]]

    
    if (length(aa) == 1) {
      alt1 <- paste0(aa[1],codon)
      alt2 <- paste0('p.', alt1)
      if (aa[1] != "Ter") {
        if (aa[1] %in% aa_dict$three_letter) {
          aa1 <- seqinr::a(aa[1])
          if (!is.na(aa1)) {
            alt3 <- paste0(seqinr::a(aa[1]), codon)
            hits <- c(alt1, alt2, alt3)
          }
        }
      }else{
        alt3 <- paste0("X", codon)
        hits <- c(alt1, alt2, alt3)

      }
    }else{
      if (length(aa) == 2) {
        alt1 <- stringr::str_replace(
          stringr::str_replace(paste0(aa[1],codon, aa[2]),"(FS|Fs)$","fs"),
          "ter","Ter")
        alt2 <- stringr::str_replace(
          stringr::str_replace(paste0('p.', alt1),"(FS|Fs)$","fs"),
          "ter","Ter")
        if (aa[1] != "Ter" & 
           aa[2] != "Ter" & 
           aa[2] != "Fs" &
           nchar(aa[1]) == 3 & 
           nchar(aa[2]) == 3 &
           aa[1] %in% aa_dict$three_letter &
           aa[2] %in% aa_dict$three_letter) {
          alt3 <- paste0(seqinr::a(aa[1]), codon, seqinr::a(aa[2]))
          alt4 <- paste0('p.', alt3)
          hits <- c(alt1, alt2, alt3, alt4)
        }
        else{
          if (aa[1] == "Ter") {
            alt3 <- paste0("X", codon, seqinr::a(aa[2]))
            alt4 <- paste0('p.', alt3)
            hits <- c(alt1, alt2, alt3, alt4)
          }else{
            if (aa[2] == "Ter") {
              if (aa[1] %in% aa_dict$three_letter) {
                alt3 <- paste0(seqinr::a(aa[1]), codon, "X")
                alt4 <- paste0('p.', alt3)
                hits <- c(alt1, alt2, alt3, alt4)
              }
            }
            if (aa[2] == "Fs") {
              if (aa[1] %in% aa_dict$three_letter) {
                alt3 <- paste0(seqinr::a(aa[1]), codon, "fs")
                alt4 <- paste0('p.', alt3)
                hits <- c(alt1, alt2, alt3, alt4)
              }
            }

          }
        }
      }
    }
  }
  if (stringr::str_detect(var,"^([A-Z]{1}[0-9]{1,}[A-Z]{0,1})$") &
     !stringr::str_detect(var, "FS|DUP|DEL|INS|fs|dup|del|ins|\\*")) {
    aa <- stringr::str_to_title(strsplit(var, "[0-9]{1,}")[[1]])
    codon <- stringr::str_match(var, "[0-9]{1,}")[[1]]

    if (length(aa) == 1) {
      alt1 <- paste0(aa[1],codon)
      alt2 <- paste0('p.', alt1)
      if (aa[1] != "X" & aa[1] %in% aa_dict$one_letter) {
        alt3 <- paste0(seqinr::aaa(aa[1]), codon)
        hits <- c(alt1, alt2, alt3)
      }else{
        alt3 <- paste0("Ter", codon)
        hits <- c(alt1, alt2, alt3)
      }
    }else{
      if (length(aa) == 2) {
        alt1 <- paste0(aa[1],codon, aa[2])
        alt2 <- paste0('p.', alt1)
        if (aa[1] != "X" & 
           aa[2] != "X" & 
           aa[1] %in% aa_dict$one_letter &
           aa[2] %in% aa_dict$one_letter) {
          alt3 <- paste0(seqinr::aaa(aa[1]),
                         codon, seqinr::aaa(aa[2]))
          alt4 <- paste0('p.', alt3)
          
          if(aa[2] == aa[1]){
            alt5 <-  paste0("p.", seqinr::aaa(aa[1]),
                            codon, "=")
            alt6 <- paste0("p.", aa[1],
                           codon, "=")
            hits <- c(alt1, alt2, alt3, alt4, alt5, alt6)
          }else{
            hits <- c(alt1, alt2, alt3, alt4)
            
          }
          
        }else{
          if (aa[1] == "X" & aa[2] %in% aa_dict$one_letter) {
            alt3 <- paste0("Ter",
                           codon, seqinr::aaa(aa[2]))
            alt4 <- paste0('p.', alt3)
            hits <- c(alt1, alt2, alt3, alt4)
          }else{
            if (aa[2] == "X" & aa[1] %in% aa_dict$one_letter) {
              alt3 <- paste0(seqinr::aaa(aa[1]),
                             codon, "Ter")
              alt4 <- paste0('p.', alt3)
              hits <- c(alt1, alt2, alt3, alt4)
            }

          }
        }

      }
    }
  }
  if (stringr::str_detect(var,"^C\\.")) {
    var <- stringr::str_replace(var,"^C\\.","c.")
    #hits <- res
  }
  if (stringr::str_detect(var,"^P\\.")) {
    var <- stringr::str_replace(var,"^P\\.","p.")
    #hits <- res
  }
  if (stringr::str_detect(var,"FS$")) {
    var <- stringr::str_replace(var,"FS$","fs")
  }
  if (stringr::str_detect(var,";") & stringr::str_detect(var,"\\(")) {
    var <- tolower(var)
  }
  if (stringr::str_detect(var,"^([A-Z]{1}[0-9]{1,}[A-Z]{1})$")) {
    res1 <- paste0('p.', var)
    if(add_codon_markers == TRUE){
      res2 <- stringr::str_match(var,"^([A-Z]{1}[0-9]{1,})")[[1]]
      hits <- c(hits,res1,res2)
    }else{
      hits <- c(hits, res1)
    }
  }
  if (stringr::str_detect(var,"FS$")) {
    res <- stringr::str_replace(var,"FS$","fs")
    hits <- res
  }
  if (stringr::str_detect(var,"(DEL|Del|INS|Ins|Dup|DUP)") &
     !stringr::str_detect(var, "INSERTION|DUPLICATION|NM_|NP_|DELETERIOUS|DELETION|MUTATION|Mutation|Deletion|")) {
    var <- stringr::str_replace(var,"(DEL|Del)","del")
    var <- stringr::str_replace(var,"(INS|Ins)","ins")
    var <- stringr::str_replace(var,"(DUP|Dup)","dup")
    if (!(stringr::str_detect(var,"^(c|p|rs|;|[0-9])"))) {
      res2 <- paste0('p.', var)
      hits <- c(var, res2)
    }
  }

  if (length(hits) == 0 & var != "") {
    if (!(stringr::str_detect(var,"^(c|rs|p|[0-9])")) &
       stringr::str_detect(var, "[0-9]{1,}") &
       !stringr::str_detect(
         var,"EXON|;|MUTATION|Mutation|-|NM_|NP_|DELETERIOUS|Deletion|INTRON|DELETION|DUPLICATION")) {

      if (stringr::str_detect(var,"FS$|DUP|DELINS|FS\\*")) {
        var <- stringr::str_replace(var, "FS","fs")
        var <- stringr::str_replace(var, "DUP","dup")
        var <- stringr::str_replace(var, "DELINS","delins")
        var <- stringr::str_replace(var, "DEL","del")
      }
      res1 <- paste0('p.',var)
      hits <- c(var, res1)
    }else{
      hits <- c(var)
    }
  }

  hits <- unique(hits)

  return(hits)
}

load_civic_biomarkers <- function(
    datestamp = '20240709',
    compound_synonyms = NULL,
    hg38_fasta = 
      "/Users/sigven/research/DB/hg38/hg38.fa",
    crossmap_chainfile_dir = 
      "/Users/sigven/research/DB/chainFiles",
    cache_dir = NA) {
  

  gene_alias <- geneOncoX::get_alias(
    cache_dir = cache_dir
  )
  civic_local_fnames <- list()
  for(e in c("VariantSummaries",
             "civic_accepted",
             "ClinicalEvidenceSummaries",
             "MolecularProfileSummaries")){
    
    
    civic_local_fnames[[e]] <- file.path(
      cache_dir,
      paste0("civic_", e, "_", datestamp,".tsv")
    )
    
    remote_url <- paste0(
      "https://civicdb.org/downloads/nightly/nightly-", e, ".tsv")
    if(e == "civic_accepted"){
      remote_url <- paste0(
        "https://civicdb.org/downloads/nightly/nightly-", e, ".vcf")
      civic_local_fnames[['vcf_grch37']] <- file.path(
        cache_dir,
        paste0("civic_", datestamp,"_grch37.vcf")
      )
      civic_local_fnames[['vcf_grch38']] <- file.path(
        cache_dir,
        paste0("civic_", datestamp,"_grch38.vcf")
      )
    }
    
    
    
    if(e == "civic_accepted"){
      if(!file.exists(civic_local_fnames[['vcf_grch37']])){
        download.file(
          url = remote_url,
          destfile = civic_local_fnames[['vcf_grch37']])
      }
    }else{
      if(!file.exists(civic_local_fnames[[e]])){
        download.file(
          url = remote_url,
          destfile = civic_local_fnames[[e]])
      }
    }

    if(e == 'civic_accepted'){
      vcfhelpR::crossmap_vcf(
        direction = "hg19Tohg38",
        chain_file_directory = crossmap_chainfile_dir,
        target_genome_file = hg38_fasta,
        source_vcf = civic_local_fnames[['vcf_grch37']],
        target_vcf = civic_local_fnames[['vcf_grch38']])
      
      system(paste0('bgzip -dc ',civic_local_fnames[['vcf_grch38']],
                    '.gz > ', 
                    civic_local_fnames[['vcf_grch38']]))
    }
  }
  
  aa_dict <- get_amino_acid_dictionary()
  
  variant_coordinates <- data.frame()
  for(assembly in c('grch37','grch38')){
    skip_lines <- 8
    if(assembly == "grch38"){
      skip_lines <- 7
    }
    
    variant_coordinates <- as.data.frame(
      dplyr::bind_rows(
        variant_coordinates,
        readr::read_tsv(
          file = civic_local_fnames[[paste0('vcf_',assembly)]], 
          skip = skip_lines, 
          show_col_types = F, 
          col_names = T, 
          guess_max = 1000) |>
          janitor::clean_names() |>
          dplyr::rename(variant_id = id) |>
          dplyr::mutate(
            vcf_coord = paste0(
              assembly, ":", 
              paste(number_chrom, pos, ref, alt, sep="_"))
          ) |>
          dplyr::select(variant_id, vcf_coord)
      )
    )
  }
  
  ## ignore these biomarkers
  biomarker_names_skip_regex <- 
    paste0(
      "methylat|phosphoryl|serum lev|mislocali|peri-therap|",
      "alternative|cytoplasmic|p-loop|reg_|alternative|polymorphism|",
      "loss of heterozygosity|philadelphia|internal|n-terminal")
  
  molecularProfileSummary <- as.data.frame(
    data.table::fread(
      civic_local_fnames[['MolecularProfileSummaries']],
      select = c(1:13), fill = T)) |>
    dplyr::rename(variant_id = variant_ids,
                  molecular_profile_name = name,
                  molecular_profile_summary = summary) |>
    dplyr::select(variant_id, 
                  molecular_profile_id, 
                  molecular_profile_name,
                  molecular_profile_summary) |>
    ## ignore profiles with aberration combinations for now
    tidyr::separate_rows(variant_id, sep = ",") |>
    ##dplyr::filter(!stringr::str_detect(variant_id, ",")) |>
    dplyr::mutate(variant_id = as.integer(variant_id)) |>
    clean_variant_naming(col = "molecular_profile_name") |>
    dplyr::filter(!stringr::str_detect(
      tolower(molecular_profile_name),
      biomarker_names_skip_regex))

  clinicalEvidenceSummary <- as.data.frame(
    readr::read_tsv(
      civic_local_fnames[['ClinicalEvidenceSummaries']],
      show_col_types = F, guess_max = 5000)) |>
    dplyr::filter(source_type == "PubMed") |>
    dplyr::mutate(citation_id = as.character(citation_id)) |>
    dplyr::select(molecular_profile_id,
                  evidence_id,
                  evidence_civic_url,
                  disease,
                  doid,
                  therapies,
                  citation_id,
                  evidence_statement,
                  evidence_type,
                  evidence_level,
                  evidence_direction,
                  significance,
                  variant_origin,
                  rating) |>
    dplyr::mutate(citation_id = dplyr::case_when(
      citation_id == "27577079" ~ "28947956",
      citation_id == "27556863" ~ "28978004",
      TRUE ~ as.character(citation_id)
    )) |>
    dplyr::mutate(disease = dplyr::if_else(
      is.na(disease),
      "Cancer",
      as.character(disease)
    )) |>
    dplyr::mutate(doid = dplyr::if_else(
      is.na(doid),
      "162",
      as.character(doid)
    )) |>
    dplyr::rename(disease_ontology_id = doid,
                  clinical_significance = significance,
                  evidence_description = evidence_statement,
                  evidence_url = evidence_civic_url,
                  cancer_type = disease) |>
    dplyr::mutate(
      clinical_significance =
        dplyr::case_when(
          clinical_significance == "Resistance" |
            clinical_significance == "Reduced Sensitivity" ~ "Resistance/Non-response",
          TRUE ~ as.character(clinical_significance))) |>
  dplyr::mutate(disease_ontology_id = paste0(
    "DOID:", disease_ontology_id)) |>
    dplyr::mutate(
      evidence_level =
        dplyr::case_when(
          evidence_level == 'A' ~ 'A: Validated',
          evidence_level == 'B' ~ 'B: Clinical evidence',
          evidence_level == 'C' ~ 'C: Case study',
          evidence_level == 'D' ~ 'D: Preclinical evidence',
          evidence_level == 'E' ~ 'E: Indirect evidence',
          TRUE ~ as.character(evidence_level)
        )
    ) |>
    dplyr::mutate(
      variant_origin =
        dplyr::if_else(
          variant_origin == "N/A" | 
            variant_origin == "Unknown",
          as.character(NA),
          as.character(
            stringr::str_replace(
              variant_origin,"Rare |Common ",""))))


  therapeutic_contexts <- as.data.frame(
    clinicalEvidenceSummary |>
      dplyr::select(evidence_id, therapies) |>
      dplyr::filter(!is.na(therapies) &
                      nchar(therapies) > 0) |>
      tidyr::separate_rows(therapies, sep = ",") |>
      dplyr::mutate(therapies = tolower(therapies)) |>
      dplyr::left_join(
        dplyr::select(compound_synonyms,
                      alias_lc,
                      drug_name,
                      molecule_chembl_id),
        by = c("therapies" = "alias_lc"), 
        multiple = "all", 
        relationship = "many-to-many"
      ) |>
      dplyr::filter(!is.na(drug_name)) |>
      dplyr::select(evidence_id, molecule_chembl_id) |>
      dplyr::distinct() |>
      dplyr::group_by(evidence_id) |>
      dplyr::summarise(molecule_chembl_id = paste(
        unique(sort(molecule_chembl_id)), collapse = "|"),
        .groups = "drop"
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(molecule_chembl_id = stringr::str_replace(
        molecule_chembl_id,"^NA\\||\\|NA$",""
      ))
  )

  clinicalEvidenceSummary <- clinicalEvidenceSummary |>
    dplyr::left_join(
      therapeutic_contexts,
      by = "evidence_id", multiple = "all",
      relationship = "many-to-many")
  
  variantSummary <- as.data.frame(
    data.table::fread(civic_local_fnames[['VariantSummaries']],
                      select = c(1:25), fill = T)) |>
    ## Ignore the "Factor" feature type (Kataegis, CK, methylation)
    dplyr::filter(feature_type != "Factor") |>
    #dplyr::rename(molecular_profile_id = variant_id) |>
    dplyr::mutate(alteration_type = "MUT") |>
    dplyr::mutate(variant = stringr::str_trim(variant)) |>
    dplyr::mutate(variant = dplyr::if_else(
      !is.na(variant) & 
        (variant == "Fusion" | '::' %in% variant),
      as.character(feature_name),
      as.character(variant)
    )) |>
    dplyr::mutate(variant_aliases = dplyr::if_else(
      gene == "BRAF" & variant == "V600",
      paste(variant_aliases, "V640", sep = ","),
      as.character(variant_aliases)
    )) |>
    dplyr::mutate(variant_aliases = dplyr::if_else(
      !is.na(hgvs_descriptions) & 
        hgvs_descriptions != "NA" & 
        nchar(hgvs_descriptions) > 0,
      paste(variant_aliases, hgvs_descriptions, sep=","),
      as.character(variant_aliases))) |>
    dplyr::select(-c("hgvs_descriptions")) |>
    dplyr::mutate(variant_aliases = stringr::str_replace(
      variant_aliases, "^,", ""
    )) |>
    tidyr::separate_rows(variant_aliases, sep = ",") |>
    tidyr::separate_rows(variant, sep=",") |>
    clean_variant_naming(col = "variant") |>
    clean_variant_naming(col = "variant_aliases") |>
    dplyr::rename(variant_consequence = variant_types,
                  molecular_profile_id = 
                    single_variant_molecular_profile_id,
                  variant_alias = variant_aliases) |>
    clean_variant_consequence() |>
    #dplyr::filter(variant != "Alu insertion") |>
    dplyr::filter(
        !stringr::str_detect(
            tolower(variant),
            biomarker_names_skip_regex)) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      feature_type == "Fusion" &
        nchar(variant_alias) == 0 | 
        stringr::str_detect(
          variant_alias,"multiple partners"),
      paste0(
        stringr::str_replace(
          feature_name,"v::|::v","")," Fusion"),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(variant = dplyr::if_else(
      nchar(variant_alias) > 0 &
        feature_type == "Fusion" &
        stringr::str_detect(variant,"^e[0-9]{1,}::"),
      stringr::str_replace(variant_alias," E[0-9]{1,}-E[0-9]{1,}$",""),
      as.character(variant))) |>
    dplyr::mutate(gene = dplyr::if_else(
      feature_type == "Fusion",
      as.character(feature_name),
      as.character(gene)
    )) |>
    set_alteration_type() |>
    dplyr::select(variant_id,
                  molecular_profile_id,
                  variant_consequence,
                  variant,
                  alteration_type,
                  gene,
                  variant_alias) |>
    dplyr::mutate(
      variant = stringr::str_replace_all(
        variant, 
        paste0(
          "((Inactivating|Activating|Deleterious|",
          "Promoter|Truncating|TKD) Mutation)|Loss-of-function"),
        "Mutation"
      )
    ) |>    
    dplyr::mutate(variant_alias = paste(
      variant_alias, variant, sep=","
    )) |>
    tidyr::separate_rows(
      variant_alias, sep=","
    ) |>
    dplyr::filter(
      !stringr::str_detect(
        tolower(variant_alias),
        biomarker_names_skip_regex)) |>
    dplyr::filter(
      (!(stringr::str_detect(variant,"H$") &
           stringr::str_detect(toupper(variant_alias),"TYR$"))) &
        (!(stringr::str_detect(variant,"Y$") &
             stringr::str_detect(toupper(variant_alias),"HIS$")))
      
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      variant_name_primary = dplyr::case_when(
        alteration_type == "CNA" |
          (alteration_type == "MUT_LOF" & !endsWith(variant,"X")) |
          stringr::str_detect(
            variant,
            "Mutation|Wildtype|Splice|Copy Number|Deletion|Amplification") |
          stringr::str_detect(alteration_type,"^EXP") ~
          paste0(gene," ",stringr::str_to_title(variant)),
        variant_consequence == "transcript_fusion" |
          (nchar(variant_consequence) == 0 & variant == "FUSION") ~
          paste0(gene, " Fusion"),
        variant == "REARRANGEMENT" ~ 
          paste0(gene, " Rearrangement"),
        TRUE ~ as.character(
          paste0(gene, " ", variant)
        ))
    ) |>
    dplyr::mutate(variant_consequence = dplyr::if_else(
      stringr::str_detect(variant,"^aa_region:"),
      as.character("missense_variant"),
      as.character(variant_consequence))) |>
    dplyr::filter(!nchar(variant_alias) == 0) |>
    dplyr::mutate(variant_name_primary = stringr::str_replace(
      variant_name_primary, "FS","fs")
    ) |>
    dplyr::mutate(variant_name_primary = stringr::str_replace(
      variant_name_primary, "DELINS","delins")
    ) |>
    dplyr::mutate(variant_name_primary = stringr::str_replace(
      variant_name_primary, "DEL","del")
    ) |>
    dplyr::mutate(variant_name_primary = stringr::str_replace(
      variant_name_primary, "INS$","ins")
    )
  
  variants_expanded <- list()
  variants_expanded[['mut']] <- data.frame()

  for (i in unique(variantSummary$variant_id)) {
    vrows <- variantSummary[variantSummary$variant_id == i,]

    variant_entries <- unique(c(vrows$variant,
                                vrows$variant_alias))
    gene <- unique(vrows$gene)
    variant_consequence <- unique(vrows$variant_consequence)
    alteration_type <- unique(vrows$alteration_type)
    variant_name_primary <- unique(vrows$variant_name_primary)
    
    all_aliases <- c()

    if (stringr::str_detect(alteration_type,"MUT")) {

      for (e in variant_entries) {
        #cat(e,'\n')
        clean_aliases <- expand_hgvs_terms(e, aa_dict = aa_dict)
        for (c in clean_aliases) {
         all_aliases <- unique(c(all_aliases, c))
        }
      }


      for (variant_alias in all_aliases) {
        df <- data.frame(
          'variant_id' = i,
          'symbol' = gene,
          'variant_name_primary' = variant_name_primary,
          'variant_alias' = variant_alias,
          'alteration_type' = alteration_type,
          'variant_consequence' = variant_consequence,
          stringsAsFactors = F)
        variants_expanded[['mut']] <- dplyr::bind_rows(
          variants_expanded[['mut']], df
        )
      }
    }

  }

  variants_expanded[['mut']] <- variants_expanded[['mut']] |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(
            variant_alias,
            "^(p\\.)?(([A-Z]{1}[0-9]{1,})|([A-Z]{1}[a-z]{2}[0-9]{1,}))$"),
          as.character("CODON"),
          as.character(alteration_type)
        )) |>
    dplyr::mutate(
      variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias,"FS"),
        stringr::str_replace(variant_alias,"FS","fs"),
        as.character(variant_alias)
      )
    ) |>
    dplyr::mutate(
      variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias,"DEL") &
          !stringr::str_detect(variant_alias, "DELETION|DELETERIOUS"),
        stringr::str_replace(variant_alias,"DEL","del"),
        as.character(variant_alias)
      )
    ) |>
    dplyr::mutate(
      variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias,"INS") &
          !stringr::str_detect(variant_alias,"INSERTION"),
        stringr::str_replace(variant_alias,"INS","ins"),
        as.character(variant_alias)
      )
    ) |>
    dplyr::mutate(
      variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias,"TER|ter") &
          !stringr::str_detect(
            tolower(variant_alias), 
            "alteration"),
        stringr::str_replace(variant_alias,"ter|TER","Ter"),
        as.character(variant_alias)
      )
    ) |>
    dplyr::mutate(
      variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias,"DUP") &
          !stringr::str_detect(variant_alias,"DUPLICATION"),
        stringr::str_replace(variant_alias,"DUP","dup"),
        as.character(variant_alias)
      )
    ) |>
    dplyr::filter(
      variant_alias != "p.=" &
        !stringr::str_detect(variant_alias, "^p\\.\\(") &
        !stringr::str_detect(variant_alias, "^p\\.Exon") &
        !stringr::str_detect(variant_alias, "ENSP|ENST") &
        !stringr::str_detect(variant_alias, "SPLICE"))
  
  
  aa_dict <- get_amino_acid_dictionary()
  
  for (i in 1:nrow(aa_dict)) {
    aa <- aa_dict[i,"three_letter_uc"]
    aa_three_lc <- aa_dict[i, "three_letter"]
    
    if (nrow(variants_expanded[['mut']][!is.na(variants_expanded[['mut']]$variant_alias) & 
                                 stringr::str_detect(variants_expanded[['mut']]$variant_alias,aa),]) > 0) {
      variants_expanded[['mut']][!is.na(variants_expanded[['mut']]$variant_alias) & 
                             stringr::str_detect(variants_expanded[['mut']]$variant_alias, aa),]$variant_alias <-
        stringr::str_replace_all(variants_expanded[['mut']][!is.na(variants_expanded[['mut']]$variant_alias) & 
                                                    stringr::str_detect(variants_expanded[['mut']]$variant_alias, aa),]$variant_alias,
                             aa,
                             aa_three_lc)
    }
  }
  
  variants_expanded[['mut']] <- variants_expanded[['mut']] |>
    dplyr::mutate(alias_type = "hgvsp") |>
    dplyr::mutate(variant_alias = stringr::str_trim(variant_alias)) |>
    dplyr::mutate(alias_type = dplyr::case_when(
      stringr::str_detect(variant_alias, "^rs[0-9]{1,}$") ~ "dbsnp",
      stringr::str_detect(variant_alias, "^aa_region") ~ "aa_region",
      stringr::str_detect(variant_alias, "^c\\.[0-9]{1,}") ~ "hgvsc",
      stringr::str_detect(variant_alias, "^Exon ") ~ "exon",
      stringr::str_detect(
        variant_alias, 
        "Mutation|Copy|Wildtype|Loss|Truncation|Alteration|Frameshift") ~ "other_gene",
      TRUE ~ as.character(alias_type)
    )) |>
    dplyr::left_join(variant_coordinates, 
                     by = "variant_id", 
                     relationship = "many-to-many") |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      !is.na(vcf_coord),
      paste(variant_alias, vcf_coord, sep=","),
      as.character(variant_alias)
    )) |>
    tidyr::separate_rows(
      variant_alias, sep = ","
    ) |>
    dplyr::select(-vcf_coord) |>
    dplyr::distinct() |>
    dplyr::mutate(alias_type = dplyr::if_else(
      stringr::str_detect(variant_alias, "^grch37"),
      "genomic_grch37",
      as.character(alias_type)
    )) |>
    dplyr::mutate(alias_type = dplyr::if_else(
      stringr::str_detect(variant_alias, "^grch38"),
      "genomic_grch38",
      as.character(alias_type)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic_grch37"),
      stringr::str_replace(variant_alias, "grch37:",""),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic_grch38"),
      stringr::str_replace(variant_alias, "grch38:",""),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(alteration_type = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic"),
      "MUT",
      as.character(alteration_type)
    )) |>
    dplyr::filter(
      symbol != "TERT" |
        (symbol == "TERT" &
           !stringr::str_detect(variant_alias, "Cys|^p\\.") &
           alteration_type != "CODON")) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      symbol == "TERT" & variant_alias == "-146C>T",
      paste(variant_alias, 
            "grch37:5_1295250_G_A,grch38:5_1295135_G_A", 
            sep=","),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      symbol == "MET" & 
        variant_alias == "Exon 14 Skipping Mutation",
      paste(
        variant_alias, 
        "c.3028+1G>T",
        "c.3028+1G>C",
        "c.3028+1G>A",
        "c.3028+2T>C",
        "c.3028+2T>G", 
        "c.3082+1G>T",
        "c.3082+1G>C",
        "c.3082+1G>A",
        "c.3082+2T>C",
        "c.3082+2T>G", 
        sep=","),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      symbol == "TERT" & 
        variant_alias == "-124C>T",
      paste(variant_alias, 
            "grch37:5_1295228_G_A,grch38:5_1295113_G_A", 
            sep=","),
      as.character(variant_alias)
    )) |>
    tidyr::separate_rows(
      variant_alias, sep = ","
    ) |>
    
    dplyr::mutate(alias_type = dplyr::if_else(
      stringr::str_detect(variant_alias, "^grch37"),
      "genomic_grch37",
      as.character(alias_type)
    )) |>
    dplyr::mutate(alias_type = dplyr::if_else(
      stringr::str_detect(variant_alias, "^grch38"),
      "genomic_grch38",
      as.character(alias_type)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic_grch37"),
      stringr::str_replace(variant_alias, "grch37:",""),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(variant_alias = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic_grch38"),
      stringr::str_replace(variant_alias, "grch38:",""),
      as.character(variant_alias)
    )) |>
    dplyr::mutate(alteration_type = dplyr::if_else(
      stringr::str_detect(alias_type, "^genomic"),
      "MUT",
      as.character(alteration_type)
    )) |>

    dplyr::mutate(alias_type = dplyr::if_else(
      stringr::str_detect(variant_alias, "^c\\."),
      "hgvsc",
      as.character(alias_type)
    )) |>
  
    dplyr::mutate(variant_exon = dplyr::if_else(
      stringr::str_detect(variant_alias, "Exon"),
      stringr::str_replace_all(
        variant_alias, 
        "Exon |( (Frameshift|((Skipping )?Mutation)|Deletion|Insertion))",""),
      as.character(NA)
    )) |>
    dplyr::mutate(variant_exon = dplyr::if_else(
      !is.na(variant_exon) & 
        variant_consequence == "exon_variant",
      stringr::str_replace_all(
        variant_exon,"18-21",
        "18|19|20|21"),
      as.character(variant_exon))
    ) |>
    dplyr::mutate(variant_exon = dplyr::if_else(
      !is.na(variant_exon) & 
        variant_consequence == "exon_variant",
      stringr::str_replace_all(
        variant_exon,
        "1-2","1|2"),
      as.character(variant_exon))
    ) |>
    dplyr::mutate(variant_exon = dplyr::if_else(
      !is.na(variant_exon) & 
        variant_consequence == "exon_variant",
      stringr::str_replace_all(
        variant_exon,
        "10 OR 21","10|21"),
      as.character(variant_exon))
    ) |>
    dplyr::mutate(variant_exon = dplyr::if_else(
      !is.na(variant_exon) & 
        variant_consequence == "exon_variant",
      stringr::str_replace_all(
        variant_exon,
        "10 AND 21","10&21"),
      as.character(variant_exon))
    ) |>
    dplyr::filter(
      is.na(variant_exon) 
      | (variant_exon != "9 AND 20" & 
           !stringr::str_detect(
             variant_exon,"Deletion"))) |>
    dplyr::distinct()
  
  
  ## Make sure alteration types for genomic coordinates
  ## are correct
  variants_nongenomic <- variants_expanded[['mut']] |>
    dplyr::filter(
      !stringr::str_detect(
        alias_type, "genomic_"
      )
    ) |>
    dplyr::group_by(
      variant_id
    ) |>
    dplyr::summarise(
      alteration_type =
        paste(unique(alteration_type), collapse=";")
    ) |>
    dplyr::mutate(
      alteration_type = dplyr::if_else(
        stringr::str_detect(alteration_type, ";"),
        stringr::str_replace(
          alteration_type, ";?CODON;?",""
        ),
        as.character(alteration_type)
      )
    )
  
  
  variants_expanded_mut1 <- variants_expanded[['mut']] |>
    dplyr::filter(!stringr::str_detect(
      alias_type, "genomic"))
  
  variants_expanded_mut2 <- variants_expanded[['mut']] |>
    dplyr::filter(stringr::str_detect(
      alias_type, "genomic"))
  variants_expanded_mut2$alteration_type <- NULL
  variants_expanded_mut2 <- variants_expanded_mut2 |>
    dplyr::left_join(variants_nongenomic, by = "variant_id",
                     relationship = "many-to-many")
  
  variants_expanded[['mut']] <- variants_expanded_mut1 |>
    dplyr::bind_rows(variants_expanded_mut2)
  
  
  ## Combine the two dataframes
  variants_expanded[['other']]  <- as.data.frame(
    variantSummary |>
      dplyr::filter(!stringr::str_detect(alteration_type,"MUT")) |>
      dplyr::select(variant_id, variant, variant_consequence,
                    variant_alias, variant_name_primary, alteration_type,
                    gene) |>
      dplyr::rename(symbol = gene) |>
      dplyr::mutate(
        alias_type = "other_gene"
      ) |>
      dplyr::distinct()
  )
  attr(variants_expanded[['other']],"groups") <- NULL
  
  gene_aliases <- gene_alias$records |>
    dplyr::filter(.data$n_primary_map == 1 & 
                    source == "NCBI") |>
    dplyr::select(alias, symbol, entrezgene) |>
    dplyr::distinct()
  
  biomarker_items <- list()
  biomarker_items[['variant']] <- 
    dplyr::bind_rows(
      variants_expanded[['mut']], 
      variants_expanded[['other']]) |>
    dplyr::select(variant_id, variant_alias, 
                  alias_type, dplyr::everything()) |>
    dplyr::rename(gene = symbol) |>
    dplyr::left_join(
      gene_aliases,
      by = c("gene" = "alias"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(variant_consequence = dplyr::case_when(
      nchar(variant_consequence) == 0 & stringr::str_detect(
        variant_alias, "delins([A-Z]{1,4})+$") ~ "inframe_insertion",
      nchar(variant_consequence) == 0 & stringr::str_detect(
        variant_alias, "ins([A-Z]{1,4})+$") ~ "inframe_insertion",
      nchar(variant_consequence) == 0 & stringr::str_detect(
        variant_alias, "[0-9]{1,}?del([A-Z]{1,4})+$") ~ "inframe_deletion",
      nchar(variant_consequence) == 0 & stringr::str_detect(
        variant_alias, "fs") ~ "frameshift_variant",
      nchar(variant_consequence) == 0 & stringr::str_detect(
        variant_alias, "^Mutation$") ~ "protein_altering_variant",
      nchar(variant_consequence) == 0 & 
        alias_type == "hgvsp" ~ "missense_variant",
      TRUE ~ as.character(variant_consequence)
    )) |> 
    dplyr::mutate(variant_alias = dplyr::case_when(
      (is.na(variant_alias) | nchar(variant_alias) == 0) &
        !is.na(gene)  &
        alias_type == "other_gene" &
        stringr::str_detect(gene,"-") ~ gene,
      (is.na(variant_alias) | nchar(variant_alias) == 0) &
        !is.na(gene)  &
        alias_type == "other_gene" &
        !stringr::str_detect(gene,"-") ~ 
        paste0(gene, " Fusion"),
      TRUE ~ as.character(variant_alias)
    )) |>
    dplyr::distinct()
  
  biomarker_items[['clinical']] <- 
    clinicalEvidenceSummary |> 
    dplyr::full_join(
      molecularProfileSummary, by = "molecular_profile_id",
      multiple = "all", relationship = "many-to-many") |> 
    dplyr::select(molecular_profile_id, molecular_profile_name, 
                  molecular_profile_summary, 
                  variant_id, dplyr::everything()) |>
    dplyr::rename(therapeutic_context = therapies) |>
    dplyr::mutate(evidence_id = paste0("EID",evidence_id)) |>
    dplyr::filter(evidence_id != "EIDNA") |>
    dplyr::mutate(biomarker_source = "civic",
                  biomarker_source_datestamp = datestamp,
                  biomarker_entity = T) |>
    dplyr::inner_join(dplyr::select(
      biomarker_items[['variant']], variant_id), 
      by = "variant_id",
      multiple = "all", relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::arrange(evidence_level, molecular_profile_id)
  
  biomarker_items[['clinical']] <- map_biomarker_phenotypes(
    biomarker_items[['clinical']], 
    cache_dir = cache_dir)
  
  biomarker_items <- get_literature_references(biomarker_items)
  if(is.data.frame(biomarker_items)){
    if(unique(colnames(biomarker_items) == c("source","source_id")) == T){
      lgr::lgr$error("literature missing ")
      return(biomarker_items)
    }
  }
  
  biomarker_items[['clinical']] <- biomarker_items[['clinical']] |>
    dplyr::mutate(molecular_profile_type = "Single") |>
    dplyr::mutate(molecular_profile_type = dplyr::case_when(
      stringr::str_detect(molecular_profile_name, " AND|and|And ") ~
      paste0(
        "Combo_", stringr::str_count(variant_id, " AND|and|And ")[1] + 1),
      stringr::str_detect(molecular_profile_name, " OR|or|Or ") ~
      paste0(
        "Combo_Opt_", stringr::str_count(variant_id, ",")[1] + 1),
      TRUE ~ molecular_profile_type
    )) |>
    dplyr::select(
      biomarker_source,
      biomarker_source_datestamp,
      biomarker_entity,
      molecular_profile_id,
      molecular_profile_name,
      molecular_profile_summary,
      molecular_profile_type,
      variant_id,
      source_id,
      evidence_id,
      variant_origin,
      primary_site,
      cancer_type,
      disease_ontology_id,
      do_name,
      efo_id,
      efo_name,
      cui,
      cui_name,
      dplyr::everything()
    )
  
  biomarker_items[['variant']] <- as.data.frame(
    biomarker_items[['variant']] |>
      dplyr::mutate(biomarker_source = "civic",
                    biomarker_source_datestamp = datestamp) |>
      dplyr::select(
        biomarker_source,
        biomarker_source_datestamp,
        variant_id,
        variant_alias,
        variant_name_primary,
        alteration_type,
        alias_type,
        variant_consequence,
        variant_exon,
        gene,
        symbol,
        entrezgene,
        dplyr::everything()
      ) |>
      dplyr::inner_join(
        dplyr::select(biomarker_items[['clinical']],
                      variant_id, molecular_profile_id),
        by = "variant_id",
        relationship = "many-to-many"
      ) |>
      dplyr::distinct() |>
      dplyr::select(-c("variant")) |>
      dplyr::select(
        biomarker_source,
        biomarker_source_datestamp,
        variant_id,
        molecular_profile_id,
        dplyr::everything()
      )
  )

  return(biomarker_items)

}


load_cgi_biomarkers <- function(compound_synonyms = NULL,
                                datestamp = '20221027',
                                cache_dir = NA) {

  aa_dict <- get_amino_acid_dictionary()
  gene_alias <- geneOncoX::get_alias(
    cache_dir = cache_dir
  )
  cancer_type_abbreviations <-
    read.table(
      file = file.path(cache_dir, "cgi_cancer_subtypes.tsv"),
      sep = "\t", header = T,
      stringsAsFactors = F, quote = "", 
      na.strings = "",fill = T) |>
    dplyr::mutate(
      disease_ontology_id = paste0(
        "DOID:",disease_ontology_id))

  ## SNVs/InDels
  ## SKIP/IGNORE mutations that are
  # 1) not explicitly mapped to protein alterations (MUT)
  # 2) marked with uncertainty, rechecking etc (Comments)
  # 3) involved in toxicity
  # 4) curated from abstracts etc (not PMIDs or guidelines)
  # 5) Non-malignant/Non-cancer conditions
  
  ## NEW
  cgi_biomarkers <- as.data.frame(
    readr::read_tsv(
      file.path(cache_dir, "cgi_biomarkers_20221017.tsv"),
      col_names = T, show_col_types = F, guess_max = 1000) |>
      janitor::clean_names() |>
      dplyr::mutate(
        biomarker = 
          stringr::str_replace_all(
            biomarker,
            "\\.",""
          )) |>
      dplyr::select(-c(curator, curation_date, targeting,
                       drug_family, drug_status, drug, tcgi_included,
                       metastatic_tumor_type, 
                       primary_tumor_type_full_name, assay_type)) |>
      dplyr::distinct() |>
      dplyr::mutate(biomarker = stringr::str_replace_all(
        biomarker, "12,13,59,61,117,146",
        "G12,G13,A59,Q61,K117,A146")) |>
      dplyr::mutate(alteration_type = stringr::str_replace(
        alteration_type, "EXPR", "EXP"
      )) |>
      dplyr::mutate(alteration_type = stringr::str_replace(
        alteration_type, "FUS", "FUSION"
      )) |>
      dplyr::mutate(biomarker = stringr::str_replace_all(
        biomarker, "\\*","X")) |>
      dplyr::mutate(biomarker = stringr::str_replace_all(
        biomarker, "618,620,634,768,791,891,918",
        "C618,C620,C634,E768,Y791,S891,M918")) |>
      dplyr::mutate(biomarker = stringr::str_replace_all(
        biomarker, "12,13","G12,G13")) |>
      dplyr::mutate(biomarker = stringr::str_replace(
        biomarker, "undexpression","Underexpression")) |>
      dplyr::mutate(biomarker = stringr::str_replace_all(
        biomarker, ",537,538,",",Y537,D538,")) |>
      dplyr::rename(
        molecular_profile = biomarker
      ) |>
      dplyr::mutate(molecular_profile = stringr::str_replace_all(
        molecular_profile, ";|(, )", ","
      )) |>
      dplyr::mutate(molecular_profile = stringr::str_replace_all(
        molecular_profile, " (\\+|-)$", ""
      )) |>
      dplyr::filter(
        !stringr::str_detect(source, "^(ENA|ESMO|ASCO|AACR)")
      ) |>
      tidyr::separate_rows(
        primary_tumor_type, sep = ";") |>
      dplyr::left_join(
        cancer_type_abbreviations,
        by = c("primary_tumor_type" = "cgi_abbreviation"),
        relationship = "many-to-many") |>
      dplyr::mutate(molecular_profile_id = dplyr::row_number()) |>
      dplyr::mutate(molecular_profile_type = dplyr::if_else(
        stringr::str_detect(molecular_profile, " \\+ "),
        "Combo",
        "Single"
      )) |>
      dplyr::mutate(molecular_profile_name = molecular_profile) |>
      tidyr::separate_longer_delim(molecular_profile, delim = " + ")
  )
      
  combo_index <- cgi_biomarkers |>
    dplyr::group_by(molecular_profile_id) |>
    dplyr::reframe(combo_index = dplyr::row_number()) |>
    dplyr::select(-molecular_profile_id)
  
  cgi_biomarkers <- as.data.frame(
    cgi_biomarkers |>
      dplyr::bind_cols(combo_index) |>                               
      dplyr::mutate(molecular_profile_type = dplyr::if_else(
        molecular_profile_type == "Combo",
        paste(
          molecular_profile_type, combo_index, sep="_"),
        as.character(molecular_profile_type)
      )) |>
      dplyr::select(-combo_index) |>
      dplyr::mutate(alteration2 = 
                      stringr::str_replace_all(
                        stringr::str_match(
                          molecular_profile, "\\([A-Z0-9]+(,[A-Z0-9]+){0,}\\)")[,1],
                        "\\)|\\(","")
      ) |>
      dplyr::mutate(alteration_type = dplyr::case_when(
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(molecular_profile, " amplification") ~ "CNA",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(molecular_profile, " deletion") &
          !stringr::str_detect(molecular_profile, "inframe") ~ "CNA",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(molecular_profile, " fusion") ~ "FUSION",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(tolower(molecular_profile), "overexpression") ~ "EXP_OVER",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(tolower(molecular_profile), "underexpression") ~ "EXP_OVER",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(molecular_profile, " expression") ~ "EXP",
        stringr::str_detect(molecular_profile_type, "Combo") &
          stringr::str_detect(alteration_type,"MUT") & 
          stringr::str_detect(tolower(molecular_profile), "oncogenic|\\(") ~ "MUT_ONC",
        TRUE ~ as.character(alteration_type)
      )) |>
      tidyr::separate(gene, c("gene1","gene2"), sep = ";") |>
      dplyr::mutate(gene = dplyr::if_else(
        !is.na(gene1) & !is.na(gene2) & 
          startsWith(molecular_profile, gene2),
        as.character(gene2),
        as.character(gene1)
      )) |>
      dplyr::select(-c(gene1,gene2)) |>
      dplyr::select(molecular_profile, gene, dplyr::everything()) |>
      dplyr::mutate(alteration_type = dplyr::case_when(
        alteration_type == "BIA" ~ "MUT_LOF_BIALLELIC",
        stringr::str_detect(
          molecular_profile, " wildtype") ~ "WT",
        stringr::str_detect(
          molecular_profile, " oncogenic") ~ "MUT_ONC",
        stringr::str_detect(
          tolower(molecular_profile), "overexpression")  ~ "EXP_OVER",
        stringr::str_detect(
          tolower(molecular_profile), "underexpression") ~ "EXP_UNDER",
        stringr::str_detect(
          molecular_profile, " fusion$") ~ "FUSION",
        TRUE ~ as.character(alteration_type)
      )) |>
      
      dplyr::mutate(variant_consequence = dplyr::case_when(
        
        alteration_type == "WT" ~ as.character(NA),
        alteration_type == "CNA" & 
          stringr::str_detect(
            alteration,":amp") ~ "transcript_amplification",
        alteration_type == "CNA" & 
          stringr::str_detect(
            alteration,":del") ~ "transcript_ablation",
        alteration_type == "EXP_OVER" ~ "transcript_amplification",
        alteration_type == "FUSION" ~ "transcript_fusion",
        alteration_type == "EXP_UNDER" ~ "transcript_ablation",
        alteration_type == "MUT_LOF_BIALLELIC" ~ "loss_of_function_variant",
        alteration_type == "MUT_ONC" ~ "protein_altering_variant",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          !stringr::str_detect(molecular_profile,"inframe|splice|insertion|deletion") ~ "missense_variant",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          stringr::str_detect(molecular_profile,"splice donor") ~ "splice_donor_variant",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          stringr::str_detect(molecular_profile,"splice acceptor") ~ "splice_acceptor_variant",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          stringr::str_detect(molecular_profile,"inframe") &
          stringr::str_detect(molecular_profile,"insertion") &
          stringr::str_detect(molecular_profile, "deletion") ~ "inframe_deletion,inframe_insertion",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          stringr::str_detect(molecular_profile,"inframe") &
          stringr::str_detect(molecular_profile, "deletion") ~ "inframe_deletion",
        molecular_profile == "FLT3-ITD" ~ "inframe_insertion",
        molecular_profile == "PDGFRA exon 18 mutations" ~ "inframe_deletion,inframe_insertion,missense_variant",
        (alteration_type == "MUT" | alteration_type == "MUT_ONC") &
          stringr::str_detect(molecular_profile,"inframe") &
          stringr::str_detect(molecular_profile, "insertion") ~ "inframe_insertion",
        (alteration_type == "MUT") &
          stringr::str_detect(molecular_profile, " exon ") &
          stringr::str_detect(alteration, "inframe_insertion") ~ "inframe_insertion",
        (alteration_type == "MUT") &
          stringr::str_detect(molecular_profile, " exon ") &
          stringr::str_detect(alteration, "inframe_deletion") ~ "inframe_deletion",
        TRUE ~ as.character(NA)
      )) |>
      dplyr::rename(
        variant_alias = alteration2
      ) |>
      dplyr::mutate(variant_alias = dplyr::case_when(
        alteration_type == "EXP" ~ "Expression",
        stringr::str_detect(molecular_profile, "wildtype") ~ "Wildtype",
        is.na(variant_alias) & 
          (alteration_type == "MUT" | 
             alteration_type == "MUT_ONC" |
             alteration_type == "MUT_LOF_BIALLELIC") ~ "Mutation",
        alteration_type == "CNA" & 
          variant_consequence == "transcript_amplification" ~ "Amplification",
        alteration_type == "CNA" & 
          variant_consequence == "transcript_ablation" ~ "Deletion",
        alteration_type == "FUSION" ~ "Fusion",
        alteration_type == "EXP_UNDER" & 
          variant_consequence == "transcript_ablation" ~ "Underexpression",
        alteration_type == "EXP_OVER" & 
          variant_consequence == "transcript_amplification" ~ "Overexpression",
        TRUE ~ as.character(variant_alias)
      )) |>
      dplyr::mutate(gene = dplyr::if_else(
        alteration_type == "FUSION" &
          stringr::str_detect(molecular_profile, "-"),
        stringr::str_replace(
          stringr::str_replace(molecular_profile, " fusion",""),
          "ABL1-BCR","BCR-ABL1"), 
        as.character(gene))
      ) |>
      dplyr::mutate(gene = dplyr::case_when(
        gene == "MLL2" ~ "KMT2D",
        TRUE ~ as.character(gene)
      )) |>
      dplyr::mutate(variant_origin = "Somatic") |>
      dplyr::mutate(variant_origin = dplyr::if_else(
        !is.na(comments) & 
          stringr::str_detect(tolower(comments),"germline"),
        "Germline",
        as.character(variant_origin)
      )) |>
      dplyr::mutate(evidence_type = 'Predictive',
                    therapeutic_context = drug_full_name,
                    clinical_significance = NA) |>
      dplyr::mutate(
        clinical_significance =
          dplyr::case_when(
            association == "Resistant" ~ "Resistance/Non-response",
            association == "No Responsive" ~ "Resistance/Non-response",
            association == "Responsive" ~ "Sensitivity/Response",
            stringr::str_detect(association, "Toxicity") ~ "Toxicity",
            TRUE ~ as.character("Other"))) |>
      dplyr::rename(evidence_level_raw = evidence_level) |>
      dplyr::mutate(mapping_rank = 1) |>
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(
            stringr::str_detect(evidence_level_raw,"guidelines"),
            "A: FDA/NCCN/ELN guidelines",
            as.character(NA))) |>
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(
            stringr::str_detect(evidence_level_raw,"Late trials"),
            "B1: Clinical evidence: late trials",
            as.character(evidence_level))) |>
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(
            stringr::str_detect(evidence_level_raw,"Early trials"),
            "B2: Clinical evidence: early trials",
            as.character(evidence_level))) |>
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(
            stringr::str_detect(evidence_level_raw,"Case report"),
            "C: Case study",
            as.character(evidence_level))) |>
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(
            stringr::str_detect(evidence_level_raw,"Pre-clinical"),
            "D: Preclinical evidence",
            as.character(evidence_level))) |>
      dplyr::filter(!is.na(evidence_level_raw)) |>
      dplyr::select(-comments) |>
      dplyr::mutate(alteration = dplyr::if_else(
        molecular_profile == "KIT mutation in exon 9,11,13,14 or 17",
        "KIT:449-514,550-592,627-664,664-714,788-828",
        as.character(alteration)
      )) |>
      dplyr::mutate(alteration = dplyr::if_else(
        stringr::str_detect(molecular_profile, " exon |-ITD$"),
        stringr::str_replace_all(
          alteration, 
          "KIT:|([A-Z0-9]{3,}::consequence::(skipping_mutation|(inframe_(insertion|variant|deletion))):)",
          ""
        ),
        as.character(alteration)
      )) |>
      dplyr::mutate(alteration = dplyr::if_else(
        stringr::str_detect(molecular_profile, " exon |-ITD$"),
        stringr::str_replace_all(
          alteration, "::(inframe_insertion|missense_variant):", ""),
        as.character(alteration)
      )) |>
      dplyr::mutate(variant_alias = dplyr::if_else(
        stringr::str_detect(molecular_profile, " exon |-ITD$") &
          !is.na(alteration),
        as.character(alteration),
        as.character(variant_alias)
      )) |>
      tidyr::separate_longer_delim(variant_alias, delim = ",") |>
      dplyr::distinct() |>
      dplyr::mutate(variant_alias = dplyr::if_else(
        stringr::str_detect(variant_alias, "-"),
        paste0("aa_region:", variant_alias),
        as.character(variant_alias)
      )) |>
      dplyr::mutate(variant_alias = dplyr::if_else(
        stringr::str_detect(
          molecular_profile, "^(IL7R inframe insertion \\(237)"),
        "aa_region:237-255",
        as.character(variant_alias)
      )) |>
      dplyr::left_join(
        dplyr::select(gene_alias$records, 
                      alias, symbol, entrezgene),
        by = c("gene" = "alias"), multiple = "all",
        relationship = "many-to-many"
      ) |>
      dplyr::mutate(alteration = stringr::str_replace_all(
        alteration,"__\\.","::v")
      ) |>
      dplyr::mutate(variant_consequence = dplyr::if_else(
        stringr::str_detect(variant_alias, "X$"),
        "stop_gained",
        as.character(variant_consequence)
      )) |>
      dplyr::mutate(alteration_type = dplyr::if_else(
        stringr::str_detect(
          variant_consequence, 
          "stop_gained|splice_donor|splice_acceptor"),
        "MUT_LOF",
        as.character(alteration_type)
      ))
  )
      
  unique_variants <- as.data.frame(
    cgi_biomarkers  |>
      dplyr::filter(!is.na(variant_alias) &
                      !is.na(gene) &
                      !is.na(alteration_type)) |>
      dplyr::select(gene, symbol, entrezgene, variant_alias,
                    alteration_type, variant_consequence) |>
      dplyr:: arrange(symbol, alteration_type, variant_alias) |>
      dplyr::distinct() |>
      dplyr::mutate(variant_id = dplyr::row_number()) |>
      dplyr::mutate(variant_name_primary = dplyr::case_when(
        !is.na(symbol) ~ paste(symbol, variant_alias),
        is.na(symbol) &
          !is.na(gene) ~ paste(gene, variant_alias),
        TRUE ~ as.character(NA)
      ))
  )

  cgi_variants <- data.frame()
  for (i in 1:nrow(unique_variants)) {

    alt <- unique_variants[i,"variant_alias"]
    terms <- expand_hgvs_terms(alt, aa_dict = aa_dict)
    for (t in terms) {
      df <- data.frame(
        'variant_id' = 
          unique_variants[i,"variant_id"],
        'variant_alias' = t,
        'variant_name_primary' = unique_variants[i,"variant_name_primary"],
        'gene' = unique_variants[i,"gene"],
        'entrezgene' = unique_variants[i,"entrezgene"],
        'symbol' = unique_variants[i,"symbol"],
        'variant_consequence' = unique_variants[i,"variant_consequence"],
        'alteration_type' = unique_variants[i,"alteration_type"],
        stringsAsFactors = F)
      cgi_variants <- cgi_variants |>
        dplyr::bind_rows(df)
    }
  }
  
  cgi_variants <- cgi_variants |>
    dplyr::mutate(alias_type = "hgvsp") |>
    dplyr::mutate(alias_type = dplyr::case_when(
      stringr::str_detect(variant_alias, "aa") ~ "aa_region",
      stringr::str_detect(variant_alias, "^(Mut|Wild|Exp|Overexp|Underexp|Fus|Ampli|Delet)") ~ "other_gene",
      TRUE ~ as.character(alias_type)
    )) |>
    dplyr::mutate(variant_exon = NA) |>
    dplyr::mutate(biomarker_source = "cgi") |>
    dplyr::mutate(biomarker_source_datestamp = datestamp)  |>
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          !is.na(alteration_type) & !is.na(variant_alias) & 
            stringr::str_detect(
              variant_alias, 
              "^(([A-Z]{1}[0-9]{1,})|([A-Z]{1}[a-z]{2}[0-9]{1,}))$"),
          as.character("CODON"),
          as.character(alteration_type)
        )) |>
    dplyr::select(
      biomarker_source,
      biomarker_source_datestamp,
      variant_id,
      variant_alias,
      variant_name_primary,
      alteration_type,
      alias_type,
      variant_consequence,
      variant_exon,
      gene,
      symbol,
      entrezgene,
    ) |>
    dplyr::distinct()
    

  cgi_clinical <- cgi_biomarkers |>

    ## get variant id's per molecular profile
    dplyr::inner_join(
      dplyr::select(
        cgi_variants, variant_id, 
        gene, alias_type, variant_alias),
      by = c("gene", "variant_alias"),
      relationship = "many-to-many"
    ) |>
    
    dplyr::select(-c(entrezgene, variant_consequence, variant_alias,
                     alias_type, gene, symbol, alteration_type,
                     alteration)) |>
    dplyr::distinct() |>
    dplyr::mutate(therapeutic_context = dplyr::if_else(
      stringr::str_detect(
        therapeutic_context,
        "(ib|ide|ole|ant|stat|ine|us|one|mab|in|el|[0-9])( |  )\\("),
      stringr::str_replace(
        therapeutic_context,
        "( |  )\\(.*\\)( )?$",""),
      as.character(therapeutic_context)
    ))
  

  therapeutic_contexts <- as.data.frame(
    cgi_clinical |>
      dplyr::select(molecular_profile_id, therapeutic_context) |>
      dplyr::filter(!is.na(therapeutic_context) &
                      nchar(therapeutic_context) > 0) |>
      tidyr::separate_rows(therapeutic_context, sep = " \\+ ") |>
      dplyr::mutate(therapeutic_context = tolower(therapeutic_context)) |>
      dplyr::left_join(
        dplyr::select(compound_synonyms,
                      alias_lc, 
                      drug_name, 
                      molecule_chembl_id),
        by = c("therapeutic_context" = "alias_lc"), 
        multiple = "all", relationship = "many-to-many"
      ) |>
      dplyr::filter(!is.na(drug_name)) |>
      dplyr::select(molecular_profile_id, molecule_chembl_id) |>
      dplyr::distinct() |>
      dplyr::group_by(molecular_profile_id) |>
      dplyr::summarise(molecule_chembl_id = paste(
        unique(sort(molecule_chembl_id)), collapse="|"),
        .groups = "drop"
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(molecule_chembl_id = stringr::str_replace(
        molecule_chembl_id,"^NA\\||\\|NA$",""
      ))
  )

  cgi_clinical <- cgi_clinical |>
    dplyr::left_join(
      therapeutic_contexts, by = "molecular_profile_id", 
      multiple = "all", relationship = "many-to-many") |>
    dplyr::mutate(evidence_url = "https://www.cancergenomeinterpreter.org/biomarkers") |>
    dplyr::mutate(biomarker_source = "cgi") |>
    dplyr::mutate(biomarker_source_datestamp = datestamp) |>
    dplyr::mutate(evidence_description = NA,
                  evidence_id = paste0("CGI", molecular_profile_id),
                  molecular_profile_summary = NA,
                  rating = NA,
                  biomarker_entity = T) |>
    dplyr::select(-c(association, evidence_level_raw,
                     primary_tumor_type, drug_full_name)) |>
    dplyr::rename(cancer_type = disease_name,
                  citation_id = source)

  cgi_clinical <- map_biomarker_phenotypes(
    biomarkers_clinical = cgi_clinical,
    cache_dir = cache_dir) |>
    dplyr::distinct()
  
  biomarker_items <- list()
  
 
  biomarker_items[['clinical']] <- cgi_clinical
  biomarker_items[['clinical']]$molecular_profile <- NULL
  
  
  biomarker_items <- get_literature_references(biomarker_items)
  biomarker_items[['clinical']] <- biomarker_items[['clinical']] |>
    dplyr::select(
      biomarker_source,
      biomarker_source_datestamp,
      biomarker_entity,
      molecular_profile_id,
      molecular_profile_name,
      molecular_profile_summary,
      molecular_profile_type,
      variant_id,
      source_id,
      evidence_id,
      variant_origin,
      primary_site,
      cancer_type,
      disease_ontology_id,
      do_name,
      efo_id,
      efo_name,
      cui,
      cui_name,
      dplyr::everything()
    ) |>
    dplyr::distinct()
  
  biomarker_items[['variant']] <- as.data.frame(
    cgi_variants |>
    dplyr::inner_join(
      dplyr::select(
        biomarker_items$clinical,
        variant_id
      ), by = "variant_id", 
      relationship = "many-to-many"
    ) |>
     dplyr::distinct() 
  )

  return(biomarker_items)
}

map_biomarker_phenotypes <- function(biomarkers_clinical = NULL,
                                     cache_dir = NA) {

  cancer_term_map <- 
    phenOncoX::get_terms(cache_dir = cache_dir)$records |>
    dplyr::select(do_id, efo_id, efo_name, do_name,
                  cui, cui_name, primary_site) |>
    dplyr::filter(!is.na(do_id) & !is.na(cui)) |>
    dplyr::rename(disease_ontology_id = do_id) |>
    dplyr::distinct()
  
  cancer_aux_pheno_maps <- 
    phenOncoX::get_aux_maps(cache_dir = cache_dir)

  umls_terms_all <- 
    cancer_aux_pheno_maps$records$umls$concept |>
    dplyr::filter(main_term == T) |>
    dplyr::select(cui, cui_name)
  
  do_terms_all <- 
    cancer_aux_pheno_maps$records$do |>
    dplyr::select(
      cui, do_id, do_name
    ) |>
    dplyr::distinct()
  
  efo_terms_all <- 
    cancer_aux_pheno_maps$records$efo$efo2xref |>
    dplyr::select(
      cui, efo_id, efo_name, primary_site) |>
    dplyr::distinct()

  biomarkers_all_phenotypes <- biomarkers_clinical |>
    dplyr::left_join(
      cancer_term_map, by = "disease_ontology_id", 
      multiple = "all", relationship = "many-to-many")

  ## if conditions in biomarkers is not found in phenOncoX terms
  ## match against individual dictionaries
  missing_do_phenotypes <- biomarkers_all_phenotypes |>
    dplyr::filter(is.na(cui)) |>
    dplyr::select(-c(cui, cui_name, primary_site,
                     efo_id, do_name, efo_name)) |>
    dplyr::inner_join(
      do_terms_all,
      by = c("disease_ontology_id" = "do_id"), 
      multiple = "all", relationship = "many-to-many") |>
    dplyr::left_join(
      umls_terms_all, by = "cui", 
      multiple = "all", relationship = "many-to-many") |>
    dplyr::left_join(
      efo_terms_all,
      by = "cui", multiple = "all",  
    relationship = "many-to-many") |>
    dplyr::distinct() |>
    dplyr::mutate(primary_site = dplyr::case_when(
      do_name == "cancer" ~ as.character(NA),
      is.na(primary_site) & stringr::str_detect(do_name,"brain") ~ "CNS/Brain",
      is.na(primary_site) & stringr::str_detect(do_name,"breast") ~ "Breast",
      is.na(primary_site) & stringr::str_detect(do_name,"colon|rectum") ~ "Colon/Rectum",
      is.na(primary_site) & stringr::str_detect(do_name,"gastric|stomach|esophag") ~ "Esophagus/Stomach",
      is.na(primary_site) & stringr::str_detect(do_name,"prostate") ~ "Prostate",
      is.na(primary_site) & stringr::str_detect(do_name,"pancrea") ~ "Pancreas",
      is.na(primary_site) & stringr::str_detect(do_name,"lung") ~ "Lung",
      is.na(primary_site) & stringr::str_detect(do_name,"myeloid") ~ "Myeloid",
      is.na(primary_site) & stringr::str_detect(do_name,"lymphoma") ~ "Lymphoid",
      is.na(primary_site) & stringr::str_detect(do_name,"ovary|ovarian") ~ "Ovary",
      TRUE ~ as.character(primary_site)
    ))

  biomarkers_phenotype_mapped <-
    dplyr::filter(biomarkers_all_phenotypes, !is.na(cui)) |>
    dplyr::bind_rows(missing_do_phenotypes)

  return(biomarkers_phenotype_mapped)
}

load_custom_fusion_db <- function() {

  other_translocations <- as.data.frame(read.table(
    file = "data-raw/biomarkers/translocation_patterns.tsv", header = T,
    sep = "\t", stringsAsFactors = F) |>
      dplyr::select(ngram_pattern, biomarker) |>
      dplyr::rename(alias = ngram_pattern,
                    variant = biomarker) |>
      dplyr::mutate(biomarker_entity = T,
                    symbol = variant,
                    alias = toupper(alias)) |>
      tidyr::separate_rows(alias, sep = "\\|") |>
      dplyr::mutate(biomarker_source_db = "fusion_customdb") |>
      dplyr::distinct()
  )

  other_translocations_lc <- as.data.frame(read.table(
    file = "data-raw/biomarkers/translocation_patterns.tsv", header = T,
    sep = "\t", stringsAsFactors = F) |>
      dplyr::select(ngram_pattern, biomarker) |>
      dplyr::rename(alias = ngram_pattern,
                    variant = biomarker) |>
      dplyr::mutate(biomarker_entity = T,
                    symbol = variant) |>
      tidyr::separate_rows(alias, sep = "\\|") |>
      dplyr::mutate(biomarker_source_db = "fusion_customdb") |>
      dplyr::distinct()
  )

  biomarker_items <- list()
  
  biomarker_items[['variant']] <- other_translocations |>
    dplyr::bind_rows(other_translocations_lc) |>
    dplyr::distinct() |>
    dplyr::rename(variant_alias = alias,
                  gene = symbol) |>
    dplyr::mutate(alteration_type = "FUSION",
                  variant_consequence = "transcript_fusion",
                  variant_id = variant,
                  variant_name_primary = variant) |>
    dplyr::mutate(alias_type = "other_gene") |>
    dplyr::select(variant_id, variant_name_primary,
                  variant_alias, alias_type,
                  gene, alteration_type, variant_consequence)

  return(biomarker_items)
}

load_depmap_fusions <- function(db_datestamp = "24Q4"){
  
  # Load DepMap fusions
  depmap_data <- list()
  depmap_data[['fusions']] <- as.data.frame(read.csv(
    file = "data-raw/depmap/OmicsFusionFiltered.csv.gz", header = T))
  
  depmap_data[['models']] <- as.data.frame(read.csv(
    file = "data-raw/depmap/Model.csv.gz", header = T)) |>
    dplyr::select(
      ModelID, CellLineName, OncotreeLineage,
      OncotreePrimaryDisease, OncotreeCode,
      Age, Sex, PrimaryOrMetastasis, SampleCollectionSite,
      SourceType
    )
  
  return(depmap_data)
}

load_mitelman_db <- function(cache_dir = NA,
                             db_datestamp = "20241015") {

  # Load Mitelman database
  # dos2unix -q -n MBCA.TXT.DATA MBCA.TXT
  
  morphology <- read.table(
    file = file.path(
      cache_dir, "mitelmandb", "KODER.TXT.DATA"),
    sep = "\t", stringsAsFactors = F, 
    header = T, fill = T) |>
    dplyr::filter(KodTyp == "MORPH" & nchar(Kod) > 0) |>
    dplyr::mutate('morphology_code' = as.character(Kod),
                  'morphology_name' = Benamning) |>
    dplyr::select(morphology_code, morphology_name) |>
    dplyr::distinct()
  
  topography <- read.table(
    file = file.path(
      cache_dir, "mitelmandb", "KODER.TXT.DATA"),
    sep = "\t", stringsAsFactors = F, 
    header = T, fill = T) |>
    dplyr::filter(KodTyp == "TOP" & nchar(Kod) > 0) |>
    dplyr::mutate('topography_code' = as.character(Kod),
                  'topography_name' = Benamning) |>
    dplyr::select(topography_code, topography_name) |>
    dplyr::distinct()

  pmid_data <- readr::read_tsv(
    file = file.path(
      cache_dir, "mitelmandb", "REF.TXT.DATA"),
    quote = "none", show_col_types = F, guess_max = 1000000) |>
    dplyr::select(RefNo, Pubmed) |>
    dplyr::filter(nchar(Pubmed) > 0) |>
    dplyr::mutate(reference_id = as.character(RefNo)) |>
    dplyr::distinct() |>
    dplyr::rename(source_id = Pubmed) |>
    dplyr::select(reference_id, source_id) |>
    dplyr::distinct()

  fusion_event_data <- as.data.frame(readr::read_tsv(
    file = file.path(
      cache_dir, "mitelmandb", "MBCA.TXT.DATA"),
    show_col_types = F, guess_max = 100000)) |>
    dplyr::filter(stringr::str_detect(GeneShort,"::")) |>
    dplyr::rename(variant = GeneShort,
                  karyotype = KaryShort,
                  reference_id = RefNo,
                  morphology_code = Morph,
                  topography_code = Topo) |>
    dplyr::mutate(reference_id = as.character(reference_id),
                  morphology_code = as.character(morphology_code)) |>
    dplyr::select(MolClin, 
                  reference_id, 
                  morphology_code, 
                  topography_code,
                  variant, karyotype) |>
    tidyr::separate_rows(variant, sep = ",") |>
    dplyr::left_join(pmid_data, by = "reference_id", 
                     multiple = "all", 
                     relationship = "many-to-many") |>
    dplyr::left_join(morphology, by = "morphology_code", 
                     multiple = "all", 
                     relationship = "many-to-many") |>
    dplyr::left_join(topography, by = "topography_code", 
                     multiple = "all", 
                     relationship = "many-to-many") |>
    dplyr::filter(stringr::str_detect(variant,"::")) |>
    dplyr::mutate(variant = stringr::str_replace(
      variant, "\\+", "")) |>
    dplyr::mutate(evidence_id = dplyr::row_number()) |>
    dplyr::mutate(evidence_id = paste0("MITDB_",evidence_id)) |>
    dplyr::filter(stringr::str_count(variant,",") ==
                    stringr::str_count(karyotype,",")) |>
    tidyr::separate_rows(c(variant,karyotype),
                         sep = ",") |>
    dplyr::mutate(topography_name2 = dplyr::case_when(
      is.na(topography_name) & 
        stringr::str_detect(
          tolower(morphology_name), "leukemia|myeloma") ~ "Myeloid",
      is.na(topography_name) & 
        stringr::str_detect(
        tolower(morphology_name), "lymphoma") ~ "Lymphoid",
      topography_name == "Ovary" ~ "Ovary/Fallopian Tube",
      topography_name == "Fallopian tube" ~ "Ovary/Fallopian Tube",
      topography_name == "Uterus, corpus" ~ "Uterus",
      stringr::str_detect(topography_name, 
        paste0(
          "(Nasal cavity/Paranasal sinuses|Nasopharynx|",
          "Oral cavity|Oro- and hypopharynx|Salivary gland|",
          "Larynx|Teeth|Tongue)")) ~ "Head and Neck",
      topography_name == "Gallbladder/Biliary system" ~ "Biliary Tract",
      topography_name == "Soft tissue" ~ "Soft Tissue",
      topography_name == "Brain" ~ "CNS/Brain",
      morphology_name == "Mature B-cell neoplasm, NOS" ~ "Lymphoid",
      topography_name == "Brain stem" ~ "CNS/Brain",
      topography_name == "Cerebellum" ~ "CNS/Brain",
      topography_name == "Pineal body" ~ "CNS/Brain",
      topography_name == "Vagina" ~ "Vulva/Vagina",
      topography_name == "Skeleton" ~ "Bone",
      stringr::str_detect(
        tolower(morphology_name),
        "myelodysplastic|myeloproliferative") ~ "Myeloid",
      stringr::str_detect(
        tolower(morphology_name),
        "miscellaneous hematopoietic/lymphoid|hodgkin disease") ~ "Lymphoid",
      topography_name == "Ureter" ~ "Kidney",
      topography_name == "Unknown site" ~ "Any",
      topography_name == "Urethra" ~ "Bladder/Urinary Tract",
      topography_name == "Spinal cord" ~ "Peripheral Nervous System",
      topography_name == "Uterus, cervix" ~ "Cervix",
      stringr::str_detect(
        topography_name,"Large intestine|Anus") ~ "Colon/Rectum",
      topography_name == "Bladder" ~ "Bladder/Urinary Tract",
      topography_name == "Adrenal" ~ "Adrenal Gland",
      stringr::str_detect(
        topography_name,
        "Intraabdominal|Oesophagus|Stomach") ~ "Esophagus/Stomach",
      TRUE ~ topography_name
    )) |>
    dplyr::filter(
      topography_name2 != "Heart" &
        topography_name2 != "Spleen" &
        topography_name2 != "Penis" &
        topography_name2 != "Parathyroid" &
        topography_name2 != "Small intestine" &
        topography_name2 != "Intrathoracal" &
        topography_name2 != "Lymph node"
    ) |>
    dplyr::rename(primary_diagnosis = morphology_name,
                  primary_site = topography_name2) |>
    dplyr::select(variant, 
                  karyotype, 
                  primary_diagnosis,
                  primary_site,
                  source_id, 
                  evidence_id) |>
    dplyr::distinct() |>
    dplyr::mutate(variant_id = paste(
      variant, karyotype, sep = " - "
    ))

  evidence_id <- nrow(fusion_event_data) + 1

  fusion_event_data$variant_alias <- paste(
    fusion_event_data$variant,
    paste(
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,2],
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,1],
      sep = "-"
    ),
    paste(
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,1],
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,2],
      sep = "/"
    ),
    paste(
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,1],
      stringr::str_split_fixed(fusion_event_data$variant,"::",2)[,2],
      sep = "-"
    ),
    fusion_event_data$karyotype,
    sep = "@@@"
  )
  
  fusion_event_data$variant_name_primary <- fusion_event_data$variant

  bcr_abl_custom <-
    data.frame('variant_id' = "BCR::ABL1 - t(9;22)(q34;q11)",
               'variant_alias' = 'BCR::ABL',
               'karyotype' = 't(9;22)(q34;q11)',
               'variant_name_primary' = 'BCR::ABL1',
               'evidence_id' = paste0('MITDB_',evidence_id),
               'primary_diagnosis' = NA,
               'primary_site' = NA,
               'source_id' = NA,
               'symbol' = 'BCR-ABL1',
               stringsAsFactors = F)

  mitelman_db <- as.data.frame(
    fusion_event_data |>
      dplyr::select(variant_alias, 
                    variant_id,
                    karyotype,
                    variant_name_primary,
                    primary_diagnosis, 
                    primary_site,
                    evidence_id, 
                    source_id) |>
      tidyr::separate_rows(variant_alias, sep = "@@@") |>

      ## Remove aliases where fusion partners are mentioned in reverse order
      #dplyr::filter(
        # !stringr::str_detect(
        #   variant_alias,
        #   paste0(
        #     "RARA::PML|ABL1::BCR|ERG::TMPRSS2|RUNX1T1::RUNX1|NSD1::NUP98|",
        #     "MYH11::CBFB|RET::KIF5B|RET::CCDC6|ROS1::CD74|CAMTA1::WWTR1|",
        #     "NPM1::MLF1|BRAF::KIAA1549|CREBBP::KAT6A|MLF1::NPM1"))) |>
      dplyr::filter(
        !stringr::str_detect(
          variant_name_primary,
          paste0(
            "RARA::PML|ABL1::BCR|ERG::TMPRSS2|RUNX1T1::RUNX1|NSD1::NUP98|",
            "MYH11::CBFB|RET::KIF5B|RET::CCDC6|ROS1::CD74|CAMTA1::WWTR1|",
            "NPM1::MLF1|BRAF::KIAA1549|CREBBP::KAT6A|MLF1::NPM1"))) |>
      
      dplyr::bind_rows(
        bcr_abl_custom
      ) |>
      dplyr::mutate(alteration_type = "FUSION",
                    variant_origin = "somatic",
                    variant_consequence = "transcript_fusion",
                    evidence_url = "https://mitelmandatabase.isb-cgc.org/",
                    biomarker_source = "mitelmandb",
                    #biomarker_entity = T,
                    #alias_type = "other_gene",
                    biomarker_source_datestamp = db_datestamp) |>
      dplyr::arrange(evidence_id) |>
      dplyr::rename(citation_id = source_id) |>
      dplyr::select(biomarker_source,
                    biomarker_source_datestamp,
                    variant_id,
                    variant_alias,
                    #alias_type,
                    variant_name_primary,
                    alteration_type,
                    dplyr::everything())

  )

  mitelman_db$symbol <- stringr::str_split_fixed(
    mitelman_db$variant_id, " - ", 2)[,1]

  
  biomarker_items <- list()
  biomarker_items[['clinical']] <- mitelman_db
  
  biomarker_items[['variant']] <- biomarker_items[['clinical']] |>
    #dplyr::select(variant_id, variant_alias, alteration_type, symbol) |>
    dplyr::mutate(alias_type = "other_gene") |>
    dplyr::rename(gene = symbol) |>
    #dplyr::rename(variant_alias = variant) |>
    dplyr::select(biomarker_source,
                  biomarker_source_datestamp,
                  variant_id, 
                  variant_alias, 
                  variant_name_primary,
                  alteration_type,
                  alias_type,
                  variant_consequence,
                  gene) |>
    dplyr::distinct()
  
  biomarker_items[['clinical']] <- biomarker_items[['clinical']] |>
    dplyr::select(-c("variant_alias", "variant_name_primary",
                     "alteration_type","variant_consequence",
                     "symbol")) |>
    dplyr::distinct()
  
  biomarker_items <- 
    get_literature_references(biomarker_items)
  
  return(biomarker_items)
}
