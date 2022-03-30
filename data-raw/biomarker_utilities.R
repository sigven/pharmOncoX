
expand_hgvs_terms <- function(var) {
  hits <- c()
  if(stringr::str_detect(toupper(var),"([A-Z]{3}[0-9]{1,}(([A-Z]{3}){0,1})([0-9]{1,})?)") &
     !stringr::str_detect(toupper(var),"(DEL|INS|DUP)|=|-|\\*")){
    aa <- stringr::str_to_title(strsplit(var, "[0-9]{1,}")[[1]])
    codon <- stringr::str_match(var, "[0-9]{1,}")[[1]]

    #cat(var,' - ', aa, ' - ', codon, '\n')
    if(length(aa) == 1){
      alt1 <- paste0(aa[1],codon)
      alt2 <- paste0('p.', alt1)
      if(aa[1] != "Ter"){
        aa1 <- seqinr::a(aa[1])
        if(!is.na(aa1)){
          alt3 <- paste0(seqinr::a(aa[1]), codon)
          hits <- c(alt1, alt2, alt3)
        }
      }else{
        alt3 <- paste0("X", codon)
        hits <- c(alt1, alt2, alt3)

      }
    }else{
      if(length(aa) == 2){
        alt1 <- stringr::str_replace(
          stringr::str_replace(paste0(aa[1],codon, aa[2]),"(FS|Fs)$","fs"),
          "ter","Ter")
        alt2 <- stringr::str_replace(
          stringr::str_replace(paste0('p.', alt1),"(FS|Fs)$","fs"),
          "ter","Ter")
        if(aa[1] != "Ter" & aa[2] != "Ter" & aa[2] != "Fs" &
           nchar(aa[1]) == 3 & nchar(aa[2]) == 3){
          alt3 <- paste0(seqinr::a(aa[1]), codon, seqinr::a(aa[2]))
          alt4 <- paste0('p.', alt3)
          hits <- c(alt1, alt2, alt3, alt4)
        }
        else{
          if(aa[1] == "Ter"){
            alt3 <- paste0("X", codon, seqinr::a(aa[2]))
            alt4 <- paste0('p.', alt3)
            hits <- c(alt1, alt2, alt3, alt4)
          }else{
            if(aa[2] == "Ter"){
              alt3 <- paste0(seqinr::a(aa[1]), codon, "X")
              alt4 <- paste0('p.', alt3)
              hits <- c(alt1, alt2, alt3, alt4)
            }
            if(aa[2] == "Fs"){
              alt3 <- paste0(seqinr::a(aa[1]), codon, "fs")
              alt4 <- paste0('p.', alt3)
              hits <- c(alt1, alt2, alt3, alt4)
            }

          }
        }
      }
    }
  }
  if(stringr::str_detect(var,"^([A-Z]{1}[0-9]{1,}[A-Z]{0,1})$") &
     !stringr::str_detect(var, "FS|DUP|DEL|INS|fs|dup|del|ins|\\*")){
    aa <- stringr::str_to_title(strsplit(var, "[0-9]{1,}")[[1]])
    codon <- stringr::str_match(var, "[0-9]{1,}")[[1]]

    if(length(aa) == 1){
      alt1 <- paste0(aa[1],codon)
      alt2 <- paste0('p.', alt1)
      if(aa[1] != "X"){
        alt3 <- paste0(seqinr::aaa(aa[1]), codon)
        hits <- c(alt1, alt2, alt3)
      }else{
        alt3 <- paste0("Ter", codon)
        hits <- c(alt1, alt2, alt3)
      }
    }else{
      if(length(aa) == 2){
        alt1 <- paste0(aa[1],codon, aa[2])
        alt2 <- paste0('p.', alt1)
        if(aa[1] != "X" & aa[2] != "X" & nchar(aa[1]) == 1 & nchar(aa[2]) == 1){
          alt3 <- paste0(seqinr::aaa(aa[1]),
                         codon, seqinr::aaa(aa[2]))
          alt4 <- paste0('p.', alt3)
          hits <- c(alt1, alt2, alt3, alt4)
        }else{
          if(aa[1] == "X"){
            alt3 <- paste0("Ter",
                           codon, seqinr::aaa(aa[2]))
            alt4 <- paste0('p.', alt3)
            hits <- c(alt1, alt2, alt3, alt4)
          }else{
            if(aa[2] == "X"){
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
  if(stringr::str_detect(var,"^C\\.")){
    var <- stringr::str_replace(var,"^C\\.","c.")
    #hits <- res
  }
  if(stringr::str_detect(var,"^P\\.")){
    var <- stringr::str_replace(var,"^P\\.","p.")
    #hits <- res
  }
  if(stringr::str_detect(var,"FS$")){
    var <- stringr::str_replace(var,"FS$","fs")
  }
  if(stringr::str_detect(var,";") & stringr::str_detect(var,"\\(")){
    var <- tolower(var)
  }
  if(stringr::str_detect(var,"^([A-Z]{1}[0-9]{1,}[A-Z]{1})$")){
    res1 <- paste0('p.', var)
    res2 <- stringr::str_match(var,"^([A-Z]{1}[0-9]{1,})")[[1]]
    hits <- c(hits,res1,res2)
  }
  if(stringr::str_detect(var,"FS$")){
    res <- stringr::str_replace(var,"FS$","fs")
    hits <- res
  }
  if(stringr::str_detect(var,"(DEL|Del|INS|Ins|Dup|DUP)") &
     !stringr::str_detect(var, "INSERTION|DUPLICATION|NM_|NP_|DELETERIOUS|DELETION|MUTATION|Mutation|Deletion|")){
    var <- stringr::str_replace(var,"(DEL|Del)","del")
    var <- stringr::str_replace(var,"(INS|Ins)","ins")
    var <- stringr::str_replace(var,"(DUP|Dup)","dup")
    if(!(stringr::str_detect(var,"^(c|p|rs|;|[0-9])"))){
      res2 <- paste0('p.', var)
      hits <- c(var, res2)
    }
  }

  if(length(hits) == 0 & var != ""){
    if(!(stringr::str_detect(var,"^(c|rs|p|[0-9])")) &
       stringr::str_detect(var, "[0-9]{1,}") &
       !stringr::str_detect(var,"EXON|;|MUTATION|Mutation|-|NM_|NP_|DELETERIOUS|Deletion|INTRON|DELETION|DUPLICATION")){

      if(stringr::str_detect(var,"FS$|DUP|DELINS|FS\\*")){
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

load_civic_biomarkers <- function(datestamp = '20211217'){


  clinicalEvidenceSummary <- as.data.frame(
    data.table::fread(paste0("data-raw/biomarkers/civic/clinical_evidence_summary_",
                             datestamp,
                             ".tsv"),
                      select = c(1:41), fill = T)) %>%
    dplyr::filter(source_type == "PubMed") %>%
    dplyr::select(variant_id,
                  evidence_id,
                  evidence_civic_url,
                  disease,
                  doid,
                  drugs,
                  citation_id,
                  evidence_statement,
                  #chromosome,
                  #start,
                  #stop,
                  #reference_bases,
                  #variant_bases,
                  evidence_type,
                  evidence_level,
                  evidence_direction,
                  clinical_significance,
                  variant_origin,
                  rating) %>%
    dplyr::rename(disease_ontology_id = doid,
                  evidence_description = evidence_statement,
                  evidence_url = evidence_civic_url,
                  cancer_type = disease) %>%
                  #chr_start = start,
                  #chr_stop = stop,
                  #refbase = reference_bases,
                  #altbase = variant_bases) %>%

  dplyr::mutate(disease_ontology_id = paste0(
    "DOID:", disease_ontology_id)
  ) %>%
    dplyr::mutate(
      evidence_level =
        dplyr::case_when(evidence_level == 'A' ~ 'A: Validated',
                         evidence_level == 'B' ~ 'B: Clinical evidence',
                         evidence_level == 'C' ~ 'C: Case study',
                         evidence_level == 'D' ~ 'D: Preclinical evidence',
                         evidence_level == 'E' ~ 'E: Indirect evidence',
                         TRUE ~ as.character(evidence_level)
        )
    ) %>%

    dplyr::mutate(
      variant_origin =
        dplyr::if_else(
          variant_origin == "N/A" | variant_origin == "Unknown",
          as.character(NA),
          as.character(stringr::str_replace(variant_origin,"Rare |Common ",""))))


  therapeutic_contexts <- as.data.frame(
    clinicalEvidenceSummary %>%
      dplyr::select(evidence_id, drugs) %>%
      dplyr::filter(!is.na(drugs) &
                      nchar(drugs) > 0) %>%
      tidyr::separate_rows(drugs, sep=",") %>%
      dplyr::mutate(drugs = tolower(drugs)) %>%
      dplyr::left_join(oncoPharmaDB::compound_synonyms,
                       by = c("drugs" = "alias")) %>%
      dplyr::filter(!is.na(nci_concept_display_name)) %>%
      dplyr::select(evidence_id, molecule_chembl_id) %>%
      dplyr::distinct() %>%
      dplyr::group_by(evidence_id) %>%
      dplyr::summarise(molecule_chembl_id = paste(
        unique(sort(molecule_chembl_id)), collapse="|"),
        .groups = "drop"
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(molecule_chembl_id = stringr::str_replace(
        molecule_chembl_id,"^NA\\||\\|NA$",""
      ))
  )

  clinicalEvidenceSummary <- clinicalEvidenceSummary %>%
    dplyr::left_join(therapeutic_contexts, by = "evidence_id")

  variantSummary <- as.data.frame(
    data.table::fread(paste0("data-raw/biomarkers/civic/variant_summary_",
                             datestamp,
                             ".tsv"),
                      select = c(1:30), fill = T) %>%
      dplyr::mutate(alteration_type = "MUT") %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant),"EXPRESSION"),
            "EXP",as.character(alteration_type))) %>%
      dplyr::mutate(
        variant_types =
          dplyr::if_else(
            stringr::str_detect(toupper(variant),"EXPRESSION") &
              variant_types == "N/A",
            "transcript_amplification",
            as.character(variant_types))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::case_when(
            alteration_type == "EXP" &
              stringr::str_detect(toupper(variant),"OVER") ~ "EXP_OVER",
            alteration_type == "EXP" &
              stringr::str_detect(toupper(variant),"UNDER") ~ "EXP_UNDER",
            TRUE ~ as.character(alteration_type)
          )
      ) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant),"METHYLAT"),
            "METHYL",as.character(alteration_type))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant),"PHOSPHORY"),
            "PHOSPHO",as.character(alteration_type))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant),"LOSS-OF-FUNCTION|LOSS|TRUNCAT|INACTIVAT") |
              stringr::str_detect(
                variant_types, "stop_gained,loss_of_function_variant|loss_of_function_variant,stop_gained") |
              (stringr::str_detect(variant_types, "loss_of_function_variant") &
                 toupper(variant) == "MUTATION") |
              stringr::str_detect(variant_types,"truncation"),
            "MUT_LOF",as.character(alteration_type))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            alteration_type == 'MUT_LOF' &
              (stringr::str_detect(toupper(variant),"FRAMESHIFT") |
                 stringr::str_detect(variant_types,"truncation")),
            "MUT_LOF_FS",as.character(alteration_type))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant), "AMPLIFICATION|DELETION|COPY") &
              stringr::str_detect(variant_types, "_amplification|_ablation"),
            "CNA",as.character(alteration_type))) %>%

      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            stringr::str_detect(toupper(variant_types),"TRANSCRIPT_FUSION|REGION_FUSION") |
              stringr::str_detect(toupper(variant),"TRANSLOCATION|FUSION|REARRANGEMENT") |
              variant_types == "" & stringr::str_detect(variant,"-"),
            "TRANSLOCATION_FUSION",as.character(alteration_type))) %>%
      dplyr::mutate(
        alteration_type =
          dplyr::if_else(
            (stringr::str_detect(toupper(variant_types),"TRANSCRIPT_FUSION") &
               stringr::str_detect(toupper(variant_types),"MISSENSE")) |
              variant_types == "missense_variant" & stringr::str_detect(variant,"BCR-ABL"),
            "TRANSLOCATION_FUSION_MUT",as.character(alteration_type))) %>%

      dplyr::select(variant_id,
                    variant_types,
                    variant,
                    alteration_type,
                    gene,
                    variant_aliases) %>%
      tidyr::separate_rows(variant_aliases, sep=",") %>%

      dplyr::mutate(
        variant =
          stringr::str_replace(variant,
                               " \\(.+\\)$","")) %>%
      dplyr::mutate(
        variant_aliases =
          stringr::str_replace(variant_aliases,
                               " \\(.+\\)$","")) %>%
      dplyr::mutate(
        variant =
          stringr::str_replace(variant,
                               "\\*$","X")) %>%
      dplyr::mutate(
        variant_aliases =
          stringr::str_replace(variant_aliases,
                               "\\*$","X")) %>%
      dplyr::mutate(variant = dplyr::if_else(
        variant_types == "N/A",
        stringr::str_replace(toupper(variant),"P16 |NUCLEAR |ISOFORM ",""),
        as.character(variant)
      )) %>%
      dplyr::mutate(
        variant =
          stringr::str_replace(
            toupper(variant)," HOMOZYGOSITY|CYTOPLASMIC |DECREASED PERI-THERAPEUTIC ","")
      ) %>%
      dplyr::mutate(
        variant =
          stringr::str_replace(
            toupper(variant)," FUSIONS| E20-E20| E6-E20| E2-E20|TYPE 1","")
      ) %>%
      dplyr::mutate(
        variant =
          stringr::str_replace(
            toupper(variant)," RECIPROCAL |RARE | AND AMPLIFICATION| FUSION","")
      ) %>%
      dplyr::mutate(
        variant_aliases =
          stringr::str_replace(
            toupper(variant_aliases),"FUSION | TYPE 1|RECIPROCAL ","")
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(variant = stringr::str_replace(
        variant, "::", "-"
      )) %>%
      dplyr::mutate(gene2 = dplyr::if_else(
        stringr::str_detect(variant,"BCR-ABL") |
          variant_types == "missense_variant,transcript_fusion" |
          variant_types == "inframe_insertion,transcript_fusion" |
          variant_types == "transcript_fusion" |
          variant_types == "transcript_regulatory_region_fusion" |
          variant_types == "gene_fusion,transcript_fusion",
        stringr::str_match(variant,"\\S{1,}-\\S{1,}")[[1]],
        as.character(gene)
      )) %>%
      dplyr::mutate(gene2 = dplyr::if_else(
        is.na(gene2) &
          (variant_types == "missense_variant,transcript_fusion" |
             variant_types == "transcript_fusion"),
        as.character(gene),
        as.character(gene2)
      )) %>%
      dplyr::mutate(variant_aliases = dplyr::if_else(
        stringr::str_detect(variant_aliases,"^RS[0-9]{5,}$"),
        stringr::str_replace(variant_aliases,"RS","rs"),
        as.character(variant_aliases))) %>%
      dplyr::mutate(variant = dplyr::if_else(
        stringr::str_detect(variant,"^RS[0-9]{5,}$"),
        stringr::str_replace(variant,"RS","rs"),
        as.character(variant))) %>%
      dplyr::mutate(variant_aliases = dplyr::if_else(
        stringr::str_detect(variant_types,"fusion"),
        stringr::str_replace(variant_aliases,"ALK(S)? |(\\S{1,}-\\S{1,}\\s?)",""),
        as.character(variant_aliases)
      )) %>%
      dplyr::mutate(variant = dplyr::if_else(
        stringr::str_detect(variant_types,"fusion"),
        stringr::str_replace(variant,"ALK(S)? |(\\S{1,}-\\S{1,}\\s?)",""),
        as.character(variant)
      )) %>%
      dplyr::mutate(variant = dplyr::if_else(
        alteration_type == "TRANSLOCATION_FUSION_MUT" |
          stringr::str_detect(variant,"BCR-ABL(1)? [A-Z]"),
        stringr::str_replace(variant,"BCR-ABL(1)? ",""),
        as.character(variant)
      )) %>%
      dplyr::mutate(variant_aliases = dplyr::if_else(
        alteration_type == "TRANSLOCATION_FUSION_MUT" |
          stringr::str_detect(variant_aliases,"BCR-ABL(1)? [A-Z]"),
        stringr::str_replace(variant_aliases,"BCR-ABL(1)? ",""),
        as.character(variant_aliases)

      )) %>%
      dplyr::mutate(variant_types = dplyr::if_else(
        stringr::str_detect(variant_types,"") &
          stringr::str_detect(variant,"EXON") &
          stringr::str_detect(variant,"DELETION"),
        "exon_loss_variant",
        as.character(variant_types)
      )) %>%
      dplyr::mutate(variant_types = dplyr::if_else(
        stringr::str_detect(variant_types,"") &
          stringr::str_detect(variant,"EXON") &
          stringr::str_detect(variant,"MUTATION"),
        "exon_variant",
        as.character(variant_types)
      )) %>%
      dplyr::mutate(variant_types = dplyr::if_else(
        stringr::str_detect(variant_types,"") &
          stringr::str_detect(variant,"^(LOSS|LOSS-OF-FUNCTION)$"),
        "loss_of_function_variant",
        as.character(variant_types)
      )) %>%
      dplyr::mutate(
        variant = stringr::str_replace_all(
          variant, "(PRPS1 |IKZF1 |EGFR |FGFR2 | E2-E11| E11-E10|MET1 |F547 SPLICE SITE )",
          ""
        )) %>%
      dplyr::mutate(
        variant = stringr::str_replace_all(
          variant, "^(Mutations|MUTATIONS|MUT|ALTERATION|Deletion And Mutation|Mutation)$",
          "MUTATION"
        )
      ) %>%
      dplyr::mutate(
        variant_aliases = dplyr::if_else(
          stringr::str_detect(variant_aliases,"^(INV|T)\\("),
          tolower(variant_aliases),
          as.character(variant_aliases)
        )
      )
  )

  expanded_variant_set <- data.frame()

  for(i in unique(variantSummary$variant_id)){
    vrows <- variantSummary[variantSummary$variant_id == i,]

    variant_entries <- unique(c(vrows$variant,
                                vrows$variant_aliases))
    gene <- unique(vrows$gene2)
    variant_type <- unique(vrows$variant_types)
    alteration_type <- unique(vrows$alteration_type)

    all_aliases <- c()

    if(stringr::str_detect(alteration_type,"MUT")){

      for(e in variant_entries){
        clean_aliases <- expand_hgvs_terms(e)
        for(c in clean_aliases){
          all_aliases <- unique(c(all_aliases, c))
        }
      }


      for(variant_alias in all_aliases){
        df <- data.frame(
          'variant_id' = i,
          'symbol' = gene,
          'alteration_type' = alteration_type,
          'variant_type' = variant_type,
          'variant' = variant_alias,
          stringsAsFactors = F)
        expanded_variant_set <- dplyr::bind_rows(
          expanded_variant_set, df
        )
      }
    }

  }

  expanded_variant_set <- expanded_variant_set %>%
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(variant,"^(p\\.)?(([A-Z]{1}[0-9]{1,})|([A-Z]{1}[a-z]{2}[0-9]{1,}))$"),
          as.character("CODON"),
          as.character(alteration_type)
        )) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"FS"),
        stringr::str_replace(variant,"FS","fs"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"DEL") &
          !stringr::str_detect(variant, "DELETION"),
        stringr::str_replace(variant,"DEL","del"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"INS") &
          !stringr::str_detect(variant,"INSERTION"),
        stringr::str_replace(variant,"INS","ins"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"TER|ter") &
          !stringr::str_detect(tolower(variant),"alternative|internal|promoter|deleterious"),
        stringr::str_replace(variant,"ter|TER","Ter"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"DUP") &
          !stringr::str_detect(variant,"DUPLICATION"),
        stringr::str_replace(variant,"DUP","dup"),
        as.character(variant)
      )
    )

  for(aa in c("PRO","ILE","HIS","GLY","GLU","LEU","THR",
              "TYR","SER","ASN","ASP","LYS","PHE")){
    aa_three_lc <- stringr::str_to_title(aa)
    if(nrow(expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]) > 0){
      expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]$variant <-
        stringr::str_replace(expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]$variant,
                             aa,
                             aa_three_lc)
    }
  }

  variant_set_other <- variantSummary %>%
    dplyr::filter(!stringr::str_detect(alteration_type,"MUT")) %>%
    dplyr::select(variant_id, variant_types,
                  variant, alteration_type,
                  gene2) %>%
    dplyr::rename(variant_type = variant_types,
                  symbol = gene2) %>%
    dplyr::filter(!stringr::str_detect(variant,"PHOSPHO")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "INSerTION","INSERTION")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "CONSerVED","CONSERVED")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "POLYMORPHisM","POLYMORPHISM")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "PHILAdelPHIA","PHILADELPHIA")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "ProMOTER","PROMOTER")) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "delETERIOUS","DELETERIOUS"))

  variant_evidence_full <-
    clinicalEvidenceSummary %>%
    dplyr::inner_join(dplyr::bind_rows(expanded_variant_set, variant_set_other)) %>%
    dplyr::rename(therapeutic_context = drugs) %>%
    dplyr::mutate(variant = stringr::str_replace_all(
      variant, "delETERIOUS","DELETERIOUS")) %>%
    dplyr::select(-c(rating)) %>%
    dplyr::mutate(evidence_id = paste0("EID",evidence_id)) %>%
    dplyr::mutate(biomarker_source_db = "civic",
                  biomarker_entity = T)

  biomarkers_all <- variant_evidence_full %>%
    dplyr::select(evidence_id,
                  evidence_url,
                  variant_id,
                  variant_type,
                  variant,
                  alteration_type,
                  symbol,
                  therapeutic_context,
                  disease_ontology_id,
                  cancer_type,
                  molecule_chembl_id,
                  dplyr::everything()
    )

  biomarkers_all <- map_biomarker_phenotypes(biomarkers_all) %>%
    dplyr::select(-do_cancer_slim) %>%


  return(biomarkers_all)

}


load_cgi_biomarkers <- function(){

  cancer_type_abbreviations <-
    read.table(file = "data-raw/biomarkers/cgi/cancer_subtypes.tsv",sep = "\t", header = T,
               stringsAsFactors = F, quote = "", na.strings = "",fill = T) %>%
    dplyr::mutate(disease_ontology_id = paste0("DOID:",disease_ontology_id))

  ## SNVs/InDels
  ## SKIP/IGNORE mutations that are
  # 1) not explicitly mapped to protein alterations (MUT)
  # 2) marked with uncertainty, rechecking etc (Comments)
  # 3) involved in toxicity
  # 4) curated from abstracts etc (not PMIDs or guidelines)
  # 5) Non-malignant/Non-cancer conditions
  cgi_biomarker_pr_variant <- as.data.frame(
    read.table(file = "data-raw/biomarkers/cgi/cgi_biomarkers_per_variant.tsv", sep = "\t",
               header = T, stringsAsFactors = F,
               na.strings = c("","[]"), fill = T, quote = "") %>%
      dplyr::select(-c(Assay.type, Curation.date, Curator,
                       TCGI.included, Targeting, region, strand)) %>%
      janitor::clean_names() %>%
      dplyr::distinct() %>%
      dplyr::select(alteration_type,
                    association, source,
                    biomarker, gene, comments,
                    drug_full_name,
                    evidence_level,
                    primary_tumor_acronym,
                    primary_tumor_type,
                    individual_mutation,
                    info, g_dna) %>%
      dplyr::mutate(variant_id = rep(1:nrow(.))) %>%
      dplyr::mutate(gene = dplyr::case_when(
        gene == "MLL" ~ "KMT2A",
        gene == "MLL2" ~ "KMT2D",
        TRUE ~ as.character(gene)
      )) %>%
      dplyr::mutate(variant_origin = "Somatic") %>%
      dplyr::mutate(variant_origin = dplyr::if_else(
        !is.na(comments) & stringr::str_detect(tolower(comments),"germline"),
        "Germline",
        as.character(variant_origin)
      )) %>%
      ## Ignore composite biomarkers for now
      dplyr::filter(!stringr::str_detect(alteration_type,";")) %>%
      dplyr::mutate(
        alteration_type = dplyr::case_when(
          alteration_type == "EXPR" &
            stringr::str_detect(toupper(biomarker),"OVER") ~ "EXP_OVER",
          alteration_type == "EXPR" &
            stringr::str_detect(toupper(biomarker),"UNDE") ~ "EXP_UNDER",
          alteration_type == "EXPR" &
            !stringr::str_detect(toupper(biomarker),"OVER|UNDE") ~ "EXP",
          alteration_type == "FUS" ~ "TRANSLOCATION_FUSION",
          TRUE ~ as.character(alteration_type)
        )) %>%

      tidyr::separate_rows(primary_tumor_acronym, sep=";") %>%
      dplyr::filter(!stringr::str_detect(association,"Toxicity")) %>%
      dplyr::filter(!stringr::str_detect(source, "^ENA|ASCO|AACR|ESMO|http")) %>%
      dplyr::filter(primary_tumor_acronym != 'SM' &
                      primary_tumor_acronym != 'AA' &
                      primary_tumor_acronym != 'HES' &
                      primary_tumor_acronym != 'WM') %>%
      dplyr::rename(citation_id = source) %>%
      dplyr::left_join(cancer_type_abbreviations,
                       by = c("primary_tumor_acronym" = "cgi_abbreviation")) %>%
      dplyr::distinct() %>%
      tidyr::separate(individual_mutation, c('symbol','alteration'),sep=":") %>%
      dplyr::mutate(symbol = dplyr::if_else(is.na(symbol) & !is.na(gene),
                                            gene,
                                            as.character(symbol))) %>%
      dplyr::mutate(
        eitem_consequence =
          dplyr::if_else(stringr::str_detect(toupper(biomarker), "AMPLIFICATION") &
                           alteration_type == "CNA",
                         "gain",
                         as.character(NA))) %>%
      dplyr::mutate(
        eitem_consequence =
          dplyr::if_else(stringr::str_detect(toupper(biomarker), "DELETION") &
                           alteration_type == "CNA",
                         "loss",
                         as.character(eitem_consequence))) %>%
      dplyr::mutate(
        eitem_consequence =
          dplyr::if_else(stringr::str_detect(biomarker,"(inframe )?insertion") &
                           alteration_type == "MUT",
                         "inframe_insertion",
                         as.character(eitem_consequence))) %>%
      dplyr::mutate(
        eitem_consequence =
          dplyr::if_else(stringr::str_detect(biomarker,"(inframe )?deletion") &
                           alteration_type == "MUT",
                         "inframe_deletion",
                         as.character(eitem_consequence))) %>%
      dplyr::mutate(
        eitem_exon =
          dplyr::if_else(stringr::str_detect(toupper(biomarker),"EXON [0-9]{1,}") &
                           alteration_type == "MUT",
                         stringr::str_match(toupper(biomarker),"EXON [0-9]{1,}"),
                         as.character(NA))) %>%
      dplyr::mutate(eitem_exon =
                      stringr::str_replace(eitem_exon,"EXON ","")) %>%
      #dplyr::mutate(biomarker)
      #dplyr::filter(!stringr::str_detect(alteration_type,";|BIA")) %>%
      dplyr::mutate(evidence_type = 'Predictive',
                    therapeutic_context = drug_full_name,
                    clinical_significance = NA) %>%
      dplyr::mutate(
        clinical_significance =
          dplyr::case_when(
            association == "Resistant" ~ "Resistance or Non-response",
            association == "No Responsive" ~ "Resistance or Non-response",
            association == "Responsive" ~ "Sensitivity",
            TRUE ~ as.character("Other"))) %>%
      dplyr::rename(evidence_level_raw = evidence_level) %>%
      dplyr::mutate(mapping_rank = 1) %>%
      dplyr::mutate(
        alteration = dplyr::if_else(
          !is.na(biomarker) & stringr::str_detect(biomarker,"oncogenic mutation"),
          "MUTATION",
          as.character(alteration)
        )
      ) %>%
      # dplyr::mutate(
      #   biomarker_mapping =
      #     dplyr::if_else(alteration_type == "MUT" & alteration != "MUTATION",
      #                    as.character("exact"),
      #                    as.character("gene"))) %>%
      # dplyr::mutate(
      #   biomarker_mapping =
      #     dplyr::if_else(!is.na(eitem_exon),
      #                    as.character("exon"),
      #                    as.character(biomarker_mapping))
      # ) %>%
    dplyr::mutate(
      evidence_level =
        dplyr::if_else(stringr::str_detect(evidence_level_raw,"guidelines"),
                       "A: FDA/NCCN/ELN guidelines",
                       as.character(NA))) %>%
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(stringr::str_detect(evidence_level_raw,"Late trials"),
                         "B1: Clinical evidence: late trials",
                         as.character(evidence_level))) %>%
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(stringr::str_detect(evidence_level_raw,"Early trials"),
                         "B2: Clinical evidence: early trials",
                         as.character(evidence_level))) %>%
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(stringr::str_detect(evidence_level_raw,"Case report"),
                         "C: Case study",
                         as.character(evidence_level))) %>%
      dplyr::mutate(
        evidence_level =
          dplyr::if_else(stringr::str_detect(evidence_level_raw,"Pre-clinical"),
                         "D: Preclinical evidence",
                         as.character(evidence_level))) %>%
      dplyr::filter(!is.na(evidence_level_raw)) %>%
      dplyr::select(-c(association, evidence_level_raw,
                       primary_tumor_acronym, drug_full_name, gene)) %>%
      dplyr::select(-c(info, g_dna)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(evidence_id = rep(1:nrow(.))) %>%
      dplyr::mutate(evidence_id = paste0("CGI_",evidence_id)) %>%
      dplyr::mutate(evidence_type = "Predictive") %>%
      #dplyr::select(-c(biomarker, biomarker_mapping)) %>%
      dplyr::distinct()
  )


  all_mut_entries <- as.data.frame(
    cgi_biomarker_pr_variant  %>%
      dplyr::filter(!is.na(alteration) &
                      alteration != "MUTATION" &
                      alteration_type == "MUT") %>%
      dplyr::mutate(
        alteration =
          stringr::str_replace(alteration,
                               "\\*$","X")) %>%
      dplyr::select(evidence_id, alteration)
  )

  all_mut_entries_expanded <- data.frame()
  for(i in 1:nrow(all_mut_entries)){

    alt <- all_mut_entries[i,"alteration"]
    terms <- expand_hgvs_terms(alt)
    for(t in terms){
      df <- data.frame('evidence_id' = all_mut_entries[i,"evidence_id"],
                       'alteration' = t,
                       stringsAsFactors = F)
      all_mut_entries_expanded <- all_mut_entries_expanded %>%
        dplyr::bind_rows(df)
    }
  }

  cgi_biomarker_other <- cgi_biomarker_pr_variant %>%
    dplyr::filter(is.na(alteration) |
                    alteration == "MUTATION" |
                    (alteration != "MUTATION" & alteration_type != "MUT"))

  cgi_biomarker_muts <-
    cgi_biomarker_pr_variant %>%
    dplyr::filter(!is.na(alteration) &
                    alteration != "MUTATION" &
                    alteration_type == "MUT") %>%
    dplyr::select(-alteration) %>%
    dplyr::inner_join(all_mut_entries_expanded,
                      by = "evidence_id")


  cgi_biomarkers_all <- cgi_biomarker_muts %>%
    dplyr::bind_rows(cgi_biomarker_other) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          !is.na(alteration_type) & !is.na(alteration) & stringr::str_detect(alteration,"^(([A-Z]{1}[0-9]{1,})|([A-Z]{1}[a-z]{2}[0-9]{1,}))$"),
          as.character("CODON"),
          as.character(alteration_type)
        )) %>%
    dplyr::rename(variant = alteration) %>%
    dplyr::select(symbol,
                  variant,
                  variant_id,
                  alteration_type,
                  evidence_type,
                  citation_id,
                  therapeutic_context,
                  clinical_significance,
                  evidence_id,
                  variant_origin,
                  disease_ontology_id,
                  primary_tumor_type,
                  evidence_level) %>%
    dplyr::mutate(evidence_direction = NA,
                  variant_type = NA) %>%
    dplyr::mutate(variant_type = dplyr::case_when(
      stringr::str_detect(alteration_type,"MUT|CODON") ~ "missense_variant",
      stringr::str_detect(alteration_type,"CNA") ~ "transcript_amplification",
      stringr::str_detect(alteration_type,"EXP") ~ "transcript_amplification",
      TRUE ~ as.character(NA)
    )) %>%
    dplyr::rename(cancer_type = primary_tumor_type)


  cgi_biomarkers_all <- cgi_biomarkers_all %>%
    dplyr::mutate(therapeutic_context = dplyr::if_else(
      stringr::str_detect(therapeutic_context,
                          "(ib|ide|ole|ant|stat|ine|us|one|mab|in|el|[0-9])( |  )\\("),
      stringr::str_replace(therapeutic_context,
                           "( |  )\\(.*\\)( )?$",""),
      as.character(therapeutic_context)
    ))


  therapeutic_contexts <- as.data.frame(
    cgi_biomarkers_all %>%
      dplyr::select(evidence_id, therapeutic_context) %>%
      dplyr::filter(!is.na(therapeutic_context) &
                      nchar(therapeutic_context) > 0) %>%
      tidyr::separate_rows(therapeutic_context, sep=" \\+ ") %>%
      dplyr::mutate(therapeutic_context = tolower(therapeutic_context)) %>%
      dplyr::left_join(oncoPharmaDB::compound_synonyms,
                       by = c("therapeutic_context" = "alias")) %>%
      dplyr::filter(!is.na(nci_concept_display_name)) %>%
      dplyr::select(evidence_id, molecule_chembl_id) %>%
      dplyr::distinct() %>%
      dplyr::group_by(evidence_id) %>%
      dplyr::summarise(molecule_chembl_id = paste(
        unique(sort(molecule_chembl_id)), collapse="|"),
        .groups = "drop"
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(molecule_chembl_id = stringr::str_replace(
        molecule_chembl_id,"^NA\\||\\|NA$",""
      ))
  )

  cgi_biomarkers_all <- cgi_biomarkers_all %>%
    dplyr::left_join(therapeutic_contexts, by = "evidence_id") %>%
    dplyr::mutate(evidence_url = "https://www.cancergenomeinterpreter.org/biomarkers") %>%
    dplyr::mutate(biomarker_source_db = "cgi") %>%
    dplyr::mutate(evidence_description = NA,
                  biomarker_entity = T)

  biomarkers_all <- cgi_biomarkers_all %>%
    dplyr::select(evidence_id,
                  evidence_url,
                  variant_id,
                  variant_type,
                  variant,
                  alteration_type,
                  symbol,
                  therapeutic_context,
                  disease_ontology_id,
                  cancer_type,
                  molecule_chembl_id,
                  dplyr::everything())

  biomarkers_all <- map_biomarker_phenotypes(biomarkers_all) %>%
    dplyr::select(-do_cancer_slim) %>%

  return(biomarkers_all)
}

map_biomarker_phenotypes <- function(biomarkers_all){

  omap <- oncoPhenoMap::oncotree_expanded_full %>%
    dplyr::select(do_id, efo_id, efo_name, do_name,
                  cui, cui_name, primary_site) %>%
    dplyr::filter(!is.na(do_id) & !is.na(cui)) %>%
    dplyr::rename(disease_ontology_id = do_id) %>%
    dplyr::distinct()

  umls_terms <- oncoPhenoMap::auxiliary_maps$umls$concept %>%
    dplyr::filter(main_term == T) %>%
    dplyr::select(cui, cui_name)

  biomarkers_all_phenotypes <- biomarkers_all %>%
    dplyr::left_join(omap, by = "disease_ontology_id")

  missing_do_phenotypes <- biomarkers_all_phenotypes %>%
    dplyr::filter(is.na(cui)) %>%
    dplyr::select(-c(cui, cui_name, primary_site,
                     efo_id, do_name, efo_name)) %>%
    dplyr::inner_join(oncoPhenoMap::auxiliary_maps$do,
                      by = c("disease_ontology_id" = "do_id")) %>%
    dplyr::left_join(umls_terms, by = "cui") %>%
    dplyr::left_join(dplyr::select(
      oncoPhenoMap::auxiliary_maps$efo$efo2xref,
      cui, efo_id, efo_name
    ), by = "cui") %>%
    dplyr::distinct() %>%
    dplyr::mutate(primary_site = dplyr::case_when(
      do_name == "cancer" ~ as.character(NA),
      stringr::str_detect(do_name,"brain") ~ "CNS/Brain",
      stringr::str_detect(do_name,"breast") ~ "Breast",
      stringr::str_detect(do_name,"colon|rectum") ~ "Colon/Rectum",
      stringr::str_detect(do_name,"gastric|stomach|esophag") ~ "Esophagus/Stomach",
      stringr::str_detect(do_name,"prostate") ~ "Prostate",
      stringr::str_detect(do_name,"pancrea") ~ "Pancreas",
      stringr::str_detect(do_name,"lung") ~ "Lung",
      stringr::str_detect(do_name,"ovary|ovarian") ~ "Ovary",


      TRUE ~ as.character(NA)
    ))

  biomarkers_phenotype_mapped <-
    dplyr::filter(biomarkers_all_phenotypes, !is.na(cui)) %>%
    dplyr::bind_rows(missing_do_phenotypes)

  return(biomarkers_phenotype_mapped)
}


load_pmkb_biomarkers <- function(){


  pmkb_biomarkers <- as.data.frame(
    readr::read_delim(file="data-raw/biomarkers/pmkb/PMKB_Interpretations_Complete_20211220-0312.csv",
                    show_col_types = F, delim = ",")
  ) %>%
    janitor::clean_names() %>%
    dplyr::mutate(gene = dplyr::if_else(
      gene == "SEPT14", "SEPTIN14", as.character(gene)
    )) %>%
    tidyr::separate_rows(variant_s, sep="\\|") %>%
    dplyr::mutate(variant_s = dplyr::if_else(
      is.na(variant_s),"Undefined",as.character(variant_s)
    )) %>%
    tidyr::separate_rows(variant_s, sep=", ") %>%
    dplyr::filter(!stringr::str_detect(variant_s,"[0-9]{1,}( any)?$")) %>%
    dplyr::mutate(variant_s = dplyr::if_else(
      is.na(variant_s),"Undefined",as.character(variant_s)
    )) %>%
    dplyr::mutate(variant_type = dplyr::case_when(
      stringr::str_detect(variant_s,"codon\\(s\\)|missense| [A-Z][0-9]{1,}[A-Z]") &
        !stringr::str_detect(variant_s, "frameshift$") ~ "missense_variant",
      stringr::str_detect(variant_s,"frameshift|fs$") ~ "frameshift_variant",
      stringr::str_detect(variant_s,"insertion|ins[A-Z]") ~ "inframe_insertion",
      stringr::str_detect(variant_s,"deletion|del$") ~ "inframe_deletion",
      stringr::str_detect(variant_s,"nonsense|\\*$") ~ "stop_gained",
      stringr::str_detect(variant_s,"any mutation") ~ "protein_altering_variant",
      stringr::str_detect(variant_s,"copy number loss") ~ "transcript_ablation",
      stringr::str_detect(variant_s,"copy number gain") ~ "transcript_amplification",
      TRUE ~ as.character(NA)
    )) %>%
    dplyr::mutate(variant = dplyr::case_when(
      stringr::str_detect(variant_s," any ") ~ "MUTATION",
      variant_type == "transcript_amplification" ~ "AMPLIFICATION",
      variant_type == "transcript_ablation" ~ "DELETION",
      TRUE ~ as.character(variant_s)
    )) %>%
    dplyr::mutate(variant = dplyr::if_else(
      stringr::str_detect(variant,"exon\\(s\\) [0-9]{1,}(-[0-9]{1,})? deletion"),
      paste0("EXON ",
             stringr::str_replace_all(variant,"\\S+ exon\\(s\\) | deletion",""),
             " DELETION"),
      as.character(variant)
    )) %>%
    dplyr::mutate(variant = dplyr::if_else(
      stringr::str_detect(variant,"exon\\(s\\) [0-9]{1,}(-[0-9]{1,})? insertion"),
      paste0("EXON ",
             stringr::str_replace_all(variant,"\\S+ exon\\(s\\) | insertion",""),
             " INSERTION"),
      as.character(variant)
    )) %>%
    dplyr::mutate(variant = dplyr::if_else(
      stringr::str_detect(variant,"exon\\(s\\) [0-9]{1,}(-[0-9]{1,})? frameshift"),
      paste0("EXON ",
             stringr::str_replace_all(variant,"\\S+ exon\\(s\\) | frameshift",""),
             " FRAMESHIFT"),
      as.character(variant)
    )) %>%
    dplyr::mutate(variant = dplyr::if_else(
      stringr::str_detect(variant,"exon\\(s\\) [0-9]{1,}(-[0-9]{1,})? (missense|nonsense)"),
      paste0("EXON ",
             stringr::str_replace_all(variant,"\\S+ exon\\(s\\) | (missense|nonsense)",""),
             " MUTATION"),
      as.character(variant)
    )) %>%
    dplyr::mutate(alteration_type = dplyr::case_when(
      variant == "AMPLIFICATION" | variant == "DELETION" ~ "CNA",
      stringr::str_detect(variant,"codon\\(s\\)") ~ "CODON",
      TRUE ~ as.character("MUT")
    )) %>%
    dplyr::mutate(gene = dplyr::if_else(
      gene == "H3F3A","H3-3A",as.character(gene)
    )) %>%
    dplyr::mutate(variant = stringr::str_replace(
      variant,"\\*$","X"
    )) %>%
    dplyr::filter(variant != "Undefined" & !stringr::str_detect(variant,"rearrangement")) %>%
    dplyr::rename(tissue = tissue_type_s,
                  evidence_description = interpretations,
                  evidence_url = pmkb_url,
                  citation_id = citations,
                  cancer_type = tumor_type_s) %>%
    dplyr::mutate(evidence_id = seq(1:nrow(.))) %>%
    dplyr::mutate(evidence_id = paste0("PMKB_",evidence_id))


  codon_mapping <- readr::read_delim("data-raw/biomarkers/pmkb/codon_mapping.tsv",
                                   show_col_types = F, delim = ";") %>%
    dplyr::mutate(variant2 = paste0(amino_acid,position)) %>%
    dplyr::select(variant, variant2) %>%
    dplyr::distinct()

  codon_hits <- pmkb_biomarkers  %>%
    dplyr::filter(stringr::str_detect(variant," codon\\(s\\)")) %>%
    dplyr::left_join(codon_mapping, by = "variant") %>%
    dplyr::filter(!is.na(variant2)) %>%
    dplyr::mutate(variant2 = paste(gene, variant2, sep=" ")) %>%
    dplyr::mutate(variant = variant2) %>%
    dplyr::select(-variant2)

  pmkb_biomarkers_final <-
    dplyr::filter(pmkb_biomarkers, alteration_type != "CODON") %>%
    dplyr::bind_rows(codon_hits) %>%
    #tidyr::separate_rows(tissue, sep="\\|") %>%
    tidyr::separate_rows(cancer_type, sep="\\|") %>%
    dplyr::mutate(evidence_level = dplyr::case_when(
      tier == 1 ~ "A",
      tier == 2 ~ "C",
      tier == 3 ~ as.character(NA),
      TRUE ~ as.character(NA))) %>%
    dplyr::mutate(therapeutic_context = NA,
                  clinical_significance = NA,
                  #citation_id = NA,
                  evidence_type = NA,
                  variant_origin = "Somatic")


  umls_aliases <- oncoPhenoMap::oncotree_expanded_full %>%
    dplyr::select(cui) %>%
    dplyr::distinct() %>%
    dplyr::left_join(oncoPhenoMap::auxiliary_maps$umls$concept, by = "cui") %>%
    dplyr::select(cui, cui_name) %>%
    dplyr::mutate(cui_name_lc = tolower(cui_name)) %>%
    dplyr::select(cui_name_lc, cui) %>%
    dplyr::distinct()

  pmkb_cancer_spec <- pmkb_biomarkers_final %>%
    dplyr::filter(!(cancer_type == "Adenocarcinoma" |
                      cancer_type == "Carcinoma" |
                      cancer_type == "Squamous Cell Carcinoma" |
                      cancer_type == "Basal Cell Carcinoma")) %>%
    dplyr::mutate(cancer_type2 = tolower(cancer_type)) %>%
    dplyr::left_join(umls_aliases, by = c("cancer_type2" = "cui_name_lc")) %>%
    #dplyr::select(-cancer_type2) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dplyr::select(
      oncoPhenoMap::oncotree_expanded_full,
      cui, cui_name, efo_id, efo_name, do_name,
      do_id, primary_site, do_cancer_slim),
      by = "cui") %>%
    dplyr::distinct() %>%
    dplyr::rename(disease_ontology_id = do_id)


  pmkb_cancer_nos <- pmkb_biomarkers_final %>%
    dplyr::filter(cancer_type == "Adenocarcinoma" |
                    cancer_type == "Carcinoma" |
                    cancer_type == "Squamous Cell Carcinoma" |
                    cancer_type == "Basal Cell Carcinoma") %>%
    tidyr::separate_rows(tissue, sep="\\|") %>%
    dplyr::mutate(tissue = stringr::str_trim(tissue)) %>%
    #dplyr::select(cancer_type, tissue) %>%
    dplyr::mutate(cancer_type2 = dplyr::case_when(
      tissue == "Colon" | tissue == "Rectum" ~ "Colorectal cancer",
      tissue == "Breast" ~ "Breast cancer",
      tissue == "Pancreas" ~ "Pancreatic adenocarcinoma",
      tissue == "Lung" ~ "Lung cancer",
      tissue == "Larynx" | tissue == "Oral Cavity" ~ "Head and neck cancer",
      # tissue == "Oral Cavity"˙~ "Head and neck cancer",
      tissue == "Stomach" ~ "Stomach cancer",
      tissue == "Prostate" ~ "Prostate cancer",
      tissue == "Bladder" ~ "Bladder cancer",
      tissue == "Ovary" ~ "Ovarian cancer",
      tissue == "Thyroid" ~ "Thyroid cancer",
      tissue == "Gallbladder" ~ "Ovarian cancer",
      tissue == "Endometrium" ~ "Endometrial cancer",
      tissue == "Gall Bladder" ~ "Biliary tract cancer",
      tissue == "Esophagus" ~ "Esophageal cancer",
      tissue == "Kidney" ~ "Renal cell carcinoma",
      tissue == "Liver" ~ "Hepatocellular carcinoma",
      tissue == "Testis" ~ "Testicular cancer",
      tissue == "Soft Tissue" ~ "Sarcoma",
      tissue == "Brain" ~ "Brain cancer"
    )) %>%
    #dplyr::mutate(cancer_type2 = tolower(cancer_type2)) %>%
    dplyr::filter(!is.na(cancer_type2)) %>%
    dplyr::mutate(cancer_type2 = tolower(cancer_type2)) %>%
    dplyr::select(-c(tissue)) %>%
    dplyr::distinct() %>%
    dplyr::left_join(umls_aliases,
                     by = c("cancer_type2" = "cui_name_lc")) %>%
    #dplyr::select(-cancer_type2) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dplyr::select(
      oncoPhenoMap::oncotree_expanded_full,
      cui, cui_name, efo_id, efo_name, do_name,
      do_id, primary_site, do_cancer_slim),
      by = "cui") %>%
    dplyr::distinct() %>%
    dplyr::rename(disease_ontology_id = do_id) %>%
    dplyr::anti_join(pmkb_cancer_spec, by = c("evidence_id","primary_site"))

  pmkb_phenotype_mapped <-
    pmkb_cancer_nos %>%
    dplyr::bind_rows(pmkb_cancer_spec) %>%
    #dplyr::filter(!stringr::str_detect())
    dplyr::arrange(evidence_id) %>%
    dplyr::select(-c(cancer_type2,
                     variant_s, tissue, do_cancer_slim)) %>%
    dplyr::filter(!(is.na(cui) & is.na(primary_site))) %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::mutate(variant = dplyr::if_else(
      !stringr::str_detect(variant,"EXON|DELETION|FRAMESHIFT|MUTATION|INSERTION") &
        !stringr::str_detect(variant, "^[0-9]{2,} "),
      stringr::str_replace(variant, "^\\S+ ",""),
      as.character(variant)
    ))


  variantSummary <- dplyr::select(pmkb_phenotype_mapped,
                                  symbol, variant_type, variant,
                                  alteration_type, evidence_id) %>%
    dplyr::filter(alteration_type != "CNA") %>%
    dplyr::distinct()

  expanded_variant_set <- data.frame()

  for(i in unique(variantSummary$evidence_id)){
    vrows <- variantSummary[variantSummary$evidence_id == i,]

    variant_entries <- unique(vrows$variant)
    gene <- unique(vrows$symbol)
    variant_type <- unique(vrows$variant_type)
    alteration_type <- unique(vrows$alteration_type)


    if(stringr::str_detect(alteration_type,"MUT|CODON")){

      all_aliases <- c()

      for(e in variant_entries){
        clean_aliases <- expand_hgvs_terms(e)
        for(c in clean_aliases){
          all_aliases <- unique(c(all_aliases, c))
        }
      }


      for(variant_alias in all_aliases){
        df <- data.frame(
          'variant' = variant_alias,
          'evidence_id' = i,
          'symbol' = gene,
          'alteration_type' = alteration_type,
          'variant_type' = variant_type,
          stringsAsFactors = F)
        expanded_variant_set <- dplyr::bind_rows(
          expanded_variant_set, df
        )
      }
    }

  }

  expanded_variant_set <- expanded_variant_set %>%
    dplyr::mutate(
      alteration_type =
        dplyr::if_else(
          stringr::str_detect(variant,"^(([A-Z]{1}[0-9]{1,})|([A-Z]{1}[a-z]{2}[0-9]{1,}))$"),
          as.character("CODON"),
          as.character(alteration_type)
        )) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"FS"),
        stringr::str_replace(variant,"FS","fs"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"DEL") &
          !stringr::str_detect(variant, "DELETION"),
        stringr::str_replace(variant,"DEL","del"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"INS") &
          !stringr::str_detect(variant,"INSERTION"),
        stringr::str_replace(variant,"INS","ins"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"TER|ter") &
          !stringr::str_detect(tolower(variant),"alternative|internal|promoter|deleterious"),
        stringr::str_replace(variant,"ter|TER","Ter"),
        as.character(variant)
      )
    ) %>%
    dplyr::mutate(
      variant = dplyr::if_else(
        stringr::str_detect(variant,"DUP") &
          !stringr::str_detect(variant,"DUPLICATION"),
        stringr::str_replace(variant,"DUP","dup"),
        as.character(variant)
      )
    )

  for(aa in c("PRO","ILE","HIS","GLY","GLU","LEU","THR",
              "TYR","SER","ASN","ASP","LYS","PHE")){
    aa_three_lc <- stringr::str_to_title(aa)
    if(nrow(expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]) > 0){
      expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]$variant <-
        stringr::str_replace(expanded_variant_set[!is.na(expanded_variant_set$variant) & stringr::str_detect(expanded_variant_set$variant,aa),]$variant,
                             aa,
                             aa_three_lc)
    }
  }



  pmkb_phenotype_mapped_expanded <- pmkb_phenotype_mapped %>%
    dplyr::filter(alteration_type != "CNA") %>%
    dplyr::select(-c(variant, alteration_type, variant_type)) %>%
    dplyr::left_join(expanded_variant_set,
                     by = c("symbol", "evidence_id")) %>%
    dplyr::distinct() %>%
    dplyr::bind_rows(
      dplyr::filter(pmkb_phenotype_mapped,
                    alteration_type == "CNA")
    ) %>%
    dplyr::filter(!(is.na(primary_site) &
                      stringr::str_detect(
                        tolower(cui_name),"myelodysplastic|leukemia"
                      ))) %>%
    dplyr::arrange(evidence_id) %>%
    dplyr::select(-tier) %>%
    dplyr::mutate(variant_id = NA,
                  evidence_direction = NA,
                  evidence_type = NA,
                  biomarker_source_db = "pmkb",
                  biomarker_entity = T,
                  molecule_chembl_id = NA) %>%
    dplyr::select(evidence_id,
                  evidence_url,
                  variant_id,
                  variant_type,
                  variant,
                  alteration_type,
                  symbol,
                  therapeutic_context,
                  disease_ontology_id,
                  cancer_type,
                  molecule_chembl_id,
                  dplyr::everything()
    )

  return(pmkb_phenotype_mapped_expanded)

#
#
#
#
#
#
#   pmkb_biomarkers <- as.data.frame(
#     openxlsx::read.xlsx(
#       "data-raw/biomarkers/pmkb/PMKB_Interpretations_Complete_20200405-0319.xlsx",
#       sheet = 1, cols = c(1,2,3,4,5), colNames = T) %>%
#       janitor::clean_names() %>%
#       dplyr::filter(tier == 1 | tier == 2) %>%
#       dplyr::rename(symbol = gene, cancer_type = tumor_type_s, evidence_level = tier) %>%
#       dplyr::mutate(evidence_level = dplyr::if_else(evidence_level == 1,"A","C")) %>%
#       dplyr::mutate(cancer_type = stringr::str_replace_all(cancer_type,", Unspecified",";Unspecified")) %>%
#       dplyr::mutate(
#         biomarker_alteration_type =
#           dplyr::case_when(
#             stringr::str_detect(variant_s, " [A-Z]{1}[0-9]{1,}[A-Z]|frameshift|mutation|codon|nonsense|missense|exon") ~ "MUT",
#             stringr::str_detect(variant_s, "copy number loss|copy number gain") ~ "CNA",
#             stringr::str_detect(variant_s, "rearrangement|translocation") ~ "TRANSLOCATION_FUSION",
#             TRUE ~ as.character("Undefined"))) %>%
#       tidyr::separate_rows(cancer_type,sep=", ") %>%
#       dplyr::group_by(symbol) %>%
#       dplyr::summarise(biomarker_cancer_type = paste(unique(cancer_type),collapse="|"),
#                        biomarker_alteration_type = paste(unique(biomarker_alteration_type),collapse="|"),
#                        biomarker_evidence_level = paste(unique(evidence_level),collapse="|"),
#                        .groups = "drop") %>%
#       dplyr::distinct() %>%
#       dplyr::mutate(biomarker_source_db = "PMKB",
#                     biomarker_clinical_significance = as.character(NA),
#                     biomarker_disease_ontology_id = as.character(NA))
#   )

  #return(pmkb_biomarkers)
}


load_custom_fusion_db <- function(){

  other_translocations <- as.data.frame(read.table(
    file="data-raw/biomarkers/translocation_patterns.tsv", header = T,
    sep = "\t", stringsAsFactors = F) %>%
      dplyr::select(ngram_pattern, biomarker) %>%
      dplyr::rename(alias = ngram_pattern,
                    variant = biomarker) %>%
      dplyr::mutate(biomarker_entity = T,
                    symbol = variant,
                    alias = toupper(alias)) %>%
      tidyr::separate_rows(alias, sep="\\|") %>%
      dplyr::mutate(biomarker_source_db = "fusion_customdb") %>%
      dplyr::distinct()
  )

  other_translocations_lc <- as.data.frame(read.table(
    file="data-raw/biomarkers/translocation_patterns.tsv", header = T,
    sep = "\t", stringsAsFactors = F) %>%
      dplyr::select(ngram_pattern, biomarker) %>%
      dplyr::rename(alias = ngram_pattern,
                    variant = biomarker) %>%
      dplyr::mutate(biomarker_entity = T,
                    symbol = variant) %>%
      tidyr::separate_rows(alias, sep="\\|") %>%
      dplyr::mutate(biomarker_source_db = "fusion_customdb") %>%
      dplyr::distinct()
  )

  custom_translocations <- other_translocations %>%
    dplyr::bind_rows(other_translocations_lc) %>%
    dplyr::distinct() %>%
    dplyr::rename(variant_id = variant,
                  variant = alias) %>%
    dplyr::mutate(variant_type = "transcript_fusion",
                  alteration_type = "TRANSLOCATION_FUSION")

}

load_mitelman_db <- function(){


  phenotypes_morph <- read.table(
    file = "data-raw/biomarkers/mitelmandb/KODER.TXT.DATA",
    sep = "\t", stringsAsFactors = F, header = T, fill = T) %>%
    dplyr::filter(KodTyp == "MORPH" & nchar(Kod) > 0)

  morph <- data.frame('Morph' = as.integer(phenotypes_morph$Kod),
                      'cancer_type' = phenotypes_morph$Benamning,
                      stringsAsFactors = F )

  pmid_data <- read.table(
    file = "data-raw/biomarkers/mitelmandb/REF.TXT.DATA",
    sep = "\t", stringsAsFactors = F, header = T, fill = T,
    quote = "") %>%
    dplyr::select(RefNo, Pubmed) %>%
    dplyr::filter(nchar(Pubmed) > 0) %>%
    dplyr::mutate(RefNo = as.integer(RefNo)) %>%
    dplyr::distinct() %>%
    dplyr::rename(citation_id = Pubmed)

  mbca_data <- read.table(
    file = "data-raw/biomarkers/mitelmandb/MBCA.TXT.DATA",
    sep = "\t", stringsAsFactors = F, header = T) %>%
    dplyr::filter(stringr::str_detect(GeneShort,"::")) %>%
    dplyr::mutate(
      variant = stringr::str_replace_all(GeneShort,"::","-")) %>%
    dplyr::rename(karyotype = KaryShort) %>%
    dplyr::select(MolClin, RefNo, InvNo, Morph, Topo,
                  variant, karyotype) %>%
    tidyr::separate_rows(variant, sep=",") %>%
    dplyr::left_join(pmid_data, by = "RefNo") %>%
    dplyr::left_join(morph, by = "Morph") %>%
    dplyr::filter(stringr::str_detect(variant,"-")) %>%
    dplyr::filter(!stringr::str_detect(variant,"(-|\\+)$")) %>%
    dplyr::mutate(evidence_id = seq(1:nrow(.))) %>%
    dplyr::mutate(evidence_id = paste0("MITDB_",evidence_id)) %>%
    dplyr::filter(stringr::str_count(variant,",") ==
                    stringr::str_count(karyotype,",")) %>%
    tidyr::separate_rows(c(variant,karyotype),
                         sep = ",") %>%
    dplyr::filter(!stringr::str_detect(variant,"\\+")) %>%
    dplyr::select(variant, karyotype, cancer_type,
                  citation_id, evidence_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(variant_id = paste(
      variant, karyotype, sep = " - "
    ))

  evidence_id <- nrow(mbca_data) + 1
  #i <- 1
  #all_translocation_variants <- data.frame()

  mbca_data$alias <- paste(
    mbca_data$variant,
    paste(
      stringr::str_split_fixed(mbca_data$variant,"-",2)[,2],
      stringr::str_split_fixed(mbca_data$variant,"-",2)[,1],
      sep="-"
    ),
    paste(
      stringr::str_split_fixed(mbca_data$variant,"-",2)[,1],
      stringr::str_split_fixed(mbca_data$variant,"-",2)[,2],
      sep="/"
    ),
    mbca_data$karyotype,
    sep="@@@"
  )

  bcr_abl_custom <-
    data.frame('variant_id' = "BCR-ABL1 - t(9;22)(q34;q11)",
               'variant' = 'BCR-ABL',
               'evidence_url' = "https://mitelmandatabase.isb-cgc.org/",
               'evidence_id' = paste0('MITDB_',evidence_id),
               'evidence_description' = NA,
               'cancer_type' = NA,
               'citation_id' = NA,
               'biomarker_entity' = T,
               'alteration_type' = 'TRANSLOCATION_FUSION',
               'biomarker_source_db' = 'mitelmandb',
               'symbol' = 'BCR-ABL1',
               stringsAsFactors = F)

  mitelman_db <- as.data.frame(
    mbca_data %>%
      dplyr::select(alias, variant_id, cancer_type,
                    evidence_id, citation_id) %>%
      tidyr::separate_rows(alias, sep="@@@") %>%
      dplyr::rename(variant = alias) %>%
      dplyr::filter(
        !stringr::str_detect(
          variant,"RARA-PML|ABL1-BCR|ERG-TMPRSS2")) %>%
      dplyr::mutate(biomarker_entity = T) %>%
      dplyr::mutate(alteration_type = "TRANSLOCATION_FUSION",
                    variant_origin = "somatic",
                    molecule_chembl_id = NA,
                    evidence_level = NA,
                    clinical_significance = NA,
                    evidence_type = NA,
                    evidence_direction = NA,
                    variant_type = "transcript_fusion",
                    therapeutic_context = NA,
                    evidence_url = "https://mitelmandatabase.isb-cgc.org/",
                    biomarker_source_db = "mitelmandb") %>%
      dplyr::bind_rows(
        bcr_abl_custom
      ) %>%
      dplyr::arrange(evidence_id) %>%
      dplyr::select(evidence_id,
                    evidence_url,
                    variant_id,
                    variant_type,
                    variant,
                    alteration_type,
                    symbol,
                    dplyr::everything())

  )

  mitelman_db$symbol <- stringr::str_split_fixed(
    mitelman_db$variant_id, " - ", 2)[,1]

  umls_aliases <- oncoPhenoMap::auxiliary_maps$umls$concept %>%
    dplyr::select(cui, cui_name) %>%
    dplyr::mutate(cui_name = tolower(cui_name)) %>%
    dplyr::distinct()

  mitelman_db <- mitelman_db %>%
    dplyr::mutate(cancer_type = dplyr::if_else(
      cancer_type == "Astrocytoma, grade I-II",
      "low-grade astrocytoma",
      as.character(cancer_type)
    )) %>%
    dplyr::mutate(cancer_type = dplyr::if_else(
      cancer_type == "Astrocytoma, grade III-IV/Glioblastoma",
      "high-grade astrocytoma",
      as.character(cancer_type)
    )) %>%
    dplyr::mutate(cancer_type_lc = tolower(cancer_type)) %>%
    dplyr::left_join(umls_aliases, by = c("cancer_type_lc" = "cui_name")) %>%
    dplyr::filter(!stringr::str_detect(cancer_type_lc,"nonneoplastic"))

  missing <- mitelman_db %>%
    dplyr::filter(is.na(cui) & !is.na(cancer_type)) %>%
    dplyr::select(-cui) %>%
    dplyr::mutate(cancer_type_lc = stringr::str_replace(
      cancer_type_lc, "leukemia/lymphoblastic ","leukemia/"
    )) %>%
    dplyr::mutate(cancer_type_lc = stringr::str_replace(
      cancer_type_lc,
      ", (nos|special type|, cutaneous type|dedifferentiated|aberrant translocation)$","")
    ) %>%
    dplyr::left_join(umls_aliases,
                     by = c("cancer_type_lc" = "cui_name"))

  mitelman_db_final <-
    dplyr::bind_rows(
      dplyr::filter(mitelman_db, !is.na(cui) & !is.na(cancer_type)),
      missing
    ) %>%
    dplyr::select(-cancer_type_lc)


  mitelman <-
    mitelman_db_final %>%
    dplyr::left_join(oncoPhenoMap::auxiliary_maps$do, by = "cui") %>%
    dplyr::select(-do_cancer_slim) %>%
    dplyr::filter(is.na(cui) | (!is.na(cui) & !is.na(do_id))) %>%
    dplyr::filter(do_name != "malignant adenoma") %>%
    dplyr::left_join(oncoPhenoMap::auxiliary_maps$efo$efo2xref, by = "cui") %>%
    dplyr::filter(is.na(cui) | (!is.na(cui) & !is.na(do_id))) %>%
    dplyr::filter(do_name != "malignant adenoma") %>%
    dplyr::rename(disease_ontology_id = do_id) %>%
    dplyr::select(-xref_source) %>%
    dplyr::distinct()  %>%
    dplyr::filter(efo_id != "MONDO:0015667")


  ## REMAIN match against oncoPhenoMap


  return(mitelman)
}
