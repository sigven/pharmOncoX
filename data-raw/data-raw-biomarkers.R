source('data-raw/biomarker_utilities.R')
source('data-raw/drug_utilities.R')
library(magrittr)

datestamp <- '20220602'

civic_variant_summary <- paste0("data-raw/biomarkers/civic/variant_summary_",
                                datestamp,
                                ".tsv")
civic_clinical_evidence <- paste0("data-raw/biomarkers/civic/clinical_evidence_summary_",
                                  datestamp,
                                  ".tsv")

if(!file.exists(civic_variant_summary)){
  download.file(url = "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv",
                destfile = civic_variant_summary)
}
if(!file.exists(civic_clinical_evidence)){
  download.file(url = "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv",
                destfile = civic_clinical_evidence)
}

gene_info <- get_gene_info_ncbi(
  path_data_raw = file.path(here::here(), "data-raw"),
  update = T)

raw_biomarkers <- list()
raw_biomarkers[['civic']] <-
  load_civic_biomarkers(
    datestamp = datestamp)
raw_biomarkers[['cgi']] <- load_cgi_biomarkers()
raw_biomarkers[['mitelmandb']] <- load_mitelman_db()
raw_biomarkers[['pmkb']] <- load_pmkb_biomarkers()
custom_fusiondb <- load_custom_fusion_db() |>
  dplyr::anti_join(raw_biomarkers[["mitelmandb"]], by = "variant")

biomarkers_all <- do.call(rbind, raw_biomarkers)
rownames(biomarkers_all) <- NULL
biomarkers_all <- biomarkers_all |>
  dplyr::bind_rows(custom_fusiondb) |>
  dplyr::mutate(therapeutic_context =
                  stringi::stri_enc_toascii(therapeutic_context)) |>
  dplyr::mutate(cancer_type =
                  stringi::stri_enc_toascii(cancer_type)) |>
  dplyr::mutate(evidence_description =
                  stringi::stri_enc_toascii(evidence_description))

compound_biomarkers <- list()
compound_biomarkers[['curated']] <- biomarkers_all
compound_biomarkers[['version']] <- list()
compound_biomarkers[['version']][['civic']] <- datestamp
compound_biomarkers[['version']][['cgi']] <- '20180117'
compound_biomarkers[['version']][['pmkb']] <- '20200405'
compound_biomarkers[['version']][['mitelman']] <- '20220418'
compound_biomarkers[['version']][['prism_depmap']] <- '19Q4_21Q4'

targets <- list()
targets[['druggable']] <- oncoPharmaDB::oncopharmadb |>
  # dplyr::mutate(target_source_db = dplyr::if_else(
  #   !is.na(opentargets_version),
  #   "NCI_thesaurus&OTPlatform",
  #   "NCI_thesaurus"
  # )) |>
  dplyr::select(target_symbol, target_ensembl_gene_id, target_genename,
                nci_concept_display_name, molecule_chembl_id) |>
  dplyr::rename(symbol = target_symbol,
                ensembl_gene_id = target_ensembl_gene_id,
                genename = target_genename) |>
  dplyr::left_join(dplyr::select(gene_info, symbol, entrezgene), by = "symbol") |>
  dplyr::group_by(symbol, ensembl_gene_id, genename, entrezgene) |>
  dplyr::summarise(nci_concept_display_name = paste(
    sort(unique(nci_concept_display_name)), collapse="|"),
    molecule_chembl_id = paste(sort(unique(molecule_chembl_id)), collapse="|"),
    # target_source_db = paste(sort(unique(target_source_db)), collapse="&"),
    .groups = "drop") |>
  dplyr::distinct()

targets[['actionable']] <- as.data.frame(biomarkers_all |>
  dplyr::filter(biomarker_source_db == "civic" |
                  biomarker_source_db == "cgi" |
                  biomarker_source_db == "pmkb") |>
  dplyr::mutate(source_db = toupper(biomarker_source_db)) |>
  dplyr::group_by(symbol, source_db) |>
  dplyr::summarise(phenotypes = paste(sort(unique(cancer_type)), collapse="|"),
                   clinical_significance = paste(sort(unique(clinical_significance)), collapse="|"),
                   .groups = "drop") |>
  dplyr::select(symbol, phenotypes, clinical_significance, source_db) |>
  dplyr::left_join(dplyr::select(
    gene_info, symbol, ensembl_gene_id,
    name, entrezgene), by = "symbol") |>
  dplyr::rename(genename = name) |>
  dplyr::select(symbol, ensembl_gene_id,
                genename, entrezgene,
                dplyr::everything()) |>
  dplyr::distinct()
)

ccle_drugs <- readr::read_delim(
  file="data-raw/biomarkers/ccle_depmap/primary-screen-replicate-treatment-info.csv",
  delim=",", show_col_types = F) |>
  dplyr::filter(perturbation_type == "experimental_treatment") |>
  dplyr::filter(!is.na(broad_id)) |>
  dplyr::select(broad_id, name, target, moa, indication) |>
  dplyr::rename(compound_id = broad_id,
                drugname = name,
                drugtarget = target,
                drug_moa = moa,
                drug_indication = indication) |>
  dplyr::mutate(drugname_lc = tolower(drugname)) |>
  dplyr::distinct()

ccle_predictive_features <- readr::read_delim(
  file="data-raw/biomarkers/ccle_depmap/Repurposing_primary_predictability_results.csv",
  delim=",", show_col_types = F) |>
  dplyr::mutate(compound_id = stringr::str_replace(gene,"BRD:","")) |>
  dplyr::semi_join(ccle_drugs, by = c("compound_id"))

all_predictive_features <- data.frame()

for(n in 1:nrow(ccle_predictive_features)){
  row <- as.data.frame(ccle_predictive_features[n,])
  compound_id <- stringr::str_replace(row$gene,"BRD:","")
  model <- row$model
  pearson <- row$pearson
  best <- row$best

  if(model != "Core_omics"){
    next
  }

  for(i in 0:9){
    col1 <- paste0("feature",i)
    col2 <- paste0("feature",i,"_importance")
    symbol <- NA
    feature_type <- NA
    feature_value <- NA
    if(stringr::str_detect(row[col1],"(CN|MutHot|MutMis|RNAseq|MutDam)$")){
      feature_properties <- stringr::str_split(
        stringr::str_replace_all(as.character(row[col1]),"\\)|\\(",""),"_")[[1]]
      if(length(feature_properties) == 3){
        symbol <- feature_properties[1]
        entrezgene <- feature_properties[2]
        feature_value <- paste(symbol, entrezgene, sep = ";")
        feature_type <- feature_properties[3]
        if(feature_type == "CN"){
          feature_type = "CopyNumber"
        }
        if(feature_type == "MutMmis"){
          feature_type = "MissenseMutation"
        }
        if(feature_type == "MutDam"){
          feature_type = "DamagingMutation"
        }
        if(feature_type == "MutHot"){
          feature_type = "HotSpotMutation"
        }
      }
    }else{
      if(stringr::str_detect(row[col1],"Lin$")){
        feature_value = stringr::str_replace(row[col1],"_Lin$","")
        feature_type = "Lineage"
      }
      else if(stringr::str_detect(row[col1],"_Fusion$")){
        fusion_candidate = stringr::str_replace(row[col1],"_Fusion$","")
        if(stringr::str_count(fusion_candidate,"_") == 1){
          feature_type = "Fusion"
          feature_value <- fusion_candidate
        }
      }
    }
    importance <- round(as.numeric(row[col2]), digits = 7)

    if(!is.na(feature_type)){
      df <- data.frame(
        'compound_id' = compound_id,
        'model' = model,
        'feature_type' = feature_type,
        'feature_value' = feature_value,
        'feature_importance' = importance,
        'prediction_accuracy' = pearson,
        stringsAsFactors = F
      )

      all_predictive_features <- all_predictive_features |>
        dplyr::bind_rows(df)
    }

  }
}

compound_biomarkers[['invitro_screen']] <- all_predictive_features |>
  dplyr::inner_join(ccle_drugs, by = "compound_id") |>
  dplyr::inner_join(oncoPharmaDB::compound_synonyms,
                    by = c("drugname_lc" = "alias")) |>
  dplyr::select(-c(drugname_lc, drugtarget, drug_moa,
                   drug_indication, drugname)) |>
  dplyr::distinct()


rm(ccle_drugs)
rm(ccle_predictive_features)
rm(all_predictive_features)

usethis::use_data(compound_biomarkers, overwrite = T)
usethis::use_data(targets, overwrite = T)

