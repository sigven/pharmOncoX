
## get metadata from metadata_pharma_oncox.xlsx
metadata <- list()
for(elem in c('compounds','biomarkers')){
  metadata[[elem]] <- as.data.frame(openxlsx::read.xlsx(
    "data-raw/metadata_pharma_oncox.xlsx", sheet = elem, colNames = T) |>
      dplyr::mutate(version = dplyr::if_else(
        is.na(version) &
          abbreviation == "civic",
        as.character(stringr::str_replace_all(Sys.Date(),"-","")),
        as.character(version)
      ))
  )
}

nci_db_release <- metadata$compounds[metadata$compounds$abbreviation == "nci", "version"]
opentargets_version <- metadata$compounds[metadata$compounds$abbreviation == "opentargets", "version"]
package_datestamp <- stringr::str_replace_all(Sys.Date(),"-","")
chembl_pubchem_datestamp <- '20220531'

## set logging layout
lgr::lgr$appenders$console$set_layout(
  lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

suppressPackageStartupMessages(source('data-raw/drug_utilities.R'))
suppressPackageStartupMessages(source('data-raw/biomarker_utilities.R'))


nci_ftp_base <- paste0("https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/",
                       nci_db_release,
                       "_Release/")
path_data_raw <-
  file.path(here::here(), "data-raw")
path_data_tmp_processed <-
  file.path(path_data_raw, "tmp_processed")


# ####--- NCBI gene xrefs----####
gene_basic <- geneOncoX::get_basic(cache_dir = path_data_tmp_processed)
gene_gencode <- geneOncoX::get_gencode(cache_dir = path_data_tmp_processed)

gene_info <- dplyr::bind_rows(
  dplyr::inner_join(
    dplyr::select(gene_basic$records, entrezgene, symbol, name, gene_biotype),
    dplyr::select(gene_gencode$records$grch38, entrezgene, ensembl_gene_id),
    by = c("entrezgene")),
  dplyr::inner_join(
    dplyr::select(gene_basic$records, entrezgene, symbol, name, gene_biotype),
    dplyr::select(gene_gencode$records$grch37, entrezgene, ensembl_gene_id),
    by = c("entrezgene"))) |>
  dplyr::filter(gene_biotype == "protein-coding") |>
  dplyr::distinct() |>
  dplyr::mutate(association_sourceID = "nci_thesaurus_custom",
                target_type = "single_protein") |>
  dplyr::filter(!is.na(ensembl_gene_id)) |>
  dplyr::rename(genename = name,
                target_entrezgene = entrezgene,
                target_ensembl_gene_id = ensembl_gene_id)




nci_thesaurus_files <- list()
nci_thesaurus_files[['flat']] <- paste0("Thesaurus_", nci_db_release,".FLAT.zip")
nci_thesaurus_files[['owl']] <- paste0("Thesaurus_", nci_db_release,".OWL.zip")
nci_thesaurus_files[['inf_owl']] <- paste0("ThesaurusInf_", nci_db_release,".OWL.zip")


####---- DailyMed drug indications----####

# drug_indications_dailymed <-
#   get_dailymed_drug_indications(update = update_dailymed,
#                                 path_data_raw = path_data_raw)

for(elem in c('flat','owl','inf_owl')){
  remote_file <- paste0(nci_ftp_base, nci_thesaurus_files[[elem]])
  local_file <- file.path(path_data_raw,"nci_thesaurus",nci_thesaurus_files[[elem]])
  if(!file.exists(local_file)){
    download.file(url = remote_file, destfile = local_file, quiet = T)
    system(paste0('unzip -d ',file.path(path_data_raw, "nci_thesaurus"), ' -o -u ',local_file))
  }
}

antineo_agents_url <-
  'https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Drug_or_Substance/Antineoplastic_Agent.txt'
antineo_agents_local <-
  file.path(path_data_raw,"nci_thesaurus","Antineoplastic_Agent.txt")
if(!file.exists(antineo_agents_local)){
  download.file(url = antineo_agents_url, destfile = antineo_agents_local, quiet = T)
}


####---- Cancer drugs: NCI + DGIdb -----####

## Get all anticancer drugs, NCI thesaurus + DGIdb
nci_antineo_all <- get_nci_drugs(
  nci_db_release = nci_db_release,
  overwrite = T,
  path_data_raw = path_data_raw,
  path_data_processed = path_data_tmp_processed)

#### -- Open Targets Platform - drugs ---####
## Get all targeted anticancer drugs from Open Targets Platform
ot_drugs <-
  get_opentargets_cancer_drugs(
    path_data_raw = path_data_raw,
    ot_version = opentargets_version)


## Merge information from Open Targets Platform and NCI targeted drugs
## 1) By molecule chembl id
## 2) By name (if molecule chembl id does not provide any cross-ref)
## 3) Remove ambiguous names/ids
##
ot_nci_drugs <- merge_nci_open_targets(
  ot_drugs = ot_drugs,
  nci_antineo_all = nci_antineo_all)


####-- Cancer drugs: NCI custom match----####

drug_df <- map_custom_nci_targets(
  gene_info = gene_info,
  path_data_raw = path_data_raw,
  drug_df = ot_nci_drugs
)

####-- Cancer drug categories ---####
drug_df <- assign_drug_category(
  drug_df = drug_df,
  path_data_raw = path_data_raw
)

drug_index_map <- clean_final_drug_list(
  drug_df = drug_df
)

####--- Cancer drug aliases ----#####
drug_aliases <- expand_drug_aliases(
  drug_index_map = drug_index_map,
  path_data_raw = path_data_raw,
  chembl_pubchem_datestamp = chembl_pubchem_datestamp
)

drug_index_map[['id2alias']] <- drug_aliases
drug_index_map$id2synonym <- NULL

compound_synonyms <- drug_index_map[['id2alias']] |>
  dplyr::mutate(drugname_lc = tolower(alias)) |>
  dplyr::left_join(drug_index_map[['id2name']], by = "drug_id") |>
  dplyr::left_join(
    dplyr::select(drug_index_map[['id2basic']],
                  drug_id, molecule_chembl_id),
    by = "drug_id")

rm(drug_df)


##################### BIOMARKERS ###############################

civic_variant_summary <- paste0("data-raw/biomarkers/civic/variant_summary_",
                                package_datestamp,
                                ".tsv")
civic_clinical_evidence <- paste0("data-raw/biomarkers/civic/clinical_evidence_summary_",
                                  package_datestamp,
                                  ".tsv")

if(!file.exists(civic_variant_summary)){
  download.file(url = "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv",
                destfile = civic_variant_summary)
}
if(!file.exists(civic_clinical_evidence)){
  download.file(url = "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv",
                destfile = civic_clinical_evidence)
}

# gene_info <- get_gene_info_ncbi(
#   path_data_raw = file.path(here::here(), "data-raw"),
#   update = T)

raw_biomarkers <- list()
raw_biomarkers[['civic']] <-
  load_civic_biomarkers(
    compound_synonyms = compound_synonyms,
    datestamp = package_datestamp)
raw_biomarkers[['cgi']] <- load_cgi_biomarkers(
  compound_synonyms = compound_synonyms)
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

rm(raw_biomarkers)
biomarkers_curated <- list()
biomarkers_curated[['records']] <- biomarkers_all
biomarkers_curated[['metadata']] <- metadata$biomarkers[1:4,]
rm(biomarkers_all)
rm(civic_clinical_evidence)
rm(civic_variant_summary)

ccle_drugs <- readr::read_csv(
  file="data-raw/biomarkers/ccle_depmap/primary-screen-replicate-treatment-info.csv.gz",
  show_col_types = F) |>
  dplyr::filter(perturbation_type == "experimental_treatment") |>
  dplyr::filter(!is.na(broad_id)) |>
  dplyr::select(broad_id, name, target,
                moa, indication) |>
  dplyr::rename(compound_id = broad_id,
                drugname = name,
                drugtarget = target,
                drug_moa = moa,
                drug_indication = indication) |>
  dplyr::mutate(drugname_lc = tolower(drugname)) |>
  dplyr::distinct()

ccle_predictive_features <- readr::read_csv(
  file="data-raw/biomarkers/ccle_depmap/Repurposing_primary_predictability_results.csv.gz",
  show_col_types = F) |>
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

rm(col1)
rm(col2)
rm(best)
rm(feature_properties)
rm(feature_type)
rm(feature_value)
rm(compound_id)
rm(importance)
rm(model)
rm(fusion_candidate)
rm(pearson)

biomarkers_invitro <- list()
biomarkers_invitro[['metadata']] <- metadata$biomarkers[5,]
biomarkers_invitro[['records']] <- all_predictive_features |>
  dplyr::inner_join(ccle_drugs, by = "compound_id") |>
  dplyr::inner_join(compound_synonyms,
                    by = "drugname_lc") |>
  dplyr::select(-c(drugname_lc, drugtarget, drug_moa,
                   drug_indication, drugname)) |>
  dplyr::distinct() |>
  dplyr::select(-alias) |>
  dplyr::distinct()


rm(ccle_drugs)
rm(ccle_predictive_features)
rm(all_predictive_features)
rm(compound_synonyms)


## upload to Google Drive
version_minor_bumped <- paste0(
  "0.",
  as.character(as.integer(substr(as.character(packageVersion("pharmaOncoX")),3,3)) + 1),
  ".0")

gd_records <- list()
db <- list()
db[['biomarkers_curated']] <- biomarkers_curated
db[['biomarkers_invitro']] <- biomarkers_invitro
db[['drug_map_name']] <- list()
db[['drug_map_name']][['records']] <- drug_index_map[['id2name']]
db[['drug_map_target']] <- list()
db[['drug_map_target']][['records']] <- drug_index_map[['id2target']]
db[['drug_map_indication']] <- list()
db[['drug_map_indication']][['records']] <- drug_index_map[['id2indication']]
db[['drug_map_basic']] <- list()
db[['drug_map_basic']][['records']] <- drug_index_map[['id2basic']]
db[['drug_map_alias']] <- list()
db[['drug_map_alias']][['records']] <- drug_index_map[['id2alias']]

googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))

for(elem in c('biomarkers_curated',
              'biomarkers_invitro',
              'drug_map_name',
              'drug_map_target',
              'drug_map_indication',
              'drug_map_basic',
              'drug_map_alias')){

  if(elem != "biomarkers_curated" & elem != "biomarkers_invitro"){
    db[[elem]][['metadata']] <- metadata$compounds
  }
  saveRDS(db[[elem]],
          file=paste0("data-raw/gd_local/",elem,"_v", version_minor_bumped,".rds"))

  (gd_records[[elem]] <- googledrive::drive_upload(
    paste0("data-raw/gd_local/", elem, "_v", version_minor_bumped,".rds"),
    paste0("pharmaOncoX/", elem, "_v", version_minor_bumped,".rds")
  ))

}

db_id_ref <- dplyr::bind_rows(
  dplyr::select(as.data.frame(gd_records$biomarkers_curated), name, id),
  dplyr::select(as.data.frame(gd_records$biomarkers_invitro), name, id),
  dplyr::select(as.data.frame(gd_records$drug_map_name), name, id),
  dplyr::select(as.data.frame(gd_records$drug_map_target), name, id),
  dplyr::select(as.data.frame(gd_records$drug_map_indication), name, id),
  dplyr::select(as.data.frame(gd_records$drug_map_basic), name, id),
  dplyr::select(as.data.frame(gd_records$drug_map_alias), name, id)) |>
  dplyr::rename(gid = id,
                filename = name) |>
   dplyr::mutate(name =
     stringr::str_replace(filename,"_v\\S+$",""),
     ) |>
  dplyr::mutate(date = Sys.Date(),
                pVersion = version_minor_bumped)
db_id_ref$md5Checksum <- NA

for(elem in c('biomarkers_curated',
              'biomarkers_invitro',
              'drug_map_name',
              'drug_map_target',
              'drug_map_indication',
              'drug_map_basic',
              'drug_map_alias')){

  db_id_ref[db_id_ref$name == elem,]$md5Checksum <-
    gd_records[[elem]]$drive_resource[[1]]$md5Checksum
}

usethis::use_data(db_id_ref, internal = T, overwrite = T)






