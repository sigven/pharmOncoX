suppressPackageStartupMessages(source('data-raw/drug_utilities.R'))
suppressPackageStartupMessages(source('data-raw/biomarker_utilities.R'))

## get metadata from metadata_pharma_oncox.xlsx
metadata <- list()
for (elem in c('compounds','biomarkers')) {
  metadata[[elem]] <- as.data.frame(openxlsx::read.xlsx(
    "data-raw/metadata_pharm_oncox.xlsx", sheet = elem, colNames = T) |>
      dplyr::mutate(version = dplyr::if_else(
        is.na(version) &
          abbreviation == "civic",
        as.character(stringr::str_replace_all(Sys.Date(),"-","")),
        as.character(version)
      ))
  )
}

nci_db_release <- 
  metadata$compounds[metadata$compounds$abbreviation == "nci", "version"]
opentargets_version <- 
  metadata$compounds[metadata$compounds$abbreviation == "opentargets", 
                     "version"]
package_datestamp <- stringr::str_replace_all(Sys.Date(),"-","")
chembl_pubchem_datestamp <- '20220906'

## set logging layout
lgr::lgr$appenders$console$set_layout(
  lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))



nci_ftp_base <- paste0("https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/",
                       nci_db_release,
                       "_Release/")
path_data_raw <-
  file.path(here::here(), "data-raw")
path_data_tmp_processed <-
  file.path(path_data_raw, "tmp_processed")


# ####--- NCBI gene xrefs----####
gene_basic <- geneOncoX::get_basic(
  cache_dir = path_data_tmp_processed)
gene_gencode <- geneOncoX::get_gencode(
  cache_dir = path_data_tmp_processed)
gene_alias <- geneOncoX::get_alias(
  cache_dir = path_data_tmp_processed
)

gene_info <- dplyr::bind_rows(
  dplyr::inner_join(
    dplyr::select(gene_basic$records, entrezgene, 
                  symbol, name, gene_biotype),
    dplyr::select(gene_gencode$records$grch38, 
                  entrezgene, ensembl_gene_id),
    by = c("entrezgene"), multiple = "all"),
  dplyr::inner_join(
    dplyr::select(gene_basic$records, entrezgene, 
                  symbol, name, gene_biotype),
    dplyr::select(gene_gencode$records$grch37, 
                  entrezgene, ensembl_gene_id),
    by = c("entrezgene"), multiple = "all")) |>
  dplyr::filter(gene_biotype == "protein-coding") |>
  dplyr::distinct() |>
  dplyr::mutate(association_sourceID = "nci_thesaurus_custom",
                target_type = "single_protein") |>
  dplyr::filter(!is.na(ensembl_gene_id)) |>
  dplyr::rename(genename = name,
                target_entrezgene = entrezgene,
                target_ensembl_gene_id = ensembl_gene_id)



nci_thesaurus_files <- list()
nci_thesaurus_files[['flat']] <- 
  paste0("Thesaurus_", nci_db_release,".FLAT.zip")
nci_thesaurus_files[['owl']] <- 
  paste0("Thesaurus_", nci_db_release,".OWL.zip")
nci_thesaurus_files[['inf_owl']] <- 
  paste0("ThesaurusInf_", nci_db_release,".OWL.zip")


options(timeout = 50000)
for (elem in c('flat','owl','inf_owl')) {
  remote_file <- paste0(nci_ftp_base, nci_thesaurus_files[[elem]])
  local_file <- file.path(path_data_raw,"nci_thesaurus", 
                          nci_thesaurus_files[[elem]])
  if (!file.exists(local_file)) {
    download.file(url = remote_file, destfile = local_file, quiet = T)
    system(paste0('unzip -d ',
                  file.path(path_data_raw, "nci_thesaurus"), 
                  ' -o -u ',local_file))
  }
}

antineo_agents_url <-
  'https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Drug_or_Substance/Antineoplastic_Agent.txt'
antineo_agents_local <-
  file.path(path_data_raw,"nci_thesaurus","Antineoplastic_Agent.txt")
if (!file.exists(antineo_agents_local)) {
  download.file(url = antineo_agents_url, 
                destfile = antineo_agents_local, quiet = T)
}


####---- Cancer drugs: NCI + DGIdb -----####

## Get all anticancer drugs, NCI thesaurus + DGIdb
nci_antineo_all <- get_nci_drugs(
  nci_db_release = nci_db_release,
  overwrite = F,
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
ot_nci_drugs <- merge_nci_opentargets(
  ot_drugs = ot_drugs,
  path_data_raw = path_data_raw,
  nci_antineo_all = nci_antineo_all)


####-- Cancer drugs: Custom/curated target match ----####
drug_df <- map_curated_targets(
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
  dplyr::mutate(alias_lc = tolower(alias)) |>
  dplyr::left_join(
    drug_index_map[['id2name']], by = "drug_id",
    multiple = "all") |>
  dplyr::left_join(
    dplyr::select(drug_index_map[['id2basic']],
                  drug_id, molecule_chembl_id),
    by = "drug_id", multiple = "all") |>
  dplyr::select(
    drug_id,
    drug_name,
    alias,
    alias_lc,
    molecule_chembl_id
  ) |>
  dplyr::filter(alias != molecule_chembl_id) |>
  dplyr::distinct()
  

rm(drug_df)


##################### BIOMARKERS ###############################

# gene_info <- get_gene_info_ncbi(
#   path_data_raw = file.path(here::here(), "data-raw"),
#   update = T)

raw_biomarkers <- list()
raw_biomarkers[['civic']] <-
  load_civic_biomarkers(
    compound_synonyms = compound_synonyms,
    datestamp = package_datestamp,
    cache_dir = file.path(path_data_raw, "biomarkers"))
raw_biomarkers[['cgi']] <- load_cgi_biomarkers(
  compound_synonyms = compound_synonyms,
  cache_dir = file.path(path_data_raw, "biomarkers"))
raw_biomarkers[['mitelmandb']] <- load_mitelman_db(
  cache_dir = file.path(path_data_raw, "biomarkers"))
# raw_biomarkers[['pmkb']] <- load_pmkb_biomarkers(
#   cache_dir = path_data_raw
# )
raw_biomarkers[['custom_fusions']] <- load_custom_fusion_db()

raw_biomarkers[['custom_fusions']]$variant <- raw_biomarkers[['custom_fusions']]$variant |>
  dplyr::anti_join(raw_biomarkers[["mitelmandb"]][['variant']], by = "variant_alias")

biomarkers <- list()
biomarkers[['data']] <- raw_biomarkers
biomarkers[['metadata']] <- rbind(metadata$biomarkers[1:2,], metadata$biomarkers[4,])
#rm(biomarkers_all)

## upload to Google Drive
version_bump <- paste0(
  substr(as.character(packageVersion("pharmOncoX")),1,4),
  as.character(as.integer(substr(as.character(packageVersion("pharmOncoX")),5,5)) + 1))
  


db <- list()
db[['biomarkers']] <- biomarkers
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

#googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))

gd_records <- list()
db_id_ref <- data.frame()

for(elem in c('biomarkers',
              'drug_map_name',
              'drug_map_target',
              'drug_map_indication',
              'drug_map_basic',
              'drug_map_alias')){

  if(elem != "biomarkers"){
    db[[elem]][['metadata']] <- metadata$compounds
  }
  
  local_rds_fpath <- file.path(
    "data-raw", "gd_local", 
    paste0(elem,"_v", version_bump, ".rds"))
  
  saveRDS(
    db[[elem]],
          file = local_rds_fpath)

  (gd_records[[elem]] <- googledrive::drive_upload(
    media = local_rds_fpath,
    path = paste0("pharmOncoX/", elem, "_v", version_bump,".rds")
  ))

  google_rec_df <-
    dplyr::select(
      as.data.frame(gd_records[[elem]]), name, id) |>
    dplyr::rename(
      gid = id,
      filename = name) |>
    dplyr::mutate(
      name =
        stringr::str_replace(filename,"_v\\S+$",""),
      date = as.character(Sys.Date()),
      pVersion = version_bump) |>
    dplyr::mutate(
      md5Checksum =
        gd_records[[elem]]$drive_resource[[1]]$md5Checksum)

  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)
  
}

usethis::use_data(db_id_ref, internal = T, overwrite = T)






