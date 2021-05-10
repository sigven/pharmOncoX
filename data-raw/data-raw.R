library(magrittr)
pharmamine_datestamp <- '20210510'
nci_db_release <- '21.04d'
chembl_db_release <- 'ChEMBL_28'
opentargets_version <- '2021.04'
uniprot_release <- '2021_02'
dgidb_db_release <- 'v2021_01'

.libPaths("/Library/Frameworks/R.framework/Resources/library")

suppressPackageStartupMessages(source('R/op_db_utilities.R'))

nci_ftp_base <- paste0("https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/",
                       nci_db_release,
                       "_Release/")
path_data_raw <- file.path(here::here(),"data-raw")
path_data_tmp_processed <- file.path(here::here(), "tmp_processed")

## Append targets for recently developed targeted drugs in NCI thesaurus
gene_info <- as.data.frame(
  read.table(gzfile(file.path(path_data_raw,
                              "gene_info","Homo_sapiens.gene_info.gz")),
             sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "") %>%
    dplyr::select(Symbol, GeneID, dbXrefs, description) %>%
    dplyr::rename(symbol = Symbol, target_entrezgene = GeneID, genename = description) %>%
    dplyr::mutate(association_sourceID = "nci_thesaurus_custom", target_type = "single_protein") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(target_ensembl_gene_id = stringr::str_match(dbXrefs,"ENSG[0-9]{1,}")[[1]]) %>%
    dplyr::select(-dbXrefs) %>%
    dplyr::filter(!is.na(target_ensembl_gene_id))
)

nci_thesaurus_files <- list()
nci_thesaurus_files[['flat']] <- paste0("Thesaurus_", nci_db_release,".FLAT.zip")
nci_thesaurus_files[['owl']] <- paste0("Thesaurus_", nci_db_release,".OWL.zip")
nci_thesaurus_files[['inf_owl']] <- paste0("ThesaurusInf_", nci_db_release,".OWL.zip")

for(elem in c('flat','owl','inf_owl')){
  remote_file <- paste0(nci_ftp_base, nci_thesaurus_files[[elem]])
  local_file <- file.path(path_data_raw,"nci_thesaurus",nci_thesaurus_files[[elem]])
  #if(!file.exists(local_file)){
  #download.file(url = remote_file, destfile = local_file, quiet = T)
  system(paste0('unzip -d ',file.path(path_data_raw, "nci_thesaurus"), ' -o -u ',local_file))
  #}
}
antineo_agents_url <-
  'https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Drug_or_Substance/Antineoplastic_Agent.txt'
antineo_agents_local <-
  file.path(path_data_raw,"nci_thesaurus","Antineoplastic_Agent.txt")
download.file(url = antineo_agents_url, destfile = antineo_agents_local, quiet = T)

## Get all anticancer drugs, NCI thesaurus + DGIdb
nci_antineo_all <- get_nci_drugs(
  nci_version = nci_db_release,
  overwrite = T,
  path_data_raw = path_data_raw,
  path_data_processed = path_data_tmp_processed)

system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.owl"))
system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.txt"))
system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.tsv"))

## NCI anticancer drugs (targeted) - with compound identifier (CHEMBL)
nci_antineo_chembl <- nci_antineo_all %>%
  dplyr::select(nci_t,
                nci_concept_definition,
                nci_concept_display_name,
                molecule_chembl_id,
                drug_name_nci,
                nci_concept_synonym_all) %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::distinct()

## NCI anticancer drugs (non-targeted) - lacking compound identifier (CHEMBL)
nci_antineo_nochembl <- nci_antineo_all %>%
  dplyr::filter(is.na(molecule_chembl_id)) %>%
  dplyr::select(nci_t,
                nci_concept_definition,
                nci_concept_display_name,
                drug_name_nci,
                nci_concept_synonym_all) %>%
  dplyr::distinct()

chembl_pubchem_xref <-
  get_chembl_pubchem_compound_xref(
    datestamp = pharmamine_datestamp,
    chembl_release = chembl_db_release,
    path_data_raw = path_data_raw)

## Get all targeted anticancer drugs from Open Targets Platform
opentargets_targeted_cancer_drugs <-
  get_opentargets_cancer_drugs(
    path_data_raw = path_data_raw,
    ot_version = opentargets_version,
    uniprot_release = uniprot_release)

## Merge information from Open Targets Platform and NCI targeted drugs
## 1) By molecule chembl id
## 2) By name (if molecule chembl id does not provide any cross-ref)

ot_nci_matched <- opentargets_targeted_cancer_drugs %>%
  dplyr::left_join(nci_antineo_chembl, by = c("molecule_chembl_id")) %>%
  dplyr::mutate(
    nci_concept_display_name =
      dplyr::if_else(is.na(nci_concept_display_name) &
                       !stringr::str_detect(drug_name,"[0-9]"),
                     Hmisc::capitalize(tolower(drug_name)),
                     nci_concept_display_name)) %>%
  dplyr::mutate(
    nci_concept_display_name =
      dplyr::if_else(is.na(nci_concept_display_name) &
                       stringr::str_detect(drug_name,"[0-9]"),
                     drug_name,nci_concept_display_name)) %>%
  dplyr::mutate(nci_version = nci_db_release) %>%
  dplyr::mutate(chembl_version = chembl_db_release) %>%
  dplyr::mutate(opentargets_version = opentargets_version)

ot_nci_set1 <- ot_nci_matched %>%
  dplyr::filter(!is.na(nci_t))

## Check for molecule chembl ID's with multiple NCI mappings

ot_nci_set2 <- ot_nci_matched %>%
  dplyr::filter(is.na(nci_t)) %>%
  dplyr::select(-c(nci_concept_display_name,
                   nci_concept_synonym_all,
                   nci_t,nci_concept_definition,
                   drug_name_nci)) %>%
  dplyr::mutate(drug_name_lc = tolower(drug_name)) %>%
  dplyr::left_join(dplyr::select(nci_antineo_chembl,
                                 -molecule_chembl_id),
                   by = c("drug_name_lc" = "drug_name_nci")) %>%
  dplyr::rename(drug_name_nci = drug_name_lc) %>%
  dplyr::filter(!is.na(nci_t)) %>%
  dplyr::anti_join(ot_nci_set1, by = c("drug_name"))


ot_nci_set3 <- ot_nci_matched %>%
  dplyr::filter(is.na(nci_t)) %>%
  dplyr::select(-c(nci_concept_display_name,
                   nci_concept_synonym_all,
                   nci_t,nci_concept_definition,
                   drug_name_nci)) %>%
  dplyr::mutate(drug_name_lc = tolower(drug_name)) %>%
  dplyr::left_join(nci_antineo_nochembl,
                   by = c("drug_name_lc" = "drug_name_nci")) %>%
  dplyr::rename(drug_name_nci = drug_name_lc) %>%
  dplyr::anti_join(ot_nci_set1, by = c("drug_name")) %>%
  dplyr::anti_join(ot_nci_set2, by = c("drug_name"))

ot_cancer_drugs <- ot_nci_set1 %>%
  dplyr::bind_rows(ot_nci_set2) %>%
  dplyr::bind_rows(ot_nci_set3) %>%
  dplyr::left_join(
    dplyr::select(gene_info, target_ensembl_gene_id,
                  target_entrezgene),
    by = "target_ensembl_gene_id")


## NCI drugs/regimens with ChEMBL identifier (not present in Open Targets)
other_nci_chembl_chemotherapies <- nci_antineo_all %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::select(nci_t, molecule_chembl_id, nci_concept_display_name,
                nci_concept_definition, drug_name_nci, nci_concept_synonym_all) %>%
  dplyr::anti_join(ot_cancer_drugs, by = c("molecule_chembl_id")) %>%
  dplyr::mutate(nci_version = nci_db_release) %>%
  dplyr::mutate(chembl_version = chembl_db_release, opentargets_version = NA,
                drug_name = toupper(nci_concept_display_name)) %>%
  dplyr::anti_join(ot_cancer_drugs, by = c("drug_name"))

## NCI drugs/regimens without ChEMBL (not present in Open Targets)
other_nci_nochembl_chemotherapies <- nci_antineo_all %>%
  dplyr::filter(is.na(molecule_chembl_id)) %>%
  dplyr::anti_join(ot_cancer_drugs, by = "nci_concept_display_name") %>%
  dplyr::select(nci_t, nci_concept_display_name, nci_concept_definition,
                drug_name_nci, nci_concept_synonym_all, molecule_chembl_id) %>%
  dplyr::mutate(nci_version = nci_db_release) %>%
  dplyr::mutate(chembl_version = chembl_db_release, opentargets_version = NA,
                drug_name = toupper(nci_concept_display_name))


all_cancer_drugs <- ot_cancer_drugs %>%
  dplyr::bind_rows(other_nci_chembl_chemotherapies) %>%
  dplyr::bind_rows(other_nci_nochembl_chemotherapies) %>%
  dplyr::mutate(molecule_chembl_id =
                  dplyr::if_else(nci_concept_display_name == "Anti-Thymocyte Globulin",
                                 as.character(NA),
                                 as.character(molecule_chembl_id))) %>%
  dplyr::select(target_genename, target_symbol, target_type,
                target_ensembl_gene_id, target_entrezgene,
                target_uniprot_id, target_chembl_id,
                dplyr::everything()) %>%
  dplyr::arrange(drug_name) %>%
  dplyr::distinct()




drug_target_patterns <-
  read.table(file = file.path(
    path_data_raw,
    "custom_drug_target_regex_nci.tsv"),
    sep = "\t", header = T, stringsAsFactors = F, quote = "") %>%
  dplyr::inner_join(gene_info) %>%
  dplyr::distinct()


all_inhibitors_no_target <- all_cancer_drugs %>%
  dplyr::filter(is.na(target_symbol)) %>%
  dplyr::filter(stringr::str_detect(nci_concept_display_name,
                                    "(I|i)nhibitor|antagonist|blocker"))

custom_nci_targeted_drugs <- data.frame()
for(i in 1:nrow(drug_target_patterns)){
  pattern <- drug_target_patterns[i, "pattern"]
  target_symbol <- drug_target_patterns[i, "symbol"]
  target_genename <- drug_target_patterns[i, "genename"]
  target_entrezgene <- drug_target_patterns[i, "target_entrezgene"]
  target_type <- drug_target_patterns[i, "target_type"]
  target_ensembl_gene_id <- drug_target_patterns[i, "target_ensembl_gene_id"]

  hits <- all_inhibitors_no_target %>%
    dplyr::filter(stringr::str_detect(nci_concept_display_name,
                                      pattern = pattern))
  if(nrow(hits) > 0){
    hits$target_symbol <- target_symbol
    hits$target_genename <- target_genename
    hits$target_type <- target_type
    hits$target_entrezgene <- target_entrezgene
    hits$target_ensembl_gene_id <- target_ensembl_gene_id
    hits$drug_clinical_source <- "nci_thesaurus_custom"
    hits$cancer_drug <- TRUE
    custom_nci_targeted_drugs <- custom_nci_targeted_drugs %>%
      dplyr::bind_rows(hits)
  }
}

all_cancer_drugs_final <-
  dplyr::anti_join(all_cancer_drugs, custom_nci_targeted_drugs,
                   by = "nci_concept_display_name") %>%
  dplyr::bind_rows(custom_nci_targeted_drugs) %>%
  dplyr::arrange(target_symbol, nci_concept_display_name) %>%
  dplyr::mutate(drug_action_type = dplyr::if_else(
    stringr::str_detect(tolower(nci_concept_display_name),"inhibitor") &
      is.na(drug_action_type),
    "INHIBITOR",
    as.character(drug_action_type))) %>%
  dplyr::mutate(cancer_drug = dplyr::if_else(
    is.na(cancer_drug) &
      (stringr::str_detect(
        tolower(nci_concept_definition),
        "anti-tumor|chemotherapy|cancer vaccine|immunothera|monoclonal antibody|antineoplastic|treatment of cancer|treatment of metastat") |
      stringr::str_detect(tolower(nci_concept_display_name)," regimen|recombinant|carcinoma|immune checkpoint|anti-programmed cell death ")),
    as.logical(TRUE),
    as.logical(cancer_drug)
  )) %>%
  dplyr::mutate(drug_action_type = dplyr::if_else(
    is.na(drug_action_type) &
      stringr::str_detect(drug_action_type,"^(SUBSTRATE|HYDROLYTIC ENZYME|RELEASING AGENT)"),
    paste0(drug_action_type,"_OTHER"),
    as.character(drug_action_type)
  ))

cancer_drug <- all_cancer_drugs_final %>%
  dplyr::distinct()

usethis::use_data(cancer_drug, internal = T, overwrite = T)

unique_chembl_pubchem <- all_cancer_drugs %>%
  dplyr::select(molecule_chembl_id, nci_concept_display_name) %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::distinct() %>%
  dplyr::left_join(chembl_pubchem_xref,by="molecule_chembl_id")

non_ambiguous_synonyms <- as.data.frame(
  all_cancer_drugs %>%
    dplyr::select(molecule_chembl_id, nci_concept_display_name, drug_name_nci) %>%
    dplyr::distinct() %>%
    dplyr::group_by(drug_name_nci) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n == 1 & nchar(drug_name_nci) >= 4)
  )

antineopharma_synonyms <- all_cancer_drugs %>%
  dplyr::select(drug_name_nci, nci_concept_display_name, molecule_chembl_id) %>%
  dplyr::inner_join(non_ambiguous_synonyms, by = "drug_name_nci") %>%
  dplyr::rename(alias = drug_name_nci) %>%
  dplyr::select(-n) %>%
  dplyr::distinct()

tmp <- antineopharma_synonyms %>%
  dplyr::select(molecule_chembl_id, nci_concept_display_name) %>%
  dplyr::mutate(alias = tolower(nci_concept_display_name)) %>%
  dplyr::distinct()

antineopharma_synonyms <- antineopharma_synonyms %>%
  dplyr::bind_rows(tmp) %>%
  dplyr::arrange(nci_concept_display_name) %>%
  dplyr::distinct()

i <- 1
#antineopharma_properties_chembl <- data.frame()
#antineopharma_moa_chembl <- data.frame()
antineopharma_properties_pubchem <- data.frame()
#antineopharma_synonyms_pubchem <- data.frame()

antineopharma_synonyms_pubchem <-
  readRDS(file = file.path(path_data_tmp_processed, "antineopharma_synonyms.rds"))
tmp2 <- antineopharma_synonyms_pubchem %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::select(nci_concept_display_name) %>%
  dplyr::distinct()
unique_chembl_pubchem <- unique_chembl_pubchem %>%
  dplyr::anti_join(tmp2, by = "nci_concept_display_name")


while(i <= nrow(unique_chembl_pubchem)){
  id_chembl <- unique_chembl_pubchem[i,"molecule_chembl_id"]
  id_pchem <- unique_chembl_pubchem[i,"pubchem_cid"]
  name <- unique_chembl_pubchem[i,"nci_concept_display_name"]
  props <- get_compound_properties(molecule_chembl_id = id_chembl,
                                   pchem_cid = id_pchem,
                                   version_chembl = chembl_db_release)

  for(p in c('properties_pubchem')){
    if(!is.null(props[[p]]) & !is.null(props[['molecule_chembl_id']])){
      df <- props[[p]]
      if('pubchem_synonyms' %in% colnames(df)){
        synonym_df <- data.frame('molecule_chembl_id' = props[['molecule_chembl_id']],
                                 'alias' = stringr::str_split(props[[p]]['pubchem_synonyms'],"@@@@")[[1]],
                                 'nci_concept_display_name' = name,
                                 stringsAsFactors = F)
        antineopharma_synonyms_pubchem <- dplyr::bind_rows(antineopharma_synonyms_pubchem, synonym_df)
        df[,'pubchem_synonyms'] <- NULL
      }
    }
  }
  cat(i,name,'\n')
  if(i %% 50 == 0){
    saveRDS(antineopharma_synonyms_pubchem, file = file.path(path_data_tmp_processed,
                                                             paste0("antineopharma_synonyms_",
                                                          pharmamine_datestamp, "_", i, ".rds")))
  }
  i <- i + 1
}

unambiguous_drug_aliases <- antineopharma_synonyms %>%
  dplyr::bind_rows(antineopharma_synonyms_pubchem) %>%
  dplyr::select(alias, nci_concept_display_name) %>%
  dplyr::distinct() %>%
  dplyr::group_by(alias) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n == 1) %>%
  dplyr::select(alias)

antineopharma_synonyms <-
  dplyr::bind_rows(antineopharma_synonyms, antineopharma_synonyms_pubchem) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(unambiguous_drug_aliases, by = "alias") %>%
  dplyr::distinct()

cancer_drug_synonyms <- antineopharma_synonyms

usethis::use_data(cancer_drug_synonyms, overwrite = T)
