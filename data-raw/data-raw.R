library(magrittr)
pharmamine_datestamp <- '20210921'
nci_db_release <- '21.08e'
chembl_db_release <- 'ChEMBL_29'
opentargets_version <- '2021.06'
uniprot_release <- '2021_03'
dgidb_db_release <- 'v2021_05'
update_dailymed <- F

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

####---UniProt---####
uniprot_map <-
  get_uniprot_map(
    basedir = here::here(),
    uniprot_release = uniprot_release)

ensembl2up <- uniprot_map$uniprot_map %>%
  dplyr::filter(!is.na(uniprot_reviewed)) %>%
  dplyr::select(ensembl_gene_id, uniprot_id) %>%
  dplyr::rename(target_ensembl_gene_id = ensembl_gene_id,
                target_uniprot_id = uniprot_id) %>%
  dplyr::distinct()

gene_info <- gene_info %>%
  dplyr::left_join(ensembl2up, by = "target_ensembl_gene_id")


nci_thesaurus_files <- list()
nci_thesaurus_files[['flat']] <- paste0("Thesaurus_", nci_db_release,".FLAT.zip")
nci_thesaurus_files[['owl']] <- paste0("Thesaurus_", nci_db_release,".OWL.zip")
nci_thesaurus_files[['inf_owl']] <- paste0("ThesaurusInf_", nci_db_release,".OWL.zip")



drug_indications_dailymed <-
  get_dailymed_drug_indications(update = update_dailymed,
                                path_data_raw = path_data_raw)

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
download.file(url = antineo_agents_url, destfile = antineo_agents_local, quiet = T)


####---- Dataset of cancer drugs -----####

## Get all anticancer drugs, NCI thesaurus + DGIdb
nci_antineo_all <- get_nci_drugs(
  nci_db_release = nci_db_release,
  overwrite = T,
  path_data_raw = path_data_raw,
  path_data_processed = path_data_tmp_processed)

#system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.owl"))
#system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.txt"))
#system(paste0('rm -f ',path_data_raw, "/nci_thesaurus/*.tsv"))

## NCI anticancer drugs (targeted) - with compound identifier (CHEMBL)
nci_antineo_chembl <- nci_antineo_all %>%
  dplyr::select(nci_t,
                nci_concept_definition,
                nci_concept_display_name,
                molecule_chembl_id,
                drug_name_nci,
                nci_concept_synonym_all) %>%
  dplyr::filter(!(nci_t == "C1806" & drug_name_nci == "gemtuzumab")) %>%
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
  dplyr::anti_join(ot_nci_set1, by = c("nci_concept_display_name")) %>%
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
  dplyr::anti_join(ot_nci_set2, by = c("drug_name")) %>%
  dplyr::anti_join(ot_nci_set1, by = c("nci_concept_display_name")) %>%
  dplyr::anti_join(ot_nci_set2, by = c("nci_concept_display_name"))

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
  dplyr::filter(stringr::str_detect(tolower(nci_concept_display_name),
                                    "inhibitor|antagonist|antibody|blocker") |
                  stringr::str_detect(tolower(nci_concept_display_name),
                                      "ib$|mab$|mab/|^anti-") |
                  (stringr::str_detect(nci_concept_definition,"KRAS") &
                     stringr::str_detect(nci_concept_definition,"inhibitor")))

custom_nci_targeted_drugs <- data.frame()
for(i in 1:nrow(drug_target_patterns)){
  pattern <- drug_target_patterns[i, "pattern"]
  target_symbol <- drug_target_patterns[i, "symbol"]
  target_genename <- drug_target_patterns[i, "genename"]
  target_entrezgene <- drug_target_patterns[i, "target_entrezgene"]
  target_type <- drug_target_patterns[i, "target_type"]
  target_ensembl_gene_id <- drug_target_patterns[i, "target_ensembl_gene_id"]
  target_uniprot_id <- drug_target_patterns[i, "target_uniprot_id"]

  hits <- all_inhibitors_no_target %>%
    dplyr::filter(stringr::str_detect(nci_concept_display_name,
                                      pattern = pattern) |
                    (stringr::str_detect(nci_concept_display_name,"ib$|mab$|mab/") &
                    stringr::str_detect(nci_concept_definition, pattern))
    )



  if(nrow(hits) > 0){

    for(n in 1:nrow(hits)){
      hit <- hits[n,]

      if(stringr::str_detect(hit$nci_concept_display_name,"mab$|monoclonal antibody")){
        hit$drug_type <- "Antibody"
      }else{
        hit$drug_type <- "Small molecule"
      }

      hit$drug_action_type <- "INHIBITOR"
      if(stringr::str_detect(
        tolower(hit$nci_concept_display_name),
        "antagonist")){
        hit$drug_action_type <- "ANTAGONIST"
      }
      if(stringr::str_detect(
        tolower(hit$nci_concept_display_name),
        "blocker")){
        hit$drug_action_type <- "BLOCKER"
      }
      hit$target_symbol <- target_symbol
      hit$target_genename <- target_genename
      hit$target_type <- target_type
      hit$target_entrezgene <- target_entrezgene
      hit$target_ensembl_gene_id <- target_ensembl_gene_id
      hit$target_uniprot_id <- target_uniprot_id
      hit$drug_clinical_source <- "nci_thesaurus_custom"
      hit$cancer_drug <- TRUE
      custom_nci_targeted_drugs <- custom_nci_targeted_drugs %>%
        dplyr::bind_rows(hit)

    }
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
  )) %>%
  dplyr::mutate(comb_regimen_indication = F)



## Add indications retrieved from DailyMed

drugs2max_ct_phase <- all_cancer_drugs_final %>%
  dplyr::filter(!is.na(drug_max_ct_phase)) %>%
  dplyr::select(drug_name_nci, drug_max_ct_phase) %>%
  dplyr::distinct()

all_drugs2targets_no_indications <- all_cancer_drugs_final %>%
  dplyr::select(-c(cui, cui_name, disease_efo_id, disease_efo_label,
                   drug_clinical_source, drug_approved_indication,
                   primary_site, drug_max_ct_phase, drug_clinical_id,
                   drug_max_phase_indication, comb_regimen_indication)) %>%
  dplyr::distinct()

supplemental_drug_indications <- drug_indications_dailymed %>%
  dplyr::select(-drugname_trade) %>%
  dplyr::inner_join(all_drugs2targets_no_indications, by = "drug_name_nci") %>%
  dplyr::select(-drug_name_nci) %>%
  dplyr::mutate(drug_name_nci = tolower(nci_concept_synonym_all)) %>%
  tidyr::separate_rows(drug_name_nci, sep="\\|") %>%
  dplyr::select(-drug_max_ct_phase) %>%
  dplyr::distinct() %>%
  dplyr::left_join(drugs2max_ct_phase, by = "drug_name_nci") %>%
  dplyr::distinct()

all_cancer_drugs_final2 <- all_cancer_drugs_final %>%
  dplyr::bind_rows(supplemental_drug_indications) %>%
  dplyr::arrange(nci_concept_display_name)

oncopharmadb <- all_cancer_drugs_final2 %>%
  dplyr::filter(!(is.na(nci_concept_display_name) &
                    !stringr::str_detect(drug_action_type,"INHIBITOR|BLOCKER"))) %>%
  dplyr::distinct() %>%
  dplyr::filter(!(is.na(molecule_chembl_id) & drug_name_nci == "sunitinib malate")) %>%
  dplyr::mutate(nci_concept_display_name = dplyr::if_else(
    is.na(nci_concept_display_name),
    stringr::str_to_title(drug_name),
    as.character(nci_concept_display_name)
  )) %>%
  dplyr::mutate(nci_concept_display_name = dplyr::if_else(
    nci_concept_display_name == "Cediranib Maleate",
    "Cediranib",
    as.character(nci_concept_display_name)
  )) %>%
  dplyr::mutate(drug_name_nci = dplyr::if_else(
    nci_concept_display_name == "Cediranib" & drug_name_nci == "azd2171 maleate",
    "cediranib",
    as.character(drug_name_nci)
  ))


oncopharmadb <- oncopharmadb %>%
  dplyr::distinct() %>%
  dplyr::mutate(antimetabolite = dplyr::if_else(
    !is.na(nci_concept_definition) &
      stringr::str_detect(
        tolower(nci_concept_definition),
        "antimetabol|anti-metabol|nucleoside analog"),TRUE,FALSE)
    ) %>%
  dplyr::mutate(topoisomerase_inhibitor = dplyr::if_else(
    (!is.na(nci_concept_definition) &
      stringr::str_detect(
        nci_concept_definition,
        "(T|t)opoisomerase II-mediated|(T|t)opoisomerase( I|II )? \\(.*\\) inhibitor|inhibit(ion|or) of (T|t)opoisomerase|(stabilizes|interrupts|binds to|interacts with|inhibits( the activity of)?)( the)?( DNA)? (t|T)opoisomerase|(T|t)opoisomerase( (I|II))? inhibitor")) |
      (!is.na(target_genename) & stringr::str_detect(target_genename,"topoisomerase")),TRUE,FALSE)
  ) %>%
  dplyr::mutate(hedgehog_antagonist = dplyr::if_else(
    (!is.na(nci_concept_definition) &
       stringr::str_detect(
         nci_concept_definition,
         "Hedgehog") & stringr::str_detect(nci_concept_display_name,"Smoothened Antagonist|(ate|ib)$")) |
      (!is.na(nci_concept_display_name) &
         stringr::str_detect(nci_concept_display_name,"Hedgehog Inhibitor|SMO Protein Inhibitor")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(hdac_inhibitor = dplyr::if_else(
    !is.na(target_symbol) &
      stringr::str_detect(
        target_symbol,
        "^HDAC") |
      !is.na(nci_concept_definition) &
      stringr::str_detect(nci_concept_definition,"inhibitor of histone deacetylase"),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(alkylating_agent = dplyr::if_else(
    is.na(drug_moa) &
    !stringr::str_detect(nci_concept_display_name,
                         "antiangiogenic") &
    !is.na(nci_concept_definition) &
      stringr::str_detect(
        tolower(nci_concept_definition),
        "alkylating agent|alkylating activities"),TRUE,FALSE)
  ) %>%
  dplyr::mutate(parp_inhibitor = dplyr::if_else(
    !is.na(target_symbol) &
      stringr::str_detect(
        target_symbol,
        "^PARP[0-9]{1}"),TRUE,FALSE)
  ) %>%
  dplyr::mutate(bet_inhibitor = dplyr::if_else(
    !is.na(target_symbol) &
      stringr::str_detect(
        target_symbol,
        "^BRD(T|[1-9]{1})") |
      (!is.na(nci_concept_display_name) &
      stringr::str_detect(nci_concept_display_name,"BET( Bromodomain)? Inhibitor")),TRUE,FALSE)
  ) %>%
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
  ) %>%
  dplyr::mutate(ar_antagonist = dplyr::if_else(
    (!is.na(target_genename) &
       stringr::str_detect(
         tolower(target_genename),
         "androgen receptor")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(kinase_inhibitor = dplyr::if_else(
    (!is.na(target_symbol) & stringr::str_detect(target_symbol,"EGFR|PTPN11|ABL1|FGFR|PDGFR|CSF1R")) |
    (((!is.na(drug_action_type) &
        stringr::str_detect(tolower(drug_action_type),"blocker|inhibitor|antagonist")) |
       stringr::str_detect(tolower(nci_concept_display_name),"ib$")) &
      (!is.na(target_genename) &
         stringr::str_detect(tolower(target_genename),"kinase|eph receptor"))) |
      (!is.na(nci_concept_definition) & stringr::str_detect(nci_concept_definition,"kinase inhibit(or|ion)")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(angiogenesis_inhibitor = dplyr::if_else(
    stringr::str_detect(tolower(drug_action_type),"blocker|inhibitor|antagonist") &
      (!is.na(nci_concept_display_name) & stringr::str_detect(tolower(nci_concept_display_name), "antiangiogenic|angiogenesis inhibitor")) |
      (!is.na(nci_concept_definition) &
         stringr::str_detect(tolower(nci_concept_definition),"antiangiogenic activities|angiogenesis inhibitor|(inhibiting|blocking)( tumor)? angiogenesis|anti(-)?angiogenic|(inhibits|((inhibition|reduction) of))( .*) angiogenesis")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(monoclonal_antibody = dplyr::if_else(
    stringr::str_detect(tolower(nci_concept_display_name),
                        "^anti-|mab |mab$|monoclonal antibody") &
      (!is.na(nci_concept_definition) &
      stringr::str_detect(nci_concept_definition,"monoclonal antibody")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(proteasome_inhibitor = dplyr::if_else(
    (stringr::str_detect(tolower(nci_concept_display_name),
                        "^proteasome") &
       !stringr::str_detect(tolower(nci_concept_display_name),"vaccine")) |
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           tolower(nci_concept_definition),"proteasome inhibitor|inhibits the proteasome|inhibition of proteasome")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(hormone_therapy = dplyr::if_else(
    stringr::str_detect(tolower(nci_concept_display_name),
                        "aromatase inhib|estrogen receptor (inhibitor|degrader|modulator)") |
      (!is.na(nci_concept_definition) &
         stringr::str_detect(
           tolower(nci_concept_definition),"inhibitor of estrogen|estrogen receptor (modulator|inhibitor|degrader)|antiestrogen|aromatase inhibit(or|ion)") &
         !stringr::str_detect(nci_concept_definition,"antiestrogen resistance")) |
      (!is.na(target_symbol) & stringr::str_detect(target_symbol,"ESR[0-9]")),
    TRUE,FALSE)
  ) %>%
  dplyr::mutate(anthracycline = dplyr::if_else(
    (!is.na(nci_concept_definition) &
       stringr::str_detect(
         tolower(nci_concept_definition),
         "anthracycline|anthracenedione")),
    TRUE, FALSE)
  ) %>%
  dplyr::mutate(immune_checkpoint_inhibitor = dplyr::if_else(
    (!is.na(nci_concept_definition) &
      stringr::str_detect(
        tolower(nci_concept_definition),
        "immune checkpoint inhib")) |
      (stringr::str_detect(nci_concept_display_name,
                          "Tremelimumab|Milatuzumab")) |
      (!is.na(target_symbol) & (target_symbol == "CD274" |
      target_symbol == "CTLA4" | target_symbol == "PDCD1")) |
      (!is.na(nci_concept_definition) & !is.na(target_symbol) &
         stringr::str_detect(nci_concept_definition,
                          "immunemodulating|immune response") &
      target_symbol == "ADORA2A"),
    TRUE,FALSE)
  )


drug_action_types <- as.data.frame(
  oncopharmadb %>%
  dplyr::select(nci_concept_display_name, drug_action_type) %>%
  dplyr::distinct() %>%
  dplyr::group_by(nci_concept_display_name) %>%
  dplyr::summarise(drug_action_type = paste(
    drug_action_type, collapse = "/"
  ))
)

oncopharmadb$drug_action_type <- NULL
oncopharmadb <- oncopharmadb %>%
  dplyr::left_join(drug_action_types,
                   by = "nci_concept_display_name") %>%
  dplyr::select(drug_name, nci_concept_display_name, drug_type,
                drug_action_type, molecule_chembl_id, drug_moa,
                drug_max_phase_indication, dplyr::everything())

usethis::use_data(oncopharmadb, overwrite = T)

rm(all_cancer_drugs)
rm(all_cancer_drugs_final)


####---- Dataset of cancer drug aliases-----####


## unique drug aliases from NCI Thesaurus
non_ambiguous_synonyms <- as.data.frame(
  oncopharmadb %>%
    dplyr::select(molecule_chembl_id, nci_concept_display_name, drug_name_nci) %>%
    dplyr::distinct() %>%
    dplyr::group_by(drug_name_nci) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n == 1 & nchar(drug_name_nci) >= 4)
  )

## Drug alias to NCI concept display name (primary name)
antineopharma_synonyms <- oncopharmadb %>%
  dplyr::select(drug_name_nci, nci_concept_display_name, molecule_chembl_id) %>%
  dplyr::inner_join(non_ambiguous_synonyms, by = "drug_name_nci") %>%
  dplyr::rename(alias = drug_name_nci) %>%
  dplyr::select(-n) %>%
  dplyr::distinct()

## include also the primary name among aliases
tmp <- antineopharma_synonyms %>%
  dplyr::select(molecule_chembl_id, nci_concept_display_name) %>%
  dplyr::mutate(alias = tolower(nci_concept_display_name)) %>%
  dplyr::distinct()

antineopharma_synonyms <- antineopharma_synonyms %>%
  dplyr::bind_rows(tmp) %>%
  dplyr::arrange(nci_concept_display_name) %>%
  dplyr::distinct()


## Extend aliases with those found in PubChem

## get drug set that contains with PubChem cross-references
unique_chembl_pubchem <- oncopharmadb %>%
  dplyr::select(molecule_chembl_id, nci_concept_display_name) %>%
  dplyr::filter(!is.na(molecule_chembl_id)) %>%
  dplyr::distinct() %>%
  dplyr::left_join(chembl_pubchem_xref,by="molecule_chembl_id") %>%
  dplyr::filter(!is.na(pubchem_cid)) %>%
  dplyr::select(-c(chembl_db_version))


## Retrieve aliases for drugs with PubChem x-refs
pubchem_synonym_files <-
  sort(list.files(path = file.path(here::here(), "data-raw","pubchem"),
                  pattern = "CID-Synonym-filtered_",
                  full.names = T))

antineopharma_synonyms_pubchem <- data.frame()
for(f in pubchem_synonym_files){
  cat(f)
  cat('\n')
  synonym_data <- as.data.frame(readr::read_tsv(
    f, col_names = c('pubchem_cid','alias'),
    col_types = "dc",
    progress = F
  ))

  pubchem_alias_df <- synonym_data %>%
    dplyr::inner_join(unique_chembl_pubchem,
                      by = "pubchem_cid")

  if(nrow(pubchem_alias_df) > 0){
    pubchem_alias_df <- pubchem_alias_df %>%
      dplyr::select(-pubchem_cid)
    antineopharma_synonyms_pubchem <-
      antineopharma_synonyms_pubchem %>%
      dplyr::bind_rows(pubchem_alias_df)
  }
  rm(synonym_data)
}

## Only include drug aliases that are unambiguous
unambiguous_drug_aliases <- antineopharma_synonyms %>%
  dplyr::bind_rows(antineopharma_synonyms_pubchem) %>%
  dplyr::filter(!(alias == "nab-paclitaxel" &
                    nci_concept_display_name == "Paclitaxel")) %>%
  dplyr::select(alias, nci_concept_display_name) %>%
  dplyr::distinct() %>%
  dplyr::group_by(alias) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n == 1) %>%
  dplyr::select(alias)

oncopharma_synonyms <-
  dplyr::bind_rows(antineopharma_synonyms, antineopharma_synonyms_pubchem) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(unambiguous_drug_aliases, by = "alias") %>%
  dplyr::distinct()  %>%
  dplyr::mutate(alias = dplyr::if_else(
    alias == "nab-paclitaxel" & nci_concept_display_name == "Paclitaxel",
    "paclitaxel",
    as.character(alias)
  ))


usethis::use_data(oncopharma_synonyms, overwrite = T)
