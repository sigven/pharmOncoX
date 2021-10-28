#' Dataset of targeted and non-targeted cancer drugs
#'
#' A dataset containing an extensive collection of anti-cancer drugs, i.e. targeted compounds, chemotherapy regimens etc. The
#' data has been harvested from the Open Targets Platform, the NCI Thesaurus, and ChEMBL. Note that each entry corresponds
#' to a "drug synonym-target-indication association"
#'
#' @format \bold{oncopharmadb} - A data frame with 418,620 rows and 49 columns:
#' \itemize{
#'   \item \emph{drug_name} - Primary drug name (upper case, NCI Thesaurus)
#'   \item \emph{nci_concept_display_name} - Primary drug name (NCI Thesaurus)
#'   \item \emph{drug_type} - type of drug molecule (Antibody, small molecule etc)
#'   \item \emph{drug_action_type} - main action elicited by drug (antagonist, inhibitor, stabiliser etc)
#'   \item \emph{molecule_chembl_id} - ChEMBL compound identifier
#'   \item \emph{drug_max_phase_indication} - maximum clinical phase for drug (given indication)
#'   \item \emph{drug_max_ct_phase} - maximum clinical phase for drug (any indication)
#'   \item \emph{target_genename} - name/description of drug target
#'   \item \emph{target_symbol} - symbol of drug target
#'   \item \emph{target_type} - type of drug target (single protein, protein family etc.)
#'   \item \emph{target_ensembl_gene_id} - Ensembl gene identifier of drug target
#'   \item \emph{target_entrezgene} - Entrez gene identifier of drug target
#'   \item \emph{target_uniprot_id} - UniProt identifier of drug target
#'   \item \emph{disease_efo_id} - EFO (Experimental Factor Ontology) identifier for drug indication
#'   \item \emph{disease_efo_label} - EFO (Experimental Factor Ontology) label for drug indication
#'   \item \emph{cui} - UMLS metathesaurus identifier (concept unique identifier) for drug indication
#'   \item \emph{cui_name} - UMLS metahesaurus name for drug indication
#'   \item \emph{primary_site} - primary tumor site/type (obtained through https://github.com/sigven/oncoPhenoMap)
#'   \item \emph{nci_concept_synonym} - drug synonym (NCI, Open Targets Platform)
#'   \item \emph{nci_concept_synonym_all} - all drug synonyms (NCI Thesurus, Open Targets Platform)
#'   \item \emph{drug_clinical_source} - underlying source for drug entry (DailyMed, clinicaltrials.gov, FDA etc.)
#'   \item \emph{drug_approved_indication} - logical indicating if drug has an approved indication
#'   \item \emph{drug_blacbox_warning} - logical indicating if drug has blackbox warning
#'   \item \emph{drug_year_first_approval} - year drug was first approved
#'   \item \emph{drug_clinical_id} - drug clinical identifier (clinicaltrials.gov, DailyMed, FDA etc.)
#'   \item \emph{nci_t} - NCI Thesaurus identifier
#'   \item \emph{nci_concept_definition} - detailed description of drug mechanism-of-action (NCI Thesaurus)
#'   \item \emph{nci_version} - version of NCI used for this version of oncoPharmaDB
#'   \item \emph{chembl_version} - version of ChEMBL used for this version of oncoPharmaDB
#'   \item \emph{opentargets_version} - version of Open Target Platform used for this version of oncoPharmaDB
#'   \item \emph{comb_regimen_indication} - type of drug molecule (Antibody, small molecule etc)
#'   \item \emph{immune_checkpoint_inhibitor} - logical indicating if drug is an immune checkpoint inhibitor
#'   \item \emph{topoisomerase_inhibitor} - logical indicating if drug is a topoisomerase inhibitor
#'   \item \emph{tubulin_inhibitor} - logical indicating if drug is a tubulin inhibitor
#'   \item \emph{kinase_inhibitor} - logical indicating if drug is a kinase inhibitor
#'   \item \emph{hdac_inhibitor} - logical indicating if drug is a HDAC inhibitor
#'   \item \emph{parp_inhibitor} - logical indicating if drug is a PARP inhibitor
#'   \item \emph{bet_inhibitor} - logical indicating if drug is a BET inhibitor
#'   \item \emph{ar_antagonist} - logical indicating if drug is an androgen receptor antagonist
#'   \item \emph{monoclonal_antibody} - logical indicating if drug is an antibody
#'   \item \emph{antimetabolite} - logical indicating if drug is an antimetabolite
#'   \item \emph{angiogenesis_inhibitor} - logical indicating if drug is an angiogenesis inhibitor
#'   \item \emph{alkylating_agent} - logical indicating if drug is an alkylating agent
#'   \item \emph{platinum_compound} - logical indicating if drug is a platinum compound
#'   \item \emph{anthracycline} - logical indicating if drug is an anthracycline
#'   \item \emph{proteasome_inhibitor} - logical indicating if drug is a proteasome inhibitor
#'   \item \emph{hormone_therapy} - logical indicating if drug is a hormone therapy
#'   \item \emph{hedgehog_antagonist} - logical indicating if drug is a hedgehog pathway antagonist
#' }
#'
#'
#' @source \url{https://targetvalidation.org/}
#' @source \url{https://www.ebi.ac.uk/chembl/}
#' @source \url{https://ncithesaurus.nci.nih.gov/ncitbrowser/}
"oncopharmadb"

#' Dataset of drug synonyms - linking known drug aliases to a primary identifier
#'
#' @format \bold{oncopharma_synonyms} - A data frame with 9,946 rows and 52 columns:
#' \itemize{
#'   \item \emph{alias} - drug alias
#'   \item \emph{nci_concept_display_name} - Primary drug name (NCI Thesaurus)
#'   \item \emph{molecule_chembl_id} - ChEMBL compound identifier
#' }
#'
#' @source \url{https://targetvalidation.org/}
#' @source \url{https://www.ebi.ac.uk/chembl/}
#' @source \url{https://ncithesaurus.nci.nih.gov/ncitbrowser/}
"oncopharma_synonyms"
