#' Datasets of compound biomarkers - omics properties that are linked to drug response
#'
#' @format \bold{curated} - A data frame with 288,910 rows and 26 columns:
#' \itemize{
#'   \item \emph{evidence_id} - evidence identifier
#'   \item \emph{evidence_url} - evidence URL
#'   \item \emph{variant_id} - variant identifier
#'   \item \emph{variant_type} - type of variant
#'   \item \emph{variant} - variant (alias)
#'   \item \emph{alteration_type} - type of genetic alteration
#'   \item \emph{symbol} - drug target (gene symbol)
#'   \item \emph{therapeutic_context} - therapeutic context (drug sensitivity/resistance)
#'   \item \emph{disease_ontology_id} - disease ontology (DO) identifier (mapped - oncoPhenoMap)
#'   \item \emph{cancer_type} - tumor type/disease (provided by source)
#'   \item \emph{molecule_chembl_id} - ChEMBL compound identifier
#'   \item \emph{citation_id} - citation identifier (PMID etc)
#'   \item \emph{evidence_description} - evidence description
#'   \item \emph{evidence_type} - evidence type
#'   \item \emph{evidence_level} - evidence level
#'   \item \emph{evidence_direction} - evidence direction
#'   \item \emph{clinical_significance} - clinical significance
#'   \item \emph{variant_origin} - variant origin (somatic/germline)
#'   \item \emph{biomarker_source_db} - type of biomarker source database
#'   \item \emph{biomarker_entity} - logical indicating biomarker entity
#'   \item \emph{efo_id} - EFO identifier (mapped - oncoPhenoMap)
#'   \item \emph{efo_name} - EFO name (mapped - oncoPhenoMap)
#'   \item \emph{do_name} - disease ontology (DO) name (mapped - oncoPhenoMap)
#'   \item \emph{cui} - disease UMLS identifier (mapped - oncoPhenoMap)
#'   \item \emph{cui_name} - disease UMLS name (mapped - oncoPhenoMap)
#'   \item \emph{primary_site} - primary tumor site (mapped - oncoPhenoMap)
#' }
#'
#'
#' @format \bold{invitro_screen} - A data frame with 23,814 rows and 8 columns:
#' \itemize{
#'   \item \emph{compound_id} - compound identifier (PRISM)
#'   \item \emph{model} - type of omics drug prediction model (DepMap/PRISM)
#'   \item \emph{feature_type} - type of omics feature
#'   \item \emph{feature_value} - value of omics feature
#'   \item \emph{feature_importance} - omics feature importance
#'   \item \emph{prediction_accuracy} - model prediction accuracy
#'   \item \emph{nci_concept_display_name} - Primary drug name (NCI Thesaurus)
#'   \item \emph{molecule_chembl_id} - ChEMBL compound identifier
#' }
#'
#'
#' @source \url{https://civicdb.org/}
#' @source \url{https://www.cancergenomeinterpreter.org/biomarkers}
#' @source \url{https://pmkb.weill.cornell.edu/}
#' @source \url{https://mitelmandatabase.isb-cgc.org/}
#' @source \url{https://depmap.org/portal/}
#"compound_biomarkers"
