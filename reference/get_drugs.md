# Get antineoplastic drugs and drug regimens

Downloads preprocessed datasets to a local cache directory and returns a
selected set of drugs based on various criteria set by the user.

The dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding drug resources used

- `records` - a data frame with drug records

## Usage

``` r
get_drugs(
  cache_dir = NA,
  force_download = FALSE,
  exclude_salt_forms = TRUE,
  exclude_adc = FALSE,
  treatment_category = c("targeted_therapy_classified", "targeted_therapy_unclassified",
    "chemo_therapy_classified", "hormone_therapy_classified",
    "immuno_suppressants_classified", "other"),
  drug_is_approved = FALSE,
  drug_target = NULL,
  drug_action_type = NULL,
  drug_indication_main = NULL,
  drug_source_opentargets = FALSE,
  drug_cancer_indication = TRUE,
  drug_classified_cancer = TRUE,
  drug_has_blackbox_warning = FALSE,
  drug_approval_year = 1939,
  drug_minimum_phase_any_indication = 0,
  output_resolution = "drug2target2indication",
  drug_action_inhibition = F
)
```

## Arguments

- cache_dir:

  local cache directory for data retrieval

- force_download:

  force download data from remote repository even if data exists in
  cache

- exclude_salt_forms:

  exclude salt forms of drugs

- exclude_adc:

  exclude antibody-drug conjugates (ADCs)

- treatment_category:

  main treatment category, classified according to ATC or not
  ('targeted_therapy_classified',
  'targeted_therapy_unclassified','chemo_therapy_classified','hormone_therapy_classified',
  'immuno_suppressants_classified','other')

- drug_is_approved:

  logical indicating if resulting drug records should contain approved
  drugs only

- drug_target:

  character vector with drug targets (gene symbols) for drug records
  included in results

- drug_action_type:

  character vector with drug action types to include in drug record
  list - possible values "INHIBITOR","AGONIST","MODULATOR","ANTAGONIST",
  "BLOCKER","ACTIVATOR","BINDING AGENT","OPENER",
  "STABILISER","CROSS-LINKING AGENT",DISRUPTING AGENT","OTHER"

- drug_indication_main:

  character vector with main tumor types for which drug(s) are
  indicated. Possible values: "Adrenal Gland","Biliary Tract",
  "Bladder/Urinary Tract","Bone","Breast","Cervix","CNS/Brain",
  "Colon/Rectum","Esophagus/Stomach","Eye","Head and Neck","Kidney",
  "Liver","Lung","Lymphoid","Myeloid","Ovary/Fallopian Tube",
  "Pancreas","Penis","Peripheral Nervous System","Peritoneum",
  "Pleura","Prostate","Skin","Soft Tissue","Testis","Thymus",
  "Thyroid","Uterus","Vulva/Vagina"

- drug_source_opentargets:

  logical indicating if resulting drug records should contain drug
  records from Open Targets Platform/ChEMBL only

- drug_cancer_indication:

  logical indicating if resulting drug records should be for those
  indicated for cancer conditions only (from approved conditions, found
  in clinical trials etc.)

- drug_classified_cancer:

  logical indicating if resulting drug records should be for those
  classified only in the "L" class of ATC ( "ANTINEOPLASTIC AND
  IMMUNOMODULATING AGENTS") only

- drug_has_blackbox_warning:

  logical indicating if resulting drug records should contain drugs with
  black box warnings only

- drug_approval_year:

  only include records for drugs approved later than this date (year)

- drug_minimum_phase_any_indication:

  only include drug records that are in a clinical phase (any
  indication) greater or equal than this phase

- output_resolution:

  dictate output record resolution
  ('drug','drug2target','drug2target2indication')

- drug_action_inhibition:

  logical indicating to only return drug records with inhibitory
  mechanism-of-action

## Value

The `records` data frame contains the following columns (only selected
columns will be shown based on the value of `output_resolution`)

- *drug_id* - drug identifier (pharmaOncoX)

- *drug_name* - primary drug name (upper case, NCI Thesaurus)

- *drug_type* - type of drug molecule (Antibody, small molecule etc)

- *molecule_chembl_id* - ChEMBL compound identifier

- *drug_action_type* - main action elicited by drug (antagonist,
  inhibitor, stabiliser etc)

- *drug_alias* - collection of unambiguous drug aliases (separated by
  '\|')

- *nci_concept_definition* - detailed description of drug
  mechanism-of-action (NCI Thesaurus)

- *opentargets* - logical - drug is found in the Open Targets Platform
  resource

- *is_salt* - logical - drug record represents a salt form (excluded by
  default)

- *is_adc* - logical - drug record represents an antibody-drug conjugate
  (ADC - excluded by default)

- *drug_blacbox_warning* - logical indicating if drug has blackbox
  warning

- *nci_t* - NCI thesaurus identifier

- *target_symbol* - gene symbol of drug target

- *target_entrezgene* - Entrez gene identifier of drug target

- *target_genename* - gene name/description of drug target

- *target_ensembl_gene_id* - Ensembl gene identifier of drug target

- *target_type* - type of drug target (single protein, protein family
  etc.)

- *drug_max_phase_indication* - maximum clinical phase for drug (given
  indication)

- *drug_approved_indication* - logical indicating if drug has an
  approved indication

- *drug_frac_cancer_indications* - fraction of drug indications that are
  for cancers

- *drug_approved_noncancer* - logical indicating if drug is approved for
  a non-cancer disease

- *drug_n_indications* - number of indications for the given drug (from
  approved indications, clinical trials etc)

- *drug_year_first_approval* - year drug was first approved

- *drug_max_ct_phase* - maximum clinical phase for drug (any indication)

- *disease_efo_id* - EFO (Experimental Factor Ontology) identifier for
  drug indication

- *disease_efo_label* - EFO (Experimental Factor Ontology) label for
  drug indication

- *primary_site* - primary tumor site/type (obtained through
  https://github.com/sigven/oncoPhenoMap)

- *drug_clinical_id* - drug clinical identifier (clinicaltrials.gov,
  DailyMed, FDA etc.)

- *drug_clinical_source* - underlying source for drug entry (DailyMed,
  clinicaltrials.gov, FDA etc.)

- *atc_code_level1* - drug identifier ATC (level 1)

- *atc_level1* - drug label ATC (level 1)

- *atc_code_level1* - drug identifier ATC (level 2)

- *atc_level2* - drug label ATC (level 2)

- *atc_code_level3* - drug identifier ATC (level 3)

- *atc_level3* - drug label ATC (level 3)

- *atc_treatment_category* - treatment category (targeted/chemo/hormone,
  cancer/other etc)
