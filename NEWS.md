# Version 0.7.2 (March 30th 2022)

* Fixed bug in installation documentation

# Version 0.7.1 (March 30th 2022)

* Data updates: NCI Thesaurus (22.03d)
* Updated biomarkers (CIViC)

# Version 0.7.0 (March 4th 2022)

* Data updates: NCI Thesaurus (22.02d), Open Targets Platform (2022.02)

# Version 0.6.9 (February 8th 2022)

* Data update: NCI Thesaurus (22.01e)
* Added [FDA pharmacologic class](https://www.fda.gov/industry/structured-product-labeling-resources/pharmacologic-class) (established pharmacological class - EPC) phrases to `oncopharmadb` data.frame

# Version 0.6.6 (January 17th 2022)

* Fixed duplicate column in `oncopharmadb` (immune checkpoint inhibitor) 
* Updated biomarkers (CIViC)

# Version 0.6.5 (January 4th 2022)

* Data update: NCI Thesaurus (21.12d)
* Renamed exported dataset __oncopharma_synonyms__ -> __compound_synonyms__
* Added datasets on compound biomarkers, both curated (CIViC/CGI/PMKB, _oncoPharmaDB::compound_biomarkers[['curated']]_), and data from invitro drug screens (DepMap/PRISM, _oncoPharmaDB::compound_biomarkers[['invitro_screen']]_)

# Version 0.6.3.900 (December 13th 2021)

* Added missing entries from Open Targets Platform (bug with non-targeted drugs)
* Preparing addition of biomarkers for forthcoming release, providing users a list of
known molecular biomarkers pr. anti-cancer drug (mutations, copy number amplificiations etc., as retrieved from e.g. CIViC), both with respect to resistance and sensitivity

# Version 0.6.3 (December 3rd 2021)

* Data update: NCI Thesaurus (21.11e)

# Version 0.6.2 (November 25th 2021)

* Data update: Open Targets Platform (2021.11)
* General code/documentation improvements - following `devtools::check()`

# Version 0.6.1 (October 28th 2021)

* Data update: NCI Thesaurus (21.10d)
* Removed some redundant/non-informative columns of main data frame (**oncopharmadb**)
  * `target_chembl_id`
  * `drug_synonyms`
  * `drug_description`
  * `drug_moa`
  * `drug_tradenames`
  * `cancer_drug`
* Documented exported data objects (**oncopharmadb** and **oncopharma_synonyms**)
* Annotated gonadotropin releasing hormone analogues with *hormone_therapy*  label
* Added *platinum_compound* as a dedicated drug category (and as a column in **oncopharmadb**) 
  * New argument `is_platinum_compound` to `get_onco_drugs()`


# Version 0.6.0 (September 30th 2021)

* Data update: Open Targets Platform (2021.09)

# Version 0.5.4 (September 29th 2021)

## Fixed

* Bug in assignment of drug categories (some drugs were assigned ambiguous categories)

## Added

* Data update: NCI Thesaurus version 21.09d

# Version 0.5.3 (September 24th 2021)

## Added

* Added multiple warnings in output when conditions result in zero drug record hits
* Added argument 'disease_indication_main', to retrieve drugs indicated for given tumor type(s)
* When *output_resolution* is **drug** or **drug2target**, each record is appended with the following info:
  * disease_indication (list of tumor subtypes within each *primary_site* that are indicated for the drug)
  * disease_indication_max_phase (maximum phase for the indications.. random order (currently not optimal))
  * disease_main_group (list of primary tumors/tissues indicated for drug, i.e. *primary_site*)


# Version 0.5.2.1 (September 23rd 2021)

## Added

* Debugging immune checkpoint inhibitor targets from NCI drugs (CTLA4, PD-1, PD-L1)

# Version 0.5.1 (September 23rd 2021)

## Added

* Improved mapping of drug targets for drugs (monoclonal antibodies, kinase inhibitors etc.) listed solely in NCI thesaurus

# Version 0.5.0 (September 22nd 2021)

## Added

* Mapped several new targets for drugs listed in NCI thesaurus

# Version 0.4.9 (September 22nd 2021)

## Fixed

* Fixed bug related to kinase drug category 
* Added uniprot identifiers to custom targets (experimental drugs, NCI)
* Added imports (tidyverse, assertthat)

# Version 0.4.8 (September 21st 2021)

## Fixed

* Fixed bug related to missing drug indications (cancer types)
* Updated oncoPhenoMap
