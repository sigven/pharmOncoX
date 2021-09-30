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
