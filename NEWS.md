# Version 1.4.1 (July 23rd 2023)

* Added missing biomarker aliases for synonymous variants

# Version 1.4.0 (July 11th 2023)

* Updated metadata
* Updated CIViC

# Version 1.3.8 (July 5th 2023)

* Cleaned biomarker data (CIViC)

# Version 1.3.7 (June 28th 2023)

* Updated NCI Thesaurus (23.06d release)
* Updated biomarker data (CIViC)
* Updated Open Targets Platform (2023.06)

# Version 1.3.6 (June 1st 2023)

* Updated NCI Thesaurus (23.05d release)
* Cleaned biomarker alteration types

# Version 1.3.5 (May 24th 2023)

* Removed general words from drug alias list
* Fixed bug in biomarker alias parsing

# Version 1.3.3 (May 15th 2023)

* Added additional drug entries (primary indication non-cancer)
* Fixed a bug in parsing of drugs from Open Targets Platform

# Version 1.3.0 (May 10th 2023)

* Updated NCI Thesaurus (23.04d release)
* Drugs are now provided with [ATC](https://www.whocc.no/atc_ddd_index/) classification labels
* Major changes/revision to input arguments in `get_drugs()`
   - all arguments for filtering towards specific inhibitors etc are now removed (e.g. `is_alkylating_agent`,`is_angiogenesis_inhibitor` etc.), such filtering can be performed by the user, e.g. by considering the drug class labels from ATC or the drug targets
   - Other naming changes to arguments for `get_drugs()`
       * `drug_is_targeted` -> `drug_targeted_agent`
       * `inhibitor_only` -> `drug_action_inhibition`
       * `drug_approved_later_than` -> `drug_approval_year`
       * `source_opentargets_only` -> `drug_source_opentargets`
   - The `drug_targeted_agent` option has been improved basd on ATC labeling
   - New arguments for `get_drugs()`:
       * `drug_cancer_indication` - logical indicating if resulting drug records 
should be for those indicated for cancer conditions only (i.e. from approved 
conditions, listed in clinical trials etc.) - defaults to TRUE
       * `drug_classified_cancer` - logical indicating if resulting drug 
records should be for those classified only in the "L" class of ATC (
"ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS") - defaults to TRUE
       
       
# Version 1.2.0 (April 26th 2023)

* Using [sigven/phenOncoX](https://github.com/sigven/phenOncoX) v0.5.8
  for cancer phenotype matching
* Updated NCI Thesaurus (23.03d release)
* New helper function to retrieve drug biomarkers (`pharmOncoX:::get_biomarkers`)

# Version 1.1.9 (March 1st 2023)

* Updated NCI Thesaurus (23.02d release)
* Updated biomarker data (CIViC)
* Updated Open Targets Platform (2023.02)

# Version 1.1.8 (February 6th 2023)

* Resolved different drug names mapped to the same ChEMBL identifier

# Version 1.1.7 (February 1st 2023)

* Using [sigven/phenOncoX](https://github.com/sigven/phenOncoX) v0.5.6
  for cancer phenotype matching
* Updated NCI Thesaurus (23.01e release)
* Updated Mitelman/CIViC data

# Version 1.1.6 (January 26th 2023)

* Updated Mitelman/CIViC data

# Version 1.1.5 (January 19th 2023)

* Rescued some missing fusion entries from Mitelman database (`pharmOncoX:::get_biomarkers()`)

# Version 1.1.4 (January 12th 2023)

* Upgraded NCI Thesaurus (v22.12d)

# Version 1.1.3 (December 19th 2022)

* Upgraded [sigven/phenOncoX](https://github.com/sigven/phenOncoX) to v0.5.5
* Excluded immunomodulatory adenosine receptor antagonists from the immune
checkpoint inhibition category

# Version 1.1.2 (December 9th 2022)

* Added drug category *iap_inhibitor* - inhibitor of IAP 
(*Inhibitor of Apoptosis*) proteins

# Version 1.1.1 (December 5th 2022)

* Removal of multiple duplicate drug entries (Open Targets Platform 
vs. NCI thesaurus)

# Version 1.1.0 (December 5th 2022)

* Updated NCI Thesaurus (2022.11d)
* Rename package to *pharmOncoX* 
* Rescue drug aliases that were erroneously discarded

# Version 1.0.0 (November 27th 2022)

* Updated NCI Thesaurus (22.10e)
* Updated Open Targets Platform (2022.11)
* Renamed `get_onco_drugs()` to `get_drugs()`

# Version 0.9.8 (October 3rd 2022)

* Fixed bug in v0.9.7, duplicate data record entries (Google Drive, `db_id_ref` in `R/sysdata.rda`)
* Updated NCI Thesaurus (2022.09d)
* Updated Open Targets Platform (2022.09)
* Updated biomarkers (CIViC)

# Version 0.9.7 (September 21st 2022)

* Upgraded [sigven/phenOncoX](https://github.com/sigven/phenOncoX) to v0.4.0
* Upgraded [sigven/geneOncoX](https://github.com/sigven/geneOncoX)
* Updated biomarkers (CIViC database)

# Version 0.9.6 (September 10th 2022)

* Added logical argument `inhibitor_only` to `get_onco_drugs` to retrieve only drugs with an
inhibitory role with respect to mechanism-of-action
* Fixed minor bug in parsing of approved indications for drugs from Open Targets Platform

# Version 0.9.5 (September 7th 2022)

* Updated NCI Thesaurus (22.08e)
* Improved classification of alkylating agents
* Improved drug salts detection
* Excluded several radio(immuno)conjugate entries

# Version 0.9.3 (September 1st 2022)

* Removed drug aliases that coincide with english words

# Version 0.9.1 (August 30th 2022)

* Added helper function to retrieve drug biomarkers

# Version 0.9.0 (August 25th 2022)

* Renamed package name (**pharmaOncoX**)
  * Complete restructring of `get_onco_drugs`, data is now longer 
    provided within the package, but must be downloaded to a cache directory
    provided by the user, utilizing the [googledrive]() package
  * Currently, `get_onco_drugs` is the only function provided, considering
    adding extended support for drug biomarkers moving forward
* Now using [lgr](https://github.com/s-fleck/lgr) for logging


