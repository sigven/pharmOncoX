# Version 1.4.4 (September 1st 2023)

* Updated NCI Thesaurus - 23.08d
* Updated CIViC data
* Updated Mitelman database (20230803)
* `get_biomarkers()` now exported as a main function

# Version 1.4.3 (August 7th 2023)

* NCI Thesaurus 23.07e

# Version 1.4.2 (August 1st 2023)

* Added biomarker aliases
* Updated CIViC data

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


