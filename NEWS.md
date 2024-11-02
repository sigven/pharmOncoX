# Version 1.8.0

* CIViC update (20241102)
* NCI Thesaurus update (24.09e)
* OTP update (2024.09)
* Clean-up of code for ChEMBL drug mappings
  - skip DGIdb cross-references - rely on OTP

# Version 1.7.3

* CIViC update (20240918)

# Version 1.7.2

* CIViC update (20240904)

# Version 1.7.1

* CIViC update (20240826)
  - added location of MET exon 14 skipping splice site (grch37 principal transcript)

# Version 1.7.0

* CIViC update (20240807)
* NCI Thesaurus update (24.07e)
* MitelmanDB update (20240715)
* New dataset: DepMap (cell line) RNA fusion events

# Version 1.6.10 

* Fixed some erroneous drug classifications
* CIViC update (20240621)
* OTP update (2024.06)

# Version 1.6.8 (June 7th 2024)

* NCI Thesaurus update (24.05d)

# Version 1.6.7 (May 23rd 2024)

* NCI Thesaurus update (24.04e)
* CIViC update (20240523)

# Version 1.6.4 (April 30th 2024)

* Improved clinical (tumor site) annotations of fusions from MitelmanDB

# Version 1.6.3 (April 26th 2024)

* CIViC update (20240426)

# Version 1.6.2 (April 12th 2024)

* NCI Thesaurus update (24.03d)

# Version 1.6.1 (March 26th 2024)

* Open Targets Platform update (v2024.03)

# Version 1.6.0 (March 9th 2024)

* CIViC/Mitelmandb updates

# Version 1.5.9 (February 29th 2024)

* Updated CIViC

# Version 1.5.8 (February 6th 2024)

* NCI Thesaurus update (24.01e)

# Version 1.5.7 (February 3rd 2024)

* Fixed buggy alias from PubChem

# Version 1.5.6 (January 31st 2024)

* Refine possible values in `treatment_category` argument to `get_drugs()` function
* Fix anti-androgen classification
* Fixed bug in alias type notation for copy numbers and expression biomarkers
* Filter noise and rank output in helper function `get_targeted_drugs`
* Include variant primary name in biomarker variants

# Version 1.5.0 (January 25th 2024)

* Considerable improvement of drug classifications (output column `atc_level3`),
considering both targeted agents and chemotherapy
* argument `drug_targeted_agent` renamed to `treatment_category`
* helper function `get_on_off_label_drugs` renamed to `get_targeted_drugs`

# Version 1.4.10 (January 17th 2024)

* Updated NCI Thesaurus (23.12d)
* Updated CIViC data (20240114)
* Revised manual drug classification (building upon ATC)

# Version 1.4.8 (December 12th 2023)

* Updated NCI Thesaurus - 23.11d
* Updated CIViC data (20231212)

# Version 1.4.8 (November 30th 2023)

* Updated OTP (v2023.12)
* Updated Mitelman database (20231016)

# Version 1.4.7 (November 3rd 2023)

* Fixed corrupt biomarker literature data frame

# Version 1.4.6 (October 18th 2023)

* Fixed metadata
* Fixed out-dated version descriptions (GH site)
* Do not exclude compounds with missing drug action type (platins etc.)

# Version 1.4.5 (October 6th 2023)

* Updated NCI Thesaurus - 23.09d
* Updated CIViC data
* Added drug class (ATC) to `get_on_off_label_drugs()`

# Version 1.4.4 (September 1st 2023)

* Updated NCI Thesaurus - 23.08d
* Updated CIViC data
* Updated Mitelman database (20230803)
* `get_biomarkers()` now exported as a main function



