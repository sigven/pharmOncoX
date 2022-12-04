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


