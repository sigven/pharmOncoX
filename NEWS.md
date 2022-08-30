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

# Version 0.8.4 (July 27th 2022)

* Updated NCI Thesaurus (22.07d)
* Updated biomarkers (CIViC/Mitelman database)

# Version 0.8.3 (July 1st 2022)

* Updated NCI Thesaurus (22.06d)
* Updated Open Targets Platform (2022.06)
* Updated biomarkers (CIViC)

# Version 0.8.2 (June 19th 2022)

* Added more drugs from NCI thesaurus
* Drug indication updates: (EFO/DO - https://github.com/sigven/oncoPhenoMap)
* Updated biomarkers (CIViC)

# Version 0.8.1 (June 2nd 2022)

* Data updates: NCI Thesaurus (22.05e)

