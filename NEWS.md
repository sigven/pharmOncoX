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

# Version 0.8.0 (May 31st 2022)

* Updated indication matching with data from https://github.com/sigven/oncoPhenoMap


# Version 0.7.6 (April 22nd 2022)

* Fixed some ambiguous synonym entries (paclitaxel, cisplatin, seribantumab, trastuzumab emtansine)
* Updated indication matching with data from https://github.com/sigven/oncoPhenoMap (v0.3.2)
* Updated biomarker collection (MitelmanDB - April 2022 release) + CIViC (April 22nd 2022)

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


