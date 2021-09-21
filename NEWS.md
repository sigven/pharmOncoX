# Version 0.4.9 (September 22nd 2021)

## Fixed

* Fixed bug related to kinase drug category 
* Added uniprot identifiers to custom targets (experimental drugs, NCI)
* Added imports (tidyverse, assertthat)

# Version 0.4.8 (September 21st 2021)

## Fixed

* Fixed bug related to missing drug indications (cancer types)
* Updated oncoPhenoMap

# Version 0.4.7 (August 23rd 2021)

## Added

* New drug classes: _anthracycline_ and _hedgehog_antagonist_
* Possibility to show extensive or narrow output annotations per drug record (argument _output_style_)
* Improved argument checks

# Version 0.4.6 (August 21st 2021)

## Added

* Altered name for main function: _get_onco_drugs_
* Possibility to show extensive or narrow output annotations per drug record (argument _output_style_)
* Improved argument checks

# Version 0.4.5 (August 19th 2021)

## Added

* Added 14 provisional drug categories, as logical variables pr. drug, e.g.
  _bet_inhibitor_ (TRUE/FALSE), _immune_checkpoint_inhibitor_ (TRUE/FALSE),
  _antimetabolite_ (TRUE/FALSE) etc. Note that categories are not necessarily
  mutually exclusive. 
