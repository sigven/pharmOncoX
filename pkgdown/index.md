&nbsp;

# pharmOncoX <a href="https://sigven.github.io/pharmOncoX/"><img src="man/figures/logo.png" align="right" height="130" width="113"/></a>

**pharmOncoX** provides access to targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. 

The data is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)). Associations retrieved from Open Targets Platform are integrated with cancer-relevant indications/conditions (as provided in [sigven/phenOncoX](https://github.com/sigven/phenOncoX)), allowing the user to retrieve drugs indicated for main tumor types (e.g. `Lung`, `Colon/Rectum` etc.) 

Drug-target associations from the Open Targets Platform have furthermore been integrated and appended with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), showing also non-targeted cancer drugs (chemotherapeutic agents etc.), and various drug regimens.

We provide anti-cancer drugs in **pharmOncoX** with drug class labeling from [Anatomical Therapeutic Chemical (ATC) Classification System](https://www.whocc.no/atc_ddd_index/), enabling a filtering of drugs according to their main mechanisms of action.

Currently (as of late October 2023), `pharmOncoX` is built upon the following 
releases of external databases:

 - Open Targets Platform (2023.09)
 - ChEMBL (v33)
 - NCI Thesaurus (23.09d)

### Getting started

* [Installation instructions](articles/pharmOncoX.html#installation)
* [Usage examples](articles/pharmOncoX.html#retrieval-of-drugs---examples)

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/pharmOncoX/blob/main/.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
