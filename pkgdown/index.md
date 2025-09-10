&nbsp;

# pharmOncoX <a href="https://sigven.github.io/pharmOncoX/"><img src="man/figures/logo.png" align="right" height="104" width="90"/></a>

**pharmOncoX** provides access to targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. It also provides access to data on actionable genomic aberrations (i.e. molecular biomarkers), including gene fusions, mutations, copy number alterations, and expression biomarkers.

The data is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)). Associations retrieved from Open Targets Platform are integrated with cancer-relevant indications/conditions (as provided in [sigven/phenOncoX](https://github.com/sigven/phenOncoX)), allowing the user to retrieve drugs indicated for main tumor types (e.g. `Lung`, `Colon/Rectum` etc.) 

Drug-target associations from the Open Targets Platform have furthermore been integrated and appended with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), showing also non-targeted cancer drugs (chemotherapeutic agents etc.), and various drug regimens.

_pharmOncoX_ provides anti-cancer drug classification through existing entries in the [Anatomical Therapeutic Chemical (ATC) Classification System](https://www.whocc.no/atc_ddd_index/), and these have been extended significantly with manual curation, also by establishing novel drug categories that are presently missing in the ATC classificiation tree (examples include _AURK inhibitors_, _MET inhibitors_, _BET inhibitors_, _AKT inhibitors_, _PLK inhibitors_, _IAP inhibitors_, _RAS inhibitors_, _BCL2 inhibitors_ etc.) enabling a filtering of drugs according to their main mechanisms of action.

Currently (as of early September 2025), `pharmOncoX` is built upon the following 
releases of external databases:

 - Open Targets Platform (2025.06)
 - ChEMBL (v34)
 - NCI Thesaurus (25.08d)
 - MitelmanDB (20250710)
 - CIViC (20250910)

### Getting started

* [Installation instructions](articles/pharmOncoX.html#installation)
* [Usage examples](articles/pharmOncoX.html#retrieval-of-drugs---examples)

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/pharmOncoX/blob/main/.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
