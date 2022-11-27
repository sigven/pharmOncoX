&nbsp;

# pharmaOncoX <a href="https://sigven.github.io/pharmaOncoX/"><img src="man/figures/logo.png" align="right" height="130" width="113"/></a>

**pharmaOncoX** provides access to targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. 

The data is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)). Associations retrieved from Open Targets Platform are limited to cancer-relevant indications only (as provided in [sigven/oncoPhenoMap](https://github.com/sigven/oncoPhenoMap)). Drug-target associations from the Open Targets Platform have furthermore been integrated with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), where we append non-targeted cancer drugs (chemotherapies etc.), and various drug regimens.
Furthermore, we provide anti-cancer drugs in **pharmaOncoX** with the following tentative drug categories/types (not necessarily mutually exclusive), indicative of their mechanism-of-action:

* Alkylating agents
* Angiogenesis inhibitors
* Anthracyclines
* Antimetabolites
* AR antagonists
* BET inhibitors
* Hedgehog antagonists
* HDAC inhibitors
* Hormone therapies
* Immune checkpoint inhibitors
* Kinase inhibitors
* Monoclonal antibodies
* PARP inhibitors
* Platinum compounds
* Proteasome inhibitors
* Topoisomerase inhibitors
* Tubulin inhibitors

Currently (as of November 27th 2022), `pharmaOncoX` is built upon the following 
releases of external databases:

 - Open Targets Platform (2022.11)
 - ChEMBL (v31)
 - NCI Thesaurus (22.10e)

### Getting started

* [Installation instructions](articles/pharmaOncoX.html#installation)
* [Usage examples](articles/pharmaOncoX.html#retrieval-of-drugs---examples)

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/pharmaOncoX/blob/main/.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
