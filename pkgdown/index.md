
<br>

## oncoPharmaDB - Targeted and non-targeted anticancer drugs and drug regimens


This R package provides a dataset and method to query targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. 

The main dataset is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)), and where we have limited the associations to cancer-relevant indications only (as provided in [sigven/oncoPhenoMap](https://github.com/sigven/oncoPhenoMap)). Drug-target associations from the Open Targets Platform have furthermore been integrated with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), where we append non-targeted cancer drugs (chemotherapies etc.) and different drug regimens. 

`oncoPharmaDB` also provides drug biomarkers, as curated in databases such as [CIViC](https://civicdb.org), and [CGI](https://cancergenomeinterpreter.org)

Anti-cancer drugs are currently provided with the following tentative drug categories/types (not necessarily mutually exclusive), indicative of their mechanism-of-action:

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

Currently (as of May 28th 2022), `oncoPharmaDB` is built upon the following 
releases of external databases:

 - Open Targets Platform (2022.04)
 - ChEMBL (v30)
 - NCI Thesaurus (22.04d)

### Getting started

* [Installation instructions](articles/installation.html)
* [Usage examples](articles/running.html)

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/oncoPharmaDB/blob/main/.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
