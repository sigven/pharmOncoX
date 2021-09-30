### oncoPharmaDB - Targeted and non-targeted anticancer drugs and drug regimens

#### Overview

This R package provides a dataset and method to query targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. The dataset is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)), and where we have limited the associations to cancer-relevant indications only (as provided in [sigven/oncoPhenoMap](https://github.com/sigven/oncoPhenoMap)). Drug-target associations from the Open Targets Platform have furthermore been integrated with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), where we append non-targeted cancer drugs (chemotherapies etc.) and different drug regimens. 

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
* Proteasome inhibitors
* Topoisomerase inhibitors
* Tubulin inhibitors

Currently (as of September 30th 2021), the following versions are used to create the mapping:

 - Open Targets Platform (2021.09)
 - ChEMBL (v29)
 - NCI Thesaurus (21.09d)


#### Installation & Usage

##### Install from GitHub

`
install.packages("devtools"); devtools::install_github("sigven/oncoPharmaDB")
`

[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)


##### Usage examples

1. Get BRAF-targeted drugs, list records per indication

	`drugs <- oncoPharmaDB::get_onco_drugs(drug_is_targeted = T,
	drug_target = c('BRAF'))`

2. Get _approved_ BRAF-targeted drugs, list records per indication

	`drugs <- oncoPharmaDB::get_onco_drugs(drug_is_targeted = T,
	drug_target = c('BRAF'), drug_is_approved = T)`

3. Get BRAF-targeted drugs, list records per indication and drug synonym

	`drugs <- oncoPharmaDB::get_onco_drugs(drug_is_targeted = T,
	drug_target = c('BRAF'), list_per_drug_synonym = T)`

4. Get BRAF-targeted drugs, Open Targets Platform only, list per drug only

	`drugs <- oncoPharmaDB::get_onco_drugs(drug_is_targeted = T,
	drug_target = c('BRAF'), source_opentargets_only = T, output_resolution = "drug" )`
	
5. Get BRAF-targeted drugs, Open Targets Platform only, list per drug only, show key annotations only

	`drugs <- oncoPharmaDB::get_onco_drugs(drug_is_targeted = T,
	drug_target = c('BRAF'), source_opentargets_only = T, output_resolution = "drug",
	output_style = 'narrow')`
	

6. Get immune checkpoint inhibitors, list per drug-target entry

   `drugs <- oncoPharmaDB::get_onco_drugs(is_immune_checkpoint_inhibitor = T,
   output_resolution = "drug2target", output_style = "narrow")`
   
7. Get immune checkpoint inhibitors indicated for tumor subtypes within "Colon/Rectum", list per drug-target entry

   `drugs <- oncoPharmaDB::get_onco_drugs(is_immune_checkpoint_inhibitor = T,
   output_resolution = "drug2target", drug_indication_main = "Colon/Rectum", output_style = "narrow")`
   
8. Get antimetabolites

   `drugs <- oncoPharmaDB::get_onco_drugs(is_antimetabolite = T,
   output_resolution = "drug")`


### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/oncoPharmaDB/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
