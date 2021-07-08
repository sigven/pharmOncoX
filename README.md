### oncoPharmaDB - dataset with targeted and non-targeted anticancer compounds and drug regimens

#### Overview

This R package provides a dataset of targeted and non-targeted cancer drugs, including comprehensive annotations per target, drug mechanism-of-action, approval dates, clinical trial phases for various indications etc. The dataset is largely based on drug-target-indication associations provided by the [Open Targets Platform](https://targetvalidation.org) ([Ochoa et al., Nucleic Acids Res., 2021](https://doi.org/10.1093/nar/gkaa1027)), and where we have limited the associations to cancer-relevant indications only (as provided in [sigven/oncoPhenoMap](https://github.com/sigven/oncoPhenoMap)). Drug-target associations from the Open Targets Platform have furthermore been integrated with drug information from [NCI Thesaurus](https://ncithesaurus.nci.nih.gov/ncitbrowser/), where we append non-targeted cancer drugs (chemotherapies etc.) and different drug regimens.

Currently (as of July 2021), the following versions are used to create the mapping:

 - Open Targets Platform (2021.06)
 - ChEMBL (v28)
 - NCI Thesaurus (21.06e)


#### Installation & Usage

##### Install from GitHub

`
install.packages("devtools"); devtools::install_github("sigven/oncoPharmaDB")
`

[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)


##### Usage examples

1. Get BRAF-targeted drugs, list records per indication

	`drugs <- oncoPharmaDB::get_drug(drug_is_targeted = T,
	drug_target = c('BRAF'))`

2. Get _approved_ BRAF-targeted drugs, list records per indication

	`drugs <- oncoPharmaDB::get_drug(drug_is_targeted = T,
	drug_target = c('BRAF'), drug_is_approved = T)`

3. Get BRAF-targeted drugs, list records per indication and drug synonym

	`drugs <- oncoPharmaDB::get_drug(drug_is_targeted = T,
	drug_target = c('BRAF'), list_per_drug_synonym = T)`

4. Get BRAF-targeted drugs, Open Targets Platform only, list per drug only

	`drugs <- oncoPharmaDB::get_drug(drug_is_targeted = T,
	drug_target = c('BRAF'), source_opentargets_only = T, list_per_drug_only = T)`

#### Contact

sigven AT ifi.uio.no
