# Getting started

  
  

## Installation

``` r

if (!("remotes" %in% installed.packages())) {
  install.packages("remotes")
}

remotes::install_github('sigven/pharmOncoX')
#> Using github PAT from envvar GITHUB_PAT. Use `gitcreds::gitcreds_set()` and unset GITHUB_PAT in .Renviron (or elsewhere) if you want to use the more secure git credential store instead.
#> Downloading GitHub repo sigven/pharmOncoX@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpxDBIag/remotes1d6937e51c16/sigven-pharmOncoX-76f98ae/DESCRIPTION’ ... OK
#> * preparing ‘pharmOncoX’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> Omitted ‘LazyData’ from DESCRIPTION
#> * building ‘pharmOncoX_2.1.3.tar.gz’
#> Installing package into '/home/runner/work/_temp/Library'
#> (as 'lib' is unspecified)

library(pharmOncoX)
cache_dir <- tempdir()
```

  
  

## Cancer drug classifications

### Targeted agents

Plotted below are key statistics with respect to drug classification
numbers for targeted and chemotherapy agents found in *pharmOncoX*.
Existing drug classifications have been retrieved from
[ATC](https://www.whocc.no/), and these have been extended with manual
addition/curation, and also through the establishment of multiple novel
levels in the ATC tree, particularly for targeted therapies. Note that
only drugs that are indicated for cancer conditions (as harvested from
the Open Targets platform) are considered in the numbers plotted below.

  

``` r
p_targeted_classifications
```

![](pharmOncoX_files/figure-html/p_targeted-1.png)

### Chemotherapy agents

``` r
p_chemo_classifications
```

![](pharmOncoX_files/figure-html/chemo_plot-1.png)

## Retrieval of drugs - examples

### Get BRAF-targeted drugs, list records per indication

``` r

drugs <- get_drugs(
  cache_dir = cache_dir,
  treatment_category =  c("targeted_therapy_classified",
                          "targeted_therapy_unclassified"),
    drug_target = c('BRAF'))

## Number of drug records
nrow(drugs$records)
#> [1] 273

## Column names of drug records
colnames(drugs$records)
#>  [1] "drug_id"                      "drug_name"                   
#>  [3] "drug_type"                    "molecule_chembl_id"          
#>  [5] "drug_action_type"             "drug_alias"                  
#>  [7] "nci_concept_definition"       "opentargets"                 
#>  [9] "drug_cancer_relevance"        "inhibition_moa"              
#> [11] "is_salt"                      "is_adc"                      
#> [13] "drug_blackbox_warning"        "nci_t"                       
#> [15] "target_symbol"                "target_entrezgene"           
#> [17] "target_genename"              "target_ensembl_gene_id"      
#> [19] "target_type"                  "drug_max_phase_indication"   
#> [21] "drug_approved_indication"     "drug_frac_cancer_indications"
#> [23] "drug_approved_noncancer"      "drug_n_indications"          
#> [25] "drug_year_first_approval"     "drug_max_ct_phase"           
#> [27] "disease_efo_id"               "disease_efo_label"           
#> [29] "primary_site"                 "drug_clinical_id"            
#> [31] "drug_clinical_source"         "atc_code_level1"             
#> [33] "atc_level1"                   "atc_code_level2"             
#> [35] "atc_level2"                   "atc_code_level3"             
#> [37] "atc_level3"                   "atc_treatment_category"
```

  
  

### Get RAS-targeted drugs, list per drug target

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category =  c("targeted_therapy_classified"),
  drug_action_inhibition = T,
  output_resolution = "drug2target")$records |>
  dplyr::filter(atc_level3 == "RAS inhibitors")

drugs <- drugs |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    ) |>
  dplyr::select(
    drug_id, drug_name, drug_type, molecule_chembl_id,
    drug_action_type, target_symbol, dplyr::everything()
  )

dt_drugtable_ras_inhibitors <- DT::datatable(
  drugs,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get MEK inhibitors, list per drug only

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category =  c("targeted_therapy_classified"),
  drug_action_inhibition = T,
  drug_source_opentargets = T, 
  output_resolution = "drug" )$records |>
  dplyr::filter(atc_level3 == "MEK inhibitors")

drugs <- drugs |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    )

dt_drugtable_mek_inhibitors <- DT::datatable(
  drugs,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get immune checkpoint inhibitors, list per drug target

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category =  c("targeted_therapy_classified"),
  drug_source_opentargets = F,
  drug_classified_cancer = T,
   output_resolution = "drug2target")

drugs$records <- drugs$records |>
  dplyr::filter(
    (!is.na(atc_level3) & 
       (atc_level3 == "PD-1/PDL-1 inhibitors" |
         atc_level3 == "Other ICIs - LAG3 inhibitors" | 
         atc_level3 == "Other ICIs - TIGIT inhibitors" | 
         atc_level3 == "Other ICIs - CTLA4 inhibitors" |
     atc_level3 == "Other ICIs")
  )) |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    ) |>
  dplyr::select(
    drug_id,
    drug_name,
    drug_type,
    target_symbol,
    target_genename,
    dplyr::everything()
  )

dt_drugtable_ici <- DT::datatable(
  drugs$records,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get immune checkpoint inhibitors indicated for lung cancer conditions, list per drug-target entry

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  output_resolution = "drug2target", 
  treatment_category = c("targeted_therapy_classified"),
  drug_source_opentargets = T,
  drug_indication_main = "Lung")

drugs$records <- drugs$records |>
  dplyr::filter(
    (!is.na(atc_level3) & 
      (atc_level3 == "PD-1/PDL-1 inhibitors" |
         atc_level3 == "Other ICIs - LAG3 inhibitors" | 
         atc_level3 == "Other ICIs - TIGIT inhibitors" | 
         atc_level3 == "Other ICIs - CTLA4 inhibitors" |
     atc_level3 == "Other ICIs")
  )) |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    ) |>
  dplyr::select(
    drug_id,
    drug_name,
    drug_type,
    target_symbol,
    target_genename,
    dplyr::everything()
  )

dt_drugtable_ici_lung <- DT::datatable(
  drugs$records,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get antimetabolite drugs

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category = c("chemo_therapy_classified"),
  output_resolution = "drug")

drugs$records <- drugs$records |>
  dplyr::filter(
    !is.na(atc_level2) &
    stringr::str_detect(
      atc_level2, "ANTIMETABOLITES"
    )
  ) |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    )

dt_drugtable_metabolites <- DT::datatable(
  drugs$records,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get taxanes

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category = "chemo_therapy_classified",
  output_resolution = "drug")

drugs$records <- drugs$records |>
  dplyr::filter(
    stringr::str_detect(
      atc_level3, "Taxanes"
    )
  ) |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    )

dt_drugtable_taxanes <- DT::datatable(
  drugs$records,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

### Get platinum compounds

``` r

drugs <- get_drugs(
  cache_dir = cache_dir, 
  treatment_category = "chemo_therapy_classified",
  output_resolution = "drug")

drugs$records <- drugs$records |>
  dplyr::filter(
    stringr::str_detect(
      atc_level3, "Platinum compounds"
    )
  ) |>
  dplyr::select(
    -c("drug_alias",
       "disease_main_group",
       "drug_clinical_id")) |>
  dplyr::mutate(
    disease_indication = stringr::str_replace_all(
      disease_indication, "\\|",", ")
    ) |>
  dplyr::select(
    drug_id, drug_name, drug_type, molecule_chembl_id,
    drug_action_type, opentargets, dplyr::everything()
  )

dt_drugtable_platins <- DT::datatable(
  drugs$records,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

## Retrieval of biomarkers

### Reported associations between BRCA1/2 alterations and drug sensitivity

- Get evidence from [CIViC](https://civicdb.org) and
  [CGI](https://www.cancergenomeinterpreter.org/biomarkers) for cancer
  drug sensitivity of BRCA1/2 alterations (somatically (tumor) or
  inherited/germline)

``` r

biomarkers <- get_biomarkers(
  cache_dir = cache_dir)

brca1_biomarkers <- list()
for(source in c('civic','cgi')){
  brca1_biomarkers[[source]] <- 
    biomarkers$data[[source]]$variant |> 
    dplyr::filter(
      !is.na(symbol) & (symbol == "BRCA1" | symbol == "BRCA2")) |>
    dplyr::group_by(variant_id, variant_name_primary, variant_consequence) |>
    dplyr::summarise(
      variant_alias = paste(variant_alias, collapse=", "), 
      .groups = "drop") |> 
    dplyr::inner_join(
      biomarkers$data[[source]]$clinical, by = "variant_id") |>
    dplyr::select(
      variant_id, variant_name_primary, therapeutic_context,
      evidence_type, evidence_level, 
      biomarker_source, biomarker_source_datestamp,
      molecular_profile_name, evidence_id, variant_origin,
      primary_site, evidence_id, source_id,
      evidence_url, 
      evidence_description, 
      clinical_significance) |> 
    dplyr::distinct() |> 
    dplyr::rename(literature_id = source_id,
                  variant_name = variant_name_primary) |>
    dplyr::filter(evidence_type == "Predictive") |>
    dplyr::select(
      variant_name, 
      primary_site,
      therapeutic_context,
      molecular_profile_name, 
      evidence_level,
      dplyr::everything()
    )
}

brca1_biomarkers_all <- 
  dplyr::bind_rows(brca1_biomarkers[['civic']],
                   brca1_biomarkers[['cgi']]) |>
    dplyr::arrange(evidence_level)
  

dt_brca1_biomarkers <- DT::datatable(
  brca1_biomarkers_all,
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip")
)
```

  

  
  

## Session Info

``` r
# set eval = FALSE if you don't want this info (useful for reproducibility) 
# to appear
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] pharmOncoX_2.1.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10        generics_0.1.4     stringi_1.8.7      digest_0.6.37     
#>  [5] magrittr_2.0.4     RColorBrewer_1.1-3 evaluate_1.0.5     grid_4.5.2        
#>  [9] fastmap_1.2.0      jsonlite_2.0.0     processx_3.8.6     pkgbuild_1.4.8    
#> [13] googledrive_2.1.2  ps_1.9.1           httr_1.4.7         purrr_1.2.0       
#> [17] crosstalk_1.2.2    scales_1.4.0       textshaping_1.0.4  jquerylib_0.1.4   
#> [21] cli_3.6.5          rlang_1.1.6        crayon_1.5.3       remotes_2.5.0     
#> [25] withr_3.0.2        cachem_1.1.0       yaml_2.3.10        tools_4.5.2       
#> [29] gargle_1.6.0       dplyr_1.1.4        ggplot2_4.0.0      DT_0.34.0         
#> [33] curl_7.0.0         assertthat_0.2.1   vctrs_0.6.5        R6_2.6.1          
#> [37] lifecycle_1.0.4    stringr_1.6.0      fs_1.6.6           htmlwidgets_1.6.4 
#> [41] ragg_1.5.0         pkgconfig_2.0.3    desc_1.4.3         callr_3.7.6       
#> [45] pkgdown_2.2.0      pillar_1.11.1      bslib_0.9.0        gtable_0.3.6      
#> [49] glue_1.8.0         lgr_0.5.0          systemfonts_1.3.1  xfun_0.54         
#> [53] tibble_3.3.0       tidyselect_1.2.1   knitr_1.50         farver_2.1.2      
#> [57] htmltools_0.5.8.1  rmarkdown_2.30     compiler_4.5.2     S7_0.2.0
```

  
  

## References

Griffith, Malachi, Nicholas C Spies, Kilannin Krysiak, Joshua F
McMichael, Adam C Coffman, Arpad M Danos, Benjamin J Ainscough, et al.
2017. “CIViC Is a Community Knowledgebase for Expert Crowdsourcing the
Clinical Interpretation of Variants in Cancer.” *Nat. Genet.* 49 (2):
170–74. <http://dx.doi.org/10.1038/ng.3774>.

Kim, Sunghwan, Jie Chen, Tiejun Cheng, Asta Gindulyte, Jia He, Siqian
He, Qingliang Li, et al. 2021. “PubChem in 2021: New Data Content and
Improved Web Interfaces.” *Nucleic Acids Res.* 49 (D1): D1388–95.
<http://dx.doi.org/10.1093/nar/gkaa971>.

Nakken, Sigve. 2023. *phenOncoX: A Phenotype Ontology Map for Cancer*.
<https://github.com/sigven/phenOncoX>.

Ochoa, David, Andrew Hercules, Miguel Carmona, Daniel Suveges, Asier
Gonzalez-Uriarte, Cinzia Malangone, Alfredo Miranda, et al. 2021. “Open
Targets Platform: Supporting Systematic Drug-Target Identification and
Prioritisation.” *Nucleic Acids Res.* 49 (D1): D1302–10.
<http://dx.doi.org/10.1093/nar/gkaa1027>.

Sioutos, Nicholas, Sherri de Coronado, Margaret W Haber, Frank W Hartel,
Wen-Ling Shaiu, and Lawrence W Wright. 2007. “NCI Thesaurus: A Semantic
Model Integrating Cancer-Related Clinical and Molecular Information.”
*J. Biomed. Inform.* 40 (1): 30–43.
<http://dx.doi.org/10.1016/j.jbi.2006.02.013>.

Tamborero, David, Carlota Rubio-Perez, Jordi Deu-Pons, Michael P
Schroeder, Ana Vivancos, Ana Rovira, Ignasi Tusquets, et al. 2018.
“Cancer Genome Interpreter Annotates the Biological and Clinical
Relevance of Tumor Alterations.” *Genome Med.* 10 (1): 25.
<http://dx.doi.org/10.1186/s13073-018-0531-8>.
