# Match a list of drug names/aliases to the dataset of drugs provided by pharmOncoX

Match a list of drug names/aliases to the dataset of drugs provided by
pharmOncoX

## Usage

``` r
match_drug_names(
  query = NULL,
  exclude_salt_forms = TRUE,
  exclude_adc = FALSE,
  cache_dir = NA,
  force_download = FALSE
)
```

## Arguments

- query:

  A character vector of drug names/aliases

- exclude_salt_forms:

  Logical indicating if salt forms should be excluded

- exclude_adc:

  Logical indicating if antibody-drug conjugates should be excluded

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should force downloaded
