# Get curated cancer biomarker datasets

Downloads preprocessed datasets to a local cache directory and returns a
curated set of genomic biomarkers from multiple sources (CIViC, CGI,
MitelmanDB)

The dataset comes as a `list` object, with three elements:

- `metadata` - a data frame with metadata regarding drug resources used

- `data` - a list with four elements
  ('civic','cgi','mitelmandb','custom_fusions')

- `fpath` - path to cache file

## Usage

``` r
get_biomarkers(cache_dir = NA, force_download = F)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should force downloaded (i.e. set to
  TRUE to re-download even if data exists in cache)

## Value

Each entry of the source-specific (e.g. 'civic') entry in the `data`
list contains a list of three data frames:

- *variant* - list of all biomarker variants, extensively populated
  according to variant aliases (identifer - column **variant_id**)

- *clinical* - cross-references between variants recorded in the
  `variant` data frame and clinical evidence items (identifier - column
  **evidence_id**) and underlying literature evidence (identifier -
  column **source_id**)

- *literature* - lists literature for all source_id's listed in the
  `clinical` data frame
