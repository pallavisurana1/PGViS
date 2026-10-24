# PGViS 0.2.0

## New score

The PGViS score was rebuilt after benchmarking 684 scoring regimes, to remove
two sources of inflated separation between tumour and reference cohorts: the
cumulative / mean-aggregated metrics followed the length of each cohort's
variant list, and the allele-frequency shift (delta AF) was both the
variant-selection criterion and a score component.

* The default score is now, for a dataset with variant list U and a sample
  carrying the variants C:

  ```
  p_fv    = #{u in U : x_fu <= x_fv} / |U|,  f in {|delta Wi|, LOR}    (1)
  s_v     = (p_|dWi|,v * p_LOR,v)^(1/2)                                 (2)
  S_i     = sqrt( mean_{v in C_i} s_v^2 )                               (3)
  PGViS_i = (S_i - mean_ref) / SD_ref                                   (4)
  ```

  with `mean_ref` and `SD_ref` taken over the reference cohort (1000 Genomes)
  scored the same way, within superpopulation.
* delta AF is no longer part of the score; it only selects variants.
* New functions: `pgvis_variant_score()` (Eq. 1-2) and `pgvis_sample_scores()`
  (Eq. 3). `pgvis_score()` now standardises the raw score `pgvis_raw`
  (Eq. 4) and accepts `reference =` to place a new cohort on a frozen scale.
  It returns the reference mean and SD as attribute `"reference"`.
* `pgvis_run()` and `pgvis_example_run()` use the new score. The run manifest
  records the method and reference statistics instead of PC1 weights.
* `pgvis_bind_panel()` stacks tables carrying `pgvis_raw` by default.

## Behaviour change

`pgvis_score(df)` now needs the `pgvis_raw` column produced by
`pgvis_sample_scores()`. The previous PCA-weighted composite remains available
as `pgvis_score(df, method = "pca")` (with `pgvis_aggregate_samples()` and
`pgvis_weights()`), so earlier results can be reproduced. Scores from the two
methods are not on the same scale and should not be compared.

# PGViS 0.1.0

First release.

## Pipeline

* A single configuration object, `pgvis_config()`, holds every path,
  threshold and tool location for a project and lays out a stable directory
  tree beneath a project root.
* Six parameterised shell stages (`inst/bash/`) take aligned reads through
  merging, GATK variant calling, joint genotyping, BED conversion with
  per-variant carrier lists, and intersection against regulatory regions and
  dbSNP. Every stage is restartable and writes atomically.
* `pgvis_load_variants()` reads the intersected output into one tidy table,
  with column layouts defined in `pgvis_schema()` and population frequencies
  resolved from the dbSnp155 track by `pgvis_resolve_population_af()`.
* Model scoring is optional and pluggable: the `deepvregulome` Python package
  (`pgvis_predict_python()`), a local DNABERT checkout
  (`pgvis_predict_shell()`), or a prediction table of your own.
* `pgvis_score()` produces the ancestry-referenced composite score, with
  PC1-derived weights that can be frozen and reused across cohorts.

## Reference data

* `pgvis_reference_urls()`, `pgvis_prepare_dbsnp()` and
  `pgvis_check_references()` document the hg38 / dbSnp155 resources and prepare
  the per-chromosome dbSNP files.

## Reporting

* `pgvis_qc()`, `pgvis_manifest()`, `pgvis_compare_groups()`,
  `pgvis_survival()` and plotting helpers for density and survival curves.

## Documentation

* Eight vignettes: an overview, a cluster guide, and a six-part step-by-step
  walkthrough. Every runnable example executes against a bundled synthetic
  fixture (`pgvis_example()`).
