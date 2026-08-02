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
