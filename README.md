# PGViS

**P**ersonal **G**ermline **V**ariant **i**mpact **S**core for non-coding
regulatory regions.

Whole-genome sequencing turns up an enormous number of germline variants that
fall outside protein-coding sequence. Most are harmless. A few land in places
that matter — the splice sites that define exon boundaries, the short motifs
where transcription factors bind — and quietly change how a gene is regulated.
PGViS is a pipeline for finding those few and turning them into a single score
per person.

It goes from aligned reads all the way to that score: call germline variants,
keep the ones inside splice donor/acceptor sites and transcription-factor
binding sites, ask a fine-tuned DNABERT model whether each one actually
disrupts the element, weigh that against how the variant's frequency shifts
between the cohort and a matched reference population, and collapse everything
a person carries into one ancestry-referenced number.

![PGViS graphical abstract: non-coding heritable variants across lung cancer genomes are called from WGS, intersected with splice donor/acceptor sites and transcription-factor binding sites, scored with DNABERT models and by shifts in allele frequency, and combined into a cumulative per-patient PGViS risk score.](man/figures/graphical-abstract.png)

```
BAM ─▶ GVCF ─▶ joint VCF ─▶ BED + carriers ─▶ ∩ regions ─▶ ∩ dbSNP
                                                              │
                                                              ▼
    PGViS ◀── per-sample ◀── metrics ◀── model probs ◀── ref/alt sequences
   (score)     aggregate     LOR, Δscore,        (DNABERT, optional)
                             ΔWi, ΔAF
```

## How it's put together

The pipeline splits cleanly in two. The heavy upstream work — read merging,
variant calling, joint genotyping, all the bedtools and bcftools plumbing —
runs as shell, because that is where those tools live and because it wants a
job scheduler and a lot of memory. Everything from "here is a table of
variants sitting inside regulatory regions" onward is R.

The deep-learning step sits in the middle and is **optional**. PGViS prepares
the sequences a model needs, hands them off, and reads the probabilities back.
You can use the [`deepvregulome`](https://pypi.org/project/deepvregulome/)
Python package, drive a local DNABERT checkout, or supply probabilities from
anywhere else. If you already have model scores, the R side never touches
Python at all.

## Install

```r
# install.packages("remotes")
remotes::install_github("pallavisurana1/PGViS")
```

For the shell stages you'll need `gatk`, `bcftools`, `samtools`, `bedtools`,
`tabix` and `awk` on your `PATH`. For the deep-learning step (only if you use
it):

```bash
pip install deepvregulome
```

The fine-tuned models — donor and acceptor splice classifiers plus several
hundred transcription-factor, RNA-binding-protein and histone-mark models, one
per ENCODE target — are pulled from Hugging Face the first time you score and
cached locally:
[duttaprat/DeepVRegulome](https://huggingface.co/duttaprat/DeepVRegulome).

## A five-minute tour, no data required

The package ships a small synthetic dataset so you can watch the whole R-side
flow run before pointing it at anything real:

```r
library(PGViS)

res <- pgvis_example_run()
head(res$scores)
attr(res$scores, "weights")   # how the three signals were weighted
```

The individual pieces, if you want to poke at them:

```r
variants <- pgvis_example("variants")     # variants sitting in regulatory regions
preds    <- pgvis_example("predictions")  # matching model probabilities
clin     <- pgvis_example("clinical")

pgvis_qc(variants)                         # a quick health check on the table
```

None of these numbers mean anything biologically — the fixture exists so the
code runs and so you can see the shape of every table before you have your
own.

## A real run

```r
cfg <- pgvis_config(
  root            = "~/projects/lung_pgvis",
  reference_fasta = "/refs/hg38.analysisSet.fa",
  dbsnp_dir       = "/refs/dbsnp155/by_chr",
  regions         = list(acceptor = "/refs/acceptor.bed",
                         donor    = "/refs/donor.bed"),
  tfbs_dir        = "/refs/encode_tfbs",     # one subdirectory per TF
  cohorts         = c("luad", "lusc"),
  scheduler       = "slurm",
  slurm           = list(partition = "cpu", mem = "48G", time = "2-00:00:00"),
  thresholds      = list(ref_prob_min = 0.7, alt_prob_max = 0.3)
)
pgvis_init(cfg)
pgvis_check_references(cfg)          # confirm the reference data is really there

# upstream: each returns a ready-to-submit script
pgvis_merge_bams(cfg,      "samples.tsv", submit = TRUE)
pgvis_call_gvcf(cfg,       "samples.tsv", submit = TRUE)
pgvis_combine_gvcf(cfg,    submit = TRUE)
pgvis_genotype_filter(cfg, submit = TRUE)
pgvis_vcf_to_bed(cfg,      submit = TRUE)
pgvis_intersect(cfg,       submit = TRUE)

# downstream, in R
variants <- pgvis_load_variants(cfg)
clin     <- read.csv("clinical.csv")   # a `sample` column plus `race` or `superpopulation`

res <- pgvis_run(cfg, variants, clinical = clin,
                 reference_group = "1000Genomes")

head(res$scores)
pgvis_plot_density(res$scores, group_col = "cohort", stratum = "EUR")
```

`submit = FALSE` (the default) writes each shell script and stops, so you can
read it before anything runs. Nothing is a black box — the exact script that
executes is a plain file you can open, edit, and resubmit.

Every stage is restartable. Outputs that already exist and are non-empty are
skipped, and each stage writes to a temporary file that only moves into place
on success. A job killed by a walltime limit is simply resubmitted; it picks
up where it stopped.

## Scoring without the model

If you have probabilities from your own runner, skip straight past the Python:

```r
seqs  <- pgvis_build_sequences(variants)
dirs  <- pgvis_write_model_input(seqs, cfg$dirs$seq)
# score `dirs` however you like -> one predictions.tsv per feature with
# columns: func_site, var_site, p_ref, p_alt
preds <- pgvis_read_predictions(cfg$dirs$preds)
mets  <- pgvis_variant_metrics(seqs, preds, cfg)
```

## What the score is made of

Each variant that sits in a regulatory element gets three numbers:

| Signal | Function | What it captures |
|---|---|---|
| Log-odds ratio | `log_odds_ratio()` | how far the model's confidence in the element drops when the variant is introduced |
| Score change | `score_change()` | the probability shift, weighted so low-confidence changes count for less |
| ΔWi | `pgvis_delta_wi()` | how far the cohort allele frequency diverges from the matched population frequency |
| ΔAF | `pgvis_delta_wi()` | the raw cohort-minus-population alternate frequency |

Those are averaged per person, standardised against a reference cohort
(median and MAD, because the distributions are heavy-tailed), and combined.
The weights come from the first principal component of the standardised
signals, so each is weighted by how much it contributes to the shared axis of
variation rather than by a number chosen by hand.

Two things are worth understanding before you read a PGViS score:

**Ancestry matching is doing real work.** How far a person's allele
frequencies sit from a reference panel depends heavily on how genetically
distant they are from it. Standardise against the whole panel and the score
mostly reflects ancestry. PGViS standardises *within* superpopulation to
strip that out. There is a `strata = NULL` escape hatch for a genuinely
homogeneous cohort, but it is the wrong default for anything else.

**Weights are refit unless you freeze them.** By default `pgvis_score()`
derives the weights from the data in front of it. If you score a discovery
cohort and then a validation cohort separately, the two scores are not on the
same scale. To put them on the same footing, carry the weights across:

```r
w   <- attr(discovery, "weights")
val <- pgvis_score(validation, weights = w)
```

There is also `thresholds$drop_freq_sources`, which defaults to `"1000G"`.
When your reference panel *is* 1000 Genomes, scoring it against 1000G-derived
population frequencies is circular, so those variants are left out.

## Where to read next

Eight vignettes, running from a first project setup to survival analysis.
The six numbered ones walk the pipeline a stage at a time; every runnable
example in them executes against the bundled fixture.

```r
browseVignettes("PGViS")
```

- **PGViS** — the condensed overview
- **1–6** — setup, variant calling, sequences, model scoring, metrics, scoring
- **pipeline** — the whole shell side gathered in one place, for cluster runs

## Function map

| Stage | Functions |
|---|---|
| Setup | `pgvis_config()` `pgvis_init()` `pgvis_panel_config()` `pgvis_sample_sheet()` |
| Reference data | `pgvis_reference_urls()` `pgvis_prepare_dbsnp()` `pgvis_check_references()` |
| Shell (BAM → BED) | `pgvis_merge_bams()` `pgvis_call_gvcf()` `pgvis_combine_gvcf()` `pgvis_genotype_filter()` `pgvis_vcf_to_bed()` `pgvis_intersect()` |
| Reading | `pgvis_load_variants()` `pgvis_read_intersect()` `pgvis_schema()` `pgvis_resolve_population_af()` |
| Sequences | `pgvis_build_sequences()` `mutate_sequence()` `pgvis_write_model_input()` |
| Model (optional) | `pgvis_predict_python()` `pgvis_predict_shell()` `pgvis_read_predictions()` |
| Metrics | `log_odds_ratio()` `score_change()` `pgvis_delta_wi()` `pgvis_variant_metrics()` `pgvis_filter_disruptive()` |
| Cohorts | `pgvis_combine_results()` `pgvis_parse_label()` `pgvis_filter_replicated()` `pgvis_summarise_by()` |
| Score | `pgvis_explode_carriers()` `pgvis_aggregate_samples()` `pgvis_bind_panel()` `pgvis_score()` `pgvis_weights()` |
| Reporting | `pgvis_plot_density()` `pgvis_compare_groups()` `pgvis_survival()` `pgvis_plot_survival()` `pgvis_qc()` `pgvis_manifest()` |
| End to end | `pgvis_run()` `pgvis_example()` `pgvis_example_run()` |

## Reference data

PGViS is built around the hg38 assembly and UCSC's dbSnp155 track. The exact
download locations, and a helper that prepares the per-chromosome dbSNP files,
are built in:

```r
pgvis_reference_urls("hg38")
pgvis_prepare_dbsnp("/refs/dbsnp155/by_chr", subset = "common", dry_run = TRUE)
```

The population-frequency parsing is tuned to the layout of UCSC's dbSnp155
track. A different build, or the raw NCBI VCF, arranges the per-project
frequencies differently, and `pgvis_qc()` will tell you if that has gone wrong
(most variants coming back with no resolved frequency is the tell).

## Citation

The models and the framework behind them:

> Dutta P, Vaska J, Surana P, Sathian R, Chao M, Zhou Z, Liu H, Davuluri RV,
> et al. *DeepVRegulome: DNABERT-based deep-learning framework for predicting
> the functional impact of short genomic variants on the human regulome.*
> arXiv:2511.09026.

## License

MIT.
