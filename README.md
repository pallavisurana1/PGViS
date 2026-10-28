# PGViS

**P**ersonal **G**ermline **V**ariant **i**mpact **S**core for
non-coding regulatory regions.

Whole-genome sequencing turns up a huge number of germline variants
outside protein-coding sequence. Most are harmless. A few land where it
matters — the splice sites that mark exon boundaries, or the short
motifs where transcription factors bind — and quietly change how a gene
is regulated. PGViS finds those few and turns them into a single score
per person.

The pipeline runs end to end: call germline variants, keep the ones
inside splice donor/acceptor sites and transcription-factor binding
sites, ask a fine-tuned DNABERT model whether each one actually breaks
the element, weigh that against how the variant's frequency shifts
between the cohort and a matched reference population, and collapse
everything a person carries into one ancestry-referenced number.

![PGViS graphical abstract: heritable non-coding variants across lung
cancer genomes are called from WGS, intersected with splice sites and
transcription-factor binding sites, scored with DNABERT models and by
allele-frequency shifts, and combined into a cumulative per-patient
PGViS risk score.](man/figures/graphical-abstract1.png)

```         
BAM ─▶ GVCF ─▶ joint VCF ─▶ BED + carriers ─▶ ∩ regions ─▶ ∩ dbSNP
                                                              │
                                                              ▼
    PGViS ◀── per-sample ◀── metrics ◀── model probs ◀── ref/alt sequences
   (score)     aggregate     LOR, Δscore,        (DNABERT)
                             ΔWi, ΔAF
```

## How it fits together

The pipeline has two halves. Everything upstream of the variant table —
read merging, variant calling, joint genotyping, the bedtools and
bcftools work — runs as shell scripts, because that is where those tools
live and because it wants a scheduler and real memory. Everything
downstream, once you have a table of variants sitting inside regulatory
regions, is R.

The DNABERT step is the heart of the method, and PGViS keeps it at arm's
length on purpose: the package builds the sequences a model needs, hands
them off, and reads the probabilities back. Run them with the
[`deepvregulome`](https://pypi.org/project/deepvregulome/) Python
package, with a local DNABERT checkout, or with any other tool — if you
already have probabilities, PGViS takes them as they are.

## Install

``` r
# install.packages("remotes")
remotes::install_github("pallavisurana1/PGViS")
```

The shell stages need `gatk`, `bcftools`, `samtools`, `bedtools`,
`tabix` and `awk` on your `PATH`. The scoring step, if you run the
models through R, needs:

``` bash
pip install deepvregulome
```

The fine-tuned models — donor and acceptor splice classifiers, plus
transcription-factor and histone-mark models built on GENCODE and ENCODE
annotations, one per target — are described and hosted here, and cached
locally the first time you score:
[duttaprat/DeepVRegulome](https://huggingface.co/duttaprat/DeepVRegulome).

## A quick walkthrough

The package ships a small synthetic dataset, so you can run the whole R
side before you have anything of your own:

``` r
library(PGViS)

res <- pgvis_example_run()
head(res$scores)
attr(res$scores, "reference") # reference mean / SD per ancestry stratum
```

Or look at the pieces one at a time:

``` r
variants <- pgvis_example("variants")     # variants inside regulatory regions
preds    <- pgvis_example("predictions")  # matching model probabilities
clin     <- pgvis_example("clinical")

pgvis_qc(variants)                         # a quick health check on the table
```

## A real run

``` r
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

# upstream — each call returns a ready-to-submit script
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

Leave `submit = FALSE` (the default) and each call writes its shell
script and stops, so you can read it before it runs. The script that
executes is a plain file you can open, edit, and resubmit — nothing is
hidden.

Stages are restartable. Anything already finished is skipped, and each
stage writes to a temporary file that only moves into place on success,
so a job killed by a walltime limit just needs resubmitting and picks up
where it left off.

## Running the shell stages

The upstream stages are shell, but you don't write shell to use them.
Each R wrapper (`pgvis_merge_bams`, `pgvis_call_gvcf`, and so on) takes
your config and writes a ready-to-run script into `<root>/logs/`, with
your paths, thread counts and — on SLURM — the `#SBATCH` header filled
in. What runs is a plain file you can read.

**The sample sheet.** The two calling stages take a tab-separated file
with three columns and no header: cohort, tissue, sample.

```         
luad    blood   PT-001
luad    solid   PT-001
lusc    blood   PT-002
```

`pgvis_sample_sheet()` validates one before you commit a cluster to it,
and warns about samples with no matching BAM.

**Inspect, then submit.** With the default `submit = FALSE`, the wrapper
writes the script and returns its path — nothing runs:

``` r
script <- pgvis_call_gvcf(cfg, "samples.tsv")   # submit = FALSE
cat(readLines(script), sep = "\n")              # read it first
```

The rendered script opens with a generated block you can edit, then the
stage body:

``` bash
#!/usr/bin/env bash
#SBATCH --job-name=pgvis
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-24
set -euo pipefail

PGVIS_ROOT='/scratch/me/lung_pgvis'
REF_FASTA='/refs/hg38.analysisSet.fa'
SAMPLE_SHEET='samples.tsv'
CHROMS='1 2 3 ... 22 X Y'
# ---- end generated block ----
# Stage 02 -- per-sample, per-chromosome GVCF calling with GATK HaplotypeCaller.
...
```

When it looks right, either set `submit = TRUE` to send it to the
scheduler, or run the file yourself:

``` bash
sbatch /scratch/me/lung_pgvis/logs/02_call_gvcf.run.sh   # SLURM
bash   /scratch/me/lung_pgvis/logs/02_call_gvcf.run.sh   # no scheduler
```

**Without R at all.** The stage scripts live in the installed package
and read their settings from environment variables, so you can run one
straight from the shell if you'd rather not touch R for the upstream
work:

``` bash
STAGE=$(Rscript -e 'cat(system.file("bash", package="PGViS"))')

export PGVIS_ROOT=/scratch/me/lung_pgvis
export REF_FASTA=/refs/hg38.analysisSet.fa
export SAMPLE_SHEET=samples.tsv
export PGVIS_BASH_LIB="$STAGE"

bash "$STAGE/02_call_gvcf.sh"
```

`pgvis_stages()` lists every stage, what it does, and where its script
is. `set PGVIS_FORCE=1` to redo work that's already finished; by default
anything complete is skipped so an interrupted job just needs
resubmitting.

## Scoring with your own probabilities

If a different tool produced the model scores, hand PGViS the sequences
it needs and rejoin afterward:

``` r
seqs  <- pgvis_build_sequences(variants)
dirs  <- pgvis_write_model_input(seqs, cfg$dirs$seq)
# score `dirs` however you like -> one predictions.tsv per feature,
# with columns: func_site, var_site, p_ref, p_alt
preds <- pgvis_read_predictions(cfg$dirs$preds)
mets  <- pgvis_variant_metrics(seqs, preds, cfg)
```

## What the score is made of

Each variant inside a regulatory element gets three signals:

| Signal         | Function           | What it captures                                                                   |
|------------------------|------------------------|------------------------|
| Log-odds ratio | `log_odds_ratio()` | how far the model's confidence in the element drops when the variant is introduced |
| Score change   | `score_change()`   | the probability shift, weighted so low-confidence changes count for less           |
| ΔWi            | `pgvis_delta_wi()` | how far the cohort allele frequency diverges from the matched population frequency |
| ΔAF            | `pgvis_delta_wi()` | the raw cohort-minus-population alternate frequency                                |

The score itself uses two of them, $|\Delta Wi|$ and the log-odds ratio.
For a dataset with variant list $U$ and a person carrying the variants
$C_i$:

$$p_{fv} = \frac{\#\{u \in U : x_{fu} \le x_{fv}\}}{|U|}, \qquad f \in \{|\Delta Wi|, \mathrm{LOR}\} \tag{1}$$

$$s_v = \left(p_{|\Delta Wi|,v}\; p_{\mathrm{LOR},v}\right)^{1/2} \tag{2}$$

$$S_i = \sqrt{\frac{1}{|C_i|}\sum_{v \in C_i} s_v^{2}} \tag{3}$$

$$\mathrm{PGViS}_i = \frac{S_i - \mathrm{mean}_{\mathrm{ref}}}{\mathrm{SD}_{\mathrm{ref}}} \tag{4}$$

Each feature is ranked within its own dataset's variant list (1), the
two ranks are combined by geometric mean (2), a person's score is the
root mean square over the variants they carry (3), and it is
standardised against the reference cohort scored the same way (4):

``` r
long   <- pgvis_variant_score(long)    # Eq. 1-2
person <- pgvis_sample_scores(long)    # Eq. 3
scores <- pgvis_score(person)          # Eq. 4
```

**ΔAF is not in the score.** It is the criterion used to *select*
variants (|ΔAF| ≥ 0.1), so tumour ΔAF is enriched by construction while
the reference panel's is near zero; scoring on it again makes the
comparison circular. The cumulative sum is avoided too, because it
tracks how many variants a cohort's list contains rather than anything
about the person. The root mean square is only weakly related to variant
count.

Because ranks are taken within each list, the score reflects *which*
variants of the list a person carries (high-ranking or low-ranking), not
absolute feature values. Tumour and reference lists are built
separately, so read the score as moderate separation, not as a
diagnostic.

The pre-0.2.0 PCA-weighted composite is still available as
`pgvis_score(df, method = "pca")`.

Two things are worth knowing before you read a score:

**Ancestry matching earns its keep.** How far someone's allele
frequencies sit from a reference panel depends heavily on how
genetically distant they are from it, so a score standardised against
the whole panel mostly reports ancestry. PGViS standardises *within*
superpopulation to remove that. `strata = NULL` turns it off for a
genuinely homogeneous cohort, but that's the wrong default anywhere
else.

**The reference is refit unless you freeze it.** `pgvis_score()` takes
the reference mean and SD from whatever reference samples you give it, so
two cohorts scored separately aren't on the same scale. Carry the
reference across to line them up:

``` r
ref <- attr(discovery, "reference")
val <- pgvis_score(validation, reference = ref)
```

`thresholds$drop_freq_sources` defaults to `"1000G"`: when the reference
panel *is* 1000 Genomes, scoring it against 1000G-derived frequencies is
circular, so those variants are dropped.

## Vignettes

Eight of them, from first setup through survival analysis. The six
numbered ones walk the pipeline a stage at a time, and every runnable
example runs against the bundled dataset.

``` r
browseVignettes("PGViS")
```

-   **PGViS** — a condensed overview
-   **1–6** — setup, variant calling, sequences, model scoring, metrics,
    scoring
-   **pipeline** — the shell side in one place, for cluster runs

## Function map

| Stage             | Functions                                                                                                                          |
|------------------------------------|------------------------------------|
| Setup             | `pgvis_config()` `pgvis_init()` `pgvis_panel_config()` `pgvis_sample_sheet()`                                                      |
| Reference data    | `pgvis_reference_urls()` `pgvis_prepare_dbsnp()` `pgvis_check_references()`                                                        |
| Shell (BAM → BED) | `pgvis_merge_bams()` `pgvis_call_gvcf()` `pgvis_combine_gvcf()` `pgvis_genotype_filter()` `pgvis_vcf_to_bed()` `pgvis_intersect()` |
| Reading           | `pgvis_load_variants()` `pgvis_read_intersect()` `pgvis_schema()` `pgvis_resolve_population_af()`                                  |
| Sequences         | `pgvis_build_sequences()` `mutate_sequence()` `pgvis_write_model_input()`                                                          |
| Model (optional)  | `pgvis_predict_python()` `pgvis_predict_shell()` `pgvis_read_predictions()`                                                        |
| Metrics           | `log_odds_ratio()` `score_change()` `pgvis_delta_wi()` `pgvis_variant_metrics()` `pgvis_filter_disruptive()`                       |
| Cohorts           | `pgvis_combine_results()` `pgvis_parse_label()` `pgvis_filter_replicated()` `pgvis_summarise_by()`                                 |
| Score             | `pgvis_explode_carriers()` `pgvis_variant_score()` `pgvis_sample_scores()` `pgvis_bind_panel()` `pgvis_score()`                    |
| Legacy score      | `pgvis_aggregate_samples()` `pgvis_weights()` (`pgvis_score(method = "pca")`)                                                      |
| Reporting         | `pgvis_plot_density()` `pgvis_compare_groups()` `pgvis_survival()` `pgvis_plot_survival()` `pgvis_qc()` `pgvis_manifest()`         |
| End to end        | `pgvis_run()` `pgvis_example()` `pgvis_example_run()`                                                                              |

## Reference data

PGViS is built around the hg38 assembly and UCSC's dbSnp155 track. The
download locations and a helper to prepare the per-chromosome dbSNP
files are built in:

-   dbSnp155 —
    <http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb>
-   hg38 — <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/>

``` r
pgvis_reference_urls("hg38")
pgvis_prepare_dbsnp("/refs/dbsnp155/by_chr", subset = "common", dry_run = TRUE)
```

The population-frequency parsing follows the layout of UCSC's dbSnp155
track. A different build, or the raw NCBI VCF, arranges the per-project
frequencies differently — `pgvis_qc()` flags it when that happens, the
tell being most variants coming back with no resolved frequency.

## Citation
Surana et. al., (2026) bioRxiv https://www.biorxiv.org/content/10.64898/2026.08.01.742250v1

## License

MIT.
