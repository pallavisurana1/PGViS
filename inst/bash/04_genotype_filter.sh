#!/usr/bin/env bash
# Stage 04 -- joint genotyping, hard filtering, normalisation, multiallelic split.
# Required: PGVIS_ROOT REF_FASTA COHORTS ; Optional: CHROMS FILTER_EXPR
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_need gatk bcftools
pgvis_require_var PGVIS_ROOT REF_FASTA COHORTS
DIR="${DIR:-$PGVIS_ROOT/03_joint}"
# GATK hard-filter recommendations for germline SNPs. Override for your data.
FILTER_EXPR="${FILTER_EXPR:-QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0}"
read -r -a COHORT_ARR <<< "$COHORTS"
read -r -a CHROM_ARR  <<< "${CHROMS:-$(seq 1 22) X Y}"

tasks=(); for c in "${COHORT_ARR[@]}"; do for k in "${CHROM_ARR[@]}"; do tasks+=("$c $k"); done; done
mapfile -t MY_TASKS < <(pgvis_array_slice tasks)

for t in "${MY_TASKS[@]}"; do
  read -r cohort chrom <<< "$t"
  [[ "$chrom" == chr* ]] || chrom="chr${chrom}"
  db="$DIR/${cohort}_${chrom}.gdb"
  final="$DIR/${cohort}_${chrom}.filtered.pass.norm.vcf.gz"
  pgvis_done "$final" && continue
  [ -d "$db" ] || { pgvis_log "no GenomicsDB at $db"; continue; }

  raw="$DIR/${cohort}_${chrom}.raw.vcf.gz"
  flt="$DIR/${cohort}_${chrom}.flt.vcf.gz"

  gatk --java-options "-Xmx${JAVA_MEM:-16g}" GenotypeGVCFs \
      -R "$REF_FASTA" -V "gendb://${db}" -L "$chrom" -O "$raw"
  gatk VariantFiltration -R "$REF_FASTA" -V "$raw" \
      --filter-name "pgvis_hard" --filter-expression "$FILTER_EXPR" -O "$flt"

  # PASS only; split multiallelics and left-align, so that one row is one
  # ALT allele. Everything downstream assumes this.
  bcftools view -f PASS "$flt" -Ou \
    | bcftools norm -m -any -f "$REF_FASTA" -Ou \
    | bcftools norm -d exact -Oz -o "${final}.tmp.gz"
  mv "${final}.tmp.gz" "$final"
  bcftools index -t "$final"
  rm -f "$raw" "$raw.tbi" "$flt" "$flt.tbi"
  pgvis_log "wrote $final"
done
pgvis_log "stage 04 complete"
