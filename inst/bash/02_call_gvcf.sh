#!/usr/bin/env bash
# Stage 02 -- per-sample, per-chromosome GVCF calling with GATK HaplotypeCaller.
#
# Required: PGVIS_ROOT REF_FASTA SAMPLE_SHEET
# Optional: CHROMS BAM_DIR JAVA_MEM THREADS PGVIS_FORCE
#
# SAMPLE_SHEET is a headerless TSV: cohort <TAB> tissue <TAB> sample_id
# The BAM is located at BAM_DIR/cohort/tissue/sample_id/*.bam

source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh" 2>/dev/null \
  || source "${PGVIS_BASH_LIB:?set PGVIS_BASH_LIB to the inst/bash directory}/lib/common.sh"

pgvis_need gatk
pgvis_require_var PGVIS_ROOT REF_FASTA SAMPLE_SHEET

BAM_DIR="${BAM_DIR:-$PGVIS_ROOT/01_bam}"
OUT_DIR="${OUT_DIR:-$PGVIS_ROOT/02_gvcf}"
JAVA_MEM="${JAVA_MEM:-8g}"
read -r -a CHROM_ARR <<< "${CHROMS:-$(seq 1 22) X Y}"

mapfile -t CHROM_LIST < <(pgvis_array_slice CHROM_ARR)

mkdir -p "$OUT_DIR"

for CHROM in "${CHROM_LIST[@]}"; do
  [[ "$CHROM" == chr* ]] || CHROM="chr${CHROM}"
  pgvis_log "chromosome $CHROM"

  while IFS=$'\t' read -r cohort tissue sample; do
    [ -z "${sample:-}" ] && continue
    [[ "$cohort" == \#* ]] && continue

    bam=$(find "$BAM_DIR/$cohort/$tissue/$sample" -type f -name '*.bam' 2>/dev/null | head -1)
    if [ -z "$bam" ]; then
      pgvis_log "no BAM for $cohort/$tissue/$sample -- skipping"
      continue
    fi

    sample_out="$OUT_DIR/$cohort/$tissue/$sample"
    mkdir -p "$sample_out"
    out_vcf="$sample_out/${CHROM}.g.vcf.gz"
    pgvis_done "$out_vcf" && continue

    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
      -R "$REF_FASTA" -I "$bam" -L "$CHROM" -O "${out_vcf}.tmp.gz" -ERC GVCF \
      && mv "${out_vcf}.tmp.gz" "$out_vcf" \
      && mv "${out_vcf}.tmp.gz.tbi" "${out_vcf}.tbi" 2>/dev/null \
      || { pgvis_log "FAILED $cohort/$tissue/$sample $CHROM"; rm -f "${out_vcf}.tmp.gz"*; continue; }

    pgvis_log "done $cohort/$tissue/$sample $CHROM"
  done < "$SAMPLE_SHEET"
done

pgvis_log "stage 02 complete"
