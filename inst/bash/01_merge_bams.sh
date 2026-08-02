#!/usr/bin/env bash
# Stage 01 -- merge replicate BAMs per patient and tissue.
# Required: PGVIS_ROOT SAMPLE_SHEET ; Optional: RAW_BAM_DIR THREADS
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_need samtools
pgvis_require_var PGVIS_ROOT SAMPLE_SHEET
RAW_BAM_DIR="${RAW_BAM_DIR:-$PGVIS_ROOT/00_raw_bam}"
OUT_DIR="${OUT_DIR:-$PGVIS_ROOT/01_bam}"
THREADS="${THREADS:-4}"

while IFS=$'\t' read -r cohort tissue sample; do
  [ -z "${sample:-}" ] && continue; [[ "$cohort" == \#* ]] && continue
  src="$RAW_BAM_DIR/$cohort/$tissue/$sample"
  dst="$OUT_DIR/$cohort/$tissue/$sample"; mkdir -p "$dst"
  out="$dst/${sample}.bam"
  pgvis_done "$out" && continue

  mapfile -t bams < <(find "$src" -type f -name '*.bam' | sort)
  case "${#bams[@]}" in
    0) pgvis_log "no BAMs under $src -- skipping"; continue ;;
    1) pgvis_log "single BAM for $sample -- linking rather than merging"
       ln -sf "${bams[0]}" "$out" ;;
    *) pgvis_log "merging ${#bams[@]} replicates for $sample"
       samtools merge -@ "$THREADS" -f "${out}.tmp" "${bams[@]}" && mv "${out}.tmp" "$out" ;;
  esac
  samtools index -@ "$THREADS" "$out"
done < "$SAMPLE_SHEET"
pgvis_log "stage 01 complete"
