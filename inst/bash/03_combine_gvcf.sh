#!/usr/bin/env bash
# Stage 03 -- combine per-sample GVCFs into a cohort GenomicsDB, per chromosome.
# Required: PGVIS_ROOT REF_FASTA COHORTS ; Optional: CHROMS BATCH_SIZE
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_need gatk tabix
pgvis_require_var PGVIS_ROOT REF_FASTA COHORTS
GVCF_DIR="${GVCF_DIR:-$PGVIS_ROOT/02_gvcf}"
OUT_DIR="${OUT_DIR:-$PGVIS_ROOT/03_joint}"
BATCH_SIZE="${BATCH_SIZE:-50}"
read -r -a COHORT_ARR <<< "$COHORTS"
read -r -a CHROM_ARR  <<< "${CHROMS:-$(seq 1 22) X Y}"
mkdir -p "$OUT_DIR"

# Build the (cohort, chromosome) task grid once; an array job takes one cell.
tasks=(); for c in "${COHORT_ARR[@]}"; do for k in "${CHROM_ARR[@]}"; do tasks+=("$c $k"); done; done
mapfile -t MY_TASKS < <(pgvis_array_slice tasks)

for t in "${MY_TASKS[@]}"; do
  read -r cohort chrom <<< "$t"
  [[ "$chrom" == chr* ]] || chrom="chr${chrom}"
  db="$OUT_DIR/${cohort}_${chrom}.gdb"
  [ -d "$db" ] && [ "${PGVIS_FORCE:-0}" != "1" ] && { pgvis_log "skip existing $db"; continue; }
  rm -rf "$db"

  # Index any GVCF whose .tbi went missing, instead of maintaining a separate
  # "rerun the ones with broken tbi" script.
  map="$OUT_DIR/${cohort}_${chrom}.sample_map"; : > "$map"
  while IFS= read -r g; do
    [ -f "${g}.tbi" ] || { pgvis_log "reindexing $(basename "$g")"; tabix -p vcf "$g"; }
    printf '%s\t%s\n' "$(basename "$(dirname "$g")")" "$g" >> "$map"
  done < <(find "$GVCF_DIR/$cohort" -type f -name "${chrom}.g.vcf.gz" | sort)

  n=$(wc -l < "$map")
  [ "$n" -eq 0 ] && { pgvis_log "no GVCFs for $cohort $chrom"; continue; }
  pgvis_log "GenomicsDBImport: $cohort $chrom ($n samples)"

  gatk --java-options "-Xmx${JAVA_MEM:-16g}" GenomicsDBImport \
    --genomicsdb-workspace-path "$db" \
    --sample-name-map "$map" -L "$chrom" \
    --batch-size "$BATCH_SIZE" --reader-threads "${THREADS:-4}" \
    || { pgvis_log "FAILED $cohort $chrom"; rm -rf "$db"; continue; }
done
pgvis_log "stage 03 complete"
