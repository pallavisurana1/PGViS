#!/usr/bin/env bash
# Stage 06 -- intersect cohort variants with regulatory regions, then dbSNP.
#
# Required: PGVIS_ROOT REGION_MANIFEST ; Optional: DBSNP_DIR VAR_BED_DIR
# REGION_MANIFEST is a headerless TSV: feature <TAB> feature_type <TAB> bed_path
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_need bedtools
pgvis_require_var PGVIS_ROOT REGION_MANIFEST
VAR_BED_DIR="${VAR_BED_DIR:-$PGVIS_ROOT/04_bed}"
OUT_DIR="${OUT_DIR:-$PGVIS_ROOT/06_dbsnp}"
MID_DIR="${MID_DIR:-$PGVIS_ROOT/05_intersect}"
mkdir -p "$OUT_DIR" "$MID_DIR"

mapfile -t features < <(grep -v '^#' "$REGION_MANIFEST" | grep -v '^[[:space:]]*$')
mapfile -t MY_FEATURES < <(pgvis_array_slice features)

for line in "${MY_FEATURES[@]}"; do
  IFS=$'\t' read -r feature ftype region_bed <<< "$line"
  [ -s "$region_bed" ] || { pgvis_log "missing region BED for $feature -- skipping"; continue; }

  for var_bed in "$VAR_BED_DIR"/*.bed; do
    [ -e "$var_bed" ] || continue
    tag=$(basename "$var_bed" .bed)

    # Chromosome is taken from the variant BED filename (e.g. luad_chr7.bed)
    # so the matching per-chromosome dbSNP file can be found.
    if [[ "$tag" =~ chr([0-9]{1,2}|X|Y|M) ]]; then chr="${BASH_REMATCH[1]}"; else chr=""; fi

    mid="$MID_DIR/$feature/${tag}.bed"
    out="$OUT_DIR/$feature/${tag}_dbsnp.bed"
    mkdir -p "$(dirname "$mid")" "$(dirname "$out")"
    pgvis_done "$out" && continue

    # -a region -b variants -wa -wb keeps the region columns first, which is
    # the layout pgvis_schema(region_first = TRUE) expects.
    if ! pgvis_done "$mid"; then
      bedtools intersect -a "$region_bed" -b "$var_bed" -wa -wb > "${mid}.tmp" \
        && mv "${mid}.tmp" "$mid"
    fi
    [ -s "$mid" ] || { pgvis_log "no overlap: $feature x $tag"; : > "$out"; continue; }

    if [ -n "${DBSNP_DIR:-}" ] && [ -n "$chr" ]; then
      dbsnp="$DBSNP_DIR/dbSnp155.chr${chr}.bed"
      [ -s "$dbsnp" ] || { pgvis_log "no dbSNP for chr${chr}; copying un-annotated"; cp "$mid" "$out"; continue; }
      bedtools intersect -a "$mid" -b "$dbsnp" -wa -wb > "${out}.tmp" && mv "${out}.tmp" "$out"
    else
      cp "$mid" "$out"
    fi
    pgvis_log "$feature x $tag -> $(wc -l < "$out") rows"
  done
done
pgvis_log "stage 06 complete"
