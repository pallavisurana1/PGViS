#!/usr/bin/env bash
# Stage 05 -- cohort VCF to BED with allele counts and per-variant carrier lists.
# Required: PGVIS_ROOT ; Optional: VCF_DIR COHORTS CHROMS VCF_GLOB
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_need bcftools awk
pgvis_require_var PGVIS_ROOT
VCF_DIR="${VCF_DIR:-$PGVIS_ROOT/03_joint}"
OUT_DIR="${OUT_DIR:-$PGVIS_ROOT/04_bed}"
AWK_SCRIPT="${AWK_SCRIPT:-${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/vcf2bed_carriers.awk}"
VCF_GLOB="${VCF_GLOB:-*.filtered.pass.norm.vcf.gz}"
mkdir -p "$OUT_DIR"

mapfile -t vcfs < <(find "$VCF_DIR" -maxdepth 1 -name "$VCF_GLOB" | sort -V)
[ "${#vcfs[@]}" -eq 0 ] && pgvis_die "no VCFs matching $VCF_GLOB under $VCF_DIR"
mapfile -t MY_VCFS < <(pgvis_array_slice vcfs)

for vcf in "${MY_VCFS[@]}"; do
  base=$(basename "$vcf"); base="${base%%.*}"
  out="$OUT_DIR/${base}.bed"
  pgvis_done "$out" && continue
  pgvis_log "$base"
  bcftools view "$vcf" | awk -f "$AWK_SCRIPT" > "${out}.tmp" && mv "${out}.tmp" "$out"
  pgvis_log "$out ($(wc -l < "$out") rows)"
done

# Sidecar documenting the layout; the BED itself stays headerless per spec.
cat > "$OUT_DIR/COLUMNS.txt" <<'LEG'
1  chrom            2  chromStart (0-based)   3  chromEnd (start + len(REF))
4  ref              5  alt (one row per ALT)  6  allele_idx
7  vtype            8  ac                     9  af
10 an               11 n_carrier              12 carrier_samples ('.' if none)
LEG
pgvis_log "stage 05 complete"
