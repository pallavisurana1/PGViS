#!/usr/bin/awk -f
# VCF -> BED, one row per ALT allele. Headerless, tab-delimited, BED
# coordinates (0-based, half-open). Column order matches
# PGViS::pgvis_schema()'s "cohort" block -- if you change one, change both.
#
#   1 chrom  2 chromStart  3 chromEnd  4 ref  5 alt  6 allele_idx  7 vtype
#   8 ac     9 af          10 an       11 n_carrier  12 carrier_samples
#
#   bcftools view f.vcf.gz | awk -f vcf2bed_carriers.awk > f.bed
BEGIN { FS = "\t"; OFS = "\t" }
/^#CHROM/ { for (i = 10; i <= NF; i++) sample[i] = $i; nsamp = NF; next }
/^#/ { next }
{
    chrom = $1; pos = $2; ref = $4; info = $8
    start = pos - 1
    end   = start + length(ref)
    nalt  = split($5, alts, ",")

    ac_all = af_all = an = "."
    n = split(info, kv, ";")
    for (k = 1; k <= n; k++) {
        split(kv[k], p, "=")
        if      (p[1] == "AC") ac_all = p[2]
        else if (p[1] == "AF") af_all = p[2]
        else if (p[1] == "AN") an     = p[2]
    }
    split(ac_all, acv, ","); split(af_all, afv, ",")

    for (j = 1; j <= nalt; j++) {
        a1 = alts[j]; rl = length(ref); al = length(a1)
        if      (rl == al && rl == 1) vtype = "snv"
        else if (rl <  al)            vtype = "ins"
        else if (rl >  al)            vtype = "del"
        else                          vtype = "mnv"

        cnt = 0; names = ""
        for (i = 10; i <= nsamp; i++) {
            gt = $i; sub(/:.*/, "", gt)
            a = gt; gsub(/[|\/]/, " ", a)
            m = split(a, gl, " "); carry = 0
            for (x = 1; x <= m; x++) if (gl[x] + 0 == j) carry = 1
            if (carry) { cnt++; names = (names == "" ? sample[i] : names "," sample[i]) }
        }
        if (names == "") names = "."

        print chrom, start, end, ref, a1, j, vtype, \
              (j in acv ? acv[j] : "."), (j in afv ? afv[j] : "."), an, cnt, names
    }
}
