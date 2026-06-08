#!/bin/bash
set -euo pipefail

# ---------------- CONFIG ----------------
THRESH=1.0                     # allele freq threshold for "fixed"
CHR="Omy17"                     # chromosome
START=36050000                   # inversion start (bp)
END=43906000                    # inversion end (bp)

FREQ_DIR="freq_averages"
FREQ_RBT="$FREQ_DIR/RBT_averaged.frq"
FREQ_WCT="$FREQ_DIR/WCT_averaged.frq"
FREQ_YCT="$FREQ_DIR/YCT_averaged.frq"

OUTDIR="diagnostic_inversion_snps"
mkdir -p "$OUTDIR"

# Output files (tab-separated)
OUT_RBT="$OUTDIR/RBT_diagnostic_snps.tsv"
OUT_WCT="$OUTDIR/WCT_diagnostic_snps.tsv"
OUT_YCT="$OUTDIR/YCT_diagnostic_snps.tsv"

# -----------------------------------------
# Basic checks
for f in "$FREQ_RBT" "$FREQ_WCT" "$FREQ_YCT"; do
    [ -s "$f" ] || { echo "Missing or empty file: $f"; exit 1; }
done

echo "Searching $CHR:$START-$END for sites fixed (>= $THRESH) in one taxon and opposite in the other two"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# ---- produce CHR:POS -> REF_FREQ tables restricted to region ----
# assumes vcftools .frq format where $1=CHROM $2=POS $5=REF allele frequency (0..1)
awk -v chr="$CHR" -v s=$START -v e=$END 'NR>1 && $1==chr && $2>=s && $2<=e {print $1":"$2 "\t" $5}' "$FREQ_RBT" | sort -k1,1 > "$TMPDIR/RBT.tab"
awk -v chr="$CHR" -v s=$START -v e=$END 'NR>1 && $1==chr && $2>=s && $2<=e {print $1":"$2 "\t" $5}' "$FREQ_WCT" | sort -k1,1 > "$TMPDIR/WCT.tab"
awk -v chr="$CHR" -v s=$START -v e=$END 'NR>1 && $1==chr && $2>=s && $2<=e {print $1":"$2 "\t" $5}' "$FREQ_YCT" | sort -k1,1 > "$TMPDIR/YCT.tab"

# ---- require presence in all three (join) ----
# join with tab delimiter to preserve fields
join -t $'\t' -1 1 -2 1 "$TMPDIR/RBT.tab" "$TMPDIR/WCT.tab" > "$TMPDIR/tmp1" || true
join -t $'\t' -1 1 -2 1 "$TMPDIR/tmp1" "$TMPDIR/YCT.tab" > "$TMPDIR/all_freqs.tab" || true

# If no sites shared across all three, exit cleanly
if [ ! -s "$TMPDIR/all_freqs.tab" ]; then
    echo "No sites present in all three taxa within region. Exiting."
    exit 0
fi

# ---- header for outputs ----
printf "CHR\tPOS\tID\tCLASS\tf_RBT\tf_WCT\tf_YCT\n" > "$OUT_RBT"
printf "CHR\tPOS\tID\tCLASS\tf_RBT\tf_WCT\tf_YCT\n" > "$OUT_WCT"
printf "CHR\tPOS\tID\tCLASS\tf_RBT\tf_WCT\tf_YCT\n" > "$OUT_YCT"

# ---- classify diagnostic SNPs (REF or ALT for focal taxon) ----
awk -v t="$THRESH" -v rfile="$OUT_RBT" -v wfile="$OUT_WCT" -v yfile="$OUT_YCT" '
BEGIN { OFS="\t"; u = 1 - t }
{
    # fields after joins: $1=id, $2=f_rbt, $3=f_wct, $4=f_yct
    id = $1
    split(id, a, ":")
    chr = a[1]; pos = a[2]
    f_r = ($2+0); f_w = ($3+0); f_y = ($4+0)

    # RBT diagnostic?
    if (f_r >= t && f_w <= u && f_y <= u) {
        print chr, pos, id, "REF", f_r, f_w, f_y >> rfile
    } else if (f_r <= u && f_w >= t && f_y >= t) {
        print chr, pos, id, "ALT", f_r, f_w, f_y >> rfile
    }

    # WCT diagnostic?
    if (f_w >= t && f_r <= u && f_y <= u) {
        print chr, pos, id, "REF", f_r, f_w, f_y >> wfile
    } else if (f_w <= u && f_r >= t && f_y >= t) {
        print chr, pos, id, "ALT", f_r, f_w, f_y >> wfile
    }

    # YCT diagnostic?
    if (f_y >= t && f_r <= u && f_w <= u) {
        print chr, pos, id, "REF", f_r, f_w, f_y >> yfile
    } else if (f_y <= u && f_r >= t && f_w >= t) {
        print chr, pos, id, "ALT", f_r, f_w, f_y >> yfile
    }
}
' "$TMPDIR/all_freqs.tab"

echo "Done. Outputs:"
ls -lh "$OUT_RBT" "$OUT_WCT" "$OUT_YCT" || true
