#!/bin/bash
set -euo pipefail

# ==== CONFIG ====
THRESH=0.99          # allele freq threshold for "diagnostic"
FST_THRESH=0.95      # FST threshold for high differentiation

# Directories
FREQ_DIR="freq_averages"
FST_DIR="fst_averages"

# Averaged freq files
FREQ_RBT="$FREQ_DIR/RBT_averaged.frq"
FREQ_WCT="$FREQ_DIR/WCT_averaged.frq"
FREQ_YCT="$FREQ_DIR/YCT_averaged.frq"

# Averaged FST files
FST_RBT_WCT="$FST_DIR/fst_RBTvsWCT.weir.fst"
FST_RBT_YCT="$FST_DIR/fst_RBTvsYCT.weir.fst"
FST_WCT_YCT="$FST_DIR/fst_WCTvsYCT.weir.fst"

# Check files exist
for f in "$FREQ_RBT" "$FREQ_WCT" "$FREQ_YCT" "$FST_RBT_WCT" "$FST_RBT_YCT" "$FST_WCT_YCT"; do
    [ -s "$f" ] || { echo "Missing or empty file: $f"; exit 1; }
done

echo "Using thresholds: THRESH=$THRESH, FST_THRESH=$FST_THRESH"

# ==== Clean FST files to remove -nan lines ====
awk '$3 != "-nan"' "$FST_RBT_WCT" > fst_RBT_WCT_cleaned.fst
awk '$3 != "-nan"' "$FST_RBT_YCT" > fst_RBT_YCT_cleaned.fst
awk '$3 != "-nan"' "$FST_WCT_YCT" > fst_WCT_YCT_cleaned.fst

# ==== Parse high-FST SNPs ====
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_RBT_WCT_cleaned.fst | sort > tmp.RBTvsWCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_RBT_YCT_cleaned.fst | sort > tmp.RBTvsYCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_WCT_YCT_cleaned.fst | sort > tmp.WCTvsYCT

# ==== Intersect high-FST SNP sets for candidate diagnostic SNPs ====
comm -12 tmp.RBTvsWCT tmp.RBTvsYCT > RBT_diagnostic_candidates.txt
comm -12 tmp.RBTvsWCT tmp.WCTvsYCT > WCT_diagnostic_candidates.txt
comm -12 tmp.RBTvsYCT tmp.WCTvsYCT > YCT_diagnostic_candidates.txt

# ==== Function to classify allele identity ====
classify_snps() {
    GROUP=$1
    CANDIDATES=$2
    OUTFILE=$3
    FREQ_RBT=$4
    FREQ_WCT=$5
    FREQ_YCT=$6

    echo "Classifying SNPs for $GROUP ..."

    # Extract REF allele frequencies (SNP_ID = CHROM:POS)
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_RBT" | sort -k1,1 > RBT.frq.tab
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_WCT" | sort -k1,1 > WCT.frq.tab
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_YCT" | sort -k1,1 > YCT.frq.tab

    # Sort candidate SNPs
    sort -k1,1 "$CANDIDATES" > candidates.sorted

    # Join frequencies across groups
    join -1 1 -2 1 candidates.sorted RBT.frq.tab > tmp1
    join -1 1 -2 1 tmp1 WCT.frq.tab > tmp2
    join -1 1 -2 1 tmp2 YCT.frq.tab > tmp3

    # Classify SNPs as REF or ALT for this group
    awk -v group=$GROUP -v t=$THRESH '
    {
        id = $1
        f_rbt = ($2=="") ? 0 : $2+0
        f_wct = ($3=="") ? 0 : $3+0
        f_yct = ($4=="") ? 0 : $4+0

        if (group=="G1") {
            if (f_rbt >= t && f_wct <= 1-t && f_yct <= 1-t) print id, "REF"
            else if (f_rbt <= 1-t && f_wct >= t && f_yct >= t) print id, "ALT"
        }
        else if (group=="G2") {
            if (f_wct >= t && f_rbt <= 1-t && f_yct <= 1-t) print id, "REF"
            else if (f_wct <= 1-t && f_rbt >= t && f_yct >= t) print id, "ALT"
        }
        else if (group=="G3") {
            if (f_yct >= t && f_rbt <= 1-t && f_wct <= 1-t) print id, "REF"
            else if (f_yct <= 1-t && f_rbt >= t && f_wct >= t) print id, "ALT"
        }
    }' tmp3 | sort -u > "$OUTFILE"

    echo "  Wrote $OUTFILE"
    # cleanup
    rm -f RBT.frq.tab WCT.frq.tab YCT.frq.tab candidates.sorted tmp1 tmp2 tmp3
}

# ==== Run classification for each group ====
classify_snps G1 RBT_diagnostic_candidates.txt RBT_diagnostic_snps.txt "$FREQ_RBT" "$FREQ_WCT" "$FREQ_YCT"
classify_snps G2 WCT_diagnostic_candidates.txt WCT_diagnostic_snps.txt "$FREQ_RBT" "$FREQ_WCT" "$FREQ_YCT"
classify_snps G3 YCT_diagnostic_candidates.txt YCT_diagnostic_snps.txt "$FREQ_RBT" "$FREQ_WCT" "$FREQ_YCT"

# ==== Final cleanup ====
rm -f fst_*_cleaned.fst tmp.* RBT_diagnostic_candidates.txt WCT_diagnostic_candidates.txt YCT_diagnostic_candidates.txt
echo "All done!"
