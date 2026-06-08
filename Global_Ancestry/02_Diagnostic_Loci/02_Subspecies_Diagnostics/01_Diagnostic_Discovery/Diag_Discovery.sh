#!/bin/bash
set -euo pipefail

# ==== CONFIG ====
THRESH=0.95          # allele freq threshold for "diagnostic"
FST_THRESH=0.90      # FST threshold for high differentiation

# Directories
FREQ_DIR="freq_averages"
FST_DIR="fst_averages"

# Averaged freq files
FREQ_IRT="$FREQ_DIR/IRT_averaged.frq"
FREQ_CRT="$FREQ_DIR/CRT_averaged.frq"
FREQ_WCT="$FREQ_DIR/WCT_averaged.frq"
FREQ_YCT="$FREQ_DIR/YCT_averaged.frq"

# Averaged FST files
FST_IRT_CRT="$FST_DIR/fst_IRTvsCRT.weir.fst"
FST_IRT_WCT="$FST_DIR/fst_IRTvsWCT.weir.fst"
FST_IRT_YCT="$FST_DIR/fst_IRTvsYCT.weir.fst"
FST_CRT_WCT="$FST_DIR/fst_CRTvsWCT.weir.fst"
FST_CRT_YCT="$FST_DIR/fst_CRTvsYCT.weir.fst"
FST_WCT_YCT="$FST_DIR/fst_WCTvsYCT.weir.fst"

# Check files exist
for f in "$FREQ_IRT" "$FREQ_CRT" "$FREQ_WCT" "$FREQ_YCT" \
         "$FST_IRT_CRT" "$FST_IRT_WCT" "$FST_IRT_YCT" \
         "$FST_CRT_WCT" "$FST_CRT_YCT" "$FST_WCT_YCT"; do
    [ -s "$f" ] || { echo "Missing or empty file: $f"; exit 1; }
done

echo "Using thresholds: THRESH=$THRESH, FST_THRESH=$FST_THRESH"

# ==== Clean FST files to remove -nan lines ====
awk '$3 != "-nan"' "$FST_IRT_CRT" > fst_IRT_CRT_cleaned.fst
awk '$3 != "-nan"' "$FST_IRT_WCT" > fst_IRT_WCT_cleaned.fst
awk '$3 != "-nan"' "$FST_IRT_YCT" > fst_IRT_YCT_cleaned.fst
awk '$3 != "-nan"' "$FST_CRT_WCT" > fst_CRT_WCT_cleaned.fst
awk '$3 != "-nan"' "$FST_CRT_YCT" > fst_CRT_YCT_cleaned.fst
awk '$3 != "-nan"' "$FST_WCT_YCT" > fst_WCT_YCT_cleaned.fst

# ==== Parse high-FST SNPs ====
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_IRT_CRT_cleaned.fst | sort > tmp.IRTvsCRT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_IRT_WCT_cleaned.fst | sort > tmp.IRTvsWCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_IRT_YCT_cleaned.fst | sort > tmp.IRTvsYCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_CRT_WCT_cleaned.fst | sort > tmp.CRTvsWCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_CRT_YCT_cleaned.fst | sort > tmp.CRTvsYCT
awk -v t=$FST_THRESH '$3 >= t {print $1":"$2}' fst_WCT_YCT_cleaned.fst | sort > tmp.WCTvsYCT

# ==== Intersect high-FST SNP sets for candidate diagnostic SNPs ====
# IRT: must be high FST in IRT vs CRT, IRT vs WCT, IRT vs YCT
comm -12 tmp.IRTvsCRT <(comm -12 tmp.IRTvsWCT tmp.IRTvsYCT) > IRT_diagnostic_candidates.txt

# CRT: high in IRTvsCRT, CRTvsWCT, CRTvsYCT
comm -12 tmp.IRTvsCRT <(comm -12 tmp.CRTvsWCT tmp.CRTvsYCT) > CRT_diagnostic_candidates.txt

# WCT: high in IRTvsWCT, CRTvsWCT, WCTvsYCT
comm -12 tmp.IRTvsWCT <(comm -12 tmp.CRTvsWCT tmp.WCTvsYCT) > WCT_diagnostic_candidates.txt

# YCT: high in IRTvsYCT, CRTvsYCT, WCTvsYCT
comm -12 tmp.IRTvsYCT <(comm -12 tmp.CRTvsYCT tmp.WCTvsYCT) > YCT_diagnostic_candidates.txt

# ==== Function to classify allele identity ====
classify_snps() {
    GROUP=$1
    CANDIDATES=$2
    OUTFILE=$3
    FREQ_IRT=$4
    FREQ_CRT=$5
    FREQ_WCT=$6
    FREQ_YCT=$7

    echo "Classifying SNPs for $GROUP ..."

    # Extract REF allele frequencies (SNP_ID = CHROM:POS)
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_IRT" | sort -k1,1 > IRT.frq.tab
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_CRT" | sort -k1,1 > CRT.frq.tab
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_WCT" | sort -k1,1 > WCT.frq.tab
    awk 'NR>1 {print $1":"$2, $5}' "$FREQ_YCT" | sort -k1,1 > YCT.frq.tab

    sort -k1,1 "$CANDIDATES" > candidates.sorted

    join -1 1 -2 1 candidates.sorted IRT.frq.tab > tmp1
    join -1 1 -2 1 tmp1 CRT.frq.tab > tmp2
    join -1 1 -2 1 tmp2 WCT.frq.tab > tmp3
    join -1 1 -2 1 tmp3 YCT.frq.tab > tmp4

    awk -v group=$GROUP -v t=$THRESH '
    {
        id = $1
        f_irt = ($2=="") ? 0 : $2+0
        f_crt = ($3=="") ? 0 : $3+0
        f_wct = ($4=="") ? 0 : $4+0
        f_yct = ($5=="") ? 0 : $5+0

        if (group=="IRT") {
            if (f_irt >= t && f_crt <= 1-t && f_wct <= 1-t && f_yct <= 1-t) print id, "REF"
            else if (f_irt <= 1-t && f_crt >= t && f_wct >= t && f_yct >= t) print id, "ALT"
        }
        else if (group=="CRT") {
            if (f_crt >= t && f_irt <= 1-t && f_wct <= 1-t && f_yct <= 1-t) print id, "REF"
            else if (f_crt <= 1-t && f_irt >= t && f_wct >= t && f_yct >= t) print id, "ALT"
        }
        else if (group=="WCT") {
            if (f_wct >= t && f_irt <= 1-t && f_crt <= 1-t && f_yct <= 1-t) print id, "REF"
            else if (f_wct <= 1-t && f_irt >= t && f_crt >= t && f_yct >= t) print id, "ALT"
        }
        else if (group=="YCT") {
            if (f_yct >= t && f_irt <= 1-t && f_crt <= 1-t && f_wct <= 1-t) print id, "REF"
            else if (f_yct <= 1-t && f_irt >= t && f_crt >= t && f_wct >= t) print id, "ALT"
        }
    }' tmp4 | sort -u > "$OUTFILE"

    echo "  Wrote $OUTFILE"
    rm -f IRT.frq.tab CRT.frq.tab WCT.frq.tab YCT.frq.tab candidates.sorted tmp1 tmp2 tmp3 tmp4
}

# ==== Run classification for each group ====
classify_snps IRT IRT_diagnostic_candidates.txt IRT_diagnostic_snps.txt "$FREQ_IRT" "$FREQ_CRT" "$FREQ_WCT" "$FREQ_YCT"
classify_snps CRT CRT_diagnostic_candidates.txt CRT_diagnostic_snps.txt "$FREQ_IRT" "$FREQ_CRT" "$FREQ_WCT" "$FREQ_YCT"
classify_snps WCT WCT_diagnostic_candidates.txt WCT_diagnostic_snps.txt "$FREQ_IRT" "$FREQ_CRT" "$FREQ_WCT" "$FREQ_YCT"
classify_snps YCT YCT_diagnostic_candidates.txt YCT_diagnostic_snps.txt "$FREQ_IRT" "$FREQ_CRT" "$FREQ_WCT" "$FREQ_YCT"

# ==== Final cleanup ====
rm -f fst_*_cleaned.fst tmp.*
echo "All done!"
