#!/usr/bin/env bash

# Usage: ./reorder_table.sh sample_order.txt input_table.tsv > ELAI_Results_Ordered.tsv

order_file="$1"
input_table="$2"

awk -v OFS="\t" '
    FNR==NR { order[NR]=$1; next }       # Read order file, store IDs in array
    FNR==1 { header=$0; next }          # Save header line
    { table[$1]=$0 }                     # Store each table line keyed by Sample_ID
    END {
        print header                     # Print header first
        for (i=1; i<=length(order); i++) {
            id=order[i]
            if (id in table) print table[id]
            else print id "\tMISSING"   # Warn if missing in table
        }
    }
' "$order_file" "$input_table"