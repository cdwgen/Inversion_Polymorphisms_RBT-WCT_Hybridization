#!/bin/bash

# Input file
MAPTXT="Gerrard_MapInfo.txt"
OUTTSV="Gerrard_linkage_map.tsv"

# Extract header index for needed columns
awk '
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == "CatID_X") name_idx = i;
    else if ($i == "LG") lg_idx = i;
    else if ($i == "ORDER") order_idx = i;
  }
  print "name", "LG", "cM"
}
NR > 1 {
  print $name_idx, $lg_idx, $order_idx
}
' "$MAPTXT" > "$OUTTSV"
