awk 'BEGIN{OFS="\t"}
NR==1 {print; next}
{
  sum = $2 + $3 + $4
  if (sum == 0) {print $1, 0, 0, 0}
  else {print $1, $2/sum, $3/sum, $4/sum}
}' input_ancestry.tsv > proportions_normalized.tsv

