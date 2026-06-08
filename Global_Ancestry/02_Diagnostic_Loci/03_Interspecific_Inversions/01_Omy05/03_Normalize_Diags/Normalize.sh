awk 'BEGIN{OFS="\t"}
NR==1 {print; next}
{
  sum = $2 + $3 + $4
  if (sum == 0) {print $1, 0, 0, 0}
  else {print $1, $2/sum, $3/sum, $4/sum}
}' Omy05_all.tsv > Omy05_all_normalized.tsv

awk 'BEGIN{OFS="\t"}
NR==1 {print; next}
{
  sum = $2 + $3 + $4
  if (sum == 0) {print $1, 0, 0, 0}
  else {print $1, $2/sum, $3/sum, $4/sum}
}' Omy05_Inv1.tsv > Omy05_Inv1_normalized.tsv


awk 'BEGIN{OFS="\t"}
NR==1 {print; next}
{
  sum = $2 + $3 + $4
  if (sum == 0) {print $1, 0, 0, 0}
  else {print $1, $2/sum, $3/sum, $4/sum}
}' Omy05_Inv2.tsv > Omy05_Inv2_normalized.tsv

