#!/bin/sh

# Funzione di aiuto
show_help() {
  echo "Usage:"
  echo "$(basename "$0") modkit_output.tsv"
  exit 0
}

# Arguments check
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
  show_help
fi

# Input arguments
MODKIT=$(readlink -f "$1")
SCRIPTDIR=$(dirname "$(readlink -f "$0")")

# Selecting and reshaping fields, collapsing strand information
RESHAPE="${MODKIT%.tsv}_reshaped.tsv"
FINAL_M="${MODKIT%.tsv}_sorted_5mC.tsv"
FINAL_TMP_M="${FINAL_M%.tsv}_tmp.tsv"
FINAL_H="${MODKIT%.tsv}_sorted_5hmC.tsv"
FINAL_TMP_H="${FINAL_H%.tsv}_tmp.tsv"

awk -v OFS="\t" '{if ($6 == "+") {print $4,$3,$3+1,$1,$11,$15,$12} else {print $4,$3-1,$3,$1,$11,$15,$12}}' "$MODKIT" > $RESHAPE


# Sorting and merging

ANNOTYPE=$(head "$RESHAPE" -n 2 | tail -n 1 | head -c 3)

if [ "$ANNOTYPE" = 'chr' ]; then
  for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
  do
    grep -w "^$c" "$RESHAPE" | bedtools sort | awk '{if($7=="m") {print >> "'$FINAL_TMP_M'" } else {print >> "'$FINAL_TMP_H'"}}'
  done
else 
  for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
  do
    grep -w "^$c" "$RESHAPE" | bedtools sort | awk '{if($7=="m") {print >> "'$FINAL_TMP_M'" } else {print >> "'$FINAL_TMP_H'"}}'
  done
fi

rm $RESHAPE

# bgzipping and indexing
if [ -e $FINAL_TMP_M ]; then 
bgzip -c $FINAL_TMP_M > ${FINAL_M}.gz
tabix -f -p bed ${FINAL_M}.gz
zcat "${FINAL_M}.gz" | PERL_HASH_SEED=0 perl "$SCRIPTDIR/ParseModkit.pl" > "${FINAL_M%.tsv}.entropy.file.tsv"

# cleaning
rm "${FINAL_M}.gz"
rm ${FINAL_M}.gz.tbi
rm "${FINAL_TMP_M}"

fi

if [ -e $FINAL_TMP_H ]; then 

bgzip -c $FINAL_TMP_H > ${FINAL_H}.gz
tabix -f -p bed ${FINAL_H}.gz
zcat "${FINAL_H}.gz" | PERL_HASH_SEED=0 perl "$SCRIPTDIR/ParseModkit.pl" > "${FINAL_H%.tsv}.entropy.file.tsv"

# cleaning
rm "${FINAL_H}.gz"
rm ${FINAL_H}.gz.tbi
rm "${FINAL_TMP_H}"

fi

sed -i '1i chrom\tpos\tentropy\tentropy_cov\tbeta\tbeta_cov' "${FINAL_H%.tsv}.entropy.file.tsv"

echo "Entropy file generated. Check the input folder" 
