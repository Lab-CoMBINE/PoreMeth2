#!/bin/sh

# helps
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
FINAL="${MODKIT%.tsv}_sorted.tsv"
FINAL_TMP="${FINAL%.tsv}_tmp.tsv"

char=$(head -n 6 $MODKIT | tail -n 5 | cut -f 12)
MCOUNTER=$(echo "$char" | grep -c "m")
HCOUNTER=$(echo "$char" | grep -c "h")

if [ $MCOUNTER -gt 1 ] && [ $HCOUNTER -gt 1 ]; then
  tail -n +2 $MODKIT  | awk ' NR % 2 == 1 {
      key = $1;
      chr = $4;
      start = $3;
      strand = $6;
      hydro = $11;
  } NR % 2 == 0 {
      methy = $11;
      if (strand == "+") {
          print chr, start, start+1, key, hydro, methy;
      } else {
          print chr, start-1, start, key, hydro, methy;
      }
  } ' OFS="\t" > $RESHAPE

elif [ $MCOUNTER -gt 1 ] && [ $HCOUNTER -eq 0 ]; then
  awk '{if ($12 == "m") {print}}' $MODKIT | awk -v OFS="\t" '{if ($6 == "+") {print $4,$3,$3+1,$1,0,$11} else {print $4,$3-1,$3,$1,0,$11}}' > $RESHAPE

fi


# Sorting and merging

ANNOTYPE=$(head "$RESHAPE" -n 2 | tail -n 1 | head -c 3)

if [ "$ANNOTYPE" = 'chr' ]; then
  for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
  do
    grep -w "^$c" "$RESHAPE" | bedtools sort >> $FINAL_TMP
  done
else 
  for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
  do
    grep -w "^$c" "$RESHAPE" | bedtools sort >> $FINAL_TMP
  done
fi

rm $RESHAPE

# bgzipping and indexing
if [ -e $FINAL_TMP ]; then 
bgzip -c $FINAL_TMP > ${FINAL}.gz
tabix -f -p bed ${FINAL}.gz
zcat "${FINAL}.gz" | PERL_HASH_SEED=0 perl "$SCRIPTDIR/ParseModkit.pl" > "${FINAL%.tsv}.entropy.file.tsv"

# cleaning
rm "${FINAL}.gz"
rm ${FINAL}.gz.tbi
rm "${FINAL_TMP}"

fi

echo "Entropy file generated" 
