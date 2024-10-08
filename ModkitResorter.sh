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
MODKITFILENAME=$(basename $"1")
SCRIPTDIR=$(dirname "$(readlink -f $"0")")

# Selecting and reshaping fields, collapsing strand information
RESMAPE="${MODKIT%.tsv}_reshaped.tsv"
FINAL_M="${MODKIT%.tsv}_sorted_5mC.tsv"
FINAL_TMP_M="${FINAL_M%.tsv}_tmp.tsv"

awk '{if ($12 == "m") {print}}' "$MODKIT" | awk -v OFS="\t" '{if ($6 == "+") {print $4,$3,$3+1,$1,$11,$15,$12} else {print $4,$3-1,$3,$1,$11,$15,$12}}' > "$RESMAPE"

# Sorting and merging
for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
    C_M="${FINAL_TMP_M%_tmp.tsv}_c${c}_tmp.tsv"
    grep -w "^$c" "$RESMAPE" | bedtools sort > "$C_M"
    cat "$C_M" >> "$FINAL_TMP_M"
    rm "$C_M"
done

rm $RESMAPE

# bgzipping and indexing
bgzip -c "$FINAL_TMP_M" > "${FINAL_M}.gz"
tabix -f -p bed "${FINAL_M}.gz"
rm "${FINAL_M%.tsv}_tmp.tsv"

zcat "${FINAL_M}.gz" | perl "$SCRIPTDIR/ParseModkit.pl" > "${FINAL_M%.tsv}.entropy.file.tsv"
rm "${FINAL_M}.gz"