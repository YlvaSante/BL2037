#!/bin/bash

ORIGINAL_FILE="text_file.txt"
TARGET_WORD="Apollo"
OUTPUT_DIR="resultat_filer"
RESULT_CSV="counts.csv"

mkdir -p $OUTPUT_DIR
rm -f $OUTPUT_DIR/*.txt
awk -v dir="$OUTPUT_DIR" '/^#/{f=dir "/" $2 ".txt"} {print > f}' "$ORIGINAL_FILE"

echo "Date,Word,Count" > $RESULT_CSV

for fil in $OUTPUT_DIR/*.txt; do
    # Extrahera bara datumet från filnamnet för snyggare output
    date=$(basename "$fil" .txt)

    # Räkna förekomsten
    number=$(grep -ci "$TARGET_WORD" "$fil")

    # Skriv till skärmen (för att se att det händer nåt)
    echo "Bearbetar $date: hittade $number st."

    # Spara till CSV
    echo "$date,$TARGET_WORD,$number" >> $RESULT_CSV
done

