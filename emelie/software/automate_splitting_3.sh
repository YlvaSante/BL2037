#!/bin/bash

ORIGINAL_FILE="text_file.txt"
# Vi skapar en lista med orden separerade med mellanslag för loopen
WORDS=("Apollo" "moon" "landing")
OUTPUT_DIR="resultat_filer"
RESULT_CSV="multi_counts.csv"

mkdir -p $OUTPUT_DIR
rm -f $OUTPUT_DIR/*.txt
awk -v dir="$OUTPUT_DIR" '/^#/{f=dir "/" $2 ".txt"} {print > f}' "$ORIGINAL_FILE"

# Skapa CSV-huvud dynamiskt baserat på orden i listan
# Detta skapar: Date,Apollo,moon,landing
header="Date"
for w in "${WORDS[@]}"; do header="$header,$w"; done
echo "$header" > $RESULT_CSV

for fil in $OUTPUT_DIR/*.txt; do
    date=$(basename "$fil" .txt)
    
    # Starta raden med datumet
    row="$date"

    for w in "${WORDS[@]}"; do
        # Räkna varje ord individuellt (-w ser till att det är exakta ord)
        count=$(grep -icw "$w" "$fil")
        # Lägg till resultatet i raden
        row="$row,$count"
    done

    # Spara hela raden till CSV-filen
    echo "$row" >> $RESULT_CSV
    echo "Bearbetat $date"
done
