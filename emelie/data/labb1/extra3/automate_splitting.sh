#!/bin/bash

ORIGINAL_FILE="text_file.txt"
TARGET_WORD="Apollo" "moon" "landing"
OUTPUT_DIR="resultat_filer"
RESULT_CSV="multi_counts.csv"
SUMMARY_REPORT="summary_report.txt"

mkdir -p $OUTPUT_DIR
rm -f $OUTPUT_DIR/*.txt
awk -v dir="$OUTPUT_DIR" '/^#/{f=dir "/" $2 ".txt"} {print > f}' "$ORIGINAL_FILE"

echo "Date,Word,Count" > $RESULT_CSV

echo "----- SUMMARY REPORT - SEARCH RESULTS -----" >> $SUMMARY_REPORT
echo "------------------------------------------" >> $SUMMARY_REPORT
printf "%-20s %-10s\n" "FILE NAME (DATE)" "COUNT" >> $SUMMARY_REPORT
echo "------------------------------------------" >> $SUMMARY_REPORT




for fil in $OUTPUT_DIR/*.txt; do
    # Extrahera bara datumet från filnamnet för snyggare output
    date=$(basename "$fil" .txt)

    # Räkna förekomsten
    number=$(grep -Eic "$TARGET_WORD" "$fil")

    # Skriv till skärmen (för att se att det händer nåt)
    echo "Processing $date: Found; $number "

    # Spara till CSV
    echo "$date,$TARGET_WORD,$number" >> $RESULT_CSV

	#prints into summary_report.text
    printf "%-20s %-10s\n" "$date" "$number" >> $SUMMARY_REPORT


done




