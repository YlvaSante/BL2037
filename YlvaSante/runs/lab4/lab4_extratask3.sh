#!/bin/bash

LIST="/home/ylva/BL2037/YlvaSante/data/lab4/list.txt"
MATCH_CSV="/home/ylva/BL2037/YlvaSante/results/lab4/MATCH.csv"
UNMATCH_CSV="/home/ylva/BL2037/YlvaSante/results/lab4/UNMATCHED.csv"
PDB_DIR="/home/ylva/BL2037/YlvaSante/data/lab4/pdb_files"

mkdir -p $PDB_DIR
# mkdir -p /home/ylva/BL2037/YlvaSante/results

# Emptying the output .csv files if already existing, to avoid appending to old results
# rm -f $MATCH_CSV $UNMATCH_CSV

# Overwriting previous .csv files with headers
echo "UniProt_ID,PDB_Len,UniProt_Len" > $MATCH_CSV
echo "UniProt_ID,PDB_Len,UniProt_Len" > $UNMATCH_CSV

# Looping through the list with 20 ID:s
while read -r UniProt_ID; do
    # Skipping empty lines
    [ -z "$UniProt_ID" ] && continue
    
    echo "--- Processing $UniProt_ID ---"

    # STEP A: Downloading .pdb file if non existing in the pdb_files directory, using the fetch_pdb.py script
    if [ ! -f "$PDB_DIR/${UniProt_ID}.pdb" ]; then
        python /home/ylva/BL2037/YlvaSante/software/lab4/fetch_pdb.py "$UniProt_ID" "$PDB_DIR"
    fi
    
    PDB_FILE="$PDB_DIR/${UniProt_ID}.pdb" # Building the expected path to the .pdb file for this UniProt ID
    if [ -f "$PDB_FILE" ]; then
        python /home/ylva/BL2037/YlvaSante/software/lab4/validate.py \
                           --id "$UniProt_ID" \
                           --pdb_file "$PDB_FILE" \
                           --match_csv "$MATCH_CSV" \
                           --unmatch_csv "$UNMATCH_CSV"
    else
        echo "Warning: No PDB file found for $UniProt_ID in $PDB_DIR"
    fi

done < "$LIST"