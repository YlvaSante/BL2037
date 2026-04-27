# Lab 1 - extra tasks

# Extra Task 1: Skapa ett automatiserat Bash-script.
# Extra Task 2: Gör sökningen "smartare" (okänslig för stora/små bokstäver och flera ord).
# Extra Task 3: Skapa en snygg rapport och sortera den

# Skapar Bash-filen
#: nano extra_tasks.sh

#!/bin/bash

# --- Settings ---
ORIGINAL_FILE="/drives/c/Users/Admin/Ylva/GitHub/Programming_for_bio/Data/Lab1/text_file.txt"
DATA_DIR="/drives/c/Users/Admin/Ylva/GitHub/Programming_for_bio/Data/Lab1/extratasks"
RES_FILE="/drives/c/Users/Admin/Ylva/GitHub/Programming_for_bio/Results/Lab1/extratasks/summary_report.txt"


# Generating folders if needed (-p = parents)
mkdir -p "$DATA_DIR"
mkdir -p "$(dirname "$RES_FILE")"

# --- Extra Task 1: Splitting the file automatization ---
echo "Splitting the original file - wait..."
awk -v dir="$DATA_DIR" '/^#/{filename=dir "/" $2 ".txt"} {print > filename}' "$ORIGINAL_FILE"

# --- Extra Task 2 & 3: Searching and generating report ---
echo "Searching for target Words, generating report file - wait..."

# Creating headers for report file (ET3)
echo "=== SUMMARY REPORT ===" > "$RES_FILE"
echo "Generated on: $(date)" >> "$RES_FILE"
echo "-----------------------" >> "$RES_FILE"
printf "%-20s %-15s %-5s\n" "Filename" "Word" "Count" >> "$RES_FILE"
echo "-----------------------" >> "$RES_FILE"

# Looping through splitted files
# Searching for target words "Apollo" OR "Planning" (ET2)
# -i = ignore case, -E = extended regex (for |)
# basename = extracting only filename, no directory address

TARGETS="Apollo|Planning"

for file in "$DATA_DIR"/*.txt; do
    # Fetching filename w/o full directory address
    base_name=$(basename "$file")
    
    # Counting queries (-E = extended search; | = AND and OR 
    count=$(grep -Ei "$TARGETS" "$file" | wc -l)
    
    # Printing outputs to report file (ET3)
    printf "%-20s %-15s %-5s\n" "$base_name" "Multiple" "$count" >> "$RES_FILE"
done

# --- Bonus Task: Sorting results by count (Extra Task 3) ---
# Saving header and sorting the lines/rows by numeric count (-n) reverse (-r) in 3d column (-k3)
echo -e "\n--- Sorted by Count (Highest First) ---" >> "$RES_FILE"
tail -n +6 "$RES_FILE" | sort -k3 -nr >> "$RES_FILE"

echo "Done! Check the report file using: cat $RES_FILE"

# Saving in Nano: Ctrl+O, Enter
# Exit Nano: Ctrl+X.

# Preparing the script for runs
#: chmod +x extra_tasks.sh

# Runnng
#: ./extra_tasks.sh
