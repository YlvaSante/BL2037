#!/bin/bash

# 1. Defining variables for input and output
# This makes the script easy to read and simple to modify later
INPUT_FILE="/home/ylva/BL2037/YlvaSante/data/gene_expression.csv"
OUTPUT_FILE="/home/ylva/BL2037/YlvaSante/results/lab4/hist.png"

# 2. Running the python script using the variables as arguments
# plus using the flags defined in the python script
# adding --input_data and --output_results in front of the variables to match the flags in the python script
# Using python3 to ensure it runs with Python 3, which is the version in use
echo "Starting analysis of $INPUT_FILE..."
python /home/ylva/BL2037/YlvaSante/software/lab4/lab4_1.py --input_data "$INPUT_FILE" --output_results "$OUTPUT_FILE"

echo "Done! Check the folder ../results/lab4/ for the results."