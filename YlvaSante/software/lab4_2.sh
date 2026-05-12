#!/bin/bash

# 1. Defining variables for input and output
# This makes the script easy to read and simple to modify later
INPUT_FILE="../data/gene_expression.csv"
OUTPUT_FILE="../results/lab4/hist.png"

# 2. Running the python script using the variables as arguments
# Using python3 to ensure it runs with Python 3, which is the version in use
echo "Starting analysis of $INPUT_FILE..."
python3 lab4_1.py $INPUT_FILE $OUTPUT_FILE

echo "Done! Check the folder ../results/lab4/ for the results."