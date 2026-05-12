#!/bin/bash

# Creating directory for results if non existent
mkdir -p /home/ylva/BL2037/YlvaSante/results/lab4

# Looping through all .csv files in the examplesdirectory
for input_file in /home/ylva/BL2037/YlvaSante/lab4/gene_expression*.csv; do
    
    # Creating unique name for the image based on the filename, so that:
    # if input file = gene_expression2.csv, then output file -> gene_expression2_hist.png
    filename=$(basename "$input_file" .csv)
    output_file="/home/ylva/BL2037/YlvaSante/results/lab4/${filename}_hist.png"
    
    echo "Processing $input_file => $output_file"
    
    # Running the Python script for each file, with the input and output file as arguments
    # and using the flags defined in the python script
    python /home/ylva/BL2037/YlvaSante/software/lab4/lab4_1.py --input_data "$input_file" --output_results "$output_file"
done

echo "✅ All files are processed and analyzed!"