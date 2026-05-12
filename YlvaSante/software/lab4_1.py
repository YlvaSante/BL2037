# LAB 4
# CORE TASK 1: Automating the script from lab 3 by using external file specifications

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import matplotlib

matplotlib.use('Agg')  # Agg is a "non-interactive backend" that only writes to files

def main():
    # 1. Initiating argparse
    parser = argparse.ArgumentParser(description="Creating histogram from gene expression data via kommandoraden.")

    # 2. Adding two arguments defining the input and output paths: args.input and args.output, resp.
    parser.add_argument("input_data", help="Pathway and name to the input file (.csv)")
    parser.add_argument("output_results", help="Pathway and name for the output file (.png)")

    # 3. Parsing arguments
    args = parser.parse_args()

### CODING FROM LAB 3 ###

    # 4. Loading data
    file_path = args.input_data
    df = pd.read_csv(file_path)

    # 5. Printing stats (same as lab 3)
    print(df.describe())

    # 6. Creating histogram for column 'TPM' 
    # (Note: Task is saying 'expression_level', though the file is named 'TPM')
    plt.figure(figsize=(10, 6))
    plt.hist(df['TPM'], bins=20, color='skyblue', edgecolor='black')

    # 7. Adding titles and labels
    plt.title("Gene Expression Distribution")
    plt.xlabel("expression_level (TPM)")
    plt.ylabel("Frequency")

    # 8. Saving the plot in the results/lab4/ folder
    # Checking if the folder exists, otherwise create it
    output_dir = os.path.dirname(args.output_results)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Saving image with high resolution (600 dpi)
    plt.savefig(args.output_results, dpi=600)
    print(f"Done✅ Graph has been saved in: {args.output_results}")

    # 9. Displaying that it is done
    # plt.show()
if __name__ == "__main__":
    main()

# Running the .py script by using the command line:
# python lab4_1.py ../data/gene_expression.csv ../results/hist.png
