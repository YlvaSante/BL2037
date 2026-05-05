import pandas as pd
import os

# 1. Defining the pathway to the .csv file (standing in software for running the script, fetching the data file from the data folder)
file_path = "../../data/gene_expression.csv"

# Checking whether the file exists before trying to load it, to avoid errors and give a helpful message if it's not found
if not os.path.exists(file_path):
    print(f"ERROR: did not find the data file in {file_path}")
    print("Solution: Running the script from the correct directory (software) and fetching the original data file (gene_expression.csv) from the data folder.")
else:
    # 2. Loading data - sep=None is making it insensitive to whether comma or semicolon is used as separator, engine='python' is needed for sep=None to work 
    df = pd.read_csv(file_path, sep=None, engine='python')

    # 3. Cleaning the column names from invisible whitespaces
    df.columns = df.columns.str.strip()

    # 4. Printing column names to check whether everything looks correct (and to know which columns we have to work with)
    print("Analyzing the following columns:", df.columns.tolist())

    # 5. Print Summary Statistics (mean, max, min)
    print("\n--- Summary Statistics ---")
    print(df.describe())

    # 6. Filter: Creating new DataFrame with genes where TPM > 100
    # Checking whether 'TPM' exists in the columns before trying to filter, to avoid errors and give a helpful message if it's not found
    if 'TPM' in df.columns:
        filtered_df = df[df['TPM'] > 100]
        print("\n--- Genes with TPM > 100 ---")
        print(filtered_df)
    else:
        print("\nERROR: Column 'TPM' is not found. Please check the spelling in the CSV file.")