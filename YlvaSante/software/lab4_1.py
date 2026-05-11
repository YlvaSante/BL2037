# LAB 4
# CORE TASK 1: Automating the script from lab 3 by using external file specifications

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def main():
    # 1. Initiating argparse
    parser = argparse.ArgumentParser(description="Creating histogram from gene expression data via kommandoraden.")

    # 2. Adding two arguments defining the position: args.input and args.output, resp.
    parser.add_argument("input", help="Sökväg till input-filen (CSV)")
    parser.add_argument("output", help="Sökväg och namn för output-filen (PNG)")

    # 3. Parsa argumenten
    args = parser.parse_args()

### CODING FROM LAB 3 ###

# 5. Reaing data from CSV file
df = pd.read_csv(args.input)

# Creating a list of amino acids
amino_acids = ["Methionine", "Leucine", "Glycine", "Alanine"]

# Looping through the list of amino acids and printing each one 
for aa in amino_acids:
    # f"..." gör att vi kan skriva in variabeln {aa} direkt i texten
    print(f"Processing amino acid: {aa}")