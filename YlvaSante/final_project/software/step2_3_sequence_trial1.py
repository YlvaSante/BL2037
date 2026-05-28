# Installing Clustal Omega via micromamba, if not already installed. This is a one-time setup step.
# Clustal Omega is a powerful tool for multiple sequence alignment
# Uncomment the line below to install Clustal Omega if needed
# micromamba install -c bioconda clustalo -y

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# --- ABSOLUTE PATHS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

FASTA_INPUT = os.path.join(PROJECT_ROOT, "data", "sequences.fasta")
ALIGNED_FASTA = os.path.join(PROJECT_ROOT, "data", "aligned_sequences.fasta")
PIM_MATRIX_OUT = os.path.join(PROJECT_ROOT, "data", "sequence_pim.mat")
PLOT_DIR = os.path.join(PROJECT_ROOT, "results", "plots")

os.makedirs(PLOT_DIR, exist_ok=True)

def run_sequence_analysis():
    print("--- Initialising Step 2 & 3: Sequence Alignment & Clustering ---")
    
    # 1. Calling Clustal Omega via subprocess
    # Building the command exactly as it would be typed in the terminal
    clustalo_cmd = [
        "clustalo",
        "-i", FASTA_INPUT, # Flagging the input FASTA file with unaligned sequences
        "-o", ALIGNED_FASTA, # Flagging the output FASTA file for aligned sequences
        "--distmat-out=" + PIM_MATRIX_OUT,
        "--full",
        "--force"
    ]
    
    print("Calling Clustal Omega via subprocess...")
    try:
        subprocess.run(clustalo_cmd, check=True)
        print("-> Clustal Omega succeeded without errors!")
    except subprocess.CalledProcessError as e:
        print(f"!! ERROR: Clustal Omega failed: {e}")
        return

    # 2. Parsing Percent Identity Matrix (PIM)
    # Clustals PIM format has protein names in the first column and numbers separated by spaces
    print("Parsing the distance matrix...")
    
    with open(PIM_MATRIX_OUT, "r") as f:
        lines = f.readlines()
    
    # The first row often contains the number of sequences 
    # If first row is a number >> start parsing from the second row, otherwise start from the first row
    start_row = 1 if len(lines[0].split()) == 1 else 0
    
    labels = []
    matrix_data = []
    
    for line in lines[start_row:]:
        parts = line.split()
        if not parts:
            continue
        labels.append(parts[0]) # First element is protein ID (e.g. 1tup_c)
        matrix_data.append([float(x) for x in parts[1:]])
        
    # Creating a DataFrame with percentage unit identity (0-100)
    # Clustal Omega often gives distance or identity in percentage unit >> normalising to 0-100%
    pim_df = pd.DataFrame(matrix_data, index=labels, columns=labels)
    
    # Saving matrix as CSV for future use in the report
    pim_df.to_csv(os.path.join(PROJECT_ROOT, "data", "sequence_identity_matrix.csv"))

    # 3. Generating a Heatmap illustrating sequence similarity
    plt.figure(figsize=(12, 10))
    sns.heatmap(pim_df, annot=False, cmap="viridis", cbar_kws={'label': 'Percent Sequence Identity'})
    plt.title("Sequence Similarity Heatmap (Clustal Omega)")
    plt.tight_layout()
    heatmap_path = os.path.join(PLOT_DIR, "sequence_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"-> Saved heatmap to: {heatmap_path}")

    # 4. Generating a Clustering Tree (Dendrogram)
    # For clustering a DISTANCE matrix is required. The higher the identity, the shorter the distance.
    # Distance = 100 - Percentage identity
    distance_matrix = 100 - pim_df.values
    
    # Ensuring the matrix is symmetric (small rounding errors can occur in Clustal)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    np.fill_diagonal(distance_matrix, 0)
    
    # Converting to condensed form required for hierarchical clustering (i.e. upper triangular values in a 1D array)
    condensed_distance = squareform(distance_matrix)
    
    # Calculating hierarchical clustering (UPGMA / average linkage)
    Z = linkage(condensed_distance, method='average')
    
    plt.figure(figsize=(12, 6))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=10)
    plt.title("Sequence-Based Clustering Tree (Dendrogram)")
    plt.ylabel("Distance (100 - % Identity)")
    plt.tight_layout()
    dendrogram_path = os.path.join(PLOT_DIR, "sequence_dendrogram.png")
    plt.savefig(dendrogram_path, dpi=300)
    plt.close()
    print(f"-> Saved clustering tree to: {dendrogram_path}")
    print("\nStep 2 & 3 Completed!")

if __name__ == "__main__":
    run_sequence_analysis()