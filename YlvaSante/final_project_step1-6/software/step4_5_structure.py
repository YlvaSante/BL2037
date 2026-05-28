#!/usr/bin/env python3
# Installing TM-align via micromamba, if not already installed. This is a one-time setup step.
# TM-align is a powerful tool for structural alignment and RMSD calculation.
# Uncomment the line below to install TM-align if needed
# micromamba install -c bioconda tmalign -y

import os
import glob
import subprocess
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# --- DEFINING ABSOLUTE PATHSWAYS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) # Defining the directory of the current script, i.e. <this> file's location. This allows us to build paths relative to the script's location, ensuring it works regardless of where it's run from.
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CLEAN_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "clean_pdb")
PLOT_DIR = os.path.join(PROJECT_ROOT, "results", "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

def run_tmalign(pdb1, pdb2):
    """Running TMalign for two PDB files and extracts RMSD."""
    try:
        # Running TMalign command
        result = subprocess.run(["TMalign", pdb1, pdb2], capture_output=True, text=True, check=True)
        output = result.stdout
        
        # Using regex to find the line: "RMSD=  2.34,     "
        rmsd_match = re.search(r"RMSD=\s*([\d\.]+)", output)
        if rmsd_match:
            return float(rmsd_match.group(1))
    except Exception as e:
        print(f"!! Error running TM-align for {os.path.basename(pdb1)} and {os.path.basename(pdb2)}: {e}")
    return None

def run_structural_analysis():
    print("--- Initiating Step 4 & 5: Structural Alignment with TM-align ---")
    
    # Fetching all cleaned PDB files and preparing labels for the matrix and plots
    pdb_files = sorted(glob.glob(os.path.join(CLEAN_PDB_DIR, "*_clean.pdb")))
    labels = [os.path.basename(f).replace("_clean.pdb", "") for f in pdb_files]
    n = len(pdb_files)
    
    if n < 2:
        print(f"!! Error: Found only {n} cleaned PDB files in {CLEAN_PDB_DIR}. Need at least 2.")
        return

    print(f"Found {n} cleaned protein structures. Calculating pairwise RMSD values...")
    
    # Inititiating empty matrix for RMSD
    rmsd_matrix = np.zeros((n, n))
    
    # Running all pairwise comparisons (Matrix filled symmetrically, diagonal = 0)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                rmsd_matrix[i, j] = 0.0
            else:
                rmsd = run_tmalign(pdb_files[i], pdb_files[j])
                # In case that TM-align fails completely, the RMSD value is set to a high "penalty value" (e.g., 25 Å)
                if rmsd is None:
                    rmsd = 25.0 
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        print(f"  -> Successfully completed pairwise comparisons for protein {i+1}/{n}: {labels[i]}")

    # Creating DataFrame for analysis and saving the RMSD matrix to CSV for reference
    rmsd_df = pd.DataFrame(rmsd_matrix, index=labels, columns=labels)
    rmsd_df.to_csv(os.path.join(PROJECT_ROOT, "data", "structural_rmsd_matrix.csv"))

    # --- PLOT 1: Structural Heatmap ---
    plt.figure(figsize=(12, 10))
    # Using an inverse color scale (dark blue = low RMSD/high similarity, light = high RMSD/different)
    sns.heatmap(rmsd_df, annot=False, cmap="YlGnBu", cbar_kws={'label': 'RMSD (Ångström)'})
    plt.title("Structural Dissimilarity Heatmap (Pairwise RMSD via TM-align)")
    plt.tight_layout()
    heatmap_path = os.path.join(PLOT_DIR, "structural_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"-> Saved structural heatmap to: {heatmap_path}")

    # --- PLOT 2: Structural Dendrogram (Clustering) ---
    print("Calculating hierarchical clustering for structural dendrogram...")
    
    # Since RMSD is already a direct measure of distance (0 = identical structure), the matrix can be used directly for clustering!
    try:
        # Ensuring diagonal and symmetry of the distance matrix   
        np.fill_diagonal(rmsd_matrix, 0)
        condensed_distance = squareform(rmsd_matrix)
        
        # Calculating hierarchical clustering (average linkage / UPGMA)
        Z = linkage(condensed_distance, method='average')
        
        plt.figure(figsize=(14, 7))
        dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=10, link_color_func=lambda k: 'C1')
        plt.title("Structure-Based Clustering Tree (Dendrogram via TM-align RMSD)")
        plt.ylabel("Structural Distance (RMSD in Å)")
        plt.xlabel("Proteins")
        plt.tight_layout()
        
        dendrogram_path = os.path.join(PLOT_DIR, "structural_dendrogram.png")
        plt.savefig(dendrogram_path, dpi=300)
        plt.close()
        print(f"-> Saved structural dendrogram to: {dendrogram_path}")
        
    except Exception as e:
        print(f"!! Error calculating structural dendrogram: {e}")

    print("\nStep 4 & 5 Complete    !")

if __name__ == "__main__":
    run_structural_analysis()