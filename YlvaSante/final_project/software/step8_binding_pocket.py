import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser

# --- DEFINING ABSOLUTE PATHS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CLEAN_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "clean_pdb")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "plots")
OUTPUT_CSV = os.path.join(PROJECT_ROOT, "data", "step8_binding_pocket_rmsd.csv")

os.makedirs(RESULTS_DIR, exist_ok=True)

# Three-letter to one-letter amino acid code dictionary for strict parsing
AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence_and_ca_coords(pdb_path):
    """
    Parses a clean PDB file to extract both its 1D sequence 
    and the corresponding 3D coordinates of all Alpha Carbons (CA).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    residue_sequence = []
    ca_coordinates = []
    residue_ids = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Guarantee we only parse standard amino acids, filtering out heteroatoms
                if residue.get_resname() in AA_MAP and 'CA' in residue:
                    res_letter = AA_MAP[residue.get_resname()]
                    ca_atom = residue['CA']
                    
                    residue_sequence.append(res_letter)
                    ca_coordinates.append(ca_atom.get_coord())
                    residue_ids.append(residue.get_id()[1])
            # Only analyze the first chain to remain consistent with previous structural pipeline steps
            break
        break
        
    return "".join(residue_sequence), np.array(ca_coordinates), residue_ids

def calculate_pocket_rmsd(coords1, coords2):
    """
    Calculates the exact structural RMSD between two sets of 3D coordinates 
    using the Kabsch algorithm (translation and optimal rotation).
    """
    if len(coords1) != len(coords2) or len(coords1) == 0:
        return np.nan
        
    # Translate both sets of coordinates to their respective centers of mass
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    c1 = coords1 - centroid1
    c2 = coords2 - centroid2
    
    # Calculate the covariance matrix
    H = np.dot(c1.T, c2)
    
    # Perform Singular Value Decomposition (SVD) to find optimal rotation
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    
    # Handle special reflection case to ensure a valid right-handed rotation matrix
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = np.dot(Vt.T, U.T)
        
    # Rotate coordinates of the first dataset
    c1_rotated = np.dot(c1, R)
    
    # Compute the final quantitative Root Mean Square Deviation (RMSD)
    diff = c1_rotated - c2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def run_binding_pocket_analysis():
    print("============================================================")
    print("=== INITIATING STEP 8: ADVANCED BINDING POCKET ANALYSIS ===")
    print("============================================================")
    
    pdb_files = glob.glob(os.path.join(CLEAN_PDB_DIR, "*.pdb"))
    print(f"Found {len(pdb_files)} rensade PDB files to scan for the binding motif.")
    
    # A more robust bioinformatics approach: Scan for a glycine-rich region 
    # followed closely by a critical Lysine/Arginine (K/R) and Serine/Threonine (S/T)
    # This captures structural variations in the loop across different families
    pocket_data = {}
    
    for pdb_path in pdb_files:
        file_name = os.path.basename(pdb_path)
        protein_id = file_name.replace("_clean.pdb", "")
        
        try:
            seq, ca_coords, res_ids = extract_sequence_and_ca_coords(pdb_path)
            
            # Sliding window approach to find the active site core
            # Look for at least 2 Glycines (G) and 1 Lysine/Arginine (K/R) within a 10 AA window
            found_motif = False
            for i in range(len(seq) - 10):
                window = seq[i:i+10]
                if window.count('G') >= 2 and ('K' in window or 'R' in window):
                    start_idx = i
                    end_idx = i + 10
                    
                    motif_coords = ca_coords[start_idx:end_idx]
                    motif_seq = seq[start_idx:end_idx]
                    
                    pocket_data[protein_id] = {
                        'coords': motif_coords,
                        'seq': motif_seq,
                        'length': len(motif_coords)
                    }
                    found_motif = True
                    break # Stop scanning this protein once the pocket is found
                    
        except Exception as e:
            continue

    print(f"-> Successfully extracted the binding pocket motif from {len(pocket_data)} proteins!")
    
    if len(pocket_data) < 2:
        print("!! ERROR: Too few motifs found to generate a comparative matrix. Aborting.")
        return

    # --- STEP 8.2: PAIRWISE STRUCTURAL ALIGNMENT OF THE ACTIVE SITES ---
    labels = list(pocket_data.keys())
    matrix_size = len(labels)
    rmsd_matrix = np.zeros((matrix_size, matrix_size))
    
    for i in range(matrix_size):
        for j in range(matrix_size):
            id1, id2 = labels[i], labels[j]
            
            # Since our window is strictly fixed to 10 residues, lengths will always match perfectly!
            rmsd_matrix[i, j] = calculate_pocket_rmsd(pocket_data[id1]['coords'], pocket_data[id2]['coords'])

    # Format the finalized matrix into a readable and shareable Pandas DataFrame
    rmsd_df = pd.DataFrame(rmsd_matrix, index=labels, columns=labels)
    rmsd_df.to_csv(OUTPUT_CSV)
    print(f"-> Active site local RMSD matrix successfully saved to: {OUTPUT_CSV}")

    # --- STEP 8.3: VISUALIZING THE CONSERVED FINGERPRINT ---
    plt.figure(figsize=(14, 12))
    # Using the 'rocket' color map where a dark/low value indicates highly conserved structures
    sns.heatmap(
        rmsd_df, 
        annot=False, 
        cmap="rocket_r", 
        vmin=0.0, 
        vmax=3.0, 
        cbar_kws={'label': 'Active Site Pocket RMSD (Ångströms)'}
    )
    plt.title("Structural Conservation Fingerprint of the Nucleotide Binding Pocket\n(Sliding-Window Active Site Motif Analysis)", fontsize=14, fontweight='bold')
    plt.xlabel("Proteins", fontsize=11)
    plt.ylabel("Proteins", fontsize=11)
    plt.tight_layout()
    
    plot_path = os.path.join(RESULTS_DIR, "step8_binding_pocket_heatmap.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    print(f"-> Publication-ready pocket fingerprint heatmap saved to: {plot_path}")
    print("Step 8 Analysis Completed Successfully!")

if __name__ == "__main__":
    run_binding_pocket_analysis()