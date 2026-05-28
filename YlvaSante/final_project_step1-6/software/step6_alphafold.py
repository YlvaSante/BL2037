import os
import glob
import json
import urllib.request
import subprocess
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# --- DEFINING ABSOLUTE PATHWAYS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CLEAN_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "clean_pdb")
AF_RAW_DIR = os.path.join(PROJECT_ROOT, "data", "alphafold_raw")
AF_CLEAN_DIR = os.path.join(PROJECT_ROOT, "data", "alphafold_clean")
PLOT_DIR = os.path.join(PROJECT_ROOT, "results", "plots")

os.makedirs(AF_RAW_DIR, exist_ok=True)
os.makedirs(AF_CLEAN_DIR, exist_ok=True)

def get_uniprot_id(pdb_id):
    """Fetching UniProt ID for the given PDB-ID via PDBe API (EBI)."""
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            # PDBe API returns data structured under the PDB-id
            if pdb_id.lower() in data:
                uniprot_dict = data[pdb_id.lower()]["UniProt"]
                # Retrieving the first UniProt ID in the list (the key in the search result)
                if uniprot_dict:
                    return list(uniprot_dict.keys())[0]
    except Exception as e:
        # If PDBe API would fail, then trying a secondary public API (Uniprot org)
        try:
            url = f"https://rest.uniprot.org/uniprotkb/search?query=pdb:{pdb_id.upper()}&fields=accession"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode())
                if "results" in data and len(data["results"]) > 0:
                    return data["results"][0]["primaryAccession"]
        except Exception as e2:
            print(f"Could not fetch UniProt for {pdb_id}: {e2}")
    return None

def download_alphafold_pdb(uniprot_id, output_path):
    """Calling AlphaFold API and downloading the predicted PDB file."""
    try:
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            if data and len(data) > 0:
                # AlphaFold returns a list, where the first entry contains the PDB URL 
                pdb_url = data[0]["pdbUrl"]
                urllib.request.urlretrieve(pdb_url, output_path)
                return True
    except Exception as e:
        print(f"AlphaFold API failed for {uniprot_id}: {e}")
    return False

def run_alphafold_pipeline():
    print("=== INITIATING STEP 6: ALPHAFOLD ANALYSIS ===")
    
    # 1. Searching for the same proteins successfully run and analyzed in Step 1
    pdb_files = sorted(glob.glob(os.path.join(CLEAN_PDB_DIR, "*_clean.pdb")))
    identifiers = [os.path.basename(f).replace("_clean.pdb", "") for f in pdb_files]
    
    af_mapping = {}
    fasta_records = []
    
    print(f"Mapping {len(identifiers)} working PDB-ids to UniProt and downloading from AlphaFold...")
    
    for identifier in identifiers:
        pdb_id = identifier.split("_")[0]
        print(f"\nProcessing {identifier}...")
        
        uniprot_id = get_uniprot_id(pdb_id)
        if not uniprot_id:
            print(f"!! Skipping {identifier}: No UniProt ID found.")
            continue
            
        af_raw_path = os.path.join(AF_RAW_DIR, f"{identifier}_af.pdb")
        af_clean_path = os.path.join(AF_CLEAN_DIR, f"{identifier}_af_clean.pdb")
        
        # Downloading from AlphaFold
        if download_alphafold_pdb(uniprot_id, af_raw_path):
            print(f"-> Downloaded AlphaFold structure for {uniprot_id}")
            
            # Cleaning the AlphaFold file (same as in Step 1, removing heteroatoms if existing)
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(uniprot_id, af_raw_path)
            io = PDBIO()
            
            class CleanStructure(Select):
                def accept_residue(self, residue):
                    return 1 if residue.id[0] == " " else 0
                    
            io.set_structure(structure)
            io.save(af_clean_path, CleanStructure())
            
            # Extracting the sequence from the cleaned AlphaFold structure
            clean_struct = parser.get_structure(uniprot_id, af_clean_path)
            sequence = ""
            for model in clean_struct:
                for chain in model:
                    for residue in chain:
                        sequence += seq1(residue.get_resname())
            
            if sequence:
                fasta_records.append(f">{identifier}_af\n{sequence}\n")
                af_mapping[identifier] = af_clean_path
                print(f"-> Saved cleaned AF sequence ({len(sequence)} aa)")
        else:
            print(f"!! Failed to fetch AlphaFold data for {uniprot_id}")

    # Saving AlphaFold FASTA
    af_fasta_out = os.path.join(PROJECT_ROOT, "data", "alphafold_sequences.fasta")
    with open(af_fasta_out, "w") as f:
        f.writelines(fasta_records)
    print(f"\n-> AlphaFold FASTA saved to {af_fasta_out}")

    # ==========================================
    # 2. APPLYING STEP 2 & 3: Clustal Omega on AlphaFold
    # ==========================================
    print("\n--- Running Clustal Omega on AlphaFold sequences ---")
    af_aligned_fasta = os.path.join(PROJECT_ROOT, "data", "alphafold_aligned.fasta")
    af_pim_out = os.path.join(PROJECT_ROOT, "data", "alphafold_pim.mat")
    
    subprocess.run([
        "clustalo", "-i", af_fasta_out, "-o", af_aligned_fasta,
        f"--distmat-out={af_pim_out}", "--full", "--force"
    ], check=True)
    
    # Parsing and drawing tree for AF sequences
    with open(af_pim_out, "r") as f:
        lines = f.readlines()
    start_row = 1 if len(lines[0].split()) == 1 else 0
    labels, matrix_data = [], []
    for line in lines[start_row:]:
        parts = line.split()
        if parts:
            labels.append(parts[0])
            matrix_data.append([float(x) for x in parts[1:]])
            
    af_pim_df = pd.DataFrame(matrix_data, index=labels, columns=labels)
    if af_pim_df.values.max() <= 1.0:
        af_pim_df = af_pim_df * 100
        
    # Drawing AF-dendrogram
    distance_matrix = 100 - af_pim_df.values
    np.fill_diagonal(distance_matrix, 0)
    Z = linkage(squareform(np.clip(distance_matrix, 0, 100)), method='average')
    
    plt.figure(figsize=(14, 7))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=10, link_color_func=lambda k: 'C2')
    plt.title("AlphaFold Sequence-Based Dendrogram")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "alphafold_sequence_dendrogram.png"), dpi=300)
    plt.close()

    # ==========================================
    # 3. APPLYING STEP 4 & 5: TM-align on AlphaFold
    # ==========================================
    print("\n--- Running TM-align on AlphaFold structures ---")
    af_keys = list(af_mapping.keys())
    n_af = len(af_keys)
    af_rmsd_matrix = np.zeros((n_af, n_af))
    
    for i in range(n_af):
        for j in range(i, n_af):
            if i == j:
                af_rmsd_matrix[i, j] = 0.0
            else:
                p1 = af_mapping[af_keys[i]]
                p2 = af_mapping[af_keys[j]]
                res = subprocess.run(["TMalign", p1, p2], capture_output=True, text=True, check=True)
                match = re.search(r"RMSD=\s*([\d\.]+)", res.stdout)
                rmsd = float(match.group(1)) if match else 25.0
                af_rmsd_matrix[i, j] = rmsd
                af_rmsd_matrix[j, i] = rmsd
                
    af_labels = [f"{k}_af" for k in af_keys]
    af_rmsd_df = pd.DataFrame(af_rmsd_matrix, index=af_labels, columns=af_labels)
    
    # Drawing AF structure dendrogram
    np.fill_diagonal(af_rmsd_matrix, 0)
    Z_struct = linkage(squareform(af_rmsd_matrix), method='average')
    
    plt.figure(figsize=(14, 7))
    dendrogram(Z_struct, labels=af_labels, leaf_rotation=90, leaf_font_size=10, link_color_func=lambda k: 'C3')
    plt.title("AlphaFold Structure-Based Dendrogram (RMSD)")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "alphafold_structural_dendrogram.png"), dpi=300)
    plt.close()
    
    print("\n=== STEG 6 KLART! Alla AlphaFold-grafer sparade ===")

if __name__ == "__main__":
    run_alphafold_pipeline()