# pip install biopython # if not already installed
# sudo apt-get install clustalo # if not already installed

import subprocess
import os
from Bio import Phylo
from Bio import AlignIO

# 1. Defining absolute paths
DATA_DIR = "/home/ylva/BL2037/YlvaSante/data/lab5/"
SOFTWARE_DIR = "/home/ylva/BL2037/YlvaSante/software/lab5/"
RESULTS_DIR = "/home/ylva/BL2037/YlvaSante/results/lab5/"

# Ensuring the results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

# Creating the list of files by inputting/defining filenames (including the reference 8JIQ_R)
pdb_filenames = ["8JIQ_R.pdb", "8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]
input_files = [os.path.join(DATA_DIR, f) for f in pdb_filenames]

# Defining output files with absolute paths in the results directory
fasta_input = os.path.join(RESULTS_DIR, "receptors.fasta")
aligned_fasta = os.path.join(RESULTS_DIR, "aligned_receptors.fasta")
tree_file = os.path.join(RESULTS_DIR, "tree_receptors.nwk")

def extract_sequence_from_pdb(pdb_path):
    """
    Extracts the amino acid sequence from ATOM records in the PDB file.
    Maps 3-letter codes to 1-letter codes.
    """
    # Mapping 3-letter amino acid codes to 1-letter codes
    d3to1 = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
    }
    
    sequence = []
    if not os.path.exists(pdb_path):
        print(f"Warning: File not found {pdb_path}")
        return ""

    with open(pdb_path, 'r') as f:
        for line in f:
            # Looking for ATOM records for the Alpha Carbon (CA) to define the sequence
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_name = line[17:20].strip()
                sequence.append(d3to1.get(res_name, 'X'))
                
    return "".join(sequence)

# 2. Generating the FASTA file
print(f"Generating FASTA file at: {fasta_input}")
with open(fasta_input, "w") as fasta_out:
    for pdb_path in input_files:
        label = os.path.basename(pdb_path).replace(".pdb", "")
        sequence = extract_sequence_from_pdb(pdb_path)
        if sequence:
            fasta_out.write(f">{label}\n{sequence}\n")

# 3. Defining and executing the Clustal Omega command
# Ensure 'clustalo' is installed in the system environment
clustal_command = [
    "clustalo", 
    "-i", fasta_input, 
    "-o", aligned_fasta, 
    "--outfmt=fasta", 
    "--guidetree-out", tree_file, 
    "--force"
]

print("Executing Clustal Omega sequence alignment...")
try:
    # Running the subprocess and checking for errors
    result = subprocess.run(clustal_command, capture_output=True, text=True, check=True)
    print(f"Alignment completed successfully. Results saved in: {RESULTS_DIR}")
except subprocess.CalledProcessError as e:
    print(f"An error occurred during Clustal Omega execution: {e.stderr}")
except FileNotFoundError:
    print("Error: 'clustalo' was not found. Please install it using 'sudo apt-get install clustalo'.")

# 4. Reading and displaying the phylogenetic tree (.nwk) using Biopython
if os.path.exists(tree_file):
    print("\n" + "="*40)
    print("PHYLOGENETIC TREE (Sequence-based)")
    print("="*40)
    tree = Phylo.read(tree_file, "newick")
    # Drawing a simple ASCII representation of the tree in the terminal
    Phylo.draw_ascii(tree)
else:
    print("Tree file generation failed or file not found.")

# --- VISUALIZATION ---
print("\n" + "="*50)
print("PART 1: ALIGNMENT PREVIEW")
print("="*50)
if os.path.exists(aligned_fasta):
    alignment = AlignIO.read(aligned_fasta, "fasta")
    # Show the first 60 residues of the alignment for all proteins
    print(alignment[:, :60])

print("\n" + "="*50)
print("PART 2: PHYLOGENETIC TREE")
print("="*50)
if os.path.exists(tree_file):
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw_ascii(tree)
else:
    print("Tree file not found.")

    # --- NEW FUNCTION TO PRINT THE TREE (Biopython) ---
def print_tree():
    """Reading the Newick file and prints an ASCII representation to the terminal."""
    if os.path.exists(tree_file):
        print("\n" + "="*50)
        print("FINAL PHYLOGENETIC TREE (From Newick file)")
        print("="*50)
        
        # Biopython reads the tree file
        tree = Phylo.read(tree_file, "newick")
        
        # Biopython prints the tree to the terminal
        Phylo.draw_ascii(tree)
        
        print("\nNote: Horizontal line length represents evolutionary distance.")
    else:
        print(f"Error: The file {tree_file} was not found.")

# --- FINAL VISUALIZATION CALLS ---

# Alignment Preview
if os.path.exists(aligned_fasta):
    print("\n" + "="*50)
    print("PART 1: ALIGNMENT PREVIEW")
    print("="*50)
    alignment = AlignIO.read(aligned_fasta, "fasta")
    print(alignment[:, :60])

# Call the function to print the tree
print_tree()