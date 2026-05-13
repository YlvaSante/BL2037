# pip install biopython # Ensure this is installed
# sudo apt-get install clustalo # Ensure Clustal Omega is installed

import subprocess
import os
from Bio import Phylo
from Bio import AlignIO

# 1. DEFINING ABSOLUTE PATHS
# All directories are defined as absolute paths to ensure reliability
DATA_DIR = "/home/ylva/BL2037/YlvaSante/data/lab5/"
SOFTWARE_DIR = "/home/ylva/BL2037/YlvaSante/software/lab5/"
RESULTS_DIR = "/home/ylva/BL2037/YlvaSante/results/lab5/"

# Ensure the results directory exists
os.makedirs(RESULTS_DIR, exist_ok=True)

# Define input filenames (Reference + Targets)
pdb_filenames = ["8JIQ_R.pdb", "8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]
input_files = [os.path.join(DATA_DIR, f) for f in pdb_filenames]

# Define output files with absolute paths in the results directory
fasta_input = os.path.join(RESULTS_DIR, "receptors.fasta")
aligned_fasta = os.path.join(RESULTS_DIR, "aligned_receptors.fasta")
tree_file = os.path.join(RESULTS_DIR, "tree_receptors.nwk")

# 2. SEQUENCE EXTRACTION LOGIC
def extract_sequence_from_pdb(pdb_path):
    """
    Extracts the amino acid sequence from ATOM records in the PDB file.
    Maps 3-letter amino acid codes to 1-letter codes.
    """
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
            # Look for Alpha Carbon (CA) atoms to extract one residue at a time
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_name = line[17:20].strip()
                sequence.append(d3to1.get(res_name, 'X'))
                
    return "".join(sequence)

# 3. GENERATING THE INPUT FASTA FILE
print(f"Step 1: Extracting sequences into FASTA format: {fasta_input}")
with open(fasta_input, "w") as fasta_out:
    for pdb_path in input_files:
        label = os.path.basename(pdb_path).replace(".pdb", "")
        sequence = extract_sequence_from_pdb(pdb_path)
        if sequence:
            fasta_out.write(f">{label}\n{sequence}\n")

# 4. EXECUTING CLUSTAL OMEGA VIA SUBPROCESS
clustal_command = [
    "clustalo", 
    "-i", fasta_input, 
    "-o", aligned_fasta, 
    "--outfmt=fasta", 
    "--guidetree-out", tree_file, 
    "--force"
]

print("Step 2: Executing Clustal Omega sequence alignment...")
try:
    # Run the external tool and capture output/errors
    subprocess.run(clustal_command, capture_output=True, text=True, check=True)
    print(f"Step 3: Alignment completed. Files saved in: {RESULTS_DIR}")
except subprocess.CalledProcessError as e:
    print(f"Error: Clustal Omega failed. {e.stderr}")
except FileNotFoundError:
    print("Error: 'clustalo' not found. Please check your installation.")

# 5. VISUALIZATION FUNCTIONS
def print_tree():
    """
    Uses Biopython to read the Newick tree file and display it as ASCII text.
    """
    if os.path.exists(tree_file):
        print("\n" + "="*50)
        print("FINAL PHYLOGENETIC TREE (Sequence-based Distance)")
        print("="*50)
        
        # Parse the Newick format and draw it in the terminal
        tree = Phylo.read(tree_file, "newick")
        Phylo.draw_ascii(tree)
        
        print("\nNote: Horizontal line length represents evolutionary distance.")
    else:
        print(f"Error: Tree file {tree_file} not found.")

# 6. FINAL OUTPUT DISPLAY
if os.path.exists(aligned_fasta):
    print("\n" + "="*50)
    print("PART 1: SEQUENCE ALIGNMENT PREVIEW (First 60 columns)")
    print("="*50)
    alignment = AlignIO.read(aligned_fasta, "fasta")
    print(alignment[:, :60])

# Call the tree visualization function
print_tree()

print(f"\nTask 2 finished. All outputs are stored in {RESULTS_DIR}")