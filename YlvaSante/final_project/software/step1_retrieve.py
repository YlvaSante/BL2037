import os
import pandas as pd
from Bio.PDB import PDBList, PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1

# Downloading the input CSV file to the project directory, if not yet downloaded
# wget -O /home/ylva/BL2037/YlvaSante/final_project/inputs_finalproject.csv "https://docs.google.com/uc?export=download&id=1Je60M38QnlohbUUjmuQDph_qLT7eh85r"

# --- AUTOMATISATION OF ABSOLUTE PATHS ---
# OS finds the absolute path to <this file>, defines the directory, and goes back one level to the project's root directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# Defining absolute paths needed for the project
CSV_PATH = os.path.join(PROJECT_ROOT, "inputs_finalproject.csv")
RAW_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "raw_pdb")
CLEAN_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "clean_pdb")
FASTA_OUTPUT = os.path.join(PROJECT_ROOT, "data", "sequences.fasta")

# Creating folders if non-existent (safe to run multiple times)
os.makedirs(RAW_PDB_DIR, exist_ok=True)
os.makedirs(CLEAN_PDB_DIR, exist_ok=True)

class CleanStructure(Select):
    """Creating a class for cleaning the PDB file: keeping only chain A, removing water, ligands, and heteroatoms."""
    def accept_chain(self, chain):
        return chain.id == 'A'

    def accept_residue(self, residue):
        # Removing water molecules and heteroatoms (id[0] is not empty for HETATM)
        if residue.id[0] != " ":
            return 0
        return 1
    

def run_step_1():
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: Could not find the input file at the absolute path:\n{CSV_PATH}")
        return

    # Reading the PDB-id list, inc. chain information from the CSV file, ensuring to strip whitespace and convert to string
    df = pd.read_csv(CSV_PATH, header=0)  # Assuming the first row is a header
    raw_identifiers = df.iloc[:, 0].astype(str).str.strip().tolist()
    
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    
    fasta_records = []
    
    print(f"Absolute root directory of the project: {PROJECT_ROOT}")
    print(f"Found {len(raw_identifiers)} proteins in the list. Starting download...\n")
    
    for identifier in raw_identifiers:
        identifier = identifier.lower()
        
        # Separating "1tup_c" into pdb_id="1tup" and target_chain="c"
        if "_" in identifier:
            pdb_id, target_chain = identifier.split("_")
        else:
            pdb_id = identifier
            target_chain = "a" # Standard if chain is not specified
            
        target_chain = target_chain.upper() # BioPython requires uppercase letters for chains
        print(f"Processing identifier: {identifier} -> PDB: {pdb_id}, Chain: {target_chain}")
        
        try:
            # Downloading the PDB file with the clean 4-character ID
            file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=RAW_PDB_DIR, file_format="pdb")
            
            # Parsing the downloaded PDB file to get the structure
            structure = parser.get_structure(pdb_id, file_path)
            
            # Creating a custom filter class for this specific chain
            class CustomClean(Select):
                def accept_chain(self, chain):
                    return chain.id == target_chain
                def accept_residue(self, residue):
                    if residue.id[0] != " ":
                        return 0
                    return 1
            
            # Saving the cleaned file (named after the full identifier to keep them separated)
            clean_file_path = os.path.join(CLEAN_PDB_DIR, f"{identifier}_clean.pdb")
            io.set_structure(structure)
            io.save(clean_file_path, CustomClean())
            
            # Extracting the amino acid sequence from the chain specified in the identifier
            clean_struct = parser.get_structure(identifier, clean_file_path)
            sequence = ""
            for model in clean_struct:
                for chain in model:
                    for residue in chain:
                        sequence += seq1(residue.get_resname())
            
            if sequence:
                fasta_records.append(f">{identifier}\n{sequence}\n")
                print(f"-> Successfully saved {identifier}_clean.pdb ({len(sequence)} aa)")
            else:
                print(f"!! WARNING: Found no amino acids in chain {target_chain} for {pdb_id}")

        except Exception as e:
            print(f"!! ERROR with {identifier}: {e}. Proceeding to next identifier...")

    # 5. Writing the FASTA file to its absolute path in the data folder
    with open(FASTA_OUTPUT, "w") as fasta_file:
        fasta_file.writelines(fasta_records)
        
    print(f"\nStep 1 completed! FASTA saved to:\n{FASTA_OUTPUT}")


if __name__ == "__main__":
    run_step_1() 