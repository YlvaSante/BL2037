# A simple script to read a protein sequence from a PDB file with argparse
#

# Step 1: Tell Python which tools we need.
# We need PDBParser to understand the PDB file format.
# We need seq1 to convert three-letter amino acid codes (like 'ALA') to one-letter codes (like 'A').
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import argparse  # Import to input the PDB file from the command line

#Define inputs
# NEW: 1. Use argparse to read the PDB filename from the command line.
parser = argparse.ArgumentParser(description="Read a protein sequences from a PDB file.")
# 2. Define the positional argument for the PDB filename
parser.add_argument('pdb_filename', help="The path to the PDB file (e.g., 1ssc.pdb)")
#3. add parser to parse the arguments
args = parser.parse_args()


# Step 2: Create our "PDB reader" tool and parse arguments.
# This line creates a reader object that knows how to read PDB files.
reader = PDBParser()

# Step 3: Read the file define above
# We tell our reader to open and understand the file.
# The result is stored in a variable called 'structure'.
pdb_filename = args.pdb_filename # define the parsed argument pdb_filename
structure_id = pdb_filename.split('.')[0] # Use the filename prefix as the structure ID
structure = reader.get_structure(structure_id, pdb_filename) # Uses the command-line pdb filename input
print(f"Input pdbid: {structure}")

# Step 4: Loop through the file to find the amino acids.
# A PDB file is organized like nested boxes: Structure -> Model -> Chain -> Residue.
# We use loops to open each box and get to the amino acids (residues).
for model in structure:
    for chain in model:
        
        # Create an empty string to build our sequence
        sequence = "" 
        
        # Look at each amino acid (residue) in the chain
        for residue in chain:
            # PDB files often include water molecules ('HOH'). We want to ignore them.
            if residue.get_resname() != "HOH":
                # Convert the 3-letter code to a 1-letter code and add it to our sequence string
                sequence += seq1(residue.get_resname())
        
        # Step 5: Print the result for this chain.
        print(f"Found sequence for Chain {chain.id}:")
        print(sequence)
        print() # Add a blank line for readability
