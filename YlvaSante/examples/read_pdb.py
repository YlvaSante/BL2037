#
# A simple script to read a protein sequence from a PDB file
#

# Step 1: Tell Python which tools we need.
# We need PDBParser to understand the PDB file format.
# We need seq1 to convert three-letter amino acid codes (like 'ALA') to one-letter codes (like 'A').
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# Step 2: Create our "PDB reader" tool.
# This line creates a reader object that knows how to read PDB files.
reader = PDBParser()

# Step 3: Read the file.
# We tell our reader to open and understand the file "1ssc.pdb".
# The result is stored in a variable called 'structure'.
#structure = reader.get_structure("1ssc", "1ssc.pdb")
structure = reader.get_structure("9ekb", "9ekb.pdb")
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
                sequence += seq1(residue.get_resname())

        # Step 5: Print the result for this chain.
        print(f"Found sequence for Chain {chain.id}:")
        print(sequence)
        print() # Add a blank line for readability