import Bio.PDB

# 1. Initiating parser
parser = Bio.PDB.PDBParser(QUIET=True)
structure = parser.get_structure("1ssc", "../../data/1ssc.pdb")

# 2. Extracting sequence from Chain A
# Using PPBuilder (Polypeptide Builder) to build polypeptides and extract sequences from the structure.
ppb = Bio.PDB.PPBuilder()
sequence_a = ""

for pp in ppb.build_peptides(structure[0]["A"]):
    sequence_a += pp.get_sequence()
print(f"AA sequence for Chain A:\n{sequence_a}")