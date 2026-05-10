from Bio import ExPASy, SwissProt
from Bio.PDB import PDBParser, PPBuilder

# 1. Fetching the PDB sequence (same as in Task 2)
parser = PDBParser(QUIET=True)
structure = parser.get_structure("1ssc", "../../data/1ssc.pdb")
ppb = PPBuilder()
pdb_seq = ""
for pp in ppb.build_peptides(structure[0]["A"]):
    pdb_seq += str(pp.get_sequence())

# 2. Fetching the official sequence from UniProt (P61823)
print("Fetching data from UniProt...")
with ExPASy.get_sprot_raw("P61823") as handle:
    record = SwissProt.read(handle)
    uniprot_seq = record.sequence

# 3.Comparing sequences
print(f"\nPDB Seq:     {pdb_seq[:50]}...")
print(f"UniProt Seq: {uniprot_seq[:50]}...")

if pdb_seq == uniprot_seq:
    print("\nOUTCOME: Match")
else:
    print("\nOUTCOME: Mismatch")
    print(f"Length PDB: {len(pdb_seq)}, Length UniProt: {len(uniprot_seq)}")