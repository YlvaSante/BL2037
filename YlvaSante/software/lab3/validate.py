import argparse
import os
import csv
from Bio import ExPASy, SwissProt
from Bio.PDB import PDBParser, PPBuilder

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", required=True, help="UniProt ID")
    parser.add_argument("--pdb_file", required=True, help="Path to the .pdb file")
    parser.add_argument("--match_csv", required=True, help="Path to MATCH.csv")
    parser.add_argument("--unmatch_csv", required=True, help="Path to UNMATCHED.csv")
    args = parser.parse_args()

# 1. Fetching the PDB sequence (same as in Task 2)
parser = PDBParser(QUIET=True) # QUIET=True to suppress warnings about missing atoms
structure = parser.get_structure(args.id, args.pdb_file)
ppb = PPBuilder()
pdb_seq = ""
for pp in ppb.build_peptides(structure[0]["A"]):
    pdb_seq += str(pp.get_sequence())

# 2. Fetching the official sequence from UniProt by using id from the command line
print("Fetching data from UniProt...")
try:
    with ExPASy.get_sprot_raw(args.id) as handle:
        record = SwissProt.read(handle)
        uniprot_seq = record.sequence
except Exception as e:
    print(f"Error fetching UniProt data {args.id}: {e}")
    return

# 3.Comparing sequences and saving results in MATCH.csv and UNMATCHED.csv, instead of printing the sequences to the terminal window.
print(f"\nPDB Seq:     {pdb_seq[:50]}...")
print(f"UniProt Seq: {uniprot_seq[:50]}...")

if pdb_seq == uniprot_seq:
    print("\nOUTCOME: Match")
    # output_file = args.match_csv
    # print(f"{args.id}: MATCH")
else:
    print("\nOUTCOME: Mismatch")
    print(f"Length PDB: {len(pdb_seq)}, Length UniProt: {len(uniprot_seq)}")
    # output_file = args.unmatch_csv
    # print(f"{args.id}: UNMATCHED")

    # Saving into .csv file (append-mode 'a' in order not to overwrite existing data (previous results)
    file_exists = os.path.isfile(output_file) # Checking whether the file already exists to know if needing to write a header
    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f)
        # Writing header for files only
        if not file_exists:
            writer.writerow(["UniProt_ID", "PDB_Len", "UniProt_Len"])
        writer.writerow(result)

if __name__ == "__main__":
    main()