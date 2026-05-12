# This script is picking a UniProt ID from the list.txt file, fetching the first available PDB ID and downloading the file from the RCSB PDB database, and saving it in the data/pdb_files/ folder. The script takes two command line arguments: the UniProt ID and the output directory for the PDB file.
import requests
import sys
import os

def get_pdb_from_uniprot(uniprot_id, output_dir):
    # 1. Sending inquiry to UniProt for connections to the PDB database via their REST API
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Could not find UniProt data for {uniprot_id}")
        return None

    # 2. Searching in the text file (i.e. response) for lines starting with 'DR   PDB;'
    pdb_id = None
    for line in response.text.split('\n'):
        if line.startswith('DR   PDB;'):
            # This is what the line looks like: DR   PDB; 1SSC; X-ray; 1.50 A.
            pdb_id = line.split(';')[1].strip()
            break # Picking the first found

    if not pdb_id:
        print(f"No PDB structure found for {uniprot_id}")
        return None

    # 3. Downloading the PDB file from RCSB PDB
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_response = requests.get(pdb_url)
    
    if pdb_response.status_code == 200:
        # Saving the file named after the UniProt ID as a .pdb file so that validate.py can find it
        file_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
        with open(file_path, 'w') as f:
            f.write(pdb_response.text)
        print(f"Downloaded {pdb_id} as {uniprot_id}.pdb")
        return file_path
    else:
        print(f"Couldn't download .pdb file {pdb_id}")
        return None

if __name__ == "__main__":
    uid = sys.argv[1]
    out = sys.argv[2]
    get_pdb_from_uniprot(uid, out)