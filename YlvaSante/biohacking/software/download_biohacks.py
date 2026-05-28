import os
import urllib.request
import json

# Defining the molecules for Biohacking Workshops
BIOHACKING_MOLECULES = {
    "caffeine": "caffeine",
    "theobromine": "theobromine",
    "anandamide": "anandamide",
    "phenethylamine": "phenethylamine", # PEA
    "serotonin": "serotonin",
    "melatonin": "melatonin",
    "dopamine": "dopamine",
    "noradrenaline": "norepinephrine",
    "adrenaline": "epinephrine",
    "gaba": "gamma-aminobutyric acid",
    "glutamate": "glutamic acid",
    "acetylcholine": "acetylcholine",
    "oxytocin": "oxytocin" 
}

# Creating an output directory for the workshop materials
OUTPUT_DIR = "biohacking_workshop_molecules"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=========================================================")
print("=== DOWNLOADING MOLECULES FOR BIOHACKING WORKSHOP ===")
print("=========================================================\n")

for name, search_term in BIOHACKING_MOLECULES.items():
    print(f"Fetching 3D structure for: {name.upper()}...")
    try:
        # Step 1: Get the PubChem Compound ID (CID) using the name
        name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{search_term.replace(' ', '%20')}/cids/JSON"
        response = urllib.request.urlopen(name_url)
        data = json.loads(response.read().decode())
        cid = data["IdentifierList"]["CID"][0]
        
        # Step 2: Download the 3D structure in SDF format (Much more reliable on PubChem)
        # We specify ?record_type=3d to get the actual spatial coordinates for PyMOL
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF/?record_type=3d"
        output_file = os.path.join(OUTPUT_DIR, f"{name}.sdf")
        
        urllib.request.urlretrieve(sdf_url, output_file)
        print(f"-> Successfully saved to {output_file}\n")
        
    except Exception as e:
        print(f"!! Failed to download {name}. Error: {e}\n")

print("=========================================================")
print(f"All downloads complete! Files are saved in: './{OUTPUT_DIR}'")
print("You can open these .sdf files directly in PyMOL just like PDBs!")
print("=========================================================")