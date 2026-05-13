import subprocess
import re
import csv
import os

# 1. Defining pathways (absolute address for the data)
DATA_DIR = "/home/ylva/BL2037/YlvaSante/data/lab5/"
OUTPUT_DIR = "/home/ylva/BL2037/YlvaSante/results/lab5/"

# Creating absolute paths to the files
reference_file = os.path.join(DATA_DIR, "8JIQ_R.pdb")

target_filenames = ["8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]
# We creating a list where each file has its full path
target_files = [os.path.join(DATA_DIR, f) for f in target_filenames]

# 2. Creating the function for parsing output
def parse_tmalign_output(output_text):
    # re.search searching for patterns in the text 
    # (\d\.\d+) fetching the digit (t.ex. 0.85) in a "capture group".
    tm_match = re.search(r"TM-score\s*=\s*(\d\.\d+)", output_text)
    rmsd_match = re.search(r"RMSD\s*=\s*(\d\.\d+)", output_text)
    
    # Returning values if found, else "N/A"
    tm_score = tm_match.group(1) if tm_match else "N/A"
    rmsd = rmsd_match.group(1) if rmsd_match else "N/A"
    
    return tm_score, rmsd

# 2. Defining files (Make sure the names match exactly with the ones in the folder)
# reference_file = "8JIQ_R.pdb"
# target_files = ["8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]

# List for collecting all results for CSV writing
all_results = []

print(f"{'PDB ID':<12} | {'TM-score':<10} | {'RMSD':<10}")
print("-" * 40)

# 4. Looping through all targets
for target in target_files:
    # Defining the command. "-o" creates a file named 'output.pdb' (optional)
    # Adding a prefix to the output in order not to overwrite existing files each time
    base_name = os.path.basename(target).split('.')[0]
    output_prefix = os.path.join(OUTPUT_DIR, f"align_{base_name}")
    command = ["TMalign", reference_file, target, "-o", output_prefix]
    
    # 5. Running the tool with subprocess.run
    # capture_output=True makes it possible to read the output in the variable 'result.stdout'
    result = subprocess.run(command, capture_output=True, text=True)
    
    # 6. Calling the parsing function to extract TM-score and RMSD
    tm_score, rmsd = parse_tmalign_output(result.stdout)
    
    # Saving the result in our list
    all_results.append({
        "PDB ID": os.path.basename(target),
        "TM-score": tm_score,
        "RMSD": rmsd
    })
    
    # Printing a neat table in the terminal while we run
    print(f"{os.path.basename(target):<12} | {tm_score:<10} | {rmsd:<10}")

# 7. Saving to a CSV file
csv_path = os.path.join(OUTPUT_DIR, 'tmalign_summary.csv')
with open(csv_path, mode='w', newline='') as csv_file:
    fieldnames = ['PDB ID', 'TM-score', 'RMSD']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in all_results:
        writer.writerow(row)

print("\nDone! Results have been saved in 'tmalign_summary.csv'.")