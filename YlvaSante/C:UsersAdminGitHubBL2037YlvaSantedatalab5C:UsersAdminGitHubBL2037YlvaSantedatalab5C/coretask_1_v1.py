import subprocess
import re
import csv

# 1. Creating the function for parsing output
def parse_tmalign_output(output_text):
    # re.search searching for patterns in the text 
    # (\d\.\d+) fetching the digit (t.ex. 0.85) in a "capture group".
    tm_match = re.search(r"TM-score\s*=\s*(\d\.\d+)", output_text)
    rmsd_match = re.search(r"RMSD\s*=\s*(\d\.\d+)", output_text)
    
    # Returning values if fiund, else "N/A"
    tm_score = tm_match.group(1) if tm_match else "N/A"
    rmsd = rmsd_match.group(1) if rmsd_match else "N/A"
    
    return tm_score, rmsd

# 2. Defining files (Make sure the names match exactly with the ones in the folder)
reference_file = "8JIQ_R.pdb"
target_files = ["8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]

# List for collecting all results for CSV writing
all_results = []

print(f"{'PDB ID':<12} | {'TM-score':<10} | {'RMSD':<10}")
print("-" * 40)

# 3. Looping through all targets
for target in target_files:
    # Defining the command. "-o" creates a file named 'output.pdb' (optional)
    # Adding a prefix to the output in order not to overwrite existing files each time
    output_prefix = f"align_{target.split('.')[0]}"
    command = ["TMalign", reference_file, target, "-o", output_prefix]
    
    # 4. Running the tool with subprocess.run
    # capture_output=True makes it possible to read the output in the variable 'result.stdout'
    result = subprocess.run(command, capture_output=True, text=True)
    
    # 5. Calling the parsing function to extract TM-score and RMSD
    tm_score, rmsd = parse_tmalign_output(result.stdout)
    
    # Saving the result in our list
    all_results.append({
        "PDB ID": target,
        "TM-score": tm_score,
        "RMSD": rmsd
    })
    
    # Printing a neat table in the terminal while we run
    print(f"{target:<12} | {tm_score:<10} | {rmsd:<10}")

# 6. Saving to a CSV file
with open('tmalign_summary.csv', mode='w', newline='') as csv_file:
    fieldnames = ['PDB ID', 'TM-score', 'RMSD']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in all_results:
        writer.writerow(row)

print("\nDone! Results have been saved in 'tmalign_summary.csv'.")