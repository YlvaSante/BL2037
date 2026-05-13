import subprocess
import re
import csv

# 1. Skapa funktionen för att parsa utdata
def parse_tmalign_output(output_text):
    # re.search letar efter mönstret i textmassan. 
    # (\d\.\d+) fångar siffran (t.ex. 0.85) i en "capture group".
    tm_match = re.search(r"TM-score\s*=\s*(\d\.\d+)", output_text)
    rmsd_match = re.search(r"RMSD\s*=\s*(\d\.\d+)", output_text)
    
    # Returnera värdena om de hittas, annars "N/A"
    tm_score = tm_match.group(1) if tm_match else "N/A"
    rmsd = rmsd_match.group(1) if rmsd_match else "N/A"
    
    return tm_score, rmsd

# 2. Definiera filer (Se till att namnen matchar exakt de du drog in i mappen)
reference_file = "8JIQ_R.pdb"
target_files = ["8E3Y_R.pdb", "6WI9_R.pdb", "6X18_R.pdb", "7VQX_R.pdb", "7YON_R.pdb"]

# Lista för att samla alla resultat inför CSV-skrivning
all_results = []

print(f"{'PDB ID':<12} | {'TM-score':<10} | {'RMSD':<10}")
print("-" * 40)

# 3. Loopa igenom alla targets
for target in target_files:
    # Definiera kommandot. "-o" skapar en fil som heter 'output.pdb' (valfritt)
    # Vi lägger till prefix på output så vi inte skriver över samma fil hela tiden
    output_prefix = f"align_{target.split('.')[0]}"
    command = ["TMalign", reference_file, target, "-o", output_prefix]
    
    # 4. Kör verktyget med subprocess.run
    # capture_output=True gör att vi kan läsa resultatet i variabeln 'result.stdout'
    result = subprocess.run(command, capture_output=True, text=True)
    
    # 5. Anropa parsningsfunktionen
    tm_score, rmsd = parse_tmalign_output(result.stdout)
    
    # Spara resultatet i vår lista
    all_results.append({
        "PDB ID": target,
        "TM-score": tm_score,
        "RMSD": rmsd
    })
    
    # Printa en snygg tabell i terminalen medan vi kör
    print(f"{target:<12} | {tm_score:<10} | {rmsd:<10}")

# 6. Spara till en CSV-fil
with open('tmalign_summary.csv', mode='w', newline='') as csv_file:
    fieldnames = ['PDB ID', 'TM-score', 'RMSD']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in all_results:
        writer.writerow(row)

print("\nKlar! Resultaten har sparats i 'tmalign_summary.csv'.")