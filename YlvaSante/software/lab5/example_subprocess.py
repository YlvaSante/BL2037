import subprocess
import re
# --- Configuration ---
PDB_FILE_1 = "/Users/david/Programming_for_Biologists/module_5/7VDE.pdb"
PDB_FILE_2 = "/Users/david/Programming_for_Biologists/module_5/3GOU.pdb"
OUTPUT_PREFIX = "/Users/david/Programming_for_Biologists/module_5/7VDE_3GOU_superposed"

def parse_tmalign_output(output_text):
    """Parses the text output from TM-align to find key values."""
    results = {}
    tm_score_match = re.search(r"TM-score\s*=\s*(\d\.\d+)", output_text)
    if tm_score_match:
        results['tm_score'] = float(tm_score_match.group(1))
    rmsd_match = re.search(r"RMSD\s*=\s*(\d\.\d+)", output_text)
    if rmsd_match:
        results['rmsd'] = float(rmsd_match.group(1))
    return results

def main():
    """Runs two TM-align commands using subprocess."""
    print("--- Running Command 1: Standard Alignment ---")
    command1 = ["TMalign", PDB_FILE_1, PDB_FILE_2]

    result1 = subprocess.run(
        command1,
        capture_output=True,
        text=True,
        check=True
    )

    parsed_results = parse_tmalign_output(result1.stdout)
    print("✅ Alignment successful. Parsed results:")
    print(f"   - TM-score: {parsed_results.get('tm_score', 'Not Found')}")
    print(f"   - RMSD: {parsed_results.get('rmsd', 'Not Found')} Å")

    print("\n--- Running Command 2: Generate Superposed File ---")
    command2 = ["TMalign", PDB_FILE_1, PDB_FILE_2, "-o", OUTPUT_PREFIX]

    subprocess.run(command2, check=True)
    print(f"✅ Command successful. Superposed files created with prefix '{OUTPUT_PREFIX}'.")

if __name__ == "__main__":
    main()