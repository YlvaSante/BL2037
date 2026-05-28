import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer

# --- DEFINING ABSOLUTE PATHWAYS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

CLEAN_PDB_DIR = os.path.join(PROJECT_ROOT, "data", "clean_pdb")
AF_CLEAN_DIR = os.path.join(PROJECT_ROOT, "data", "alphafold_clean")
PLOT_DIR = os.path.join(PROJECT_ROOT, "results", "plots", "local_rmsd")

os.makedirs(PLOT_DIR, exist_ok=True)

def calculate_local_rmsd(pdb_path, af_path, identifier):
    """Superimposing AlphaFold structure on the PDB structure and calculating RMSD per residue."""
    parser = PDBParser(QUIET=True)
    
    try:
        ref_struct = parser.get_structure("PDB", pdb_path)
        sample_struct = parser.get_structure("AF", af_path)
    except Exception as e:
        print(f"!! Error parsing {identifier}: {e}")
        return None
    
    # Extracting all CA (Carbon-Alpha) atoms from both structures
    ref_atoms = [atom for atom in ref_struct.get_atoms() if atom.get_name() == 'CA']
    sample_atoms = [atom for atom in sample_struct.get_atoms() if atom.get_name() == 'CA']
    
    # Safeguard alignment based on index sequence to account for differing database numbering systems
    min_len = min(len(ref_atoms), len(sample_atoms))
    
    if min_len < 10: 
        print(f"!! Too few common residues for {identifier} (found only {min_len}).")
        return None
        
    # Truncating atom arrays to matching operational sizes in structural sequential order
    ref_atoms_aligned = ref_atoms[:min_len]
    sample_atoms_aligned = sample_atoms[:min_len]
    
    # 3D Structural Superimposition
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms_aligned, sample_atoms_aligned)
    super_imposer.apply(sample_atoms_aligned)
    
    residue_positions = []
    distances = []
    
    # Computing exact Euclidean distance per matched residue index
    for idx in range(min_len):
        a1 = ref_atoms_aligned[idx]
        a2 = sample_atoms_aligned[idx]
        distance = a1 - a2 
        
        # Generating a relative index position for continuous x-axis graphing
        residue_positions.append(idx + 1)
        distances.append(distance)
        
    # --- PLOTTING ---
    plt.figure(figsize=(10, 4))
    plt.plot(residue_positions, distances, color='#d95f02', linewidth=1.5, label='Local RMSD')
    plt.axhline(y=2.0, color='gray', linestyle='--', alpha=0.7, label='Significant Deviation (2.0 Å)')
    
    plt.title(f"Local Structural Deviation: PDB vs AlphaFold ({identifier.upper()})")
    plt.xlabel("Residue Position (Sequence Index)")
    plt.ylabel("CA Distance (Ångströms)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_path = os.path.join(PLOT_DIR, f"{identifier.lower()}_local_rmsd.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    avg_rmsd = np.mean(distances)
    high_error_regions = [residue_positions[i] for i, d in enumerate(distances) if d > 4.0]
    
    return avg_rmsd, len(high_error_regions)

def run_step7_pipeline():
    print("\n" + "="*60)
    print("=== INITIATING STEP 7: LOCAL RMSD PER RESIDUE (Amino Acid Alpha Carbons) ===")
    print("="*60)
    
    # Searching for files matching the new format: e.g., 1B9V_A_clean.pdb
    pdb_files = sorted(glob.glob(os.path.join(CLEAN_PDB_DIR, "*_clean.pdb")))
    print(f"Found {len(pdb_files)} cleaned PDB files (experimental) to analyze and compare to the AlphaFold structure (bioinformatically predicted).")
    
    results_summary = []
    
    for pdb_path in pdb_files:
        filename = os.path.basename(pdb_path)
        # Extracting identifier (e.g., 1B9V_A)
        identifier = filename.replace("_clean.pdb", "")
        
        # Path to the corresponding AlphaFold file (e.g., 1B9V_A_af_clean.pdb)
        af_path = os.path.join(AF_CLEAN_DIR, f"{identifier}_af_clean.pdb")
        
        if os.path.exists(af_path):
            res = calculate_local_rmsd(pdb_path, af_path, identifier)
            if res:
                avg_rmsd, num_high_errors = res
                results_summary.append({
                    "Protein_Chain": identifier,
                    "Avg_Local_RMSD_Å": round(avg_rmsd, 2),
                    "High_Error_Residues_Count": num_high_errors
                })
        
    # Saving results and presenting the most interesting cases for the report
    if results_summary:
        df = pd.DataFrame(results_summary)
        summary_path = os.path.join(PROJECT_ROOT, "data", "step7_local_rmsd_summary.csv")
        df.to_csv(summary_path, index=False)
        print(f"\n-> Step 7 complete! Summary saved to: {summary_path}")
        print(f"-> All individual profile plots saved to: {PLOT_DIR}")
        
        # Presenting automatically the 5 most interesting cases for the report.
        print("\n" + "-"*50)
        print(" TOP 5 CANDIDATES FOR PYMOL VISUALIZATION ")
        print(" (These have the highest structural discrepancies to discuss) ")
        print("-"*50)
        top5 = df.sort_values(by="Avg_Local_RMSD_Å", ascending=False).head(5)
        print(top5.to_string(index=False))
        print("-"*50)
    else:
        print("!! Could not find any matching PDB and AlphaFold pairs. Please check your directories.")

if __name__ == "__main__":
    run_step7_pipeline()