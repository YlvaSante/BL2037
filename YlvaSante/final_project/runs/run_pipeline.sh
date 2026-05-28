#!/bin/bash

# Run this script by executing:
# bash /home/ylva/BL2037/YlvaSante/final_project/runs/run_pipeline.sh

# Retrieving the absolute path to where this shell script is located (should be the "runs" directory)
RUNS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$RUNS_DIR" )"

echo "Starting pipeline from absolute environment..."
echo "Root directory: $PROJECT_ROOT"

# Running Step 1 (by pointing to the absolute path of the python file):
# Retrieving and cleaning sequences
python3 "$PROJECT_ROOT/software/step1_retrieve.py"

# Running Step 2 & 3: Sequence Alignment & Clustering
python3 "$PROJECT_ROOT/software/step2_3_sequence.py"

# Running Step 4 & 5: Structural Alignment (TM-align RMSD)
python3 "$PROJECT_ROOT/software/step4_5_structure.py"

# Running Step 6: Fetching and analyzing AlphaFold predictions
python3 "$PROJECT_ROOT/software/step6_alphafold.py"

# Running Step 7: Local RMSD analysis - PDB and AlphaFold Comparison
python3 "$PROJECT_ROOT/software/step7_local_rmsd.py"

# Running Step 8: Motif Detection and Binding Pocket Analysis - Convergent Evolution
python3 "$PROJECT_ROOT/software/step8_binding_pocket.py"