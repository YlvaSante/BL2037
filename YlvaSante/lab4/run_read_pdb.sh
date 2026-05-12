#!/bin/bash

# define arguments here
data_dir="./data/PDBS" # data directory containing PDB input files
PDB_input="$data_dir/1SSC.pdb" # pattern to match all PDB input files in the data directory

echo "Procesing PDB input file path: $PDB_input" # print the PDB input file path that is currently being processed
python ./software/read_pdb_argparsed.py \
    $PDB_input 
