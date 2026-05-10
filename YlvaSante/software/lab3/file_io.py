# Defining paths (provided standing in current directory software/lab3)
input_file = "../../data/raw_sequences.txt"
output_file = "../../data/clean_sequences.txt"

# Opening raw_sequences.txt in Read mode ('r')
with open(input_file, 'r') as infile:
    lines = infile.readlines()

# Creating/opening clean_sequences.txt in Write mode ('w')
with open(output_file, 'w') as outfile:
    for line in lines:
        sequence = line.strip() # Removing invisible line breaks
        if sequence: # Checking that the line is not empty
            outfile.write(f"Processed: {sequence}\n")

print("File has been processed and saved as clean_sequences.txt in directory /YlvaSante/data/")