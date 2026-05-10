import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. Loading data
# Exiting folder software (..) going into folder data
file_path = "../../data/gene_expression.csv"
df = pd.read_csv(file_path)

# 2. Creating histogram for column 'TPM' 
# (Note: Task is saying 'expression_level', though the file is named 'TPM')
plt.figure(figsize=(10, 6))
plt.hist(df['TPM'], bins=20, color='skyblue', edgecolor='black')

# 3. Adding titles and labels
plt.title("Gene Expression Distribution")
plt.xlabel("expression_level (TPM)")
plt.ylabel("Frequency")

# 4. Saving the plot in the results folder
# Exiting folder software (../) going into folder results
output_path = "../../results/hist.png"

# Checking if the folder exists, otherwise create it
os.makedirs("../../results", exist_ok=True)

# Saving with high resolution (600 dpi)
plt.savefig(output_path, dpi=600)
print(f"Graph has been saved in: {output_path}")

# Displayingpython that it is done
plt.show()