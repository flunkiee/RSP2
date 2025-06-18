import pandas as pd

# Read the CSV file
df = pd.read_csv('helix_residues.csv')  # or your file name

# Ensure 'Residue' is numeric (in case it's read as string)
df['Residue'] = pd.to_numeric(df['Residue'], errors='coerce')

# Group by 'Chain' (to handle multi-chain proteins)
df['Residue_diff'] = df.groupby('Chain')['Residue'].diff()  # Calculate difference between consecutive residues

# Detect jumps (where the difference > 1)
jumps = df[df['Residue_diff'] > 1]

print("Detected jumps in residue numbering:")
print(jumps)

# Save jumps to a new CSV (optional)
jumps.to_csv('residue_jumps.csv', index=False)