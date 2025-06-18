import csv

stored = {'helix_residues': []}
cmd.dss()
cmd.iterate('ss H', 'stored["helix_residues"].append((chain, resi))', space=globals())

# Write to CSV file
with open('helix_residues.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Chain', 'Residue'])  # Write header
    writer.writerows(stored['helix_residues'])  # Write data

print("Helix residues written to helix_residues.csv")