import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
from Bio.Data import CodonTable
import csv

def codon_diffs(seq1, seq2, table):
    synonymous = 0
    nonsynonymous = 0
    length = min(len(seq1), len(seq2))
    for i in range(0, length, 3):
        codon1 = seq1[i:i+3].upper()
        codon2 = seq2[i:i+3].upper()
        if len(codon1) < 3 or len(codon2) < 3:
            continue
        if '-' in codon1 or '-' in codon2:
            continue
        if codon1 != codon2:
            aa1 = table.forward_table.get(codon1, '*')
            aa2 = table.forward_table.get(codon2, '*')
            if aa1 == aa2:
                synonymous += 1
            else:
                nonsynonymous += 1
    return synonymous, nonsynonymous

def pairwise_syn_nonsyn(fasta_file, output_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    ids = [seq.id for seq in sequences]

    n = len(sequences)
    syn_matrix = [[0]*n for _ in range(n)]
    nonsyn_matrix = [[0]*n for _ in range(n)]
    omega_matrix = [[""]*n for _ in range(n)]

    total_syn = 0
    total_nonsyn = 0

    for i in range(n):
        for j in range(i+1, n):
            syn, nonsyn = codon_diffs(str(sequences[i].seq), str(sequences[j].seq), table)
            syn_matrix[i][j] = syn_matrix[j][i] = syn
            nonsyn_matrix[i][j] = nonsyn_matrix[j][i] = nonsyn
            total_syn += syn
            total_nonsyn += nonsyn

            # Calculate omega
            if syn == 0:
                omega = "Inf" if nonsyn > 0 else 0
            else:
                omega = round(nonsyn / syn, 3)
            omega_matrix[i][j] = omega_matrix[j][i] = omega


    # Write CSV output
    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Synonymous matrix
        writer.writerow(["Synonymous Differences"])
        writer.writerow([""] + ids)
        for i, row in enumerate(syn_matrix):
            writer.writerow([ids[i]] + row)

        writer.writerow([])

        # Nonsynonymous matrix
        writer.writerow(["Nonsynonymous Differences"])
        writer.writerow([""] + ids)
        for i, row in enumerate(nonsyn_matrix):
            writer.writerow([ids[i]] + row)

        writer.writerow([])

        # Omega matrix
        writer.writerow(["Omega (dN/dS)"])
        writer.writerow([""] + ids)
        for i, row in enumerate(omega_matrix):
            writer.writerow([ids[i]] + row)

        writer.writerow([])

def main():
    root = tk.Tk()
    root.withdraw()

    fasta_file = filedialog.askopenfilename(
        title="Select input FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )
    if not fasta_file:
        messagebox.showinfo("Info", "No FASTA file selected. Exiting.")
        return

    output_file = filedialog.asksaveasfilename(
        title="Select output CSV file",
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not output_file:
        messagebox.showinfo("Info", "No output file selected. Exiting.")
        return

    pairwise_syn_nonsyn(fasta_file, output_file)
    messagebox.showinfo("Done", f"Results written to {output_file}")

if __name__ == "__main__":
    main()
