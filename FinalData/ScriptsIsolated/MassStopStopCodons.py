from Bio import SeqIO
import tkinter as tk
from tkinter import filedialog
import os

def remove_stop_codons(seq):
    seq = seq.upper()
    new_seq = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in ['TAA', 'TAG', 'TGA']:
            new_seq += '---'  # blank codon
        else:
            new_seq += codon
    return new_seq

# Set up GUI
root = tk.Tk()
root.withdraw()  # Hide main window

# Ask for multiple input FASTA files
input_fastas = filedialog.askopenfilenames(
    title="Select FASTA File(s)",
    filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")]
)

if not input_fastas:
    print("No files selected. Exiting.")
    exit()

# Process each file
for input_fasta in input_fastas:
    base, ext = os.path.splitext(input_fasta)
    output_fasta = base + "_NS.fasta"

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            cleaned_seq = remove_stop_codons(str(record.seq))
            out_handle.write(f">{record.id}\n{cleaned_seq}\n")

    print(f"Processed: {os.path.basename(input_fasta)} → {os.path.basename(output_fasta)}")

print("All selected files have been processed.")
