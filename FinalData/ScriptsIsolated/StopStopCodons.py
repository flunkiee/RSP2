###Entirely optional script, used for removal of stop codons.

from Bio import SeqIO
import tkinter as tk
from tkinter import filedialog

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

# Ask for input FASTA file
input_fasta = filedialog.askopenfilename(
    title="Select FASTA File",
    filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")]
)

if not input_fasta:
    print("No file selected. Exiting.")
    exit()

# Ask for output file name
output_fasta = filedialog.asksaveasfilename(
    title="Save Output FASTA As",
    defaultextension=".fasta",
    filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")]
)

if not output_fasta:
    print("No output file specified. Exiting.")
    exit()

# Process sequences
with open(output_fasta, "w") as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        cleaned_seq = remove_stop_codons(str(record.seq))
        out_handle.write(f">{record.id}\n{cleaned_seq}\n")

print(f"Done. Cleaned sequences saved to: {output_fasta}")
