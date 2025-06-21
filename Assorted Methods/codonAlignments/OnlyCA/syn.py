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
            aa1 = table.forward_table.get(codon1, '*')  # '*' is stop codon
            aa2 = table.forward_table.get(codon2, '*')
            if aa1 == aa2:
                synonymous += 1
            else:
                nonsynonymous += 1
    return synonymous, nonsynonymous

def pairwise_syn_nonsyn(fasta_file, output_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Seq1", "Seq2", "Synonymous_diff", "Nonsynonymous_diff"])
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                syn, nonsyn = codon_diffs(str(sequences[i].seq), str(sequences[j].seq), table)
                writer.writerow([sequences[i].id, sequences[j].id, syn, nonsyn])

if __name__ == "__main__":
    fasta = "codon_aligned.fasta"  # Change to your fasta file
    output = "pairwise_syn_nonsyn.csv"  # Output CSV filename
    pairwise_syn_nonsyn(fasta, output)
    print(f"Results written to {output}")
