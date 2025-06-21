from Bio import AlignIO

alignment = AlignIO.read("codon_aligned.fasta", "fasta")

seq1 = alignment[0].seq
seq2 = alignment[1].seq

synonymous = 0
nonsynonymous = 0

codon_table = { ... }  # Standard genetic code

for i in range(0, len(seq1), 3):
    codon1 = str(seq1[i:i+3])
    codon2 = str(seq2[i:i+3])
    if '-' in codon1 or '-' in codon2:
        continue
    if codon1 != codon2:
        aa1 = codon_table.get(codon1.upper(), 'X')
        aa2 = codon_table.get(codon2.upper(), 'X')
        if aa1 == aa2:
            synonymous += 1
        else:
            nonsynonymous += 1

print(f"Synonymous differences: {synonymous}")
print(f"Nonsynonymous differences: {nonsynonymous}")
