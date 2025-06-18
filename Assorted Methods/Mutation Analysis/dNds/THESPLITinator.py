import os
from Bio import SeqIO

def split_fasta_biopython(input_fasta, output_dir="split_fastas"):
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        count += 1
        filename = os.path.join(output_dir, f"{record.id}.fasta")
        print(f"Writing sequence {record.id} to {filename}")
        SeqIO.write(record, filename, "fasta")
    print(f"Done! Processed {count} sequences.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python split_fasta.py your_input.fasta")
        sys.exit(1)
    input_fasta = sys.argv[1]
    split_fasta_biopython(input_fasta)
