from Bio import SeqIO
from Bio.Seq import Seq
import re

# Load all contigs from the genome(s)
genomes = SeqIO.to_dict(SeqIO.parse("all_genomes.fasta", "fasta"))

# Output FASTA
out_fasta = open("extracted_genes.fasta", "w")

# Read your GFF
with open("custom_annotations.gff") as gff_file:
    for line in gff_file:
        if line.startswith("#") or line.strip() == "":
            continue

        fields = line.strip().split("\t")
        full_seqid = fields[0]  # e.g., CP061164.1:1663228-1681566
        start = int(fields[3]) - 1  # convert to 0-based
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]

        # Parse contig and offset from full_seqid
        match = re.match(r"([^:]+):(\d+)-(\d+)", full_seqid)
        if not match:
            print(f"Skipping invalid line: {line}")
            continue

        contig = match.group(1)
        region_start = int(match.group(2))
        region_end = int(match.group(3))

        # Check if contig exists
        if contig not in genomes:
            print(f"Contig {contig} not found in FASTA")
            continue

        # Calculate actual coordinates on genome
        abs_start = region_start + start
        abs_end = region_start + end

        # Extract sequence
        full_seq = genomes[contig].seq
        if abs_end > len(full_seq):
            print(f"Skipping out-of-bounds match on {contig}: {abs_start}-{abs_end}")
            continue

        subseq = full_seq[abs_start:abs_end]
        if strand == "-":
            subseq = subseq.reverse_complement()

        # Get gene ID for header
        gene_id = re.search(r'ID=([^;]+)', attrs)
        gene_name = gene_id.group(1) if gene_id else "unknown"

        # Write to output
        out_fasta.write(f">{gene_name}|{contig}:{abs_start + 1}-{abs_end}({strand})\n{subseq}\n")

out_fasta.close()
print("Extraction complete! Results in 'extracted_genes.fasta'")
