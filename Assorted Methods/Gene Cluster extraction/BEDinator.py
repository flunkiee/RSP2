import pandas as pd

# Load BLAST hits
df = pd.read_csv("CA_hits.tsv", sep="\t", header=None, names=[
    "qseqid", "sseqid", "pident", "length", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
])

# Normalize coordinates
df["start"] = df[["sstart", "send"]].min(axis=1)
df["end"] = df[["sstart", "send"]].max(axis=1)

# Flanking region size in bp
FLANK = 5000

# Output BED entries
bed_entries = []

# Group by each contig/genome ID
for sseqid, group in df.groupby("sseqid"):
    min_start = max(group["start"].min() - FLANK, 0)
    max_end = group["end"].max() + FLANK
    bed_entries.append([sseqid, int(min_start), int(max_end), f"cluster_{sseqid}"])

# Save to BED file
bed_df = pd.DataFrame(bed_entries, columns=["chrom", "start", "end", "name"])
bed_df.to_csv("CA_clusters.bed", sep="\t", header=False, index=False)
