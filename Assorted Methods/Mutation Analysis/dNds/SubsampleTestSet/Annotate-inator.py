import argparse
from pathlib import Path

def convert_to_gff(input_tsv: str, output_gff: str):
    """
    Converts CAsub.tsv (custom BLAST-like format) to GFF3.
    Column order: PUL_name, subject_id, pident, length, ?, ?, start, end, evalue, bitscore
    """
    with open(input_tsv, "r") as f_in, open(output_gff, "w") as f_out:
        f_out.write("##gff-version 3\n")
        
        for line_num, line in enumerate(f_in, 1):
            if line.strip().startswith("#") or not line.strip():
                continue
            
            fields = line.strip().split("\t")
            try:
                # Extract columns (adjust indices if needed)
                pul_name = fields[0]                # e.g., "PUL1"
                subject_id = fields[1]             # e.g., "CP061164.1:1663228-1681566"
                pident = fields[2]                 # e.g., "100.000"
                start = int(fields[6])             # BLAST hit start on subject
                end = int(fields[7])               # BLAST hit end on subject
                evalue = fields[8]                 # e.g., "0.0"
                bitscore = fields[9]               # e.g., "2405"

                # Clean subject_id (remove ":start-end" if present)
                seqid = subject_id.split(":")[0]   # e.g., "CP061164.1"

                # Determine strand
                strand = "+" if start < end else "-"
                start, end = sorted((start, end))  # Ensure start < end

                # Write GFF3 line
                f_out.write(
                    f"{seqid}\tBLAST\tmatch\t{start}\t{end}\t{bitscore}\t{strand}\t.\t"
                    f"ID={pul_name}_hit_{line_num};Name={pul_name};Identity={pident};Evalue={evalue}\n"
                )
            except (IndexError, ValueError) as e:
                print(f"Error in line {line_num}: {e}. Skipping: {line.strip()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CAsub.tsv to GFF3.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file (e.g., CAsub.tsv)")
    parser.add_argument("-o", "--output", help="Output GFF3 file (default: <input>.gff)")
    args = parser.parse_args()

    output_gff = args.output if args.output else f"{Path(args.input).stem}.gff"
    convert_to_gff(args.input, output_gff)
    print(f"GFF3 saved to: {output_gff}")
