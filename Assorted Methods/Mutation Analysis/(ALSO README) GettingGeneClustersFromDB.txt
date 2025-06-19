blastn -query CAQuery.fasta -db CAblastDB \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  -out CA_hits.tsv

put this output into BEDinator.py

then put bed file from THAT into the following. ALONGISIDE all of the genomes IN ONE FASTA FILE

bedtools getfasta -fi <COMBINED GENOME FASTA> -bed <BED OUTPUT> -fo <FASTA CLUSTER OUTPUT>

be careful with the padding set in BEDinator, since if its out of bounds, it skips it entirely, MP had to be set to 0 for this exact reason.

after this you gotta split the big fasta full of extracted clusters. (MAKE SURE THE EXTENSION IS .FASTA)

