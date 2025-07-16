awk '
function revcomp(seq,    rev, i, c) {
    rev = ""
    for (i = length(seq); i > 0; i--) {
        c = substr(seq, i, 1)
        if (c == "A") rev = rev "T"
        else if (c == "T") rev = rev "A"
        else if (c == "G") rev = rev "C"
        else if (c == "C") rev = rev "G"
        else rev = rev c  # Preserve N or ambiguous bases
    }
    return rev
}

/^>/ {
    if (seq != "") {
        gsub(/[ \t\r\n]/, "", seq)
        if (header ~ /yvmC/ || header ~ /cypX/) {
            seq = revcomp(seq)
        }
        print header
        for (i = 1; i <= length(seq); i += 60)
            print substr(seq, i, 60)
    }
    header = $0
    seq = ""
    next
}

{
    gsub(/[ \t\r\n]/, "", $0)
    seq = seq $0
}

END {
    if (seq != "") {
        gsub(/[ \t\r\n]/, "", seq)
        if (header ~ /yvmC/ || header ~ /cypX/) {
            seq = revcomp(seq)
        }
        print header
        for (i = 1; i <= length(seq); i += 60)
            print substr(seq, i, 60)
    }
}
' combined.fasta > inv-combined.fasta

