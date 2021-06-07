#json to fasta

import json, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


input_file = sys.argv[1]
output_file = sys.argv[2]

fasta = []
for line in open(input_file, 'r'):
    sequence = json.loads(line.strip())
    header = []

    for h in ["f1.organism", "f1.name", "f1.qualifiers_locus_tag"]:
        if h in sequence:
            header.append(sequence[h])

    header_str = "|".join(header)

    if "f1.qualifiers_translation" in sequence and len(header_str) > 0:
        s = SeqRecord(Seq(sequence["f1.qualifiers_translation"]), id=header_str, annotations={"molecule_type": "Protein"}, description = "")
        fasta.append(s)

SeqIO.write(fasta, output_file, "fasta")
