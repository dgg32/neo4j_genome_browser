#json to embl

import json, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


new_input_json = sys.argv[1]

old_embl = sys.argv[2]

output_file = sys.argv[3]



dna_sequence = ""

for seq_record in SeqIO.parse(old_embl, "embl"):
    dna_sequence =  seq_record.seq

seq = {}

for line in open(new_input_json, "r"):
    data = json.loads(line.strip())


    node = data[list(data.keys())[0]]

    #print (node)
    if "Contig" in node["labels"]:
        seq["seq"] = dna_sequence
        seq["id"] = node["properties"]["id"]
        seq["annotations"] = {}
    
        for p in node["properties"]:
            if p.startswith("annotations."):
                seq["annotations"][p.replace("annotations.", "")] = node["properties"][p]

        seq["features"] = []
    elif "Lead" in node["labels"] or "Gene" in node["labels"]:
        qualifiers = {}
        for k in node["properties"]:
            if k.startswith("qualifiers"):
                qualifiers[k.replace("qualifiers.", "")] = node["properties"][k]


        seq["features"].append(SeqFeature(type=node["properties"]["type"], location=FeatureLocation(node["properties"]["location_start"], node["properties"]["location_end"], strand=node["properties"]["location_strand"]), qualifiers=qualifiers))


#print (seq["seq"])

seq_record = SeqRecord(seq["seq"], id=seq["id"], annotations=seq["annotations"], features=seq["features"])

SeqIO.write(seq_record, output_file, "embl")