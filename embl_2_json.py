#embl to json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import json
import sys

input_embl = sys.argv[1]
output_json = sys.argv[2]

#{"id":"0","type":"relationship","label":"KNOWS","properties":{"since":1993,"bffSince":"P5M1DT12H"},"start":{"id":"0","labels":["User"]},"end":{"id":"1","labels":["User"]}}

rid = 0

content = ""

for seq_record in SeqIO.parse(input_embl, "embl"):
    del seq_record.annotations['references']
    contig_node = {"type": "node", "id": f"{seq_record.id}", "labels": ["Contig"], "properties": {"id": str(seq_record.id), "name": seq_record.annotations["organism"], "organism": seq_record.annotations["organism"],  "annotations": seq_record.annotations
    #"seq": str(seq_record.seq)
    }}
    #print (json.dumps(contig_node))
    content += json.dumps(contig_node) + "\n"
    
    for index, f in enumerate(seq_record.features):
        
        #print (f)
        q = {k[0]: k[1] for k in list(f.qualifiers.items())}
        #print (f.location)
        #print (q)
        node = {"type": "node", "id": f"{seq_record.id}_{index}", "labels": [], "properties": {"type": f.type, "organism": seq_record.annotations["organism"], "location_start": f.location.start, "location_end": f.location.end, "location_strand": f.location.strand}}

        if "product" in q:
            node["properties"]["name"] = q["product"][0]
        else:
            node["properties"]["name"] = "no_product_name"

        if index == 0:
            node["labels"].append("Lead")
            node["properties"]["name"] = "source"
        else:
            node["labels"].append("Gene")

        for k in q:
            #if k != "translation":
            if len(q[k]) == 1:
                node["properties"][f'qualifiers_{k}'] = q[k][0]
            else:
                node["properties"][f"qualifiers_{k}"] = q[k]



        #print (json.dumps(node))
        content += json.dumps(node) + "\n"

        if index == 0:
            relation = {"id": f"{seq_record.id}_r_{rid}", "type": "relationship", "label": "OWNS", "start": {"id": f"{seq_record.id}", "labels": ["Contig"]}, "end": {"id": f"{seq_record.id}_{index}", "labels": ["Lead"]}}
        elif index == 1:
            relation = {"id": f"{seq_record.id}_r_{rid}", "type": "relationship", "label": "LEADS", "start": {"id": f"{seq_record.id}_{index-1}", "labels": ["Lead"]}, "end": {"id": f"{seq_record.id}_{index}", "labels": ["Gene"]}}
        else:
            relation = {"id": f"{seq_record.id}_r_{rid}", "type": "relationship", "label": "NEXT", "start": {"id": f"{seq_record.id}_{index-1}", "labels": ["Gene"]}, "end": {"id": f"{seq_record.id}_{index}", "labels": ["Gene"]}}
        
        rid += 1
        content += json.dumps(relation) + "\n"



with open(output_json, "w") as output:
    output.write(content)