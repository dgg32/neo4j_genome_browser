

# Introduction


This repository contains code and data for my article "Turn Neo4j into a Genome Browser".

1. The three Python scripts are for data processing

  
2. The embl folder contains the original embl files.


3. The json folder contains the formated json files ready for import into Neo4j.

  

# Prerequisite

Neo4j Desktop

  

# Run

1. Convert original embl file into Neo4j-ready json file:

```
python embl_2_json.py [embl_file] [output_json_file] 
```

2. Convert the exported json file into a fasta file:

```
python json_2_fasta.py [exported_json_file] [output_file]
```

3. Convert the exported json file into an embl file:

```
python json_2_embl.py [exported_json_file] [original_embl_file]  [output_embl_file]
```

## Authors

  

*  **Sixing Huang** - *Concept and Coding*

  

## License

  

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
