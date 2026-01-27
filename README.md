Universal Plasmid Maker  
BBL-434 Assignment 1

This repository contains a Python-based tool to generate a plasmid DNA sequence from a given DNA insert and a plasmid design specification file. The script reads biological input files, adds essential plasmid backbone components and produces a final plasmid sequence in FASTA format.
The project demonstrates basic bioinformatics scripting, file parsing, and sequence handling using Biopython.

--------------------------------------------------

Problem Description

The task is to design a universal plasmid generator that takes:
1. A DNA sequence file (Input.fa)
2. A design file (Design.txt) specifying multiple cloning sites and antibiotic markers

Based on the design file and the concepts discussed in the reference paper, the plasmid must include genes required for replication by default. The final output should be a complete plasmid DNA sequence written in FASTA format.

--------------------------------------------------

Methodology

The script performs the following steps:

1. Reads the DNA insert sequence from a FASTA file.
2. Parses the design file to identify restriction enzymes and antibiotic resistance markers.
3. Adds a predefined replication module representing the plasmid backbone.
4. Constructs the plasmid by concatenating the replication module, multiple cloning site, insert DNA sequence, and antibiotic marker sequences.
5. Writes the final plasmid sequence to Output.fa in FASTA format.

--------------------------------------------------

Repository Structure

BBL-434-ASSIGNMENT1/

├── Plasmid_maker.py    Main Python script 

├── Input.fa            Input DNA sequence (FASTA format)

├── Design.txt          Plasmid design specification file

├── Output.fa           Generated plasmid sequence

└── README.md           Project documentation

--------------------------------------------------

Requirements

- Python 3
- Biopython
