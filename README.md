Universal Plasmid Maker
BBL-434 Assignment 1

This repository contains a Python tool that generates a plasmid DNA sequence from a given DNA input and a design file provided by the user. The tool constructs a plasmid by combining a replication module, multiple cloning sites, selected antibiotic markers, and the inserted sequence, then outputs the final plasmid sequence in FASTA format.

--------------------------------------------------

Problem Overview

You will be given the sequence of an unknown organism as Input (for example, E. coli). You must identify the origin of replication (ORI) within the organism using the tools and methods taught in class. Then, using this newly identified ORI and a plasmid design file provided by the user, you must output a plasmid that works in the unknown organism.

A list of valid markers that can be incorporated in the design file is provided (markers.tab). It is acceptable to assume certain markers in your solution, but the tool must handle the non-existence of a marker in its dictionary without crashing.

The expected solution should also include a testing file.

For test cases, use the plasmid sequence pUC19.fa as input. For the design file, use Design_pUC19.txt, which deletes the EcoRI site present in the original sequence. The output should not contain this EcoRI site.

---

Repository Structure

BBL-434-ASSIGNMENT1/
|
|-- data/
|   |-- Input.fa            Input DNA sequence in FASTA format
|   |-- Design.txt          Plasmid design specification file
|   |-- markers.tab         List of restriction sites and markers
|
|-- tests/
|   |-- pUC19.fa            Test plasmid input sequence (pUC19)
|   |-- Design_pUC19.txt    Test design file for pUC19
|
|-- plasmid_maker.py        Main Python script to construct plasmid
|-- README.md               Project documentation

---

Requirements

- Python 3
- Biopython library

Biopython can be installed using:

pip install biopython

If using a conda environment:

conda activate assignment1

---

Usage

From the repository root, run the plasmid generator with:

python plasmid_maker.py data/Input.fa data/Design.txt

This command reads the input DNA and design, constructs the plasmid sequence, and writes the result to Output.fa in FASTA format.

---

Test Case

To run the test case using the provided pUC19 sequence:

python plasmid_maker.py tests/pUC19.fa tests/Design_pUC19.txt

In this test case, the design file instructs the removal of the EcoRI site. After execution, the resulting Output.fa should not contain the EcoRI recognition sequence (GAATTC).

---

Input File Format

Input.fa:
A standard FASTA file containing the DNA sequence for the unknown organism.

Design.txt:
A plain text file listing features of the plasmid to be included (multiple cloning sites and antibiotic markers). Each line contains a key and value separated by a comma.

Example:

Multiple_Cloning_Site1, EcoRI
Multiple_Cloning_Site2, BamHI
Antibiotic_marker1, Ampicillin
Antibiotic_marker2, Kanamycin

---

Markers File

markers.tab contains marker definitions for restriction enzymes and antibiotic resistance genes. The script uses this file to map design names to actual sequence motifs.

---

Output

The output file Output.fa contains the complete plasmid sequence in FASTA format with a descriptive header:

>Designed_Plasmid Auto-generated plasmid with backbone, MCS, insert, and markers

Followed by the constructed DNA sequence.

---

Notes

- The replication module and markers are defined or loaded as needed.
- The tool safely skips unknown markers specified in the design without crashing.
- EcoRI deletion is handled correctly for the specified test case.

---

Author

Sanskruti K  
BBL-434 Bioinformatics 
