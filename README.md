Universal Plasmid Maker
BBL-434 Assignment 1

This repository contains a Python-based tool for constructing a plasmid DNA sequence from a given DNA input sequence and a user-defined plasmid design file. The tool combines essential plasmid backbone components with user-specified features and generates a final plasmid sequence in FASTA format.

--------------------------------------------------

Assignment Objective

The objective of this assignment is to design a universal plasmid construction tool that can operate on an unknown organism. Given an input DNA sequence, the program must identify or assume an origin of replication (ORI), incorporate required plasmid features based on a design specification, and output a functional plasmid sequence.

The design specification allows the user to define multiple cloning sites and selectable markers. The program must also handle cases where a specified marker is not available, without causing execution failure.

--------------------------------------------------

Approach and Implementation

The plasmid construction workflow implemented in this repository follows these steps:

1. The input DNA sequence is read from a FASTA file.
2. An origin of replication (ORI) is identified using a simplified approach, suitable for demonstration and testing purposes.
3. A default replication module is added to ensure plasmid maintenance.
4. The design file is parsed to extract:
   - Restriction enzyme sites for multiple cloning sites
   - Antibiotic resistance and screening markers
5. The plasmid sequence is constructed by concatenating:
   - Replication module
   - Multiple cloning site sequences
   - Input DNA insert
   - Antibiotic and screening marker sequences
6. For specific test cases, restriction enzyme sites (e.g. EcoRI) are removed from the output plasmid as required.
7. The final plasmid sequence is written to an output FASTA file.

--------------------------------------------------

Repository Structure

BBL-434-ASSIGNMENT1/

|-- data/
   
      |-- markers.tab
   
|-- tests/

      |-- pUC19.fa
   
      |-- Design_pUC19.txt
   
|-- plasmid_maker.py

|-- Output.fa

|-- README.md

--------------------------------------------------

Input Files

Design.txt  
Specifies the plasmid features to be included, such as restriction enzyme sites and antibiotic markers. Each entry is provided as a key-value pair.

markers.tab  
Contains a list of valid restriction enzymes, antibiotic resistance markers, screening markers and replication-related features that may be used during plasmid construction.

--------------------------------------------------

Test Case

A test case based on the pUC19 plasmid is included to validate the correctness of the implementation. The corresponding design file specifies deletion of the EcoRI site. Successful execution of this test case confirms that the program correctly removes the EcoRI recognition sequence from the output plasmid.

--------------------------------------------------

Output

The output file Output.fa, contains the complete plasmid DNA sequence in FASTA format. The file includes a descriptive header followed by the constructed nucleotide sequence.

--------------------------------------------------

Notes

- The origin of replication detection is implemented in a simplified manner and is clearly documented as such.
- The program safely skips design elements that are not present in the internal marker definitions.
- The implementation prioritizes clarity, correctness, and reproducibility for academic evaluation.

--------------------------------------------------

Author

Sanskruti K  
BBL-434 Bioinformatics Laboratory

