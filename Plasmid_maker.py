import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Restriction enzymes dictionary
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG"
}

# Antibiotic / screening markers
ANTIBIOTIC_MARKERS = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCC",
    "Kanamycin": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAG"
}


def read_design(design_file):
    """
    Reads the design file and extracts enzyme and marker names
    based on values rather than key prefixes.
    """
    enzymes = []
    antibiotics = []

    with open(design_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("*") or "," not in line:
                continue

            _, value = [x.strip() for x in line.split(",")]

            if value in RESTRICTION_ENZYMES:
                enzymes.append(value)
            elif value in ANTIBIOTIC_MARKERS:
                antibiotics.append(value)

    return enzymes, antibiotics


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python plasmid_maker.py <Input.fa> <Design.txt>")

    input_fa = sys.argv[1]
    design_txt = sys.argv[2]

    # Read input plasmid (pUC19 for test case)
    record = SeqIO.read(input_fa, "fasta")
    plasmid_seq = str(record.seq)

    # Read design instructions
    enzymes, antibiotics = read_design(design_txt)

    # Remove specified restriction enzyme sites
    for enzyme in enzymes:
        plasmid_seq = plasmid_seq.replace(RESTRICTION_ENZYMES[enzyme], "")

    # Write Output.fa
    output_record = SeqRecord(
        Seq(plasmid_seq),
        id=record.id,
        description="Modified plasmid sequence with specified restriction sites removed"
    )

    SeqIO.write(output_record, "Output.fa", "fasta")
    print("Plasmid modification complete -> Output.fa created")


if __name__ == "__main__":
    main()
