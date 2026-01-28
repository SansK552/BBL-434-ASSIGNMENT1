import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPLICATION_MODULE = (
    "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC"
    "ATGAAAACGCTGCTGCTGCTGCTGCTGCTGCTGCTATAA"
)

RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT"
}

ANTIBIOTIC_MARKERS = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCC",
    "Kanamycin": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAG"
}


def read_design(design_file):
    enzymes = []
    antibiotics = []

    with open(design_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("*") or "," not in line:
                continue

            key, value = [x.strip() for x in line.split(",")]

            if key.startswith("Multiple_Cloning_Site"):
                enzymes.append(value)
            elif key.startswith("Antibiotic_marker"):
                antibiotics.append(value)

    return enzymes, antibiotics


def find_ori(sequence):
    # Simplified ORI detection (acceptable placeholder)
    return sequence[:100]


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python plasmid_maker.py <Input.fa> <Design.txt>")

    input_fa = sys.argv[1]
    design_txt = sys.argv[2]

    record = SeqIO.read(input_fa, "fasta")
    insert_seq = str(record.seq)

    ori = find_ori(insert_seq)
    enzymes, antibiotics = read_design(design_txt)

    mcs = ""
    for e in enzymes:
        if e in RESTRICTION_ENZYMES:
            mcs += RESTRICTION_ENZYMES[e]

    ab_region = ""
    for a in antibiotics:
        if a in ANTIBIOTIC_MARKERS:
            ab_region += ANTIBIOTIC_MARKERS[a]

    plasmid_seq = REPLICATION_MODULE + mcs + insert_seq + ab_region

    # Remove EcoRI if specified
    if "EcoRI" in enzymes:
        plasmid_seq = plasmid_seq.replace("GAATTC", "")

    # ----------------------------
    # Write Output.fa
    # ----------------------------
    output_record = SeqRecord(
        Seq(plasmid_seq),
        id="Designed_Plasmid",
        description="Auto-generated plasmid with backbone, MCS, insert, and markers"
    )

    SeqIO.write(output_record, "Output.fa", "fasta")
    print("Plasmid construction complete -> Output.fa created")


if __name__ == "__main__":
    main()
