import subprocess
import os
from Bio import SeqIO
import pysam

def map_reads_to_constant2(reads_fastq, constant2_fasta, constant2_index, output_sam, threads=8):
    """
    Build a minimap2 index for the Constant2 sequence if it doesn't exist.
    Map reads to the Constant2 reference.
    Save output to a SAM file.
    
    Args:
        reads_fastq (str): Path to input FASTQ reads.
        constant2_fasta (str): Path to Constant2 FASTA file.
        constant2_index (str): Path where Constant2 index (.mmi) will be saved.
        output_sam (str): Path to save the mapping results (SAM file).
        threads (int): Number of threads for minimap2.
    """
    # 1. Build the index if it does not exist
    if not os.path.exists(constant2_index):
        print(f"Building minimap2 index for Constant2: {constant2_fasta}")
        index_cmd = [
            "minimap2",
            "-d", constant2_index,
            constant2_fasta
        ]
        subprocess.run(index_cmd, check=True)
    else:
        print(f"Found existing Constant2 index: {constant2_index}")

    # 2. Map reads to Constant2
    print(f"Mapping reads to Constant2...")
    map_cmd = [
        "minimap2",
        "-x", "map-pb",             # preset for long noisy reads
        "-a",
        "-t", str(threads),           # number of threads
        "-I64G", # change these for different computers
        "-K8G",  # change these for different computers
        constant2_index,              # indexed Constant2 file
        reads_fastq                   # reads to map
    ]

    # Open the output SAM file and write
    with open(output_sam, "w") as out_fh:
        subprocess.run(map_cmd, check=True, stdout=out_fh)

    print(f"Mapping complete. Output saved to {output_sam}")

