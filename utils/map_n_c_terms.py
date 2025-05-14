import subprocess
import os

def map_n_c_terms(
    n_term_fastq,
    c_term_fastq,
    n_candidates_fasta,
    c_candidates_fasta,
    n_index,
    c_index,
    n_term_sam_out,
    c_term_sam_out,
    threads=8
):
    """
    Map N-term and C-term sequences to their respective candidate libraries using minimap2.
    Save output as SAM files.
    
    Args:
        n_term_fastq (str): FASTQ file of N-terminal split sequences.
        c_term_fastq (str): FASTQ file of C-terminal split sequences.
        n_candidates_fasta (str): FASTA file of candidate N-terminal sequences.
        c_candidates_fasta (str): FASTA file of candidate C-terminal sequences.
        n_index (str): Path to minimap2 index for N candidates.
        c_index (str): Path to minimap2 index for C candidates.
        n_term_sam_out (str): Output SAM file for N-term mappings.
        c_term_sam_out (str): Output SAM file for C-term mappings.
        threads (int): Number of threads to use for minimap2.
    """

    # 1. Build N-candidates index if not exist
    if not os.path.exists(n_index):
        print(f"Building minimap2 index for N candidates: {n_candidates_fasta}")
        subprocess.run(
            ["minimap2", "-d", n_index, n_candidates_fasta],
            check=True
        )
    else:
        print(f"Found existing N-candidate index: {n_index}")

    # 2. Build C-candidates index if not exist
    if not os.path.exists(c_index):
        print(f"Building minimap2 index for C candidates: {c_candidates_fasta}")
        subprocess.run(
            ["minimap2", "-d", c_index, c_candidates_fasta],
            check=True
        )
    else:
        print(f"Found existing C-candidate index: {c_index}")

    # 3. Map N-term reads
    print(f"Mapping N-term reads to N candidates...")
    n_map_cmd = [
        "minimap2",
        "-x", "map-pb",
        "-a",
        "-t", str(threads),
        "-I64G", # change these for different computers
        "-K8G",  # change these for different computers
        "--end-bonus", "100",
        n_index,
        n_term_fastq
    ]
    with open(n_term_sam_out, "w") as out_fh:
        subprocess.run(n_map_cmd, check=True, stdout=out_fh)
    print(f"N-term mapping complete. Output: {n_term_sam_out}")

    # 4. Map C-term reads
    print(f"Mapping C-term reads to C candidates...")
    c_map_cmd = [
        "minimap2",
        "-x", "map-pb",
        "-a",
        "-t", str(threads),
        "-I64G", # change these for different computers
        "-K8G",  # change these for different computers
        "--end-bonus", "100",
        c_index,
        c_term_fastq
    ]

    with open(c_term_sam_out, "w") as out_fh:
        subprocess.run(c_map_cmd, check=True, stdout=out_fh)
    print(f"C-term mapping complete. Output: {c_term_sam_out}")
