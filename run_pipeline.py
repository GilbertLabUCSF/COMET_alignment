import os
import argparse

from utils.find_constant2 import map_reads_to_constant2
from utils.split_reads import split_reads_at_constant2
from utils.split_read_parallel import split_reads_at_constant2_parallel
from utils.map_n_c_terms import map_n_c_terms
from utils.summarize_results import summarize_n_c_matches

def main(args):
    # Setup paths
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Index paths
    constant2_index = os.path.join(output_dir, "constant2.mmi")
    n_index = os.path.join(output_dir, "n_candidates.mmi")
    c_index = os.path.join(output_dir, "c_candidates.mmi")

    # Intermediate files
    constant2_mapping_sam = os.path.join(output_dir, "reads_vs_constant2.sam")
    n_term_fastq = os.path.join(output_dir, "n_term.fastq")
    c_term_fastq = os.path.join(output_dir, "c_term.fastq")
    bad_reads_fastq = os.path.join(output_dir, "bad_reads.fastq")
    n_term_mapping_sam = os.path.join(output_dir, "n_term.sam")
    c_term_mapping_sam = os.path.join(output_dir, "c_term.sam")
    summary_tsv = os.path.join(output_dir, "n_c_summary.tsv")

    # 1. Map reads to Constant2
    map_reads_to_constant2(
        args.reads,
        args.constant2,
        constant2_index,
        constant2_mapping_sam,
        threads=args.threads
    )

    # 2. Split reads at Constant2
    if args.parallel:
        split_reads_at_constant2_parallel(
            constant2_mapping_sam,
            n_term_fastq,
            c_term_fastq,
            bad_reads_fastq
        )
    else:
        split_reads_at_constant2(
            constant2_mapping_sam,
            n_term_fastq,
            c_term_fastq,
            bad_reads_fastq
        )

    # 3. Map N-term and C-term sequences to candidates
    map_n_c_terms(
        n_term_fastq,
        c_term_fastq,
        args.n_candidates,
        args.c_candidates,
        n_index,
        c_index,
        n_term_mapping_sam,
        c_term_mapping_sam,
        threads=args.threads
    )

    # 4. Summarize best matches
    summarize_n_c_matches(
        n_term_mapping_sam,
        c_term_mapping_sam,
        summary_tsv,
        n_index,
        c_index
    )

    print("\n✅ Pipeline complete!")
    print(f"Results written to: {summary_tsv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline to process long-read sequencing data for N- and C-terminal identification.")

    parser.add_argument("--reads", required=True, help="Input FASTQ file with reads.")
    parser.add_argument("--constant2", required=True, help="FASTA file of Constant2 (e.g., dCas9 sequence).")
    parser.add_argument("--n_candidates", required=True, help="FASTA file with candidate N-terminal sequences.")
    parser.add_argument("--c_candidates", required=True, help="FASTA file with candidate C-terminal sequences.")
    parser.add_argument("--output_dir", default="output", help="Directory to save all outputs. (default: output)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for minimap2. (default: 8)")
    parser.add_argument("--parallel", action="store_true", help="Enable multithreading in split_reads") # i don't think this makes anything faster (might even make code slower)

    args = parser.parse_args()
    main(args)

