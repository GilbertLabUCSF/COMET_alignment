import os
import argparse
import yaml
from utils.find_constant2 import map_reads_to_constant2
from utils.split_reads import split_reads_at_constant2
from utils.split_read_parallel import split_reads_at_constant2_parallel
from utils.map_n_c_terms import map_n_c_terms
from utils.summarize_results import summarize_n_c_matches

def process_dataset(dataset_name, dataset_info, output_dir, args):
    # Setup paths
    os.makedirs(output_dir, exist_ok=True)

    # Index paths
    constant2_index = os.path.join(output_dir, "constant2.mmi")
    n_index = os.path.join(output_dir, "n_candidates.mmi")
    c_index = os.path.join(output_dir, "c_candidates.mmi")

    # Intermediate files
    constant2_mapping_sam = os.path.join(output_dir, f"reads_vs_constant2_{dataset_name}.sam")
    n_term_fastq = os.path.join(output_dir, f"n_term_{dataset_name}.fastq")
    c_term_fastq = os.path.join(output_dir, f"c_term_{dataset_name}.fastq")
    bad_reads_fastq = os.path.join(output_dir, f"bad_reads_{dataset_name}.fastq")
    n_term_mapping_sam = os.path.join(output_dir, f"n_term_{dataset_name}.sam")
    c_term_mapping_sam = os.path.join(output_dir, f"c_term_{dataset_name}.sam")
    summary_tsv = os.path.join(output_dir, f"read_calls_{dataset_name}.tsv")
    counts_tsv = os.path.join(output_dir, f"counts_summary_{dataset_name}.tsv")

    # 1. Map reads to Constant2
    map_reads_to_constant2(
        dataset_info['CCS_fastq_fn'],
        dataset_info['dCas9'],
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
        dataset_info['n_candidates'],
        dataset_info['c_candidates'],
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
        counts_tsv
    )

    # Remove intermediate files after they are done being used
    if not args.save_intermediates:
        os.remove(constant2_mapping_sam)
        os.remove(n_term_fastq)
        os.remove(c_term_fastq)
        os.remove(bad_reads_fastq)
        os.remove(n_term_mapping_sam)
        os.remove(c_term_mapping_sam)
        os.remove(constant2_index)
        os.remove(n_index)
        os.remove(c_index)

    print("\n✅ Pipeline complete!")



def main(args):
    # Read the YAML file with datasets
    with open(args.yaml_file, 'r') as file:
        datasets = yaml.safe_load(file)

    # Process each dataset specified in the YAML file
    for dataset_name, dataset_info in datasets.items():
        output_dir = os.path.join(args.output_dir, dataset_name)
        process_dataset(dataset_name, dataset_info, output_dir, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline to process long-read sequencing data for N- and C-terminal identification.")

    parser.add_argument("--yaml_file", required=True, help="YAML configuration file with datasets.")
    parser.add_argument("--output_dir", default="output", help="Directory to save all outputs. (default: output)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for minimap2. (default: 8)")
    parser.add_argument("--parallel", action="store_true", help="Enable multithreading in split_reads") # i don't think this makes anything faster (might even make code slower)
    parser.add_argument("--save_intermediates", action="store_true", help="Enable to save the sam and fastq files that are created throughout the pipeline")

    args = parser.parse_args()
    main(args)

