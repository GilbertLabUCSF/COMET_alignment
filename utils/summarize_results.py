import pysam
from Bio import SeqIO
import csv
from tqdm import tqdm

def summarize_n_c_matches(n_term_sam, c_term_sam, summary_tsv_out, n_index, c_index, mapq_threshold=20):
    """
    Summarize best N and C matches for each read based on SAM alignments.
    
    Args:
        n_term_sam (str): SAM file of N-term mappings.
        c_term_sam (str): SAM file of C-term mappings.
        summary_tsv_out (str): Output TSV summarizing best N and C matches.
        mapq_threshold (int): Minimum MAPQ to consider a mapping valid.
    """

    # ref_n_lengths = {
    #     record.id: 0.8 * len(record.seq)
    #     for record in SeqIO.parse(n_index, "fasta")
    # }

    # ref_c_lengths = {
    #     record.id: 0.8 * len(record.seq)
    #     for record in SeqIO.parse(c_index, "fasta")
    # }

    


    print(f"Parsing N-term mappings from {n_term_sam}")
    n_best_hits = {}
    n_samfile = pysam.AlignmentFile(n_term_sam, "r")

    ref_n_lengths = {
        ref: length * 0.95
            for ref, length in zip(n_samfile.references, n_samfile.lengths)
        }


    for read in tqdm(n_samfile.fetch(until_eof=True)):
        if read.is_unmapped:
            continue
        if read.mapping_quality < mapq_threshold:
            continue


        ref_name = read.reference_name

        # if less than 80% of the read is aligning throw it away
        if read.query_alignment_end - read.query_alignment_start < ref_n_lengths[ref_name]:
            continue
        
        read_id = read.query_name
        mapq = read.mapping_quality

        # Keep only the best match (highest MAPQ)
        if read_id not in n_best_hits or mapq > n_best_hits[read_id][1]:
            n_best_hits[read_id] = (ref_name, mapq)

    n_samfile.close()

    print(f"Parsing C-term mappings from {c_term_sam}")
    c_best_hits = {}
    c_samfile = pysam.AlignmentFile(c_term_sam, "r")

    ref_c_lengths = {
        ref: length * 0.95
            for ref, length in zip(c_samfile.references, c_samfile.lengths)
        }

    for read in tqdm(c_samfile.fetch(until_eof=True)):
        if read.is_unmapped:
            continue
        if read.mapping_quality < mapq_threshold:
            continue

        read_id = read.query_name
        ref_name = read.reference_name
        mapq = read.mapping_quality

        if read.query_alignment_end - read.query_alignment_start < ref_c_lengths[ref_name]:
            continue

        # Keep only the best match (highest MAPQ)
        if read_id not in c_best_hits or mapq > c_best_hits[read_id][1]:
            c_best_hits[read_id] = (ref_name, mapq)
    c_samfile.close()

    # Find all unique reads
    all_read_ids = set(list(n_best_hits.keys()) + list(c_best_hits.keys()))

    print(f"Summarizing {len(all_read_ids)} reads...")

    # Write TSV output
    with open(summary_tsv_out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["read_id", "n_term_match", "n_term_mapq", "c_term_match", "c_term_mapq"])
        
        for read_id in sorted(all_read_ids):
            n_info = n_best_hits.get(read_id, ("None", 0))
            c_info = c_best_hits.get(read_id, ("None", 0))
            writer.writerow([read_id, n_info[0], n_info[1], c_info[0], c_info[1]])

    print(f"Summary written to {summary_tsv_out}")
