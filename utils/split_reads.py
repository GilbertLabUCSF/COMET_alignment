from Bio import SeqIO
import pysam
import gzip
from tqdm import tqdm

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def split_reads_at_constant2(reads_fastq, constant2_mapping_sam, n_term_out, c_term_out, bad_reads_out, mapq_threshold=20):
    """
    Splits reads based on Constant2 mapping.
    
    Args:
        reads_fastq (str): Input FASTQ file with all reads.
        constant2_mapping_sam (str): SAM file with reads mapped to Constant2.
        n_term_out (str): Output FASTQ file for N-terminal sequences.
        c_term_out (str): Output FASTQ file for C-terminal sequences.
        bad_reads_out (str): Output FASTQ file for reads with no good Constant2 alignment.
        mapq_threshold (int): Minimum mapping quality to consider a read as aligned.
    """
    print(f"Parsing Constant2 mapping file: {constant2_mapping_sam}")
    # First, parse SAM to find valid alignments
    samfile = pysam.AlignmentFile(constant2_mapping_sam, "r")

    # Map from read_id -> (alignment_start, alignment_end)
    good_alignments = {}

    no_dCas9 = []
    low_quality = []
    
    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            no_dCas9.append(read.query_name)
            continue
        
        if read.mapping_quality < mapq_threshold:
            low_quality.append((read.query_name, read.mapping_quality))
            continue

        
        read_id = read.query_name
        # is_reverse = read.is_reverse

        aln_start = read.reference_start   # 0-based start on reference (Constant2)
        aln_end = read.reference_end       # 0-based end on reference
        
        # if is_reverse and aln_start > aln_end:
        #     aln_start, aln_end = aln_end, aln_start

        # Save alignment positions
        good_alignments[read_id] = (read.query_alignment_start, read.query_alignment_end)

    samfile.close()

    # Save lists
    with open("no_dCas9_reads.txt", "w") as f:
        for read_id in no_dCas9:
            f.write(f"{read_id}\n")

    with open("low_quality_reads.txt", "w") as f:
        for read_id, quality in low_quality:
            f.write(f"{read_id} {quality}\n")


    print(f"Found {len(good_alignments)} good Constant2 alignments.")

    # Now, process reads
    n_term_records = []
    c_term_records = []
    bad_records = []

    with open_maybe_gzip(reads_fastq, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq")):
            read_id = record.id
            
            if read_id in good_alignments:
                aln_start, aln_end = good_alignments[read_id]

                print(f"[DEBUG] Splitting read {read_id}: query start={aln_start}, query end={aln_end}, total len={len(record.seq)}")

                n_seq = record.seq[:aln_start]
                c_seq = record.seq[aln_end:]

                if len(n_seq) > 0:
                    # Make a copy of the record
                    n_record = record[:]

                    # Save original annotations temporarily
                    letter_annots = n_record.letter_annotations.copy()

                    # Clear annotations
                    n_record.letter_annotations = {}

                    # Modify sequence
                    n_record.seq = n_seq

                    # Rebuild annotations
                    # Slice the quality scores or any per-letter annotations
                    for key, value in letter_annots.items():
                        n_record.letter_annotations[key] = value[:aln_start]

                    n_term_records.append(n_record)

                    # n_record = record[:]
                    # n_record.seq = n_seq
                    # # Slice the phred quality scores too
                    # n_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][:aln_start]
                    # n_term_records.append(n_record)


                if len(c_seq) > 0:
                    c_record = record[:]
                    letter_annots = c_record.letter_annotations.copy()
                    c_record.letter_annotations = {}
                    c_record.seq = c_seq
                    for key, value in letter_annots.items():
                        c_record.letter_annotations[key] = value[aln_end:]
                    c_term_records.append(c_record)

                    # c_record = record[:]
                    # c_record.seq = c_seq
                    # c_record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][aln_end:]
                    # c_term_records.append(c_record)

            else:
                bad_records.append(record)


    print(f"Writing {len(n_term_records)} N-term reads to {n_term_out}")
    print(f"Writing {len(c_term_records)} C-term reads to {c_term_out}")
    print(f"Writing {len(bad_records)} bad reads to {bad_reads_out}")

    # Write outputs
    SeqIO.write(n_term_records, n_term_out, "fastq")
    SeqIO.write(c_term_records, c_term_out, "fastq")
    SeqIO.write(bad_records, bad_reads_out, "fastq")
