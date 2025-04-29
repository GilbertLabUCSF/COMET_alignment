from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

    n_term_records = []
    c_term_records = []
    bad_dCas9_records = []


    for read in tqdm(samfile.fetch(until_eof=True), desc="Processing dCas9 alignment sam file"):
        aln_start = read.query_alignment_start   # 0-based start on reference (Constant2)
        aln_end = read.query_alignment_end
        seq = read.query_sequence or ""
        record = SeqRecord(
            Seq(seq),
            id=read.query_name,
            description=f"converted from SAM"
        )
        record.letter_annotations["phred_quality"] = read.query_qualities or []

        # Add core SAM fields as annotations
        record.annotations = {
            "SAM_flag": read.flag,
            "is_reverse": read.is_reverse,
            "is_unmapped": read.is_unmapped,
            "contig": read.reference_name or "*",
            "position": read.reference_start + 1 if read.reference_start is not None else None,  # 1-based
            "MAPQ": read.mapping_quality,
            "CIGAR": read.cigarstring,
            "mate_contig": read.next_reference_name or "*",
            "mate_position": read.next_reference_start + 1 if read.next_reference_start != -1 else None,
            "template_length": read.template_length
        }

        # Add optional tags (e.g., NM:i:0)
        for tag, value in read.tags:
            record.annotations[f"SAM_tag_{tag}"] = value

        # if read does not map, is too short, or is bad quality it is a bad read
        if read.is_unmapped or aln_end - aln_start < 3000 or read.mapping_quality < mapq_threshold:
            bad_dCas9_records.append(record)
        else:
            n_seq = record.seq[:aln_start]
            c_seq = record.seq[aln_end:]

            if len(n_seq) > 0:
                # Make a copy of the record
                n_record = record[:]

                # Save original annotations temporarily
                letter_annots = n_record.letter_annotations.copy()

                # # Rebuild annotations
                n_record.letter_annotations = {}
                n_record.seq = n_seq

                # Slice the quality scores or any per-letter annotations
                n_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][:aln_start]

                # Modify sequence
                n_term_records.append(n_record)
            

            if len(c_seq) > 0:
                c_record = record[:]
                letter_annots = c_record.letter_annotations.copy()
                c_record.letter_annotations = {}
                c_record.seq = c_seq
                c_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][aln_end:]
                c_term_records.append(c_record)


        
        # read_id = read.query_name
        # is_reverse = read.is_reverse
        # # print("is reverse:",is_reverse)

        # aln_start = read.reference_start   # 0-based start on reference (Constant2)
        # aln_end = read.reference_end       # 0-based end on reference
        
        # # Save alignment positions
        # good_alignments[read_id] = (read.query_alignment_start, read.query_alignment_end)

                

                
    samfile.close()


    print(f"Writing {len(n_term_records)} N-term reads to {n_term_out}")
    print(f"Writing {len(c_term_records)} C-term reads to {c_term_out}")
    print(f"Writing {len(bad_dCas9_records)} bad reads to {bad_reads_out}")

    # Write outputs
    SeqIO.write(n_term_records, n_term_out, "fastq")
    SeqIO.write(c_term_records, c_term_out, "fastq")
    SeqIO.write(bad_dCas9_records, bad_reads_out, "fastq")
