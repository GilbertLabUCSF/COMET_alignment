import multiprocessing
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import gzip
from tqdm import tqdm

CHUNK_SIZE = 10000  # Adjust depending on memory/cores

def chunked_iterator(iterable, size):
    iterator = iter(iterable)
    for first in iterator:
        yield [first] + list(islice(iterator, size - 1))

def process_read_chunk(reads_chunk, reference_length_cutoff, mapq_threshold):
    n_term_records = []
    c_term_records = []
    bad_dCas9_records = []

    for read in reads_chunk:
        aln_start = read['query_alignment_start']
        aln_end = read['query_alignment_end']
        seq = read['query_sequence'] or ""
        qual = read['query_qualities'] or []

        record = SeqRecord(
            Seq(seq),
            id=read['query_name'],
            description="converted from SAM"
        )
        record.letter_annotations["phred_quality"] = qual

        # Add core fields as annotations
        record.annotations = {
            "SAM_flag": read['flag'],
            "is_reverse": read['is_reverse'],
            "is_unmapped": read['is_unmapped'],
            "contig": read['reference_name'] or "*",
            "position": read['reference_start'] + 1 if read['reference_start'] is not None else None,
            "MAPQ": read['mapping_quality'],
            "CIGAR": read['cigarstring'],
            "mate_contig": read['next_reference_name'] or "*",
            "mate_position": read['next_reference_start'] + 1 if read['next_reference_start'] != -1 else None,
            "template_length": read['template_length']
        }

        for tag, value in read['tags']:
            record.annotations[f"SAM_tag_{tag}"] = value

        # Filter based on criteria
        if read['is_unmapped'] or read['reference_length'] < reference_length_cutoff or read['mapping_quality'] < mapq_threshold:
            bad_dCas9_records.append(record)
            continue

        # Extract N- and C-terminal fragments
        n_seq = record.seq[:aln_start]
        c_seq = record.seq[aln_end:]

        if len(n_seq) > 0:
            n_record = record[:]
            letter_annots = n_record.letter_annotations.copy()
            n_record.letter_annotations = {}
            n_record.seq = n_seq
            n_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][:aln_start]
            n_term_records.append(n_record)

        if len(c_seq) > 0:
            c_record = record[:]
            letter_annots = c_record.letter_annotations.copy()
            c_record.letter_annotations = {}
            c_record.seq = c_seq
            c_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][aln_end:]
            c_term_records.append(c_record)

    return n_term_records, c_term_records, bad_dCas9_records


# def process_read_chunk(reads_chunk, reference_length_cutoff, mapq_threshold):
#     n_term_records = []
#     c_term_records = []
#     bad_dCas9_records = []

#     for read in reads_chunk:
#         aln_start = read.query_alignment_start
#         aln_end = read.query_alignment_end
#         seq = read.query_sequence or ""
#         record = SeqRecord(
#             Seq(seq),
#             id=read.query_name,
#             description="converted from SAM"
#         )
#         record.letter_annotations["phred_quality"] = read.query_qualities or []

#         record.annotations = {
#             "SAM_flag": read.flag,
#             "is_reverse": read.is_reverse,
#             "is_unmapped": read.is_unmapped,
#             "contig": read.reference_name or "*",
#             "position": read.reference_start + 1 if read.reference_start is not None else None,
#             "MAPQ": read.mapping_quality,
#             "CIGAR": read.cigarstring,
#             "mate_contig": read.next_reference_name or "*",
#             "mate_position": read.next_reference_start + 1 if read.next_reference_start != -1 else None,
#             "template_length": read.template_length
#         }

#         for tag, value in read.tags:
#             record.annotations[f"SAM_tag_{tag}"] = value

#         if read.is_unmapped or read.reference_length < reference_length_cutoff or read.mapping_quality < mapq_threshold:
#             bad_dCas9_records.append(record)
#         else:
#             n_seq = record.seq[:aln_start]
#             c_seq = record.seq[aln_end:]

#             if len(n_seq) > 0:
#                 n_record = record[:]
#                 letter_annots = n_record.letter_annotations.copy()
#                 n_record.letter_annotations = {}
#                 n_record.seq = n_seq
#                 n_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][:aln_start]
#                 n_term_records.append(n_record)

#             if len(c_seq) > 0:
#                 c_record = record[:]
#                 letter_annots = c_record.letter_annotations.copy()
#                 c_record.letter_annotations = {}
#                 c_record.seq = c_seq
#                 c_record.letter_annotations["phred_quality"] = letter_annots["phred_quality"][aln_end:]
#                 c_term_records.append(c_record)

#     return (n_term_records, c_term_records, bad_dCas9_records)

def convert_read(read):
    return {
        'query_name': read.query_name,
        'query_sequence': read.query_sequence,
        'query_qualities': read.query_qualities,
        'flag': read.flag,
        'is_reverse': read.is_reverse,
        'is_unmapped': read.is_unmapped,
        'reference_name': read.reference_name,
        'reference_start': read.reference_start,
        'mapping_quality': read.mapping_quality,
        'cigarstring': read.cigarstring,
        'next_reference_name': read.next_reference_name,
        'next_reference_start': read.next_reference_start,
        'template_length': read.template_length,
        'query_alignment_start': read.query_alignment_start,
        'query_alignment_end': read.query_alignment_end,
        'reference_length': read.reference_length,
        'tags': read.tags
    }


def split_reads_at_constant2_parallel(constant2_mapping_sam, n_term_out, c_term_out, bad_reads_out, mapq_threshold=60):
    print(f"Parsing Constant2 mapping file: {constant2_mapping_sam}")
    samfile = pysam.AlignmentFile(constant2_mapping_sam, "r")
    reference_length_cutoff = samfile.lengths[0] * 0.8

    # STEP 1: Convert reads into dictionaries so they are serializable

    # STEP 2: Chunk the reads into serializable chunks
    chunks = chunked_iterator(
        (convert_read(read) for read in samfile.fetch(until_eof=True)),
        CHUNK_SIZE
    )

    # STEP 3: Run multiprocessing on chunks
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.starmap(
            process_read_chunk,
            [(chunk, reference_length_cutoff, mapq_threshold) for chunk in tqdm(chunks, desc="Processing chunks")]
        )

    samfile.close()

    # STEP 4: Merge results
    n_term_records, c_term_records, bad_dCas9_records = [], [], []
    for n_chunk, c_chunk, bad_chunk in results:
        n_term_records.extend(n_chunk)
        c_term_records.extend(c_chunk)
        bad_dCas9_records.extend(bad_chunk)

    print(f"Writing {len(n_term_records)} N-term reads to {n_term_out}")
    print(f"Writing {len(c_term_records)} C-term reads to {c_term_out}")
    print(f"Writing {len(bad_dCas9_records)} bad reads to {bad_reads_out}")

    # STEP 5: Write output FASTQs
    SeqIO.write(n_term_records, n_term_out, "fastq")
    SeqIO.write(c_term_records, c_term_out, "fastq")
    SeqIO.write(bad_dCas9_records, bad_reads_out, "fastq")


# def split_reads_at_constant2_parallel(constant2_mapping_sam, n_term_out, c_term_out, bad_reads_out, mapq_threshold=60):
#     print(f"Parsing Constant2 mapping file: {constant2_mapping_sam}")
#     samfile = pysam.AlignmentFile(constant2_mapping_sam, "r")
#     reference_length_cutoff = samfile.lengths[0] * 0.8

#     n_term_records = []
#     c_term_records = []
#     bad_dCas9_records = []

#     with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
#         results = pool.starmap(
#             process_read_chunk,
#             [(chunk, reference_length_cutoff, mapq_threshold)
#              for chunk in tqdm(chunked_iterator(samfile.fetch(until_eof=True), CHUNK_SIZE),
#                                desc="Processing chunks")]
#         )

#     # Flatten results
#     for n_chunk, c_chunk, bad_chunk in results:
#         n_term_records.extend(n_chunk)
#         c_term_records.extend(c_chunk)
#         bad_dCas9_records.extend(bad_chunk)

#     samfile.close()

#     print(f"Writing {len(n_term_records)} N-term reads to {n_term_out}")
#     print(f"Writing {len(c_term_records)} C-term reads to {c_term_out}")
#     print(f"Writing {len(bad_dCas9_records)} bad reads to {bad_reads_out}")

#     SeqIO.write(n_term_records, n_term_out, "fastq")
#     SeqIO.write(c_term_records, c_term_out, "fastq")
#     SeqIO.write(bad_dCas9_records, bad_reads_out, "fastq")
