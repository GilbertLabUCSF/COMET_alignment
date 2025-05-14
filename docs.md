this is my work on creating an alignment pipeline for COMET effectors.

this expects data in the format:
    Constant1 - VariableN - Constant2 - VariableC - Constant3

where Constant1 and Constant3 are up and downstream constant regions that were amplified from gDNA
Constant 2 is Xten-dCas9-Xten
VariableN and VariableC are variable sequences that will be used for calling domain combinations.



Input:
Long-read sequencing data (reads.fastq)
Known Constant2 sequence (constant2.fasta)
Candidate N-terminal sequences (N_candidates.fasta)
Candidate C-terminal sequences (C_candidates.fasta)

Pipeline Steps:
1. Mapping Reads to Constant2
Each sequencing read is aligned to the provided Constant2 sequence using minimap2. This identifies whether the Constant2 sequence is present within the read.

Reads that fail to align are saved to a separate file (bad_reads.fastq).

Reads that successfully align are annotated with the precise start and end positions of the Constant2 region.

Rationale: Constant2 is ~4 kb, so full-length or partial (~80%+) alignments are required to confidently detect its presence.


2. Splitting Reads into N- and C-Terminal Regions
Reads containing Constant2 are split into two fragments:

N-terminal region: the sequence preceding the Constant2 alignment.

C-terminal region: the sequence following the Constant2 alignment.

Each fragment is output as a new FASTQ file (n_term.fastq, c_term.fastq).


3. Mapping N- and C-Terminal Regions to Candidate Sequences
Each N-terminal and C-terminal fragment is independently mapped against the provided candidate databases:

N-terminal fragments are aligned to the N candidates (N_candidates.mmi).

C-terminal fragments are aligned to the C candidates (C_candidates.mmi).

Minimap2 is used again for alignment, with parameters optimized for long, noisy reads (e.g., ONT, PacBio CLR). Only high-confidence alignments (based on mapping quality and alignment coverage) are retained.


4. Assignment and Classification
For each read:

The best N-terminal match and best C-terminal match are recorded.

Optional: Reads with ambiguous mappings or low-confidence alignments can be flagged or discarded.

A summary table is generated, listing for each read:

| Read ID | Best N-terminal match | Best C-terminal match | Additional metrics (e.g., alignment identity, coverage) |


Output
bad_reads.fastq — Reads without detectable Constant2.
good_reads_n_term.fastq — Extracted N-terminal sequences.
good_reads_c_term.fastq — Extracted C-terminal sequences.
n_term.sam — Alignments of N-terms to N-candidates.
c_term.sam — Alignments of C-terms to C-candidates.
n_c_summary.tsv — Table mapping each read to its N and C candidates.




<!-- 316 and 290 are the same sequence so I will manually delete 316 from the reference -->



usage: run_pipeline.py [-h] --reads READS --constant2 CONSTANT2 --n_candidates N_CANDIDATES --c_candidates C_CANDIDATES [--output_dir OUTPUT_DIR] [--threads THREADS]

Pipeline to process long-read sequencing data for N- and C-terminal identification.

optional arguments:
  -h, --help            show this help message and exit
  --reads READS         Input FASTQ file with reads.
  --constant2 CONSTANT2 FASTA file of Constant2 (e.g., dCas9 sequence).
  --n_candidates N_CANDIDATES FASTA file with candidate N-terminal sequences.
  --c_candidates C_CANDIDATES FASTA file with candidate C-terminal sequences.
  --output_dir OUTPUT_DIR Directory to save all outputs. (default: output)
  --threads THREADS     Number of threads for minimap2. (default: 8)
