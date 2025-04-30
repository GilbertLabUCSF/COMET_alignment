create environment with following command: 
conda env create --name COMET  --file=environment.yml



example run script:

python run_pipeline.py \
    --reads path_to_fastq.gz \
    --constant2 Tests/references/dCas9.fasta \
    --n_candidates Tests/references/trim_nterm.fasta \
    --c_candidates Tests/references/trim_cterm.fasta \
    --output_dir Tests/outputs \
    --threads 128

    
