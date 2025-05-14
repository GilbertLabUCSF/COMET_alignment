create environment with following command: 
conda env create --name COMET  --file=environment.yml


example run script:

python run_pipeline.py \
    --yaml_file /home/nicholass/COMET/Tests/inputs/test.yml \
    --output_dir Tests/outputs \
    --threads 128
