#!/bin/bash
#SBATCH --job-name=Mdm4_colab       # Give your job a name
#SBATCH --partition=norm                # Specify the correct GPU partition
#SBATCH --mem=200G                      # Request 256GB of memory (adjust as necessary)
#SBATCH --cpus-per-task=16              # Request 32 CPU cores (adjust as necessary)
#SBATCH --gres=lscratch:50             # Request 100GB of local scratch space               
#SBATCH --time=3-24:00:00                   # Set a time limit for the job (6 hours)
#SBATCH --output=output/Mdm4_colab.log    # Unique output log for each job
#SBATCH --error=error/Mdm4_colab.log      # Unique error log for each job
#SBATCH --mail-type=END,FAIL   # Send email on job completion (END) or failure (FAIL)
#SBATCH --mail-user=alex.castroverde@nih.gov  # Replace with your email address

set -e
set -x  # Enable debugging

# Load necessary modules
module load colabfold alphapulldown/0.30.7

# Specify your input and output directories
INPUT_DIR="/data/CBLCCBR/af_new/Mdm4"   # Set input directory dynamically
OUTPUT_DIR="/data/CBLCCBR/af_new/Mdm4"      # Set output directory dynamically

# Create individual features
create_individual_features.py --fasta_paths=$INPUT_DIR/full.fasta,$INPUT_DIR/tiled.fasta --output_dir=$OUTPUT_DIR/template_pulldown_cf_msas --use_precomputed_msas=True --max_template_date=2023-01-01 --use_mmseqs2=True --skip_existing=True

# Run colabfold search
colabfold_search --threads $SLURM_CPUS_PER_TASK $OUTPUT_DIR/for_colabfold.fasta $COLABFOLD_DB $OUTPUT_DIR/no_pulldown_cf_msas

# Create individual features
create_individual_features.py --fasta_paths=$INPUT_DIR/full.fasta,$INPUT_DIR/tiled.fasta --output_dir=$OUTPUT_DIR/no_pulldown_cf_msas --use_precomputed_msas=False --max_template_date=1960-01-01 --use_mmseqs2=True --skip_existing=True
