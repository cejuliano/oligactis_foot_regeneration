#!/bin/bash -l
#SBATCH --job-name=dovetail_ref
#SBATCH -t 4-0
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --error=OLI_ref.err
#SBATCH --output=OLI_ref.out

module load rsem/1.2.31

module load bowtie2

./prep_reference.sh
