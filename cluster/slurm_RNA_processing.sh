#!/bin/bash -l
#SBATCH --job-name=mapping_PE
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --array=0-18
#SBATCH --mem=100G
#SBATCH --error=mapping%a.err
#SBATCH --output=mapping%a.out
module load bowtie2
module load fastqc
module load rsem
module load jdk
Module load trimmomatic

array=(SC_1 \
SC_2 \
SC_3 \
SC_4 \
SC_5 \
SC_6 \
SC_7 \
SC_8 \
SC_9 \
SC_11 \
SC_12 \
SC_13 \
SC_14 \
SC_16 \
SC_17 \
SC_19 \
SC_20 \
SC_21 \
SC_22 \
SC_24 \
SC_25 \
SC_26 \
SC_28 \
SC_29 \
SC_30 \
SC_31 \
SC_33 \
SC_34 \
SC_36 \
SC_37 \
SC_39 \
SC_40 \
SC_41 \
SC_42 \
SC_43 \
SC_44 \
SC_45 \
SC_46 \
SC_47 \
SC_50 \
SC_51 \
SC_52 \
AL_24_1 \
AL_24_2 \
AL_24_3 \
FR_0H_1 \
FR_0H_2 \
FR_0H_3 \
FR_24_1 \
FR_24_2 \
FR_24_3 \
Foot_1 \
Foot_2 \
HR_0H_1 \
HR_0H_2 \
HR_0H_3 \
HR_24_1 \
HR_24_2 \
HR_24_3 \
Head_1 \
Head_2
)

. RNA_processing.sh ${array[$SLURM_ARRAY_TASK_ID]}
