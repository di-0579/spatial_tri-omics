#!/bin/bash
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=di.zhang@yale.edu
#SBATCH --job-name=2.CITE-seq
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem-per-cpu=16g 
#SBATCH --time=24:00:00


module load miniconda

source activate DBITpl

sample=P5_Q43
tmp_data="/vast/palmer/scratch/tmp_data"
outs=$tmp_data/${sample}_Prot

CITE_outs="./${sample}_citeouts"
mkdir -p $CITE_outs


# FASTQ reads
FW=$outs/${sample}_Prot_S1_L001_R2_001.fastq.gz
RV=$outs/${sample}_Prot_S1_L001_R1_001.fastq.gz

# Output folder and experiment name
OUTPUT=$CITE_outs

CITE-seq-Count -R1 $FW -R2 $RV -t TAG_LIST_Mouse_ECM_all.csv -trim 21 -cbf 1 -cbl 16 -umif 17 -umil 26 \
-cells 10000 -wl WHITELIST.csv -o $OUTPUT -T 20
