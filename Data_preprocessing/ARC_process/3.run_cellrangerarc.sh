#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=3.cellrangerARC
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem=64g
#SBATCH --time=120:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=di.zhang@yale.edu
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

module load miniconda
source activate 10XARC

# define variables
sample=sample
tmp_data="./tmp_data"
cellranger_dir="./10X_ARC/cellranger-arc-2.0.2/cellranger-arc"
ref_dir="./10X_ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0/"

file1=$(dirname $tmp_data/${sample}_gex/${sample}_gex_S1_L001_R1_001.fastq)
file2=$(dirname $tmp_data/${sample}_atac/${sample}_atac_S1_L001_R1_001.fastq)


output="../${sample}_libraries.csv"
file1_abs=$(readlink -f $file1)
file2_abs=$(readlink -f $file2)

# create libraries file
echo "fastqs,sample,library_type" > $output
echo "${file1_abs}"/",${sample}_gex,Gene Expression" >> $output
echo "${file2_abs}"/",${sample}_atac,Chromatin Accessibility" >> $output


# run cellranger-arc count
${cellranger_dir} count \
    --id=${sample}_ARC \
    --reference=${ref_dir} \
    --libraries=$output \
    --localcores=16 \
    --localmem=64
