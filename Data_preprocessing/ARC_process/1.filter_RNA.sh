#!/bin/bash
#SBATCH --partition=general 
#SBATCH --job-name=1.filter_RNA
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem=64g
#SBATCH --time=120:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=di.zhang@yale.edu
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

module load miniconda
source activate DBITpl

# Define the variables used in the code
sample=sample1
tmp_data="./tmp_data"
bbmap_path="./useful/bbmap"
fastq_process_path="./ARC_process/RNA_BC_process.py"


tmp=$tmp_data/${sample}
mkdir -p $tmp
outs=$tmp_data/${sample}_gex
mkdir -p $outs


# Run the commands one after another
${bbmap_path}/bbduk.sh \
in1="./usftp21.novogene.com/raw_data/DBiT_2/DBiT_2_CKDL220007951-1a_HK7HLDSX3_L2_1.fq.gz" \
in2="./usftp21.novogene.com/raw_data/DBiT_2/DBiT_2_CKDL220007951-1a_HK7HLDSX3_L2_2.fq.gz" \
outm1=$tmp/${sample}_gex_raw_qc_primer_R1.fastq.gz \
outm2=$tmp/${sample}_gex_raw_qc_primer_R2.fastq.gz \
k=22 mm=f rcomp=f restrictleft=30 skipr1=t \
hdist=2 \
stats=stats.gex_primer.txt \
threads=20 \
literal=CAAGCGTTGGCTTCTCGCATCT

${bbmap_path}/bbduk.sh \
in1=$tmp/${sample}_gex_raw_qc_primer_R1.fastq.gz \
in2=$tmp/${sample}_gex_raw_qc_primer_R2.fastq.gz \
outm1=$tmp/${sample}_gex_raw_qc_linker1_R1.fastq.gz \
outm2=$tmp/${sample}_gex_raw_qc_linker1_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=98 skipr1=t \
hdist=3 \
stats=stats.gex_linker1.txt \
threads=20 \
literal=GTGGCCGATGTTTCGCATCGGCGTACGACT

${bbmap_path}/bbduk.sh \
in1=$tmp/${sample}_gex_raw_qc_linker1_R1.fastq.gz \
in2=$tmp/${sample}_gex_raw_qc_linker1_R2.fastq.gz \
outm1=$tmp/${sample}_gex_raw_qc_R1.fastq.gz \
outm2=$tmp/${sample}_gex_raw_qc_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=60 skipr1=t \
hdist=3 \
stats=stats.gex_linker2.txt \
threads=20 \
literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG

input=$tmp/${sample}_gex_raw_qc_R2.fastq.gz
output=$outs/${sample}_gex_S1_L001_R1_001.fastq

python ${fastq_process_path} --input $input --output $output


input_file=$tmp/${sample}_gex_raw_qc_R1.fastq.gz
output_file=$outs/${sample}_gex_S1_L001_R2_001.fastq.gz

cp $input_file $output_file


