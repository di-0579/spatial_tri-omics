#!/bin/bash
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=di.zhang@yale.edu
#SBATCH --job-name=1.filter_ADT
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --mem-per-cpu=16g 
#SBATCH --time=24:00:00

module load Java/17
module load miniconda
source activate DBITpl

# Define the variables used in the code
sample=sample
tmp_data="./tmp_data"
bbmap_path="./uesful/bbmap"
fastq_process_path="./CITE_BC_process.py"


tmp=$tmp_data/${sample}
mkdir -p $tmp
outs=$tmp_data/${sample}_Prot
mkdir -p $outs


# Run the commands one after another
${bbmap_path}/bbduk.sh \
in1="./CITE_43c_CKDL240010125-1A_H5TJJDSXC_L4_1.fq.gz" \
in2="./CITE_43c_CKDL240010125-1A_H5TJJDSXC_L4_2.fq.gz" \
outm1=$tmp/${sample}_Prot_raw_qc_primer_R1.fastq.gz \
outm2=$tmp/${sample}_Prot_raw_qc_primer_R2.fastq.gz \
k=22 mm=f rcomp=f restrictleft=30 skipr1=t \
hdist=2 \
stats=stats.Prot_primer.txt \
threads=20 \
literal=CAAGCGTTGGCTTCTCGCATCT

${bbmap_path}/bbduk.sh \
in1=$tmp/${sample}_Prot_raw_qc_primer_R1.fastq.gz \
in2=$tmp/${sample}_Prot_raw_qc_primer_R2.fastq.gz \
outm1=$tmp/${sample}_Prot_raw_qc_linker1_R1.fastq.gz \
outm2=$tmp/${sample}_Prot_raw_qc_linker1_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=98 skipr1=t \
hdist=3 \
stats=stats.Prot_linker1.txt \
threads=20 \
literal=GTGGCCGATGTTTCGCATCGGCGTACGACT

${bbmap_path}/bbduk.sh \
in1=$tmp/${sample}_Prot_raw_qc_linker1_R1.fastq.gz \
in2=$tmp/${sample}_Prot_raw_qc_linker1_R2.fastq.gz \
outm1=$tmp/${sample}_Prot_raw_qc_R1.fastq.gz \
outm2=$tmp/${sample}_Prot_raw_qc_R2.fastq.gz \
k=30 mm=f rcomp=f restrictleft=60 skipr1=t \
hdist=3 \
stats=stats.Prot_linker2.txt \
threads=20 \
literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG

input=$tmp/${sample}_Prot_raw_qc_R2.fastq.gz
output=$outs/${sample}_Prot_S1_L001_R2_001.fastq

python ${fastq_process_path} --input $input --output $output


input_file=$tmp/${sample}_Prot_raw_qc_R1.fastq.gz
output_file=$outs/${sample}_Prot_S1_L001_R1_001.fastq.gz

cp $input_file $output_file


input_file=$outs/${sample}_Prot_S1_L001_R2_001.fastq
output_file=$outs/${sample}_Prot_S1_L001_R2_001.fastq.gz

gzip $input_file $output_file
