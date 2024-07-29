# Spatial_epigenome-transcriptome-proteome_co-sequencing

#### Introduction
This repository aims to share the raw data processing and visualization codes used in the spatial tri-omics project.

## Spatial dynamics of mammalian brain development and neuroinflammation by multimodal tri-omics mapping


![Fig 1-0716](https://github.com/user-attachments/assets/1ac0fb86-8a1c-44a4-ab5b-9419f5cc8c50)

### Data processing
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq 6000 sequencer (pair-end 150 bp mode).
#### Spatial_ATAC-RNA-seq
##### 1.Raw Fastq data
Read 1: contains the spatial Barcode A and Barcode B
Read 2: contains the genome sequences
##### 2. Reformat raw Fastq file to Cell Ranger ARC format (10x Genomics)
**Raw read 1 -> New Read 1 + New Read 2**
- New Read 1: contains the genome sequences
- New Read 2: contains the spatial Barcode A and Barcode B

**Raw read 2 -> New Read 3**

Reformatting raw data was implemented by BC_process.py in the Data_preprocessing folder.


##### 3. Sequence alignment and generation of fragments file
The reformated data was processed using Cell Ranger ARC v2.0.2 with the following references:
Mouse reference (mm10):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz

Human reference (GRCh38):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz




#### 2. Identify useful pixels (pixel on tissue) from microscope image using Matlab
Useful pixels were generated from the Matlab script. Basically, it divide the real tissue microscope image into 100x100 small sqaures which match with DBiT-seq pixels. Then, the intensity inside each pixel was calculated and only pixels have signals above a threashold will be selected.

There two steps: To run the Matlab script "01_pixel100.m"

Use Photoshop or other photo editing software to crop the microscope image into exactly the size of the DBiT-seq covering area. 

Use threashold function under Image->adjustment menu to adjust the image, so that your tissue is black and background is compeletely white.
Invert the color of the image. 

Run the matlab script and a postion.txt file will be generated, which contains only the useful pixels.

Then run 02_tissue_positions_list.R to generate the information needed for the spatial file.







