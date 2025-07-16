#!/bin/bash
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=di.zhang@yale.edu
#SBATCH --job-name=merge_spots
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --mem-per-cpu=8g 
#SBATCH --time=0:30:00

module load Java/17
module load miniconda

python merge_spots_squaredgrid.py --spatial_position_file ./spatial/tissue_positions_list.csv  --n_merge_row 2 --n_merge_col 2 --fragments_file ./outs/atac_fragments.tsv.gz --rna_feature_bc_matrix_folder ./outs/raw_feature_bc_matrix/ --out_directory ./merge_spots/out_directory/
