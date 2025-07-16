import sys
import argparse
import numpy as np
import scipy
import pandas as pd
from pathlib import Path
import gzip
import subprocess
import scanpy as sc


def read_spatial_coordinates(filename):
    """
    Read spatial coordinates from a file.
    """
    if filename.endswith('tissue_positions_list.csv'):
        coords = pd.read_csv(filename, header=None, index_col=None, sep=',', names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres'])
    elif filename.endswith('tissue_positions.csv'):
        coords = pd.read_csv(filename, header=1, index_col=None, sep=',')
    else:
        raise ValueError(f'Non-standard spatial file: {filename}! File name should be either tissue_positions_list.csv (without header) or tissue_positions.csv (with heaer)')
    
    coords = coords[coords.in_tissue == 1]

    return coords


def merge_spots(coords, n_merge_row, n_merge_col):
    """
    Merge spots in a squared grid. Each merged spot contains n_merge_row spots in a row and n_merge_col spots in a column.

    Attributes
    ----------
    coords: pandas.DataFrame
        Spatial coordinates of spots. It should contain columns 'barcode, 'array_row' and 'array_col'.
    n_merge_row: int
        Number of spots to be merged in a row.
    n_merge_col: int
        Number of spots to be merged in a column.

    Returns
    -------
    merged_coords: pandas.DataFrame
        Merged spatial coordinates of spots. It contains columns 'barcode', 'array_row' and 'array_col', and other columns in coords. The 'barcode' is generated after merging.
    merging_table: pandas.DataFrame
        A table that records the merged spot barcode corresponding to each original barcode.
    """
    merging_table = coords.copy()
    merging_table['merged_row_idx'] = merging_table['array_row'] // n_merge_row # floor division
    merging_table['merged_col_idx'] = merging_table['array_col'] // n_merge_col # floor division
    merging_table['merged_combine_idx'] = merging_table['merged_row_idx'].astype(str) + ',' + merging_table['merged_col_idx'].astype(str)
    
    # generate new barcode for each merged spot: the new barcode is one of the original barcode in the merged spot + suffix '_merged'
    # create a dictionary that maps the merged_combine_idx to the first 'barcode' in the group
    barcode_idx_mapper = merging_table.groupby('merged_combine_idx')['barcode'].first().to_dict()
    merging_table['merged_barcode'] = merging_table['merged_combine_idx'].map(barcode_idx_mapper) + '_merged'

    # create a table of merged_coords by groupping the merging_table and aggregate the 'array_row' by min, 'array_col' by min, and other columns by first
    agg_dict = {col: 'first' for col in merging_table.columns if col not in ['array_row', 'array_col']}
    agg_dict['array_row'] = 'min'
    agg_dict['array_col'] = 'min'
    merged_coords = merging_table.groupby('merged_barcode').agg(agg_dict)
    # drop the columns 'merged_row_idx', 'merged_col_idx', 'merged_combine_idx', 'barcode' in merged_coords
    merged_coords = merged_coords.drop(columns=['merged_row_idx', 'merged_col_idx', 'merged_combine_idx', 'barcode'])
    # rename the columns 'merged_barcode' to 'barcode'
    merged_coords = merged_coords.rename(columns={'merged_barcode': 'barcode'})
    # re-order the columns to be the same as coords
    merged_coords = merged_coords[coords.columns]

    # for merging_table, only retain the columns 'barcode' and 'merged_barcode'
    merging_table = merging_table[['barcode', 'merged_barcode']]

    return merged_coords, merging_table


def update_fragments_using_merging(merging_table, fragments_file, out_dir):
    """
    Update the fragments using the merging table.

    Attributes
    ----------
    merging_table: pandas.DataFrame
        A table that records the merged spot barcode corresponding to each original barcode.
    fragments: str
        Path to the fragments.tsv.gz file. It is the standard output from cellranger multiome ATAC pipeline.
    out_directory: Path object
        Path to the output directory to store a new fragments.tsv.gz file. The output file contains the merged barcode instead of the orignal ones.
    """
    print('Converting the new fragments.tsv file...')

    # convert merging_table to a dictionary
    merging_dict = merging_table.set_index('barcode')['merged_barcode'].to_dict()

    # open the gzipped fragments_file and read it line by line
    with gzip.open(fragments_file, 'rt') as f:
        with open(out_dir / 'fragments.tsv', 'w') as out_f:
            for line in f:
                if line.startswith('#'):
                    out_f.write(line)
                else:
                    strs = line.strip().split('\t')
                    # strs[3] corresponds to barcode
                    if strs[3].split('-')[0] in merging_dict:
                        strs[3] = strs[3].split('-')[0]
                    if strs[3] in merging_dict:
                        strs[3] = merging_dict[strs[3]]
                        out_f.write('\t'.join(strs) + '\n')

    # use bgzip to compress the output fragments.tsv file
    p = subprocess.Popen(['bgzip', str(out_dir / 'fragments.tsv')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        # print log out
        print(f"bgzip error! Please run bgzip {str(out_dir / 'fragments.tsv')} manually.")
        print(f"After bgzip, please run tabix -p bed {str(out_dir / 'fragments.tsv')}.gz to index the file.")
    else:
        # run tabix to index the fragments.tsv.gz file
        p = subprocess.Popen(['tabix', '-p', 'bed', str(out_dir / 'fragments.tsv.gz')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode != 0:
            # print log out
            print(f"tabix error! Please run tabix -p bed {str(out_dir / 'fragments.tsv')}.gz manually.")


def update_RNA_using_merging(merging_table, ordered_new_barcodes, rna_feature_bc_matrix_folder, out_dir):
    """
    Aggregate transcript counts for the merged spots. It makes a new filtered_feature_bc_matrix under out_dir and generates new barcodes.tsv.gz file and matrix.mtx.gz file.

    Attributes
    ----------
    merging_table: pandas.DataFrame
        A table that records the merged spot barcode corresponding to each original barcode.

    ordered_new_barcodes: list
        A list of new barcodes in the order of the new spatial tissue_positions_list.csv file
    
    rna_feature_bc_matrix_folder: str
        Path to the *_feature_bc_matrix folder. It contains barcodes.tsv(.gz), features.tsv(.gz) or genes.tsv(.gz), and matrix.mtx(.gz) files.
    """
    print('Converting the new RNA count matrix...')

    # convert merging_table to a dictionary
    merging_dict = merging_table.set_index('barcode')['merged_barcode'].to_dict()

    # load RNA count matrix
    adata = sc.read_10x_mtx(rna_feature_bc_matrix_folder)

    # create a multiplier matrix (sparse matrix) of dimension (n_old_barcodes, n_new_barcodes), each entry is 0 or 1 to indicate whether the old barcode is merged to the new barcode
    # old barcodes are ordered according to adata.obs.index, new barcodes are ordered according to ordered_new_barcodes
    new_barcodes_order_dict = {x:i for i,x in enumerate(ordered_new_barcodes)}
    rowidx = np.arange(adata.shape[0])
    colidx = np.array([new_barcodes_order_dict.get(merging_dict.get(x, x), -1) for x in adata.obs.index])
    rowidx = rowidx[colidx != -1]
    colidx = colidx[colidx != -1]
    mult_mat = scipy.sparse.csr_matrix((np.ones(len(rowidx)), (rowidx, colidx)), shape=(adata.shape[0], len(ordered_new_barcodes)))

    # multiply the multiplier matrix with the count matrix
    new_counts = mult_mat.T @ adata.X
    
    # write matrix.mtx.gz file and barcodes.tsv.gz file
    (out_dir / 'filtered_feature_bc_matrix').mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(out_dir / 'filtered_feature_bc_matrix' / 'matrix.mtx', new_counts)
    p = subprocess.Popen(['gzip', str(out_dir / 'matrix.mtx')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    with gzip.open(out_dir / 'filtered_feature_bc_matrix' / 'barcodes.tsv.gz', 'wt') as f:
        f.write('\n'.join(ordered_new_barcodes) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge spots in a squared grid and output a new tissue_positions_list.csv and fragments.tsv.gz file.')
    parser.add_argument('--spatial_position_file', type=str, help='Path to the spatial coordinates file.')
    parser.add_argument('--n_merge_row', type=int, default=2, help='Number of spots to be merged in a row.')
    parser.add_argument('--n_merge_col', type=int, default=2, help='Number of spots to be merged in a column.')
    parser.add_argument('--fragments_file', type=str, help='Path to the fragments.tsv.gz file.')
    parser.add_argument('--rna_feature_bc_matrix_folder', type=str, help='Path to the *_feature_bc_matrix folder.')
    parser.add_argument('--out_directory', type=str, help='Path to the output fragments.tsv.gz file.')

    args = parser.parse_args()

    coords = read_spatial_coordinates(args.spatial_position_file)
    merged_coords, merging_table = merge_spots(coords, args.n_merge_row, args.n_merge_col)

    # create output directory
    out_dir = Path(args.out_directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'spatial').mkdir(parents=True, exist_ok=True)
    merged_coords.to_csv(out_dir / 'spatial'/ 'tissue_positions_list.csv', index=False, header=False, sep=',')
    merging_table.to_csv(out_dir / 'merging_table.csv', index=False, sep=',')

    update_fragments_using_merging(merging_table, args.fragments_file, out_dir)

    update_RNA_using_merging(merging_table, merged_coords.barcode.values, args.rna_feature_bc_matrix_folder, out_dir)