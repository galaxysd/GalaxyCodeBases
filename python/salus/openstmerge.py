#!/usr/bin/env python3

import sys
import os
import glob
from xopen import xopen
import fast_matrix_market as fmm
import scipy.sparse
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd

import logging
logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)-8s|%(message)s\t|%(module)s:%(funcName)s:%(lineno)d|%(process)d:%(thread)d|%(pathname)s',
    datefmt='%Y-%m-%dT%H:%M:%S',level=logging.INFO)
logger = logging.getLogger(__name__)

from concurrent.futures import ThreadPoolExecutor, as_completed
def export_adata(adata, path):
    """
    Export AnnData object to 10x Genomics format with gzip compression using xopen.
    Parameters:
    adata : AnnData
        The AnnData object to export.
    path : str
        The directory where the files should be saved.
    """
    # Create destination directory if it doesn't exist
    os.makedirs(path, exist_ok=True)
    # Step 1: Append spatial coordinates directly to obs
    if 'spatial' in adata.obsm:
        adata.obs['m_spatial_x'] = adata.obsm['spatial'][:, 0]  # X coordinates
        adata.obs['m_spatial_y'] = adata.obsm['spatial'][:, 1]  # Y coordinates
    # Prepare file paths
    genes_path = os.path.join(path, "features.tsv.gz")
    barcodes_path = os.path.join(path, "barcodes.tsv.gz")
    spatial_path = os.path.join(path, "spatial.txt.gz")
    matrix_path = os.path.join(path, "matrix.mtx.gz")
    # Prepare dataframes for writing
    genes_df = pd.DataFrame({
        'gene_name': adata.var.index,
        'feature_types': "Gene Expression"
    }, index=['g' + str(i + 1) for i in range(len(adata.var.index))])  # Set gene_id as index
    barcodes_with_spatial = adata.obs[['m_spatial_x', 'm_spatial_y']]
    obs_df = adata.obs
    # Ensure that adata.X is a sparse matrix
    if not isinstance(adata.X, scipy.sparse.spmatrix):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    # Check equality of float and int representations
    float_matrix = adata.X
    int_matrix = float_matrix.astype(int)
    # If the two matrices are equal, convert to integer
    if (float_matrix != int_matrix).nnz == 0:
        adata.X = int_matrix  # Set X to its integer representation

    compression_parameters = {'method': 'gzip', 'mtime': 1}
    def write_genes():
        genes_df.to_csv(
            genes_path,
            sep="\t",
            index=True,  # Write index as well, which will now be gene_id
            header=False,
            compression=compression_parameters
        )
    def write_barcodes():
        # Write only the index (barcodes) to the file, without a header
        barcodes_with_spatial.index.to_series().to_csv(
            barcodes_path,
            sep="\t",
            index=True,
            header=False,
            compression=compression_parameters
        )
    def write_spatial():
        barcodes_with_spatial.to_csv(
            spatial_path,
            sep=" ",
            index=True,
            header=False,
            compression={**compression_parameters, 'compresslevel': 1}
        )
    def write_matrix():
        with xopen(matrix_path, 'wb') as f:
            fmm.mmwrite(f, adata.X.T)

    # Use ThreadPoolExecutor to write files concurrently
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(write_genes): "genes",
            executor.submit(write_barcodes): "barcodes",
            executor.submit(write_spatial): "spatial",
            executor.submit(write_matrix): "matrix"
        }
        for future in as_completed(futures):
            task_name = futures[future]
            future.result()  # Get result to raise exceptions if any
# export_adata(combined_adata, 'path/to/output_directory')

def cs_read(incsv, dgepath):
    csdf = pd.read_csv(incsv)
    csdf['h5ad_path'] = pd.NA
    pattern = 'dge.all.polyA_adapter_trimmed.mm_included.spatial_beads_*.h5ad'
    search_files = glob.glob(os.path.join(dgepath, pattern))
    for file_path in search_files:
        base_name = os.path.basename(file_path)
        puck_id = base_name.split('.')[-2].replace('spatial_beads_', '')
        if puck_id in csdf['puck_id'].values:
            # Get the corresponding row
            row_index = csdf[csdf['puck_id'] == puck_id].index[0]
            csdf.at[row_index, 'h5ad_path'] = file_path
            x_offset = csdf.at[row_index, 'x_offset']
            y_offset = csdf.at[row_index, 'y_offset']
            #print(f"[{puck_id}, {x_offset}, {y_offset}] {file_path}")
    csdfused = csdf[csdf['h5ad_path'].notna()].drop(columns=['z_offset'])
    return csdfused

def sp_concat(csdf):
    adatas = []
    for index, row in csdf.iterrows():
        h5ad_path = row['h5ad_path']
        adata = sc.read_h5ad(h5ad_path)
        logger.info(f'Load anndata {row["puck_id"]}: {adata.shape}')
        # Apply x_offset and y_offset
        x_offset = row['x_offset']
        y_offset = row['y_offset']
        # Shift the spatial coordinates
        if 'spatial' in adata.obsm.keys():
            # Create a 1D array for offsets
            offsets = np.array([x_offset, y_offset])
            adata.obsm['spatial'] += offsets  # Adjust the spatial coordinates
            logger.info(f'Shift anndata {row["puck_id"]}: x+={x_offset}, y+={y_offset}')
        # Append the modified AnnData object to the list
        adatas.append(adata)
    # Step 2: Use anndata.concat to merge all AnnData objects with outer join on .var
    combined_adata = ad.concat(adatas, join='outer', label='source')
    return combined_adata

def main() -> None:
    if len(sys.argv) < 4:
        print(f'Usage: {sys.argv[0]} <coordinate_system.csv> <spacemake dge path> <outpath>', file=sys.stderr, flush=True)
        exit(0);
    elif len(sys.argv) >= 4:
        incsv = sys.argv[1]
        dgepath = sys.argv[2]
        outpath = sys.argv[3]
        logger.info(f'Coordinate:[{incsv}],DGE:[{dgepath}] => [{outpath}]/xxx')
    csdf = cs_read(incsv,dgepath)
    logger.info(f'Load Coordinate_system: {len(csdf)}')
    print(csdf)
    os.makedirs(outpath, exist_ok=True)
    combined_adata = sp_concat(csdf)
    logger.info(f'Merged anndata: {combined_adata.shape}')
    h5adfilename = os.path.join(outpath, 'raw_puck_collection.h5ad')
    combined_adata.write_h5ad(h5adfilename,compression='lzf')
    logger.info(f'Saved anndata: {h5adfilename}')
    mtxpath = os.path.join(outpath, 'raw')
    export_adata(combined_adata, mtxpath)
    logger.info(f'Exported anndata: {mtxpath}')

if __name__ == "__main__":
    main()
