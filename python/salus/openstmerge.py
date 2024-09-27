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

def main() -> None:
    if len(sys.argv) < 4:
        print(f'Usage: {sys.argv[0]} <coordinate_system.csv> <spacemake dge path> <outpath>', file=sys.stderr, flush=True)
        exit(0);
    elif len(sys.argv) >= 4:
        incsv = sys.argv[1]
        dgepath = sys.argv[2]
        outpath = sys.argv[3]
        logger.info(f'Coordinate:[{incsv}],DGE:[{dgepath}] => [{outpath}]/xxx')
    csdf = pd.cs_read(incsv,dgepath)
    logger.info(f'Load Coordinate_system: {len(csdf)}')
    print(csdf)

if __name__ == "__main__":
    main()
