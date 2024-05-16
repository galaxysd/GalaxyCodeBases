#!/usr/bin/env python3

import os
import sys
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
os.environ['NCCL_P2P_LEVEL'] = 'SYS'
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
import copy

import torch
torch.set_float32_matmul_precision('high')

import scanpy as sc
sc.logging.print_header()
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mpl.rcdefaults()
mpl.rc('ps', fonttype=42, papersize='figure')
mpl.rc('pdf', fonttype=42, compression=9)   #pdf.fonttype: 3  # Output Type 3 (Type3) or Type 42 (TrueType)
mpl.rc('figure', figsize=(129/25.4, 129/25.4), dpi=600) # autolayout=True
mpl.rc('savefig', dpi='figure') # bbox='tight'

import pandas as pd
from scipy.sparse import csr_matrix
import pathlib

import scvi
scvi.settings.seed = 42
scvi.settings.dl_num_workers = 16
scvi.settings.num_threads = 8
import cell2location

fn = 'adata31'
tp = '/share/result/spatial/test_huxs/prj/cmpmethod/testicles2/c2ltest'
prefix = '/share/result/spatial/test_huxs/prj/humlung/refdat/GSE112393'

adata_file = f'{prefix}/GSE112393_model.sc.h5ad'
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f'{prefix}/GSE112393_model', adata_ref)
mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator':'gpu'}
)

if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
#inf_aver.iloc[0:5, 0:5]

# Define a function to hijack plt.show()
def hijacked_show(*args, **kwargs):
    pdf.savefig()
    plt.close()
    plt.figure()

sqrt2val = 1.41422
mynacolor = (1,1,1,0)
mycmap = mpl.colormaps.get_cmap('nipy_spectral')  # viridis is the default colormap for imshow
maxArray = 724041728

import copy
mod0 = copy.copy(mod)

spf = f'{tp}/{fn}.h5ad'
adata_vis = sc.read_h5ad(spf)
adata_ori = adata_vis.copy()

adata_vis.obs['sample'] = 'Testicles40'
adata_vis.var['SYMBOL'] = adata_vis.var_names

sc.preprocessing.filter_cells(adata_vis,min_genes=9)
sc.preprocessing.filter_genes(adata_vis,min_cells=3)
adata_vis.X = csr_matrix(adata_vis.X)
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=200
)
mod.view_anndata_setup()

pdf = PdfPages(f"{tp}/{fn}_cell2loc.pdf")
plt.figure()
ax = plt.gca()
ax.set_rasterized(True)
sc.pl.spatial(adata_vis,color='total_counts',spot_size=40,scale_factor=1,title="nCount_Spatial",ax=ax)
pdf.savefig()
plt.close()
plt.figure()
ax = plt.gca()
ax.set_rasterized(True)
sc.pl.spatial(adata_vis,color='n_genes_by_counts',spot_size=40,scale_factor=1,title="nFeature_Spatial",ax=ax)
pdf.savefig()
plt.close()

eee = 30000
#eee = 1000
mod.train(max_epochs=eee,
    # train using full data (batch_size=None)
    batch_size=None,
    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1, accelerator="gpu"
)

plt.figure()
mod.plot_history(1000)
plt.legend(labels=['full data training']);
pdf.savefig()
plt.close()
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator':'gpu'}
)
mod.save(f'{tp}/{fn}_cell2loc_model', overwrite=True)
adata_file = f"{tp}/{fn}_cell2loc_model_sp.h5ad"
adata_vis.write(adata_file)
plt._original_show = plt.show
plt.show = hijacked_show
plt.figure()
mod.plot_QC()
plt.close()
plt.show = plt._original_show
mod.adata.uns['spatial'] = {'Testicles40':{
    'images':{'hires':np.array([0,0,0],ndmin=3)},
}}
#fig = mod.plot_spatial_QC_across_batches()
#pdf.savefig()
#plt.close()
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
torch.cuda.empty_cache()
print('[i]Ploting CellTypes:', end=' ', flush=True)
for onefactor in adata_vis.uns['mod']['factor_names']:
    mytype = ''
    print(f'[{onefactor}]', end=' ', flush=True)
    with mpl.rc_context({'axes.facecolor': 'black'}):
        plt.figure()
        ax = plt.gca()
        ax.set_rasterized(True)
        sc.pl.spatial(adata_vis, #cmap='magma',
            color=onefactor,
            #ncols=int(math.ceil(math.sqrt(len(adata_vis.uns['mod']['factor_names'])))),
            img_key = None,
            spot_size=90,
            cmap=mycmap,na_color=mynacolor,size=sqrt2val, scale_factor=1,
            vmin=0, vmax='p99.2',
            title=f'{onefactor}',
            ax=ax,
        )
        pdf.savefig()
        plt.close()
pdf.close()
print(f'.\n[i]Done.', flush=True)


