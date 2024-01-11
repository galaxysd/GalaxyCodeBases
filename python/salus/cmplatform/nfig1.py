#!/usr/bin/env python3
import sys
import os
from typing import NamedTuple

PlatformTuple = ('Illumina', 'Salus')
SamplesDict = {
    'mbrain': {
        'sid' : 'mbrain',
        'sub' : 'Mouse Brain Sptial',
        'type': 'visium',
        'fltPct'   : 99.5,
        'prefix'   : '/share/result/spatial/data/BoAo_sp',
        'suffixOut': dict.fromkeys(PlatformTuple,"outs"),
        'suffixMtx': 'filtered_feature_bc_matrix',
        'platforms': {PlatformTuple[0]:'illumina', PlatformTuple[1]: 'salus'},
        'pattern': ('prefix', 'platformV', 'sid', 'suffixOutV', 'suffixMtx')
    },
    'mkidney': {
        'sid' : 'mkidney',
        'sub' : 'Mouse Kindey Sptial',
        'type': 'visium',
        'fltPct'   : 99.5,
        'prefix'   : '/share/result/spatial/data/BoAo_sp',
        'suffixOut': dict.fromkeys(PlatformTuple,"outs"),
        'suffixMtx': 'filtered_feature_bc_matrix',
        'platforms': {PlatformTuple[0]:'illumina', PlatformTuple[1]: 'salus'},
        'pattern': ('prefix', 'platformV', 'sid', 'suffixOutV', 'suffixMtx')
    },
    'human': {
        'sid' : 'human',
        'sub' : 'Human Single Cell',
        'type': 'mobivision',
        'fltPct'   : 85,
        'prefix'   : '/share/result/spatial/data/MoZhuo_sc/FX20230913',
        'suffixOut': {PlatformTuple[0]: 'out/R22045213-220914-LYY-S11-R03-220914-LYY-S11-R03_combined_outs',
                      PlatformTuple[1]: 'out_subset/20221124-LYY-S09-R03_AGGCAGAA_fastq_outs'},
        'suffixMtx': 'filtered_cell_gene_matrix',
        'platforms': {PlatformTuple[0]:'illumina', PlatformTuple[1]: 'sailu'},
        'pattern': ('prefix', 'platformV', 'suffixOutV', 'suffixMtx')
    }
}

thisID = 'mbrain'
if __name__ == "__main__":
    if len(sys.argv) > 1:
        thisID = sys.argv[1]
        if thisID not in SamplesDict:
            print(f"[x]sid can only be {SamplesDict.keys()}", file=sys.stderr)
            exit(1)
    print(sys.argv, file=sys.stderr)
    print(f"[i]{thisID}")
    sys.stdout.flush()
nfoDict = SamplesDict[thisID]

import matplotlib; matplotlib.use("module://mplcairo.base")
from matplotlib import pyplot as plt
import mplcairo

plt.rcParams['figure.figsize'] = (6.0, 6.0) # set default size of plots
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams["savefig.transparent"] = True
font = {'family' : 'STIX Two Text',
        #'size'   : 22,
        'weight' : 'normal'}
matplotlib.rc('font', **font)

import numpy as np
import pandas as pd
import fast_matrix_market
import anndata as ad
import leidenalg
import scanpy as sc
sc._settings.ScanpyConfig.n_jobs = -1
#import squidpy as sq
import seaborn as sns
import scipy
import pynndescent

import warnings
warnings.filterwarnings('ignore')
from copy import deepcopy

class scDatItem(NamedTuple):
    name: str
    rawDat: ad.AnnData
    annDat: ad.AnnData
    def __repr__(self) -> str:
        return f'[sc:{self.name}, BC*Gene: Raw={self.rawDat.shape}, Filtered={self.annDat.shape}]'

def main() -> None:
    scDat = []
    #nfoDict = SamplesDict[thisID]
    print("[i]Start.", file=sys.stderr)
    for platform in PlatformTuple:
        nfoDict['platformK']  = platform
        nfoDict['platformV']  = nfoDict['platforms'][platform]
        nfoDict['suffixOutV'] = nfoDict['suffixOut'][platform]
        mtxPath = os.path.join( *[nfoDict[v] for v in nfoDict['pattern']] )
        print(f"[i]Reading {mtxPath}", file=sys.stderr)
        adata=sc.read_10x_mtx(mtxPath, var_names='gene_symbols', make_unique=True, gex_only=True)
        adata.var_names_make_unique()
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
        rdata = deepcopy(adata)
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
        sc.pp.filter_cells(adata, min_genes=1)
        sc.pp.filter_genes(adata, min_cells=1)
        scDat.append(scDatItem(platform,rdata,adata))
        rdata.write_h5ad(f"{nfoDict['sid']}_{platform}.raw.h5ad",compression='lzf')
    print("\n".join(map(str,scDat)))
    with pd.option_context("mode.copy_on_write", True):
        obsmbi = scDat[0].annDat.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
        obsmbs = scDat[1].annDat.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
        p1df = pd.concat([obsmbi.assign(Platform=scDat[0].name), obsmbs.assign(Platform=scDat[1].name)], ignore_index=True).replace([np.inf, -np.inf, 0], np.nan).dropna()
        p2df = obsmbi.join(obsmbs,lsuffix='_'+scDat[0].name,rsuffix='_'+scDat[1].name,how='inner').replace([np.inf, -np.inf, 0], np.nan).dropna()
        p3tuple = (frozenset(scDat[0].annDat.var_names), frozenset(scDat[1].annDat.var_names))

    metapdf={'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'}
    print("[i]Begin fig A. 1D", file=sys.stderr)
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params, font="STIX Two Text")
    figA=sns.JointGrid(data=p1df, x="total_counts", y="n_genes_by_counts", hue='Platform', dropna=True)
    #figA.plot(sns.scatterplot, sns.histplot, alpha=.7, edgecolor=".2", linewidth=.5)
    figA.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
    figA.plot_marginals(sns.histplot, kde=False, alpha=.618)
    figA.figure.suptitle(f"Gene to UMI plot - {nfoDict['sub']}")
    figA.set_axis_labels(xlabel='UMIs per Barcode', ylabel='Genes per Barcode')
    figA.savefig(f"1D_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'Gene to UMI plot'})

    print("[i]Begin fig B. 1E", file=sys.stderr)
    figB=sns.JointGrid(data=p2df, x="total_counts_Illumina", y="total_counts_Salus", dropna=True)
    figB.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
    figB.plot_marginals(sns.histplot, kde=True, alpha=.618)
    figB.figure.suptitle(f"UMI per Barcode Counts Comparing - {nfoDict['sub']}")
    figB.set_axis_labels(xlabel='UMI Counts from Illumina', ylabel='UMI Counts from Salus')
    figB.savefig(f"1E_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'UMI per Barcode Counts Comparing'})

    print("[i]Begin fig . 1G", file=sys.stderr)
    from matplotlib_venn import venn2
    plt.figure(figsize=(4,4))
    plt.title(f"Genes Venn diagram - {nfoDict['sub']}")
    p3intersection = p3tuple[0] & p3tuple[1]
    p3veen = (p3tuple[0]-p3intersection, p3tuple[1]-p3intersection, p3intersection)
    GenesA = scDat[0].annDat.var.loc[p3veen[0]-p3veen[2]]
    GenesB = scDat[1].annDat.var.loc[p3veen[1]-p3veen[2]]
    GenesC = scDat[0].annDat.var.loc[p3veen[2]]
    p3vd=venn2(subsets=tuple(map(len,p3veen)), set_labels=(scDat[0].name, scDat[1].name))
    plt.savefig(f"1G_Genes_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'Veen of Genes'})
    GenesA.to_csv(f"1G_Genes_{nfoDict['sid']}_{scDat[0].name}_only.csv",encoding='utf-8')
    GenesB.to_csv(f"1G_Genes_{nfoDict['sid']}_{scDat[1].name}_only.csv",encoding='utf-8')
    GenesC.to_csv(f"1G_Genes_{nfoDict['sid']}_intersection.csv.zst",encoding='utf-8',compression={'method': 'zstd', 'level': 9, 'write_checksum': True})

    print("[i]Begin fig C. 2A", file=sys.stderr)
    # https://www.kaggle.com/code/lizabogdan/top-correlated-genes?scriptVersionId=109838203&cellId=21
    p4xdf = scDat[0].annDat.to_df()
    p4ydf = scDat[1].annDat.to_df()
    p4corraw = p4xdf.corrwith(p4ydf,axis=1)
    p4corr = p4corraw.dropna()
    plt.figure(figsize=(6,4))
    plt.title(f"Pearson correlation - {nfoDict['sub']}")
    figC=sns.histplot(p4corr,stat='count',binwidth=0.01)
    plt.savefig(f"2A_Correlation_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'Pearson correlation'})
    '''
    print("[i]Begin fig D. 2B", file=sys.stderr)
    var_names = scDat[0].annDat.var_names.intersection(scDat[1].annDat.var_names)
    xadata = scDat[0].annDat[:, var_names]
    yadata = scDat[1].annDat[:, var_names]
    xdf=getOBSMdf(xadata)
    ydf=getOBSMdf(yadata)
    #p4df = xdf.assign(Platform=scDat[0].name).join(ydf.assign(Platform=scDat[1].name),lsuffix='_'+scDat[0].name,rsuffix='_'+scDat[1].name,how='inner')
    p4df = pd.concat([xdf.assign(Platform=scDat[0].name), ydf.assign(Platform=scDat[1].name)], ignore_index=True).replace([np.inf, -np.inf, 0], np.nan).dropna()
    figD=sns.JointGrid(data=p4df, x="P1", y="P2", hue='Platform', dropna=True)
    figD.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
    figD.plot_marginals(sns.histplot, kde=True, alpha=.618)
    figD.figure.suptitle(f"PCA - {nfoDict['sub']}")
    figD.set_axis_labels(xlabel='PC1', ylabel='PC2')
    figD.savefig(f"2B_rawPCA_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'PCA'})
    '''
    import scvi
    for IDlist in ([0],[1],[0,1]):
        rawList = [scDat[i].rawDat for i in IDlist]
        dataIDs = [scDat[i].name for i in IDlist]
        if len(rawList) == 1:
            adata = rawList[0]
            dataID = dataIDs[0]
        elif len(rawList) == 2:
            adata=ad.concat(rawList, label='Platform', keys=PlatformTuple, index_unique='-')
            dataID = 'Both'
        print(f"[i]Begin Tab 1. 1F Dropout rates - {dataID}. With scvi {scvi.__version__}", file=sys.stderr)
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
        if dataID == 'Both':
            scvi.data.poisson_gene_selection(adata,n_top_genes=8000,n_samples=10000,batch_key='Platform')
        else:
            scvi.data.poisson_gene_selection(adata,n_top_genes=8000,n_samples=10000)
        doDropOutPlot(dataID,adata)
        adata = None

def doDropOutPlot(dataID,adata) -> None:
    adata.var['mean_'] = np.array(adata.X.mean(0))[0]
    GenesM = adata.var.sort_values(by='prob_zero_enrichment_rank', ascending=False)
    GenesM.to_csv(f"1F_GenesDropout_{nfoDict['sid']}_{dataID}_PlatformAsBatch.csv.zst",encoding='utf-8',compression={'method': 'zstd', 'level': 9, 'write_checksum': True})
    print(f"[i]Begin Fig 1. 1F GenesM3DropSelected (added) - {dataID}", file=sys.stderr)
    highly_variable_df = adata.var.query('highly_variable')
    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))
    # Create the scatter plot for the main points with color bar
    scatter = sns.scatterplot(x='mean_', y='observed_fraction_zeros', hue='prob_zero_enrichment', data=adata.var, palette='viridis', legend='brief')
    # Create the line plot for expected_fraction_zeros
    sns.lineplot(x='mean_', y='expected_fraction_zeros', data=adata.var, color='r', label='Expected Fraction Zeros')
    # Highlight highly variable points
    sns.scatterplot(x='mean_', y='observed_fraction_zeros', data=highly_variable_df, color='pink', marker='.', s=5, alpha=0.5)
    box_coords = adata.var.query('highly_variable').agg({'mean_': ['min', 'max'], 'observed_fraction_zeros': ['min', 'max']})
    # Draw a rectangle to cover highly variable points
    rect = plt.Rectangle(box_coords.loc['min'],
                         box_coords['mean_'].diff()['max'], box_coords['observed_fraction_zeros'].diff()['max'],
                         fill=None, edgecolor='blue', linewidth=2, alpha=0.5)
    ax.add_patch(rect)
    # Annotate right-top and left-bottom points
    fmt = '.4f'
    for mean_val, obs_frac_val in zip(box_coords['mean_'], box_coords['observed_fraction_zeros']):
        label = f'({mean_val:{fmt}},{obs_frac_val:{fmt}})'
        # Add padding to avoid overlapping with the rectangle
        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="white", lw=1, alpha=0.62)
        ax.text(mean_val, obs_frac_val, label, bbox=bbox_props)
    # Set x-axis to log scale
    ax.set_xscale('log')
    # Set plot title
    ax.set_title(f'Mean vs Observed Fraction Zeros - {nfoDict["sub"]} {dataID}')
    # Create a color bar for Prob Zero Enrichment
    cbar = fig.colorbar(scatter.get_children()[0], ax=ax, orientation='vertical', pad=0.1)
    cbar.set_label('Prob Zero Enrichment')
    plt.savefig(f"1F_GenesM3DropSelected_{nfoDict['sid']}_{dataID}_PlatformAsBatch.pdf", metadata={'Title': 'scvi.data.poisson_gene_selection', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.close('all')
    print(f"[i]Begin Fig 1. 1F GenesDropoutHist (added) - {dataID}", file=sys.stderr)
    plt.figure(figsize=(6,4))
    plt.title(f"Gene DropRatio Histogram - {nfoDict['sub']} {dataID}")
    histplot = sns.histplot(adata.var, x='observed_fraction_zeros', bins=30, kde=False, hue='highly_variable', multiple="dodge", shrink=.8)
    bars_heights = [p.get_height() for p in histplot.patches if p.get_facecolor()[:3] == sns.color_palette()[1]]
    plt.ylim(0, max(bars_heights)*1.1)  # Adjust the margin as needed
    plt.savefig(f"1F_GenesDropoutHist_{nfoDict['sid']}_{dataID}_PlatformAsBatch.pdf", metadata={'Title': 'Gene DropRatio Histogram', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.close('all')

def getOBSMdf(anndata, obsmkey='X_pca') -> pd.DataFrame:
    if not obsmkey in anndata.obsm:
        if obsmkey=='X_pca':
            sc.tl.pca(anndata,zero_center=True)
    data=anndata.obsm[obsmkey][0:,0:2]
    df=pd.DataFrame(data=data[0:,0:], index=[anndata.obs_names[i] for i in range(data.shape[0])], columns=['P'+str(1+i) for i in range(data.shape[1])])
    return df

if __name__ == "__main__":
    main()  # time (./nfig1.py human; ./nfig1.py mbrain ; ./nfig1.py mkidney ) | tee nplot.log

# pip install -U --force-reinstall lightning
"""
pip3 install git+https://github.com/matplotlib/mplcairo
pip3 install matplotlib_venn
"""
