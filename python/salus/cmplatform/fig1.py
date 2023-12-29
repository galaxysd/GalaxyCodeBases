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
        'prefix'   : '/share/result/spatial/data/MoZhuo_sc/FX20230913',
        'suffixOut': {PlatformTuple[0]: 'out/R22045213-220914-LYY-S11-R03-220914-LYY-S11-R03_combined_outs',
                      PlatformTuple[1]: 'out_subset/20221124-LYY-S09-R03_AGGCAGAA_fastq_outs'},
        'suffixMtx': 'filtered_cell_gene_matrix',
        'platforms': {PlatformTuple[0]:'illumina', PlatformTuple[1]: 'sailu'},
        'pattern': ('prefix', 'platformV', 'suffixOutV', 'suffixMtx')
    }
}

def checkModules() -> None:
    import importlib.metadata
    from packaging import version
    pkgname = "squidpy"
    min_ver = "1.2.3"
    got_ver = importlib.metadata.version(pkgname)
    if version.parse(got_ver) < version.parse(min_ver):
        raise importlib.VersionConflict(f"{pkgname}>={min_ver} is needed, but found {pkgname}=={got_ver}")

def main() -> None:
    import matplotlib; matplotlib.use("module://mplcairo.base")
    from matplotlib import pyplot as plt
    import mplcairo

    plt.rcParams['figure.figsize'] = (6.0, 6.0) # set default size of plots
    font = {'family' : 'STIX Two Text',
            #'size'   : 22,
            'weight' : 'normal'}
    matplotlib.rc('font', **font)

    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import squidpy as sq
    import seaborn as sns

    class scDatItem(NamedTuple):
        name: str
        bgRaw: tuple[int,int]
        bgFlt: tuple[int,int]
        annDat: ad.AnnData

        def __repr__(self) -> str:
            return f'[sc:{self.name}, Raw_BC*Gene={self.bgRaw[0]}x{self.bgRaw[1]}, NonZero_BC*Gene={self.bgFlt[0]}x{self.bgFlt[1]} ({self.annDat.n_obs}x{self.annDat.n_vars})]'

    scDat = []
    nfoDict = SamplesDict[thisID]
    for platform in PlatformTuple:
        nfoDict['platformK']  = platform
        nfoDict['platformV']  = nfoDict['platforms'][platform]
        nfoDict['suffixOutV'] = nfoDict['suffixOut'][platform]
        mtxPath = os.path.join( *[nfoDict[v] for v in nfoDict['pattern']] )
        #print(mtxPath)
        adata=sc.read_10x_mtx(mtxPath, var_names='gene_symbols', make_unique=True, gex_only=True)
        adata.var_names_make_unique()  # this is necessary if using `var_names='gene_symbols'` in `sc.read_10x_mtx`
        nnRaw = (adata.n_obs,adata.n_vars)
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
        sc.pp.filter_cells(adata, min_genes=1)
        sc.pp.filter_genes(adata, min_cells=1)
        nnFlt = (adata.n_obs,adata.n_vars)
        scDat.append(scDatItem(platform,nnRaw,nnFlt,adata))

    print("\n".join(map(str,scDat)))

    with pd.option_context("mode.copy_on_write", True):
        obsmbi = scDat[0].annDat.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
        obsmbs = scDat[1].annDat.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
        p1df = pd.concat([obsmbi.assign(Platform=scDat[0].name), obsmbs.assign(Platform=scDat[1].name)], ignore_index=True).replace([np.inf, -np.inf, 0], np.nan).dropna()
        p2df = obsmbi.join(obsmbs,lsuffix='_'+scDat[0].name,rsuffix='_'+scDat[1].name,how='inner').replace([np.inf, -np.inf, 0], np.nan).dropna()
        p3tuple = (frozenset(scDat[0].annDat.var.index), frozenset(scDat[1].annDat.var.index))

    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params, font="STIX Two Text")
    figA=sns.JointGrid(data=p1df, x="total_counts", y="n_genes_by_counts", hue='Platform', dropna=True)
    #figA.plot(sns.scatterplot, sns.histplot, alpha=.7, edgecolor=".2", linewidth=.5)
    figA.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
    figA.plot_marginals(sns.histplot, kde=True, alpha=.618)
    figA.figure.suptitle(f"Gene to UMI plot - {nfoDict['sub']}")
    figA.set_axis_labels(xlabel='UMIs per Barcode', ylabel='Genes per Barcode')
    figA.savefig(f"1D_{nfoDict['sid']}.pdf", transparent=True, dpi=300, metadata={'Title': 'Gene to UMI plot', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    figB=sns.JointGrid(data=p2df, x="total_counts_Illumina", y="total_counts_Salus", dropna=True)
    figB.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
    figB.plot_marginals(sns.histplot, kde=True, alpha=.618)
    figB.figure.suptitle(f"UMI per Barcode Counts Comparing - {nfoDict['sub']}")
    figB.set_axis_labels(xlabel='UMI Counts from Illumina', ylabel='UMI Counts from Salus')
    figB.savefig(f"1E_{nfoDict['sid']}.pdf", transparent=True, dpi=300, metadata={'Title': 'UMI per Barcode Counts Comparing', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    from matplotlib_venn import venn2
    plt.figure(figsize=(4,4))
    plt.title("Sample Venn diagram")
    p3intersection = p3tuple[0] & p3tuple[1]
    p3veen = (p3tuple[0]-p3intersection, p3tuple[1]-p3intersection, p3intersection)
    GenesA = scDat[0].annDat.var.loc[p3veen[0]-p3veen[2]]
    GenesB = scDat[1].annDat.var.loc[p3veen[1]-p3veen[2]]
    GenesC = scDat[0].annDat.var.loc[p3veen[2]]
    p3vd=venn2(subsets=tuple(map(len,p3veen)), set_labels=(scDat[0].name, scDat[1].name))
    plt.savefig(f"1F_Genes_{nfoDict['sid']}.pdf", transparent=True, dpi=300, metadata={'Title': 'Veen of Genes', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    GenesA.to_csv(f"1F_Genes_{nfoDict['sid']}_{scDat[0].name}_only.csv",encoding='utf-8')
    GenesB.to_csv(f"1F_Genes_{nfoDict['sid']}_{scDat[1].name}_only.csv",encoding='utf-8')
    GenesC.to_csv(f"1F_Genes_{nfoDict['sid']}_intersection.csv.zst",encoding='utf-8',compression={'method': 'zstd', 'level': 9, 'write_checksum': True})

if __name__ == "__main__":
    if len(sys.argv) > 1:
        thisID = sys.argv[1]
        if thisID not in SamplesDict:
            print(f"[x]sid can only be {SamplesDict.keys()}", file=sys.stderr)
            exit(1)
    else:
        thisID = 'mbrain'
    print(sys.argv, file=sys.stderr)
    print(f"[i]{thisID}")
    sys.stdout.flush()
    checkModules()
    main()

'''
x1 = np.random.randn(1000)
y1 = np.random.randn(1000)
x2 = np.random.randn(1000) * 5
y2 = np.random.randn(1000)
fig, ax = plt.subplots()
# The figure and axes background must be made transparent.
fig.patch.set(alpha=0)
ax.patch.set(alpha=0)
pc1 = ax.scatter(x1, y1, c='b', edgecolors='none')
pc2 = ax.scatter(x2, y2, c='r', edgecolors='none')
mplcairo.operator_t.ADD.patch_artist(pc2)  # Use additive blending.
plt.show()
'''
