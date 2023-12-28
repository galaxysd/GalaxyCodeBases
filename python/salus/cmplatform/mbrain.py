#!/usr/bin/env python3

import importlib.metadata
from packaging import version
pkgname = "squidpy"
min_ver = "1.2.3"
got_ver = importlib.metadata.version(pkgname)
if version.parse(got_ver) < version.parse(min_ver):
    raise importlib.VersionConflict(f"{pkgname}>={min_ver} is needed, but found {pkgname}=={got_ver}")

import matplotlib; matplotlib.use("module://mplcairo.base")
from matplotlib import pyplot as plt
import mplcairo

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq
import seaborn as sns

mbi=sq.read.visium('/share/result/spatial/data/BoAo_sp/illumina/mbrain/outs/')
mbi.var_names_make_unique()
#mbic=sc.read_visium('/share/result/spatial/data/BoAo_sp/illumina/mbrain/outs/')
mbs=sq.read.visium('/share/result/spatial/data/BoAo_sp/salus/mbrain/outs/')
mbs.var_names_make_unique()

mbi.var["mt"] = mbi.var_names.str.startswith("mt-")  # mouse: mt-Nd1 mt-Nd2 mt-Co1 mt-Co2 mt-Atp8 mt-Atp6 mt-Co3 mt-Nd3 mt-Nd4l mt-Nd4 mt-Nd5 mt-Nd6 mt-Cytb
sc.pp.calculate_qc_metrics(mbi, qc_vars=["mt"], inplace=True)
mbs.var["mt"] = mbs.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(mbs, qc_vars=["mt"], inplace=True)

with pd.option_context("mode.copy_on_write", True):
    obsmbi = mbi.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
    obsmbs = mbs.obs[['n_genes_by_counts', 'total_counts']].copy(deep=False)
    p1df = pd.concat([obsmbi.assign(Platform='Illumina'), obsmbs.assign(Platform='Salus')], ignore_index=True).replace([np.inf, -np.inf, 0], np.nan).dropna()
    p2df = obsmbi.join(obsmbs,lsuffix='_Illumina',rsuffix='_Salus',how='inner').replace([np.inf, -np.inf, 0], np.nan).dropna()

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params, font="STIX Two Text")
figA=sns.JointGrid(data=p1df, x="total_counts", y="n_genes_by_counts", hue='Platform', dropna=True)
#figA.plot(sns.scatterplot, sns.histplot, alpha=.7, edgecolor=".2", linewidth=.5)
figA.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
figA.plot_marginals(sns.histplot, kde=True, alpha=.618)
figA.figure.suptitle('Gene to UMI plot - Mouse Brain')
figA.set_axis_labels(xlabel='UMIs per Barcode', ylabel='Genes per Barcode')
figA.savefig('Dmbrain.pdf', transparent=True, dpi=300, metadata={'Title': 'Gene to UMI plot', 'Subject': 'Mouse Brain Data', 'Author': 'HU Xuesong'})

figB=sns.JointGrid(data=p2df, x="total_counts_Illumina", y="total_counts_Salus", dropna=True)
figB.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
figB.plot_marginals(sns.histplot, kde=True, alpha=.618)
figB.figure.suptitle('UMI per Barcode Counts Comparing - Mouse Brain')
figB.set_axis_labels(xlabel='UMI Counts from Illumina', ylabel='UMI Counts from Salus')
figB.savefig('Embrain.pdf', transparent=True, dpi=300, metadata={'Title': 'UMI per Barcode Counts Comparing', 'Subject': 'Mouse Brain Data', 'Author': 'HU Xuesong'})


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
