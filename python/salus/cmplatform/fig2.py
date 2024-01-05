#!/usr/bin/env python3

from fig1 import *
#print(SamplesDict)

def main(thisID) -> None:
    scDat = {}
    nfoDict = SamplesDict[thisID]
    print("[i]Start.", file=sys.stderr)
    for platform in PlatformTuple:
        nfoDict['platformK']  = platform
        nfoDict['platformV']  = nfoDict['platforms'][platform]
        nfoDict['suffixOutV'] = nfoDict['suffixOut'][platform]
        h5Path = f"{nfoDict['sid']}_{platform}.h5ad"
        print(f"[i]Reading {h5Path}", file=sys.stderr)
        adata = ad.read_h5ad(h5Path)
        #adata.layers["raw"] = adata.X.copy()
        #adata.layers["prnorm"] = adata.X.copy()
        #sc.experimental.pp.normalize_pearson_residuals(adata,layer='prnorm')
        #sc.pp.normalize_total(adata,target_sum=1e6,key_added='CPMnormFactor')
        #adata.layers["norm"] = adata.X.copy()
        scDat[platform] = adata
        print(f"[i]Read {thisID}.{platform}: {adata.raw.shape} -> {adata.shape}", file=sys.stderr)
    #print(scDat)
    rawList=[scDat[v].raw.to_adata() for v in PlatformTuple]
    adata=ad.concat(rawList, label='Platform', keys=PlatformTuple, index_unique='-')
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    adata.obs['sqrt_inv_total_counts'] = 1 / np.sqrt(adata.obs['total_counts'])
    p995 = np.percentile(adata.obs['sqrt_inv_total_counts'].values, 99.5)
    c995 = 1/(p995*p995)
    print(f"[i] {nfoDict['sub']} sqrt_inv_total_counts: {p995} ({round(c995,3)})", file=sys.stderr)
    plt.figure(1)
    ax=sc.pl.violin(adata,['sqrt_inv_total_counts'],jitter=0.4, stripplot=True,show=False)
    ax.set_title(f"sqrt_inv_total_counts Violin - {nfoDict['sub']} @ {round(p995,9)}({round(c995,3)})")
    ax.axhline(y=p995, color='red', linestyle='dotted', label=f'p995={p995}')
    plt.savefig(f"2C_umiEstd_{nfoDict['sid']}_{round(c995)}.pdf", metadata={'Title': 'sqrt_inv_total_counts Violin', 'Subject': f"{nfoDict['sub']} Data @ {round(c995)}", 'Author': 'HU Xuesong'})
    plt.figure(2)
    ax=sc.pl.violin(adata,['total_counts'],jitter=0.4, stripplot=True,show=False)
    ax.set_title(f"total_counts Violin - {nfoDict['sub']} @ {round(c995)}")
    ax.axhline(y=round(c995), color='red', linestyle='dotted', label=f'c995={round(c995)}')
    plt.savefig(f"2C_umiEcnt_{nfoDict['sid']}_{round(c995)}.pdf", metadata={'Title': 'total_counts Violin', 'Subject': f"{nfoDict['sub']} Data @ {round(c995)}", 'Author': 'HU Xuesong'})
    adata.raw = adata.copy()
    sc.pp.filter_cells(adata, min_counts=round(c995))   # sqrt_inv_total_counts < p995 按照样品均值的标准差考虑。写作round(c995)。
    sc.pp.filter_genes(adata, min_cells=1)  # adata.var[adata.var['n_cells']<2].sort_values(by='sqrt_inv_total_counts') 有800个，就不过滤了。
    print(f"[i]Filtered: {adata.raw.shape} -> {adata.shape}", file=sys.stderr)

    #adata.layers["prnorm"] = adata.X.copy()
    #sc.experimental.pp.normalize_pearson_residuals(adata,layer='prnorm')
    sc.pp.normalize_total(adata,target_sum=1e6, key_added='CPMnormFactor')
    adata.layers["norm"] = adata.X.copy()
    sc.pp.log1p(adata)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata,random_state=369)
    sc.tl.draw_graph(adata,random_state=369)
    sc.tl.tsne(adata,random_state=369)
    sc.tl.leiden(adata,random_state=0)
    '''
    fig1, ax1 = plt.subplots()
    ax1.plot(x, y)
    ax1.set_title("Axis 1 title")
    ax1.set_xlabel("X-label for axis 1")
    '''
    print("[i]Begin fig E. 2Ca", file=sys.stderr)
    plt.figure(3)
    ax=sc.pl.pca(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}", annotate_var_explained=True)
    plt.savefig(f"2C_PCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'PCA', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(4)
    ax=sc.pl.umap(adata,color='Platform', show=False, title=f"UMAP - {nfoDict['sub']}")
    plt.savefig(f"2C_UMAP_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'UMAP', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(5)
    ax=sc.pl.tsne(adata, color='Platform', show=False, title=f"t-SNE - {nfoDict['sub']}")
    plt.savefig(f"2C_tSNE_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 't-SNE', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(6)
    ax=sc.pl.draw_graph(adata, color='Platform', show=False, title=f"ForceAtlas2 - {nfoDict['sub']}")
    plt.savefig(f"2C_ForceAtlas2_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'ForceAtlas2', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    fig, ax = plt.subplots()
    fig.patch.set(alpha=0)
    ax.patch.set(alpha=0)
    palette = ['#CC0000','#00CC00','#CCCC00']
    sc.pl.pca(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}", ax=ax, palette=palette, annotate_var_explained=True)
    arts=ax.findobj(matplotlib.collections.PathCollection)
    for art in arts:
        mplcairo.operator_t.ADD.patch_artist(art)
    newlabels = adata.obs['Platform'].unique().tolist() + ['Both']
    legend_elements = [plt.scatter([],[],linewidths=0, marker='o', label=label, color=color) for label, color in zip(newlabels, palette)]
    ax.legend(handles=legend_elements,bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0.2)
    plt.savefig(f"2C_mPCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'PCA', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    fig, ax = plt.subplots()
    fig.patch.set(alpha=0)
    ax.patch.set(alpha=0)
    palette = ['#CC0000','#00CC00','#CCCC00']
    sc.pl.tsne(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}", ax=ax, palette=palette)
    arts=ax.findobj(matplotlib.collections.PathCollection)
    for art in arts:
        mplcairo.operator_t.ADD.patch_artist(art)
    newlabels = adata.obs['Platform'].unique().tolist() + ['Both']
    legend_elements = [plt.scatter([],[],linewidths=0, marker='o', label=label, color=color) for label, color in zip(newlabels, palette)]
    ax.legend(handles=legend_elements,bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0.2)
    plt.savefig(f"2C_mtSNE_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 't-SNE', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    print("[i]Begin fig E. 2Cb", file=sys.stderr)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plt.subplots_adjust(wspace=0.1)
    sc.pl.umap(adata[adata.obs['Platform']=='Illumina'], color='leiden', ax=axes[0], title=f'UMAP - Illumina')
    sc.pl.umap(adata[adata.obs['Platform']=='Salus'], color='leiden', ax=axes[1], title=f'UMAP - Salus')
    axes[0].legend().set_visible(False)
    fig.suptitle(f"Clusters by Leiden - {nfoDict['sub']}")
    fig.savefig(f"2C_leidenUMAP_{nfoDict['sid']}.pdf", metadata={'Title': 'Cluster UMAP', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(figsize=(6,4))
    plt.title(f"Cluster Size Histogram - {nfoDict['sub']}")
    axB = sns.histplot(adata.obs,x='leiden',hue='Platform',multiple="dodge",shrink=.66)
    axB.set_xlabel('leiden Cluster NO.')
    axB.set_ylabel('Cluster Size')
    plt.savefig(f"2C_leidenHist_{nfoDict['sid']}.pdf", metadata={'Title': 'Cluster Size Histogram', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    print("[i]Begin fig E. 2D", file=sys.stderr)
    import pymn
    adata.obs['cell.cluster'] = adata.obs['leiden'].astype(str)
    adata.obs['study_id'] = adata.obs['Platform'].astype(str)
    pymn.variableGenes(adata, study_col='Platform')
    pymn.MetaNeighborUS(adata, study_col='study_id', ct_col='cell.cluster', fast_version=True)
    plt.figure(figsize=(6,4))
    plt.title(f"MetaNeighborUS - {nfoDict['sub']}")
    cm=pymn.plotMetaNeighborUS(adata, figsize=(10, 10), cmap='coolwarm',show=False)
    cm.savefig(f"2D_MetaNeighborUS_{nfoDict['sid']}.pdf", metadata={'Title': 'MetaNeighborUS', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    pymn.topHits(adata, threshold=0)
    mndf = adata.uns['MetaNeighborUS_topHits']
    mndf['ClusterID'] = mndf['Study_ID|Celltype_1'].str.split('|').str[1].astype(int)
    pndf=mndf[mndf['Match_type']=='Reciprocal_top_hit']
    plt.figure(figsize=(6,4))
    plt.title(f"Mean_AUROC Between Platforms - {nfoDict['sub']}")
    axC = sns.barplot(pndf,x='ClusterID',y='Mean_AUROC')
    axC.set_xlabel('leiden Cluster NO.')
    plt.savefig(f"2D_AUROC_{nfoDict['sid']}.pdf", metadata={'Title': 'AUROC', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    #adata = None

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
    main(thisID)

'''
./fig1.py human; ./fig2.py human ; ./fig1.py mbrain; ./fig2.py mbrain ; ./fig1.py mkidney; ./fig2.py mkidney
time (./fig1.py human; ./fig1.py mbrain ; ./fig1.py mkidney ) | tee plot.log
time (./fig2.py human; ./fig2.py mbrain ; ./fig2.py mkidney )
'''