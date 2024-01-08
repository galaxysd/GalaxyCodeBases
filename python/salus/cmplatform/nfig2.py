#!/usr/bin/env python3

def checkModules() -> None:
    import importlib.metadata
    from packaging import version
    pkgname = "squidpy"
    min_ver = "1.2.3"
    got_ver = importlib.metadata.version(pkgname)
    if version.parse(got_ver) < version.parse(min_ver):
        raise Exception(f"{pkgname}>={min_ver} is needed, but found {pkgname}=={got_ver}")

if __name__ == "__main__":
    checkModules()

from nfig1 import *
#print(SamplesDict)

def main(thisID) -> None:
    import squidpy as sq
    scDat = {}
    nfoDict = SamplesDict[thisID]
    print("[i]Start.", file=sys.stderr)
    for platform in PlatformTuple:
        nfoDict['platformK']  = platform
        nfoDict['platformV']  = nfoDict['platforms'][platform]
        nfoDict['suffixOutV'] = nfoDict['suffixOut'][platform]
        if nfoDict['type'] == 'mobivision':
            h5Path = f"{nfoDict['sid']}_{platform}.raw.h5ad"
            if os.path.exists(h5Path) and os.access(h5Path, os.R_OK) and os.path.getsize(h5Path) > 0:
                print(f"[i]Reading {h5Path}", file=sys.stderr)
                adata = ad.read_h5ad(h5Path)
            else:
                mtxPath = os.path.join( *[nfoDict[v] for v in nfoDict['pattern']] )
                print(f"[i]Reading {mtxPath}", file=sys.stderr)
                adata=sc.read_10x_mtx(mtxPath, var_names='gene_symbols', make_unique=True, gex_only=True)
                print(f"[i]Saving {h5Path}", file=sys.stderr)
                adata.write_h5ad(h5Path,compression='lzf')
        elif nfoDict['type'] == 'visium':
            h5Path = f"{nfoDict['sid']}_{platform}.rsp.h5ad"
            if os.path.exists(h5Path) and os.access(h5Path, os.R_OK) and os.path.getsize(h5Path) > 0:
                print(f"[i]Reading {h5Path}", file=sys.stderr)
                adata = ad.read_h5ad(h5Path)
            else:
                visiumPath = os.path.join( *[nfoDict[v] for v in nfoDict['pattern'][:-1] ] )
                print(f"[i]Reading {visiumPath}", file=sys.stderr)
                adata=sq.read.visium(visiumPath, library_id=platform)
                adata.var_names_make_unique()
                print(f"[i]Saving {h5Path}", file=sys.stderr)
                adata.write_h5ad(h5Path,compression='lzf')
        else:
            print(f"[x]Unknow Type {nfoDict['type']}", file=sys.stderr)
            exit(1)
        scDat[platform] = adata
        print(f"[i]Read {thisID}.{platform}.raw: {adata.shape}", file=sys.stderr)
    metapdf={'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'}
    rawList=[scDat[v] for v in PlatformTuple]
    adata=ad.concat(rawList, label='Platform', keys=PlatformTuple, index_unique='-')
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    adata.obs['sqrt_inv_total_counts'] = 1 / np.sqrt(adata.obs['total_counts'])
    p995 = np.percentile(adata.obs['sqrt_inv_total_counts'].values, nfoDict['fltPct'])
    c995 = 1/(p995*p995)
    print(f"[i] {nfoDict['sub']} sqrt_inv_total_counts: {p995} ({round(c995,3)})", file=sys.stderr)
    plt.figure(1)
    ax=sc.pl.violin(adata,['sqrt_inv_total_counts'],jitter=0.4, stripplot=True,show=False)
    ax.set_title(f"sqrt_inv_total_counts Violin - {nfoDict['sub']} @ {round(p995,9)}({round(c995,3)})")
    ax.axhline(y=p995, color='red', linestyle='dotted', label=f'p995={p995}')
    plt.savefig(f"2C_umiEstd_{nfoDict['sid']}_{round(c995)}.pdf", metadata={**metapdf, 'Title': 'sqrt_inv_total_counts Violin', 'Subject': f"{nfoDict['sub']} Data @ {round(c995)}"})
    plt.figure(2)
    ax=sc.pl.violin(adata,['total_counts'],jitter=0.4, stripplot=True,show=False)
    ax.set_title(f"total_counts Violin - {nfoDict['sub']} @ {round(c995)}")
    ax.axhline(y=round(c995), color='red', linestyle='dotted', label=f'c995={round(c995)}')
    plt.savefig(f"2C_umiEcnt_{nfoDict['sid']}_{round(c995)}.pdf", metadata={**metapdf, 'Title': 'total_counts Violin', 'Subject': f"{nfoDict['sub']} Data @ {round(c995)}"})
    adata.raw = adata.copy()
    sc.pp.filter_cells(adata, min_counts=round(c995))   # sqrt_inv_total_counts < p995 按照样品均值的标准差考虑。写作round(c995)。
    sc.pp.filter_genes(adata, min_cells=1)  # adata.var[adata.var['n_cells']<2].sort_values(by='sqrt_inv_total_counts') 有800个，就不过滤了。
    print(f"[i]Filtered: {adata.raw.shape} -> {adata.shape}", file=sys.stderr)
    plt.close('all')
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata,target_sum=1e6, key_added='CPMnormFactor')
    adata.layers["norm"] = adata.X.copy()
    sc.pp.log1p(adata)
    if nfoDict['type'] == 'mobivision':
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.tsne(adata,random_state=369)
    elif nfoDict['type'] == 'visium':
        import torch
        import STAGATE_pyG
        STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
        STAGATE_pyG.Stats_Spatial_Net(adata)
        adata = STAGATE_pyG.train_STAGATE(adata, device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu'))
        sc.pp.neighbors(adata, use_rep='STAGATE')
        sc.tl.tsne(adata, use_rep='STAGATE', random_state=369)
        sc.pp.pca(adata)
        data=adata.obsm['STAGATE'][0:,0:2]
        df=pd.DataFrame(data=data[0:,0:], index=[adata.obs_names[i] for i in range(data.shape[0])], columns=['STA'+str(1+i) for i in range(data.shape[1])])
        df['Platform'] = adata.obs['Platform']
        figD=sns.JointGrid(data=df, x="STA1", y="STA2", hue='Platform', dropna=True)
        figD.plot_joint(sns.scatterplot, s=12.7, alpha=.6)
        figD.plot_marginals(sns.histplot, kde=True, alpha=.618)
        figD.figure.suptitle(f"STAGATE - {nfoDict['sub']}")
        figD.savefig(f"2C_STAGATE_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'STAGATE'})
    else:
        print(f"[x]Unknow Type {nfoDict['type']}", file=sys.stderr)
        exit(1)
    sc.tl.umap(adata,random_state=369)
    sc.tl.draw_graph(adata,random_state=369)
    sc.tl.leiden(adata,random_state=0)
    '''
    fig1, ax1 = plt.subplots()
    ax1.plot(x, y)
    ax1.set_title("Axis 1 title")
    ax1.set_xlabel("X-label for axis 1")
    '''
    print("[i]Begin fig E. 2Ca", file=sys.stderr)
    plt.figure()
    ax=sc.pl.pca(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}", annotate_var_explained=True)
    plt.savefig(f"2C_PCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 'PCA'})
    plt.figure()
    ax=sc.pl.umap(adata,color='Platform', show=False, title=f"UMAP - {nfoDict['sub']}")
    plt.savefig(f"2C_UMAP_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 'UMAP'})
    plt.figure()
    ax=sc.pl.tsne(adata, color='Platform', show=False, title=f"t-SNE - {nfoDict['sub']}")
    plt.savefig(f"2C_tSNE_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 't-SNE'})
    plt.figure()
    ax=sc.pl.draw_graph(adata, color='Platform', show=False, title=f"ForceAtlas2 - {nfoDict['sub']}")
    plt.savefig(f"2C_ForceAtlas2_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 'ForceAtlas2'})
    plt.close('all')

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
    plt.savefig(f"2C_mPCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 'PCA'})
    plt.close('all')

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
    plt.savefig(f"2C_mtSNE_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={**metapdf, 'Title': 't-SNE'})
    plt.close('all')

    print("[i]Begin fig E. 2Cb", file=sys.stderr)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plt.subplots_adjust(wspace=0.1)
    sc.pl.umap(adata[adata.obs['Platform']=='Illumina'], color='leiden', ax=axes[0], title=f'UMAP - Illumina')
    sc.pl.umap(adata[adata.obs['Platform']=='Salus'], color='leiden', ax=axes[1], title=f'UMAP - Salus')
    axes[0].legend().set_visible(False)
    fig.suptitle(f"Clusters by Leiden - {nfoDict['sub']}")
    fig.savefig(f"2C_leidenUMAP_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'Cluster UMAP'})
    plt.figure(figsize=(6,4))
    plt.title(f"Cluster Size Histogram - {nfoDict['sub']}")
    axB = sns.histplot(adata.obs,x='leiden',hue='Platform',multiple="dodge",shrink=.66)
    axB.set_xlabel('leiden Cluster NO.')
    axB.set_ylabel('Cluster Size')
    plt.savefig(f"2C_leidenHist_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'Cluster Size Histogram'})
    plt.close('all')

    print("[i]Begin fig E. 2D", file=sys.stderr)
    plt.switch_backend('pdf')
    import pymn
    adata.obs['cell.cluster'] = adata.obs['leiden'].astype(str)
    adata.obs['study_id'] = adata.obs['Platform'].astype(str)
    pymn.variableGenes(adata, study_col='Platform')
    pymn.MetaNeighborUS(adata, study_col='study_id', ct_col='cell.cluster', fast_version=True)
    plt.figure(figsize=(16,16))
    plt.title(f"MetaNeighborUS - {nfoDict['sub']}")
    cm=pymn.plotMetaNeighborUS(adata, cmap='coolwarm',show=False)
    cm.savefig(f"2D_MetaNeighborUS_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'MetaNeighborUS'})
    pymn.topHits(adata, threshold=0)
    mndf = adata.uns['MetaNeighborUS_topHits']
    mndf['ClusterID'] = mndf['Study_ID|Celltype_1'].str.split('|').str[1].astype(int)
    pndf=mndf[mndf['Match_type']=='Reciprocal_top_hit']
    plt.figure(figsize=(6,4))
    plt.title(f"Mean_AUROC Between Platforms - {nfoDict['sub']}")
    axC = sns.barplot(pndf,x='ClusterID',y='Mean_AUROC')
    axC.set_xlabel('leiden Cluster NO.')
    plt.savefig(f"2D_AUROC_{nfoDict['sid']}.pdf", metadata={**metapdf, 'Title': 'AUROC'})
    plt.close('all')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        thisID = sys.argv[1]
        if thisID not in SamplesDict:
            print(f"[x]sid can only be {SamplesDict.keys()}", file=sys.stderr)
            exit(1)
    print(sys.argv, file=sys.stderr)
    print(f"[i]{thisID}")
    sys.stdout.flush()
    main(thisID)

'''
./nfig1.py human; ./nfig2.py human ; ./nfig1.py mbrain; ./nfig2.py mbrain ; ./nfig1.py mkidney; ./nfig2.py mkidney
time (./nfig1.py human; ./nfig1.py mbrain ; ./nfig1.py mkidney ) | tee nplot.log
time (./nfig2.py human; ./nfig2.py mbrain ; ./nfig2.py mkidney )
'''