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
        adata.layers["raw"] = adata.X.copy()
        adata.layers["prnorm"] = adata.X.copy()
        sc.experimental.pp.normalize_pearson_residuals(adata,layer='prnorm')
        sc.pp.normalize_total(adata,target_sum=1e6,key_added='CPMnormFactor')
        adata.layers["norm"] = adata.X.copy()
        scDat[platform] = adata
        print(f"[i]Read {thisID}.{platform}: {adata.raw.shape} -> {adata.shape}", file=sys.stderr)
    #print(scDat)
    rawList=[scDat[v].raw.to_adata() for v in PlatformTuple]
    adata=ad.concat(rawList, label='Platform', keys=PlatformTuple, index_unique='-')
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
    adata.raw = adata
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)
    print(f"[i]Filtered: {adata.raw.shape} -> {adata.shape}", file=sys.stderr)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata,random_state=369)
    sc.tl.draw_graph(adata,random_state=369)
    sc.tl.tsne(adata,random_state=369)
    '''
    fig1, ax1 = plt.subplots()
    ax1.plot(x, y)
    ax1.set_title("Axis 1 title")
    ax1.set_xlabel("X-label for axis 1")
    '''
    print("[i]Begin fig E. 2Ca", file=sys.stderr)
    plt.figure(1)
    ax=sc.pl.pca(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}")
    plt.savefig(f"2C_PCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'PCA', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(2)
    ax=sc.pl.umap(adata,color='Platform', show=False, title=f"UMAP - {nfoDict['sub']}")
    plt.savefig(f"2C_UMAP_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'UMAP', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(3)
    ax=sc.pl.tsne(adata, color='Platform', show=False, title=f"t-SNE - {nfoDict['sub']}")
    plt.savefig(f"2C_tsne_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 't-SNE', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(4)
    ax=sc.pl.draw_graph(adata, color='Platform', show=False, title=f"ForceAtlas2 - {nfoDict['sub']}")
    plt.savefig(f"2C_ForceAtlas2_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'ForceAtlas2', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    fig, ax = plt.subplots()
    fig.patch.set(alpha=0)
    ax.patch.set(alpha=0)
    sc.pl.pca(adata, color='Platform', show=False, title=f"PCA - {nfoDict['sub']}", ax=ax)
    arts=ax.findobj()
    for art in arts:
        mplcairo.operator_t.ADD.patch_artist(art)
    plt.savefig(f"2C_mPCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'PCA', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

    adata = None
    print("[i]Begin fig E. 2Cb", file=sys.stderr)
    for platform in PlatformTuple:
        adata = scDat[platform]

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