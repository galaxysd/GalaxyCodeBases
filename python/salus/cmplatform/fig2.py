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
    print("[i]Begin fig E. 2C", file=sys.stderr)
    plt.figure(1)
    plt.title(f"PCA - {nfoDict['sub']}")
    ax=sc.pl.pca(adata, color='Platform', show=False)
    plt.savefig(f"2C_PCA_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'PCA', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(2)
    plt.title(f"UMAP - {nfoDict['sub']}")
    ax=sc.pl.umap(adata,color='Platform', show=False)
    plt.savefig(f"2C_UMAP_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'UMAP', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(3)
    plt.title(f"t-SNE - {nfoDict['sub']}")
    ax=sc.pl.tsne(adata, color='Platform', show=False)
    plt.savefig(f"2C_tsne_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 't-SNE', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})
    plt.figure(4)
    plt.title(f"ForceAtlas2 - {nfoDict['sub']}")
    ax=sc.pl.draw_graph(adata, color='Platform', show=False)
    plt.savefig(f"2C_ForceAtlas2_{nfoDict['sid']}.pdf", bbox_extra_artists=(ax.get_legend(),), metadata={'Title': 'ForceAtlas2', 'Subject': f"{nfoDict['sub']} Data", 'Author': 'HU Xuesong'})

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