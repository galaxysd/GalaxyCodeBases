#!/usr/bin/env python3

import squidpy as sq
import os

prefix='/share/result/spatial/data/BoAo_sp/sub2/'
outpath='/share/result/spatial/data/BoAo_sp/sub2/h5ad/'
items=['srp_Illumina_mbrain', 'srp_Illumina_mkidney', 'srp_Salus_mbrain', 'srp_Salus_mkidney']

for ids in items:
    visiumPath = os.path.join( prefix,ids,'outs' )
    outfile = os.path.join( outpath,ids+'.rsp.h5ad' )
    print(f"{visiumPath}, {outfile}")
    adata=sq.read.visium(visiumPath, library_id=ids)
    adata.var_names_make_unique()
    adata.write_h5ad(outfile,compression='lzf')
