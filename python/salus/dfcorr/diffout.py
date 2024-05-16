#!/usr/bin/env python3

import pandas as pd

def qw(s):
    return tuple(s.split())

prefix = '/share/result/spatial/test_huxs/prj/cmpmethod/testicles2/c2ltest/'

diffout = pd.read_csv(f'{prefix}diffout.csv.zst', index_col=0)
diffcorr=diffout.corr()

# micromamba install zstandard
factors = qw('Elongating Endothelial InnateLymph Leydig Macrophage Myoid SPG STids Scytes Sertoli Unknown')

for onefactor in factors:
    colids = []
    colids.append(f"{onefactor}_1")
    colids.append(f"{onefactor}_2")
    onecorr = diffcorr[colids[0]][colids[1]]
    print(f'{onefactor}: {onecorr}')
