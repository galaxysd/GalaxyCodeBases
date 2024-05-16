#!/usr/bin/env python3

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
import seaborn as sns

def qw(s):
    return tuple(s.split())

prefix = '/share/result/spatial/test_huxs/prj/cmpmethod/testicles2/c2ltest/'

diffout = pd.read_csv(f'{prefix}diffout.csv.zst', index_col=0)
diffcorr=diffout.corr()

# micromamba install zstandard
factors = qw('Elongating Endothelial InnateLymph Leydig Macrophage Myoid SPG STids Scytes Sertoli Unknown')

pdf = PdfPages(f'{prefix}diffout.pdf')
for onefactor in factors:
    colids = []
    colids.append(f"{onefactor}_1")
    colids.append(f"{onefactor}_2")
    onecorr = diffcorr[colids[0]][colids[1]]
    print(f'{onefactor}: {onecorr}')
    tmpDF=pd.DataFrame()
    tmpDF['rank1']=diffout[colids[0]].rank()
    tmpDF['rank2']=diffout[colids[1]].rank()
    tmpDF[colids[0]]=diffout[colids[0]]
    tmpDF[colids[1]]=diffout[colids[1]]
    plt.figure()
    sns.jointplot(data=tmpDF, x=colids[0], y=colids[1],kind="scatter",marginal_ticks=True,s=1)
    fig.suptitle(f'Corr: {onecorr}')
    pdf.savefig()
    plt.close()
pdf.close()
print(f'.\n[i]Done.', flush=True)
