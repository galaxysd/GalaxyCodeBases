#!/usr/bin/env python3

import sys
import os
import tifffile
import numpy as np
from xopen import xopen
import fast_matrix_market as fmm
import scipy.sparse
import matplotlib as mpl
import matplotlib.cm
import matplotlib.pyplot as plt
import math
import logging
logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)-8s|%(message)s\t|%(module)s:%(funcName)s:%(lineno)d|%(process)d:%(thread)d|%(pathname)s',
    datefmt='%Y-%m-%dT%H:%M:%S',level=logging.INFO)
logger = logging.getLogger(__name__)

dpi = (7112000, 69)

# https://github.com/matplotlib/matplotlib/blob/a254b687df97cda8c6affa37a1dfcf213f8e6c3a/lib/matplotlib/axes/_axes.py#L4919-L5317
def hexbin_preprocess(x,y,gridsize=100,xscale='linear', yscale='linear', extent=None):
    import math
    import numpy as np
    import matplotlib.transforms as mtransforms
    # Set the size of the hexagon grid
    if np.iterable(gridsize):
        nx, ny = gridsize
    else:
        nx = gridsize
        ny = int(nx / math.sqrt(3))
    # Count the number of data in each hexagon
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    # Will be log()'d if necessary, and then rescaled.
    tx = x
    ty = y
    if xscale == 'log':
        if np.any(x <= 0.0):
            raise ValueError("x contains non-positive values, so cannot be log-scaled")
        tx = np.log10(tx)
    if yscale == 'log':
        if np.any(y <= 0.0):
            raise ValueError("y contains non-positive values, so cannot be log-scaled")
        ty = np.log10(ty)
    if extent is not None:
        xmin, xmax, ymin, ymax = extent
        if xmin > xmax:
            raise ValueError("In extent, xmax must be greater than xmin")
        if ymin > ymax:
            raise ValueError("In extent, ymax must be greater than ymin")
    else:
        xmin, xmax = (tx.min(), tx.max()) if len(x) else (0, 1)
        ymin, ymax = (ty.min(), ty.max()) if len(y) else (0, 1)
        # to avoid issues with singular data, expand the min/max pairs
        xmin, xmax = mtransforms.nonsingular(xmin, xmax, expander=0.1)
        ymin, ymax = mtransforms.nonsingular(ymin, ymax, expander=0.1)
    nx1 = nx + 1
    ny1 = ny + 1
    nx2 = nx
    ny2 = ny
    n = nx1 * ny1 + nx2 * ny2
    # In the x-direction, the hexagons exactly cover the region from xmin to xmax. Need some padding to avoid roundoff errors.
    padding = 1.e-9 * (xmax - xmin)
    xmin -= padding
    xmax += padding
    sx = (xmax - xmin) / nx
    sy = (ymax - ymin) / ny
    return padding, nx, ny, sx, sy, xmin, ymin

def hex_dither(arr2d, rows, cols, rh=5/12, rv=1/12):
    arr2d[rows, cols] = (rh * arr2d[rows-1, cols] + rh * arr2d[rows+1, cols] +
                         rv * arr2d[rows, cols-1] + rv * arr2d[rows, cols+1])
    return

def hexbinnd(x, y, ax, C=None, gridsize=100, bins=None, xscale='linear', yscale='linear', extent=None, **kwargs):
    import numpy as np
    hb = ax.hexbin(x, y, C, gridsize, bins, xscale, yscale, extent, **kwargs)
    hb_offsets = hb.get_offsets()
    hb_array = hb.get_array()   # masked_array

    padding, nx, ny, sx, sy, xmin, ymin = hexbin_preprocess(x,y,gridsize,xscale, yscale, extent)
    nx1 = nx + 1
    ny1 = ny + 1
    nx2 = nx
    ny2 = ny
    n = nx1 * ny1 + nx2 * ny2

    hb_newXY = hb_offsets.copy()
    hb_newXY[:, 0] -= xmin
    hb_newXY[:, 1] -= ymin
    hb_newXY[:, 0] /= sx
    hb_newXY[:, 1] /= sy
    hb_intXY = np.rint(hb_newXY * 2).astype(int)
    hb_img = np.zeros((1 + 2*nx, 1 + 2*ny), dtype=int)
    hb_img[hb_intXY[:, 0], hb_intXY[:, 1]] = hb_array
    #even_rows, odd_cols = np.ogrid[  :1+2*nx:2, 1:1+2*ny:2 ]   # 1:1+2*ny:2 = 1:2*ny:2
    #odd_rows, even_cols = np.ogrid[ 1:1+2*nx:2,  :1+2*ny:2 ]   # 1:1+2*nx:2 = 1:2*nx:2
    inner_even_rows, odd_cols = np.ogrid[ 2:2*nx:2, 1:2*ny:2 ]
    odd_rows, inner_even_cols = np.ogrid[ 1:2*nx:2, 2:2*ny:2 ]
    border_even_rows, border_even_cols = np.ogrid[ :1+2*nx:2*nx, :1+2*ny:2*ny ]
    hb_img[border_even_rows, border_even_cols] *=2  # Four corners are doubled, all borders are counts from 1/2 hexbin now.
    hb_pad = np.pad(hb_img, pad_width=1, mode='constant', constant_values=0).astype(float)
    hex_dither(hb_pad, border_even_rows+1, odd_cols+1)
    hex_dither(hb_pad, odd_rows+1, border_even_cols+1)
    hb_ret = hb_pad[1:-1, 1:-1]
    hb_ret[border_even_rows, 1:-1] *=2   # Four corners only need to be doubled once, thus 0 and -1 are excluded.
    hb_ret[:, border_even_cols] *=2      # Four borders are all count from 1/1 hexbin now.
    hex_dither(hb_ret, inner_even_rows, odd_cols)
    hex_dither(hb_ret, odd_rows, inner_even_cols)
    hb_ret = np.rint(hb_ret).astype(int)
    return hb_ret

def dohexbin(coomtx, bin=40):
    import math
    #dpis3 = (142240000000, 2390229)
    from fractions import Fraction
    accuracy = 400000   # Full Chip (= 2 lanes) is (40141, 224640)
    sqrt3z = int(round(math.sqrt(3),int(math.log10(accuracy)))*accuracy)
    real = math.sqrt(3)*accuracy
    diff = (real - sqrt3z)/sqrt3z
    Fraction(14*25400*100*accuracy, 345*sqrt3z).as_integer_ratio(),binsize, real, diff

    hexhorizontal = bin * math.sqrt(2) / math.sqrt(math.sqrt(3))
    hexverticalhalf = hexhorizontal / math.sqrt(3)
    gs0 = (coomtx.shape[1]/hexhorizon, coomtx.shape[0]/hexverticalhalf)
    gs = np.rint(gs0).astype(int)
    fig, ax = plt.subplots()
    myhbdat = hexbinnd(coomtx.col,coomtx.row,gridsize=gs,ax=ax,cmap='turbo')
    
    return

def main() -> None:
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <image.mtx.gz> <outprefix>', file=sys.stderr, flush=True)
        exit(0);
    elif len(sys.argv) >= 3:
        inmtx = sys.argv[1]
        outprefix = sys.argv[2]
        logger.info(f'[{inmtx}] => [{outprefix}].{hex,rect}bin.tif')
    logger.info(f'Written:[{outprefix}.fstBC.tif].')

if __name__ == "__main__":
    main()
