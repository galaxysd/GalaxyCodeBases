#!/usr/bin/env python3
# coding=utf-8

#### import Moduules
import sys, os
import math
import csv
import gzip
import itertools
#import tkinter
import zipfile
import re
#from tkinter import messagebox

import xopen
from isal import igzip, isal_zlib

import multiprocessing
from multiprocessing import Pool, cpu_count

###### Version and Date
PROG_VERSION = '0.1.0'
PROG_DATE = '2022-06-16'

###### Usage
usage = '''
     version %s
     Usage: %s <fqFile> <outputPath> <1 or 1.25> >STDOUT
''' % (PROG_VERSION,  os.path.basename(sys.argv[0]))

######  global variable  #####
p = re.compile(r'(R\d{3}(:)?C\d{3})')

# ### 解决tkinter 两个弹框的问题
# root = tkinter.Tk()
# root.withdraw()

CentRm1 = 8 - 1
CentCm1 = 66 - 1

def transCorrd(pos_x, pos_y, FovR, FovC, imageHeight, imageWidth):
    '''
    transfer the coord
    : pos_x - The original x coordinates
    : pos_y - The original y coordinates
    : FOVR - The row number of fov
    : FOVC - The col number of fov
    '''
    new_x = (FovC - CentCm1) * imageHeight - pos_y - 1
    new_y = pos_x + (FovR - CentRm1 - 1) * imageWidth
    return (math.floor(math.nextafter(100*new_x,math.copysign(math.inf,new_x)))/100, math.floor(math.nextafter(100*new_y,math.copysign(math.inf,new_y)))/100)
    #return (math.floor(100*new_x)/100, math.floor(100*new_y)/100)
    # >>> str(math.floor(78450.15*100)/100)
    #'78450.14'
    #return (new_x, new_y)

def coordTransfer(unZoomRate, fqFile, output, imageHeight, imageWidth, ratioHeight, ratioWidth, ratioHeightStart, ratioWidthStart, show = False):
    '''
    Convert Fov coordinates to chip coordinates
    '''
    fqbase = os.path.splitext(os.path.basename(fqFile))[0]
    outPutPath = os.path.join(output, "newCoord_" + fqbase + ".gz" ) ## output Fastq File
    outFastq =  xopen.xopen(outPutPath, "wb", threads=0, compresslevel=3)

    if ratioHeightStart > ratioHeight:
        print("ratioHeightStart : %f is bigger than  ratioHeight : %f" % (ratioHeightStart, ratioHeight), file=sys.stderr)
        return 0

    if ratioWidthStart > ratioWidth:
        print("ratioWidthStart : %f is bigger than  ratioWidth : %f" % (ratioWidthStart, ratioWidth), file=sys.stderr)
        return 0

    ## Get central area information
    if ratioHeightStart < 0:
        ratioHeightStart = ratioHeight / 2

    if ratioWidthStart < 0:
        ratioWidthStart = ratioWidth / 2

    x_range = [round(imageWidth * ratioWidthStart,2), round(imageWidth * (1 - ratioWidth + ratioWidthStart),2)]
    y_range = [round(imageHeight * ratioHeightStart,2), round(imageHeight * (1 - ratioHeight + ratioHeightStart),2)]
    print(f'X∈{x_range}, Y∈{y_range}. Ref:[1024,3904.34],[208.35,1951.65].', file=sys.stderr)

    #imageWidth_new = imageWidth * (1 - ratioWidth)
    #imageHeight_new = imageHeight * (1 - ratioHeight)
    imageWidth_new = 709.8 * 14.0 / 3.45
    imageHeight_new = 429.6 * 14.0 / 3.45

    if show:
        corrdRecord = []

    pre_row = 0
    #unZoomRate = 1.25
    with xopen.xopen(fqFile, threads=0) as pf:
        for line in pf:
            ### read infor for xopen.xopen
            title = line[:-1].split()[0]
            read_seq = pf.readline()[:-1]
            Links = pf.readline()
            Q_value = pf.readline()[:-1]

            pres = p.findall(title) ## get the fov num of this read

            min_x = x_range[0]
            min_y = y_range[0]
            if len(pres) != 0:
                fov = pres[0][0].replace(':','')

                if pres[0][1] == '':
                    splitSet = title.split("_")
                    THEunZoomRate = 1
                elif pres[0][1] == ':':
                    splitSet = title.split(":")
                    THEunZoomRate = unZoomRate
                else:
                    continue
                pos_y = float(splitSet[-1]) / THEunZoomRate
                pos_x = float(splitSet[-2]) / THEunZoomRate

                '''if int(fov[1:4]) == 1:
                    min_x = 0
                else:
                    min_x = x_range[0]

                if int(fov[5:8]) == 1:
                    min_y = 0
                else:
                    min_y = y_range[0]'''

                if (min_x < pos_x <= x_range[1]) and (min_y < pos_y <= y_range[1]): ##  Judge whether the coordinate is in the center of the image

                    pos_x -= min_x
                    pos_y -= min_y
                    if int(fov[1:4]) != pre_row:
                        pre_row = int(fov[1:4])
                    new_pos = transCorrd(pos_x, pos_y, int(fov[1:4]), int(fov[5:8]), imageHeight_new, imageWidth_new)
                    if show:
                        corrdRecord.append(np.array(new_pos))

                    #splitSet[-1] = str(new_pos[1])
                    #splitSet[-2] = str(new_pos[0])
                    #new_title = "_".join(splitSet)
                    new_title = f"{title} {new_pos[0]:.2f} {new_pos[1]:.2f}\037{fov}"

                    outFastq.write((new_title+ '\n' + read_seq + '\n' + Links + Q_value + '\n').encode(encoding="utf-8"))

    outFastq.close()

    if show:
        import numpy as np
        import matplotlib.pyplot as plt
        corrdRecord = np.array(corrdRecord)
        corrdRecord = np.round(corrdRecord / 100).astype(int)
        corrdRecord = corrdRecord.T

        max_x = np.max(corrdRecord[0]) + 1
        max_y = np.max(corrdRecord[1]) + 1
        min_x = np.min(corrdRecord[0])
        min_y = np.min(corrdRecord[1])
        rangeYX = (max_y - min_y, max_x - min_x)
        print(f'Y:[{min_y},{max_y}),X:[{min_x},{max_x}) -> YX:{rangeYX}', file=sys.stderr)

        corrdRecord = corrdRecord.T

        heatmap = np.zeros(rangeYX)

        for readNum in range(len(corrdRecord)):
            coordTemp = corrdRecord[readNum]
            heatmap[coordTemp[1]-min_y][coordTemp[0]-min_x] += 1
        plt.matshow(heatmap)
        outTiffPath = os.path.join(output, "newCoord_" + fqbase + '.heatMap.tif' )
        plt.savefig(outTiffPath)
        #cv2.imwrite("bins.tif", np.array(heatmap).astype(np.uint16))

def main():
    ######################### Phrase parameters #########################
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage)
    ArgParser.add_argument("--version", action="version", version=PROG_VERSION)
    ArgParser.add_argument("-r1", "--ratioWidth", action="store", dest="ratioWidth", type=float, default=0.0, metavar="FLOAT", help="Overlap ratio of width. [%(default)s]")
    ArgParser.add_argument("-r1s", "--ratioWidthStart", action="store", dest="ratioWidthStart", type=float, default=-1.0, metavar="FLOAT", help="start Overlap ratio of width. [%(default)s]")
    ArgParser.add_argument("-r2", "--ratioHeight", action="store", dest="ratioHeight", type=float, default=0.0, metavar="FLOAT", help="Overlap ratio of heigh [%(default)s]")
    ArgParser.add_argument("-r2s", "--ratioHeightStart", action="store", dest="ratioHeightStart", type=float, default=-1.0, metavar="FLOAT", help="start Overlap ratio of heigh [%(default)s]")
    ArgParser.add_argument("-s", "--showHeatmap", action="store_true", dest="showHeatmap", default=False, help="Display heatmap. [%(default)s]")

    (para, args) = ArgParser.parse_known_args()

    if len(args) < 3:
        ArgParser.print_help()
        print("\nERROR: The parameters number is not correct!", file=sys.stderr)
        sys.exit(1)
    else:
        (fqFile, output, unZoomRate) = args

    imageHight = 2160;
    imageWidth = 4096;
    para.ratioWidth = 0.296791;
    para.ratioHeight = 0.192915;
    para.ratioWidthStart = 0.25;
    coordTransfer(float(unZoomRate), fqFile, output, int(imageHight), int(imageWidth), float(para.ratioHeight), float(para.ratioWidth), float(para.ratioHeightStart), float(para.ratioWidthStart), bool(para.showHeatmap))

if __name__ == "__main__":
    main()  # 2160 4096 -r1 0.3 -r2 0.2 -r1s 0.25; 2160 4096 -r1 0.296791 -r2 0.192915 -r1s 0.25
