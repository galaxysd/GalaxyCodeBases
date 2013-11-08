#!/bin/bash

# argument 1 is path to input
# argument 2 is range

gnuplot << PLOTCMD
set term png font "/bak/seqdata/sperm/sssC2/arial.ttf" 48 size 2102,2102 truecolor linewidth 3
#set terminal png
set output "$1.png"
#set size ratio 0.5
set title "Heatmap for [$1]"

set xlabel "File 2"
set ylabel "File 1"

set tic scale 0

#set palette rgbformulae 22,13,10
set palette rgbformula 34,35,36
#set palette negative

#set cbrange [0:10000]
set cbrange [0:"$2"]
#unset cbtics

#set xrange [-1:41]
#set yrange [-1:41]

set view map

splot '$1' matrix with image
PLOTCMD

echo << EOF
sh plotkmercnt.sh t100k.40.plot 10000
eg. find out/*.plot |while read a;do sh plotkmercnt.sh $a 100000;done
EOF
