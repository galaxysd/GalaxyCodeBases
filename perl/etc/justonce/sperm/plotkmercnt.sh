#!/bin/bash

# argument 1 is path to input
# argument 2 is range

arrFN=(${1//[-.]/ })
if [[ X"" = X"${arrFN[2]}" ]]; then
	arrFN[1]="File 1";
	arrFN[2]="File 2";
fi
echo -e "[$1] -> ${arrFN[0]} ${arrFN[1]} ${arrFN[2]}\ncb_max: $2"

gnuplot << PLOTCMD
set term png font "/bak/seqdata/sperm/sssC2/arial.ttf" 48 size 2102,2102 truecolor linewidth 3
#set terminal png
set output "$1.png"
#set size ratio 0.5
set title "Heatmap for [$1]"

#set xlabel "File 2"
#set ylabel "File 1"
set ylabel "${arrFN[1]}"
set xlabel "${arrFN[2]}"

set tic scale 0

set palette rgbformulae 22,13,10 #color
set palette rgbformulae 3,3,3 #gray
#set palette rgbformula 34,35,36	#heat
#set palette negative #just for fun

#set cbrange [0:10000]
set cbrange [0:"$2"]
#unset cbtics

set tics nomirror
set xrange [-0.5:40.5]
set yrange [-0.5:40.5]
set xtics ("0" 0, "0.2" 8, "0.4" 16, "0.6" 24, "0.8" 32, "0.9" 36, ">1" 40)
set ytics ("0" 0, "0.2" 8, "0.4" 16, "0.6" 24, "0.8" 32, "0.9" 36, ">1" 40)
set tics out
set tics scale 0.7

#set cbrange [1:*]
set format cb "%3.0e"
#set cbtics in
set cbtics font "/bak/seqdata/sperm/sssC2/arial.ttf, 40"

set view map

splot '$1' matrix with image
PLOTCMD

echo << XEOF
sh plotkmercnt.sh t100k.40.plot 10000
eg. find out/*.plot |while read a;do sh plotkmercnt.sh $a 100000;done
XEOF
