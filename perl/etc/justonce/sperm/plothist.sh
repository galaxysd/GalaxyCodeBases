#!/bin/bash

# argument 1 is path to input
# argument 2 is range
arrFN=(${1//[-.]/ })
MAIN="${arrFN[0]}"
MAIN=${MAIN##*/}
echo -e "[$1] -> $MAIN ${arrFN[0]} ${arrFN[1]} ${arrFN[2]}"

gnuplot << PLOTCMD
set term png font "/bak/seqdata/sperm/sssC2/arial.ttf" 48 size 2102,2102 notruecolor linewidth 3
set output "$1.png"
set title "KmerFreq for [$1]"

Min = 1 # where binning starts
Max = 201 # where binning ends
n = 50 # the number of bins
width = (Max-Min)/n # binwidth; evaluates to 1.0
bin(x) = width*(floor( (x-Min)/width)+0.5 ) + Min

set xlabel "K-mer count"
set ylabel "Ratio"
set style fill solid 0.6 border

set xrange [0:60]
set yrange [0:0.1]
plot '<(head -n67 $1 |tail -n59)' using 1:3 with boxes title "${MAIN}" lc rgb 'black'
PLOTCMD
# plot '$1' using (bin(\$1)):(1.0) smooth freq with boxes

echo << EOF
sh plothist.sh Donor.k25.lz4.hist
eg. find hist/*.hist |while read a;do sh plothist.sh $a;done
EOF
