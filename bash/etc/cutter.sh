#!/bin/bash

pieceHeight=5408
pieceSlide=$[$pieceHeight-130]
targetdir=./work
infile="$1"
inputname=${infile%%.*}

mkdir -p ${targetdir}

function getwidth() { magick identify -format '%[height]' "$1"; }
height=$(getwidth "${infile}")

echo "[${height}]"

top="0"
fno="0"
while [ $top -lt $height ]
do
echo "$top $height"
convert -crop x${pieceHeight}+0+${top} +repage "${infile}" "${targetdir}/${inputname}_${fno}.png"
top=$[$top+$pieceSlide]
fno=$[$fno+1]
done
