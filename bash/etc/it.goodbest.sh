#!/bin/bash

for((j=1000;j<5000;j++))
do
 wget "https://itunes.apple.com/WebObjects/MZStore.woa/wa/viewGrouping?cc=jp&id=$j" --user-agent="iTunes/11.2.1 (Macintosh; OS X 10.9.2) AppleWebKit/537.74.9" -O $j.out
 sleep 5
done

