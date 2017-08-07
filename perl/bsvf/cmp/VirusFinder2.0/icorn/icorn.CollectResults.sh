#!/bin/bash

base=$1;
reads=$2;


if [ -f "$base.10.General.stats" ] 
then
	paste $base.?.General.stats $base.??.General.stats > Stats.2.Correction.csv;
	paste $ICORN_HOME/column.csv $base?/$reads*.txt $base??/$reads*.txt> Stats.Mapping.csv
else
	paste $base.?.General.stats > Stats.2.Correction.csv;
	paste $ICORN_HOME/column.csv $base?/$reads*.txt    > Stats.Mapping.csv
fi
awk 'BEGIN {sum=0;what="notset"}{ if (! ($2 ~ /[0-9]/)) {print what" "sum;what=$1;sum=0}else{sum+=$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36+$38+$40} } END {print what" "sum}' Stats.2.Correction.csv

awk '{ print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50 }' Stats.2.Correction.csv > Stats.Correction.csv

