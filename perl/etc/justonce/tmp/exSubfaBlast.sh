#!/bin/bash

samtools faidx Obimaculoides_280.fa $1 >$1.fa
blastn -task blastn -query ir813.fa -subject $1.fa -word_size 4 -outfmt 3 >$1.b3.txt
zcat annotation/Obimaculoides_280.gene.gff3.gz|egrep "^$1\s" >$1.gff

# cat matches.lst|while read a;do ./t.sh $a;done
# rar a -m5 -rr1% -mct -mcc -mt1 res.rar ir813.fa ir813.e10w4.txt matches.lst Scaffold* ir813.e1000w4.txt
# $ cat matches.lst
# Scaffold17571
# Scaffold73065
# Scaffold54938

