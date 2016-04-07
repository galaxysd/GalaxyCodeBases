## bwa2hg18.variants	86

; cp /share/users/huxs/work/bsvir/bsI/meerkat/Meerkat/merged/simout.variants bwa2hg18.variants

````bash
scripts/pre_process.pl -s 20 -k 1500 -q 0 -t 8 -b merged/simout.bam -I hg18/HomoGRCh38 -A hg18/HomoGRCh38.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 3 -t 8 -b merged/simout.bam -F hg18/ -S samtools-0.1.20/
scripts/mechanism.pl -b merged/simout.bam -R hg18rmsk/rmsk-hg18.txt
````



## meth2gx.variants	23

; cp /share/users/huxs/work/bsvir/bsI/meerkat/Meerkat/simvir/simVir4.variants meth2gx.variants

````bash
samtools merge -p simvir/simVir4.bam simVir4_aln/Fsimout_m*.bam
scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b simvir/simVir4.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 3 -t 8 -b simvir/simVir4.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b simvir/simVir4.bam -R hg18rmsk/rmsk-hg18.txt
````



## bwa2gx.variants	8

; cp /share/users/huxs/work/bsvir/bsI/meerkat/Meerkat/mergedvir/simout.variants bwa2gx.variants

````bash
 # max_depth = 7971, 12000(?)
scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b mergedvir/simout.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 3 -t 8 -b mergedvir/simout.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b mergedvir/simout.bam -R hg18rmsk/rmsk-hg18.txt
````

## New

````bash
 #scripts/pre_process.pl -s 20 -k 1500 -q 0 -t 8 -b merged/simout.bam -I hg18/HomoGRCh38 -A hg18/HomoGRCh38.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 3 -s 20 -p 5 -o 0 -t 8 -b merged/simout.bam -F hg18/ -S samtools-0.1.20/
scripts/mechanism.pl -b merged/simout.bam -R hg18rmsk/rmsk-hg18.txt

 #scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b simvir/simVir4.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 3 -s 20 -p 5 -o 0 -t 8 -b simvir/simVir4.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b simvir/simVir4.bam -R hg18rmsk/rmsk-hg18.txt

 #scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b mergedvir/simout.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 3 -s 20 -p 5 -o 0 -t 8 -b mergedvir/simout.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b mergedvir/simout.bam -R hg18rmsk/rmsk-hg18.txt

````

## Appendex

All `.blacklist.gz` is empty.

### simout.variants	93

````bash
scripts/pre_process.pl -s 20 -k 1500 -q 0 -t 8 -b merged/simout.bam -I Ref/GX.fa -A Ref/GX.fa.fai
 # fix merged/simout.isinfo to ins=100
scripts/pre_process.pl -s 20 -k 1500 -q 0 -t 8 -b merged/simout.bam -I Ref/GX.fa -A Ref/GX.fa.fai -P cl
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 3 -t 8 -b merged/simout.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b merged/simout.bam -R hg18rmsk/rmsk-hg18.txt
````
