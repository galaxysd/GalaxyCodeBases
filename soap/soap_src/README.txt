SOAP: Short Oligonucleotide Alignment Program

Current release: 1.0
2007-11-10

Author: Ruiqiang Li 
lirq@genomics.org.cn 
Bioinformatics Department
Beijing Genomics Institute, Beijing 101300, China

lirq@bmb.sdu.dk 
Department of Biochemistry and Molecular Biology 
University of Southern Denmark Campusvej 55, DK-5230 Odense M, Denmark

1. Introduction:

SOAP is an efficient program for short oligonucleotide sequence alignment and mapping. It's specially designed for the new generation DNA sequencing technology (Illumina-Solexa; ABI-SOLiD), and provides plenty of options for different applications.

a) Single-end sequencing
It will allow either a certain number of mismatches or one gap. Normally, only the best hit which have minimal number of mismatches or smaller gap will be reported. For gapped alignment, both small insertion and deletion will be identified. But to avoid ambiguous alignment for very short sequences, only one continuous gap with size ranged from 1 to 3 bp is accepted. As the connatural deficiency of the sequencing technology, errors will be accumulated, reads always have much higher sequencing errors at the 3'-end. SOAP can try to iteratively trim several basepairs at 3'-end and redo alignment if no hits were detected.

b) Pair-end sequencing
Pair-end sequencing will be able to significantly improve the accuracy of resequencing mapping. SOAP will align both ends of a pair simultaneously. A pair of reads will be aligned when two reads are mapped on the right relative orientation relationship with proper distance. A certain number of mismatches are allowed in either or both reads of a pair. If no hit exist, will try to redo alignment for trimmed reads, or try gapped alignment, depending on the options user choosed. For gapped alignment, gap can only appear on one end, and the other end should be identically aligned.

c) mRNA tag sequencing
On SAGE sequencing, there are two types of enzymes: 1)DpnII, which will recognize site "GATC" and cut 16bp after the site; 2)NlaIII, which will specially recognize site "CATG" and cut 17bp afterwards. SOAP will trim the 3'-end of read up to the proper tag length before alignment. The matched hits are required to have identical enzyme site and a certain number of mismatch in the tag region.

d) Small RNA sequencing
The new generation sequencing technology offers a powerful solution to discover small RNA genes. Small RNAs always have a size between 18 to 26bp. According to the protocol, the 3'-end of RNA sequence will be flanked by adapter sequence. SOAP will try to filter candidate adapter sequence at the 3'-end, then align the remaining sequence. Considering sequencing errors, one or two mismatches are allowed in either adapter sequence or the candidate RNA region, according to user settings.

In principle, the program works for any length of queries, but for long queries (>=60bp) or low similarity search, we suggest standard Smith-Waterman alignment, or using approximate program like Blast to have proper sensitivity.

2. Download and Installation:

SOAP is distributed under GNU Public License (GPL). The source code is freely available from http://soap.genomics.org.cn. If you find bugs or have constructive suggestions to the program, please feel free to send e-mail to the author.

SOAP is written in standard C++ language. It could be compiled and run under any type of linux or unix environment. For big reference sequences like human genome, the program has to be compiled into 64 bits, while 32 bits compiling will also works pretty well for small reference. You can simply type "make" under the source directory to compile the program.

The RAM required for the program can be roughly calculated as:
a) 64 bits
RAM=L/3+(4*3+8*6)*(4^S)+(4+1)*3*L/4+4*(2^24).

Where L is the total length of reference sequences; S is seed size. So for small reference like yeast, L=12Mb, select seed size S=10, about 200Mb RAM is needed; but for human genome, L=3Gb, select seed size S=12, about 15Gb will be needed.

b) 32 bits
RAM=L/3+(4*3+4*6)*(4^S)+(4+1)*3*L/4+4*(2^24).

The program supports multithreaded parallel computing, be sure to add "-DTHREAD" in the compiling options if you want to use the parallel version.

Compiling options:

1) "-DREAD_36":	Maximum length of reads on the query file, There are three choices now: -DREAD_36, -DREAD_48, and -DREAD_60. Smaller one will runs faster, but larger ones are compatible to smaller one.
2) "-DMAXHITS":	For extremely high frequent repeat hits, the program will stop after getting 100 hits, so in the output, 100 means there are too many equal best hits. User can set their own threshold as required, by changing the setting of "-DMAXHITS" in makefile and throughly recompile the program.
3) "-DMAXGAP":	Maximum size of gap, then "-g" option during running should not exceed this definition.
4) "-DTHREAD":	Multi-thread compiling.
5) "-DDB_CHR":	Type of reference. There are 4 choice now:
	-DDB_CHR	# of sequences <256, length of each <16Gb;
	-DDB_CONTIG	# of sequences <65,536, length of each <16Gb;
	-DDB_SHORT	# of sequences <4G, length of each <1Kb;
	-DDB_HUGE	# of sequences <4G, length of each <16Gb.
	Set proper definition will save RAM while fit your job better.

3. Usage:
The program can be run either for a single query file, or for batch query files.
For single query, you should type command as:
>soap [options]

For batch queries, you should type command as:
>soap [options] -d <ref.fa> <param.arg>

Where param.arg is a file storing the queue of commands, one command on a line.

a) Common options:
        -a      <string>        query a file, *.fq or *.fa format
        -d      <string>        reference sequences file, *.fa format
        -o      <string>        output align file
        -s      <int>   seed size, default=10, seed size*2+3 should be smaller than the size of minimal read allowed in the query file, so that will not lose sensitivity, considering the RAM, please do not make seed size exceed 12bp.
        -v      <int>   maximum number of snps allowed on a read, <=2. default=2bp. If you only want to detect identical hits or hits with at most 1 mismatch, please set the exact value, so that the program can be significantly speed up.
        -g      <int>   maximum gap size allowed on a read, default=0bp. It means that we will not try gapped alignment, if you set value 0; Or gapped alignment will be performed if no hits exist by ungapped alignment. 
        -e      <int>   will not allow gap exist inside n-bp edge of a read, default=3bp
        -q      <int>   quality threshold of snp, 0-40, default=20
        -z      <char>  initial quality, default='!'
        -c      <int>   how to trim low-quality at 3'-end? The raw sequences may have low-qualities at the ends. For example, the solexa sequencing always has high error rate at the first bp, and the last several bps, It will be able to get more reads aligned by trimming the low-quality ends if the raw sequence can not be aligned. Following are different strategies to trim reads.  
                0:      don't trim;
                1-10:   trim n-bps at 3-end;
                11-20:  trim first bp and (n-10)-bp at 3-end;
                21-30:  trim (n-20)-bp at 3-end if no hit for original read;
                31-40:  trim first bp and (n-30)-bp at 3-end if no hit for original read;
                41-50:  iteratively trim (n-40)-bp at 3-end until get hits;
                51-60:  same as 40-49, but trim first bp at beginning;
                default:        0
        -f      <int>   filter reads containing >n Ns, default=5
        -r      [0,1,2] how to report repeat hits, 0=none; 1=random one; 2=all, default=1
        -t      output read index, default output read id
        -n      <int>   do alignment for which chain? 0:both; 1:direct only; -1:complementary only. default=0
        -p      <int>   number of parallel processors, default=1
        -h      help

b) Options for pair-end alignment:
        -b      <string>        query b file
        -m      <int>   minimal insert size of pair-end, default=400
        -x      <int>   maximal insert size of pair-end, default=600

c) Options for mRNA tag alignment:
        -T      <int>   type of tag. There are two type of enzymes: a)DpnII, which will recognize site "GATC" and cut 16bp after the site; b)NlaIII, which will specially recognize site "CATG" and cut 17bp afterwards. 0:DpnII, GATC+16; 1:NlaIII, CATG+17. default=-1.

d) Options for miRNA alignment:
        -A      <string>        3-end adapter sequence, default=NULL. miRNA always has a size between 18~26bp, we will get into the 3'-end adapter if sequence a 27 or longer read. It's important to filter the adapter region before alignment.
        -S      <int>   number of mismatch allowed in adapter match, default=0
        -M      <int>   minimum length of miRNA, default=17
        -X      <int>   maximum length of miRNA, default=26

4. Format of output
One line for One hit. The columns are separated by '\t'.
1)id:	id of read;
2)seq:	full sequence of read. the read will be converted to the complementary sequence if mapped on the reverse chain of reference;
3)qual:	quality of sequence. corresponding to sequence, to be consistent with seq, it will be converted too if mapped on reverse chain;
4)number of hits:	number of equal best hits. the reads with no hits will be ignored;
5)a/b:	flag only meaningful for pair-end alignment, to distinguish which file the read is belonging to;
6)length:	length of the read, if aligned after trimming, it will report the information of trimmed read;
7)+/-:	alignment on the direct(+) or reverse(-) chain of the reference;
8)chr:	id of reference sequence;
9)location: location of first bp on the reference, counted from 1;
10)type:	type of hits.
	"0":	exact match
	"1~100	RefAllele->OffsetQueryAlleleQual":	number of mismatches, followed by detailed mutation sites and switch of allele types. Offset is relative to the initial location on reference. 'OffsetAlleleQual': offset, allele, and quality. Example: "2	A->10T30	C->13A32" means there are two mismatches, one on location+10 of reference, and the other on location+13 of reference. The allele on reference is A and C respectively, while query allele type and its quality is T,30 and A,32.
	"100+n	Offset":	n-bp insertion on read. Example: "101 15" means 1-bp insertion on read, start after location+15 on reference.
	"200+n Offset":	n-bp deletion on read. Example: "202 16" means 2-bp deletion on query, start after 16bp on reference.


5. Algorithm:
The program will load reference sequences into RAM, create hash tables for seed index. Then for each query, search for seeded hits, do alignment and output the result.

1) Load in reference sequences
On the contrary to Eland and Maq, which loaded query reads into RAM. SOAP stored the reference sequences in RAM. Two bits for each base, so one byte can store 4 bps. In theory, it will need L/4 bytes for reference with total sequence size L.

2) Create seed index tables
Suppose a read is splitted into 4 parts-a,b,c,d, two mismatches will be distributed on at most two of the 4 parts at the same time. So if use the combination of two parts as seed, and check for mismatches in the remainning parts, it will be able to get all hits with up to 2 mismatches. There are six combinations - ab,ac,ad,bc,bd,cd, and essentially 3 types of seeds-ab,ac,ad. So we build 3 index tables. To save RAM, we set a skip of 3-bp on the reference. The strategy is essentially the same as used in the Eland and Maq program.

3) look up table
We used look up table to judge how many mismatches between reference and read. To have best efficiency, the table used 3 bytes to check a fragment of 12-bp on a time. The table occupied 2^24=16Mb RAM.

4) search for hits
We will search for identical hits first, if no hits, then 1-mismatch hits will be picked up, then 2-mismatch hits, then gapped hits.

6. Evaluation:
Tested on a real dataset of 9,914,527 query reads (length 35bp) against a 5Mb human region, we compared the speed and sensitivity among Blastn, Eland, and SOAP. See following results:

	blastn (-F F)	megablast (-F F)	Eland	SOAP	SOAP iterative	SOAP iterative+gapped
Time consumed	41h	31h	190s	181s	258s	1402s
% reads aligned	80.7%	80.7%	88.3%	85.8%	90.6%	90.9%

SOAP and Eland is almost 600~800 times faster than megablast and blastn, while has much higher sensitivity. In principle, Eland and SOAP should have same sensitivity. Here, the result show difference is because Eland can look only the first 32bp. SOAP interative will be able to obviously improve sensitivity, and gapped alignment will further identify the gapped hits. We have to mention that SOAP and Eland have similar performance here, SOAP will be slower than Eland if running against huge reference like whole human genome, but for high quality data which have fewer errors, SOAP will be much faster.

7. Citation:
If you find SOAP did help your research, please cite the paper:
Ruiqiang Li, et. al. SOAP: short oligonucleotide alignment program. unpublished.
