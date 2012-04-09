1. Compile instruction

Enter into the soure directory.
Type "make", this will generate two runnable programs find_adapter and find_adapter_simple. 
Type "make install", this will move the runnable programs into ../bin/
Type "make clean", this will remove the runnable programs.


2. Program instruction

Currently, there is two useful programs:
find_adapter          只做正链的两序列比对，找出各种insert size下含有adapter成分的reads
find_adapter_simple   是find_adapter的简版，只找出某一给定insert size下含有adapter成分的reads

比对参数说明：
min_align_len       最小比对长度，默认为5bp，即reads与adapter的overlap长度不能低于5bp
max_mismatch_rate   最大错配率，默认为0.1，即10个碱基最多允许1个错配。

Memory requirment：
The program parse one reads at a time, include reading into memory and align to adapter.
So the memory requirement is quite small, you can almost ignore it.

Time consuming：
For a solexa lane with 4,047,480 reads,  align to a 33bp adapter by program find_adapter, the time usage is less than 1 minute. 


3. Adapter format

Each adapter has a uniq type and sequence, seperated by ":", such as the 3'-end adapter for genomic DNA, "gDNA-3:GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG". If using other adapters, you must also make to this format.

Following are some reference adapter types and sequences:
gDNA-3:GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
gDNA-5:ACACTCTTTCCCTACACGACGCTCTTCCGATCT
sRNA-3:TCGTATGCCGTCTTCTGCTTG
sRNA-5:GTTCAGAGTTCTACAGTCCGACGATC
mRNA-3:GTCGGACTGTAGAACTCTGAAC
mRNA-5:ACAGGTTCAGAGTTCTACAGTCCGACATG


4. Usage example

(1) run with default gDNA-3 adapter 
  ../bin/find_adapter  ./s_1_sequence.txt
(2) run with user defined adapters
  ../bin/find_adapter  -a sRNA-3:TCGTATGCCGTCTTCTGCTTG  ./s_1_sequence.txt


5. Input and output format

The input file is in standard fastaq (fq) format.

Three output files will be generated. The .list file contains all reads with adapter contents.
The .stat file contains summary information for each insert sizes, such as following:

total reads   :  4047480      the number of all reads
adapter reads :  76274        the number of reads with adapter content for all insert sizes
adapter rate  :  0.018845     = adapter reads / total reads
empty reads   :  36757        the number of reads with adapter and insert size is 0
empty rate    :  0.009081     = empty reads / total reads

There is also a png figure drawn by gnuplot, to show the adapter rate for each insert sizes.

The meaning for different insert sizes:
0: no insert fragment between solexa 5' and 3' sequencing adapters.
positive: the size of insert fragment between solexa 5' and 3' sequencing adapters
negative: deletion of bases from the start postion of adapter, and no insert frament

