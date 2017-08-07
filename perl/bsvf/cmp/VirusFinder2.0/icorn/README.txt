iCORN manual

iCORN (iterative correction of reference nucleotides) can correct genome sequences with short reads. Reads are mapped iteratively against the genome sequences, so far by SSAHA. Discrepancies between the multiple alignments of the mapping reads and reference are corrected, if by the correction the amount of perfect mapping reads doesn't decrease.

1.  Downloads and install
Download iCORN from Webpage. Uncompact the program: tar xzf icorn-0.95.tar.gz 
If you type <cd iCORN && pwd>, you will know, where you installed iCORN. Please set the global bash variable
ICORN_HOME=<result of pwd command> # bash
set ENV ICORN_HOME=<result of pwd command> # other shells
It might be good to do:
chmod 755 $ICORN_HOME


Download SSAHA_pileup.
wget ftp://ftp.sanger.ac.uk/pub/zn1/ssaha_pileup/pileup_v0.5.tgz
tar xzf pileup_v0.5.tgz
cd pileup_v0.5
./install.csh
If you type <pwd>, you will know, where you installed the SSAHA_pileup. Please set the global bash variable
PILEUP_HOME=<result of pwd command> # bash
set ENV PILEUP_HOME=<result of pwd command> # other shells

Download SNP-O-Matic.
svn co https://snpomatic.svn.sourceforge.net/svnroot/snpomatic
cd snpomatic
make
the program is call findknownsnps.
If you type <pwd>, you will know, where you installed the snpomatic (findknownsnps). Please set the global bash variable 
SNPOMATIC_HOME=<result of pwd command> # bash
set ENV SNPOMATIC_HOME=<result of pwd command> # other shells


2. Set path
First you have to set to include the path of iCORN to the PATH variable. Assuming you install iCORN in $HOME/bin/icorn. Please create the variable ICORN_HOME with
ICORN_HOME=$HOME/bin/icorn. Also export the variables, before running (export ICORN_HOME;

For the third party software, please set following both variables:
SNPOMATIC_HOME =< where you installed SNP-O-Matic_ >,
PILEUP_HOME          =< where you installed SSAHA_pileup >

export ICORN_HOME SNPOMATIC_HOME PILEUP_HOME

Best is to write those commands in the .bashrc.

3. What you need to run iCORN.
To run iCORN you need:
Reference sequence that will be corrected as single or multi fasta file

Mate pairs reads Fastq files of forward and reverse reads in sanger format
The mean of the fragment size
The range of the fragment size

[
3a RUNNING on a farm: If you have a farm with a LSF (bsub) system,
you can use the parallel version of iCORN. Just set the LSF varialbe
to 1 in the icorn.start.sh file (line 35)
]

4. Run iCORN
To start, best create a new working directory (later referred to main working directory). Reference is used to refer to the reference sequence and reads are referred to the sequence file of the short read. 
Start the program
>>$ICORN_HOME/icorn.starts.sh Reference <iteration start, mostly 1>  <iteration stop, i. e. 5> <fastq forward> <fastq reverse> <fragment size range: i.e. 100,400> <mean fragment size 300>

By just runnig icorn.start.sh 

5. What is iCORN doing:
In each iteration a directory with the name Reference.x is created, where Reference is the name of the reference and x the iteration. Reads are mapped; errors are called and double-checked. In the main working directory are written: 
Reference.<i+1>  - the new corrected reference sequence
Reference.<i+1>.<contignames>.gff - has all the corrections on in the i-th iteration, projected on the Reference.1 genome (The original file)!
Reference.<i+1>.General.stats - report of corrections done in the i-th iteration
Stats.Mapping.cvs - Mapping statistics of all iterations done so far
Stats.Correction.cvs - Correction statistics of all corrections done so far

In the directory Reference.i/plot are generated the coverage plots for each sequence/contigs of the ith iteration.

We differ 7 type of "correction tags", used in the stats and the Reference.<i+1>.<contignames>.gff files:
SNP - corrected one nucleotide error
INS - corrected 1-3bp wrongly shown deletions
DEL - corrected 1-3bp insertion
HETERO - At least 2 alleles with the ratio of 0.15 are called. The ref.<i+1> is changed to the most abundant allele.
Rej.SNP/Rej.INS/Rej.DEL - rejected correction, as the coverage of perfect mapping reads decreased. 

In the Reference.<i+1>.<contignames>.gff file you find further information about the coverage of mapping reads mapped by SSAHA. The coverage of perfect mapping reads before and after the correction, and a pileup representation of three rows upstream and downstream of the called correction. It is a representation of multiple alignments over this region. An example of an entry is:
<p align="center"><img src="Example.gfffile.tiff" width="500"/>

6. When iCORN finishes
iCORN stops after the set amount of iteration, or if it does find no more errors. 

After the last iteration, all the Reference.X.<contigname>.gff files are joined to All.<contigname>.gff. 
The files Stats.Mapping.csv and Stats.Correction.csv are generated.
Beside the ssaha coverage plots, also perfect mapping plots are generated in the directory plots/. The plots are generated for the original sequence and the "best corrected" reference (Final.corrected.fa to symbolic link to Reference.X). 

7. How to read the iCORN results
First the Stats.Mapping.csv file should be consulted. Depending on the quality of the fastq file, around 90% of the reads should map. Depending on the amount of repeats the amount of uniquely mapping reads (read pairs) should be 60-80%. If mapping to a draft genome, this number might be around 40% as most mate pairs are on different contigs.

In Stats.Correction. csv for each "correction tag" the amount of found events are reported. If the reference is a multiple fasta file, this file is difficult to read.

Next the results can be loaded in Artemis:
* Start Artemis
* Load the sequence or EMBL file. (<File> + <Read An Entry...>)
* Load feature file All.Reference.gff file. (If your reference was in multi-fasta file, you will need to open each sequence on his own. The sequence will be in Sequence.1/ (original file) and in Sequence.best/ (for the best correction). 
* Load the graph from the plot directory. (<Graph> + <Read User Plot ..>)
* Then you can analyse the correction. 
* A good tutorial for Artemis can be found <a href=http://pseudomonas-syringae.org/artemis_tutorial.htm>here</a>.



7. Optional transformation of results
Transform results to Gap4 - sorry program still beta version:
You will further need the padded consensus of the "reference" sequence. The contigs/chromosomes must have the same name in the padded and unpadded version.
Call: icorn.ToGap4.pl <padded.consensus> <Resultname> <  <(cat All.*.gff)>
A TAG file with the name Resultname.gap4 is generated that can be loaded in GAP4 or GAP5. The information are equivalent to the All.<contigname>.gff files.


