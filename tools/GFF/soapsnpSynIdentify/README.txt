	SoapSNP Syn/nonSyn Identification Tool
					by Hu Xuesong (galaxy001@gmail.com)

Requirements:
*DBD::SQLite
*threads
*Time::HiRes
*Getopt::Std
*Term::ANSIColor

Input DATA:
*GFF annotation file, with gene id being the first item in groups and the item
 name must be either 'ID' or 'Parent'. The frame doesnot need to be stranded
 as the BioPerl required as this can be fixed by fixgff.pl.
*Genome FASTA file. the sequence must be the same as the GFF file or at least
 be the same after triming out the beginning "chr(omosome)".
*SoapSNP file filtered by a Q value.All chromosome merged into a single file.
 Also the chromosome name must be the same.

Running Order:
1.	loadgff.pl
2.	checkgff.pl	(optional)
3.	fixgff.pl	(if checkgff.pl find the GFF is not with strandard frame)
4.	parsesnp.pl
5.	update_aa.pl	(with -f if you want to use the unfixed BGI frames)
6.	output.pl

Each script come with help message on --help or no arguments. For batch run,
 remember to add -b to suppress the pause for confirming the arguments.

Since some of the scripts use virtual memory for speed, DO NOT kill a script or you are going to leave something in the memory of the running node. The left files can be removed mannually or exist until that compute node get rebooted.

The help message are colored for black background, which is classical in CLI.
You may find it too bright if you are not so classical. ^_^
								Galaxy
								20090731
