/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "anytag_aux.h"

static const char *version = "2.0";

int usage(){
	printf(
			"\n"
			"Program: anytag, convert paired short reads of spread insert-sizes into long sequences\n"
			"Author : Jue Ruan <ruanjue@gmail.com>\n"
			"Version: %s\n"
			"Input  : A serial of paired-end short reads files in fq/fa format (supporting .gz). \n"
			"         Each pair of files (s1 and s2) contain one range of spread insert-size (YYY).\n"
			"         Anytag requires that the insert size of those paired-end short reads must be continuous\n"
			"Output : \n"
			"         <output_prefix>.fasta  long sequences\n"
			"         <output_prefix>.sam.gz short read alignments in long sequences\n"
			"         <output_prefix>.msa.gz multiple alignments of long sequences\n"
			"Usage: anytag <all|aln|asm|pair>\n"
			" all    run aln and asm at once, it is more fast and efficient than separately executions\n"
			" aln    *Obsolete*, align paired short reads (supporting reads SR) against anchoring reads (AR)\n"
			"        and group SRs by each AR, then output localized groups\n"
			" asm    *Obsolete*, build graph that reads as node, and edge representing overlap\n"
			"        and search a connective path start from read1 to read2 of AR, then call\n"
			"        consensus sequence\n"
			" lnk    mimic mate-paired information for FIS, requiring one mate-paired library of short reads\n"
			, version
	);
	return 1;
}

int usage_all(){
	printf(
			"Usage  : anytag all [options] <output_prefix> <inputXXX.insYYY.s1.fq/fa[.gz]> <inputXXX.insYYY.s2.fq/fa[.gz]> ...\n"
			"Options: \n"
			"#global\n"
			" -t <int>    Number of threads [1]*\n"
			" -c <int>    Max number of reads in one local assembly [400]\n"
			" -f <int>    Min number of tri-mers in one sequence, filtering low complexity [12]\n"
			" -r <int>    All reads have the same length <-r>, set to 0 if reads' length is various, [0]\n"
			"#aln\n"
			" -1          Begin parse parameters for aln (only used in options parsing)\n"
			" -X <int>    Min insert-size library to be selected as AR [auto]*\n"
			" -n <int>    Number of adject kmer in seed, 2-6, please use the default value [4]\n"
			" -k <int>    Kmer size 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap, [30]\n"
			" -s <float>  Min similiarity [0.97]\n"
			" -m <int>    Max mismatch [4]\n"
			"#asm\n"
			" -2          Begin parse parameters for asm\n"
			" -x <int>    Min length of long sequence (min insert size), [500]*\n"
			" -y <int>    Max length of long sequence (max insert size), [700]*\n"
			"             -x and -y give the range of long sequnece's length\n"
			" -n <int>    Number of adject kmer in seed, 2-6, please use the default value [4]\n"
			" -k <int>    Kmer size 3 - 16 bp, [4]\n"
			" -l <int>    Min overlap, [8]\n"
			" -s <float>  Min similiarity [0.90]\n"
			" -m <int>    Max mismatch [6]\n"
			" -g <int>    Threshold (number of reads) to trun off gap alignment [100]\n"
			"\n"
			"Example1: \n"
			"$> anytag -t 32 -x 500 -y 700 -r 80 test1 \\ \n"
			"    test1.ins150.s1.fq.gz test1.ins150.s2.fq.gz\\ \n"
			"    test1.ins200.s1.fq.gz test1.ins200.s2.fq.gz\\ \n"
			"    test1.ins250.s1.fq.gz test1.ins250.s2.fq.gz\\ \n"
			"    test1.ins300.s1.fq.gz test1.ins300.s2.fq.gz\\ \n"
			"    test1.ins350.s1.fq.gz test1.ins350.s2.fq.gz\\ \n"
			"    test1.ins400.s1.fq.gz test1.ins400.s2.fq.gz\\ \n"
			"    test1.ins450.s1.fq.gz test1.ins450.s2.fq.gz\\ \n"
			"    test1.ins500.s1.fq.gz test1.ins500.s2.fq.gz\\ \n"
			"    test1.ins550.s1.fq.gz test1.ins550.s2.fq.gz\\ \n"
			"    test1.ins600.s1.fq.gz test1.ins600.s2.fq.gz \n"
			"# use *ins600.s1* and *ins600.s2* as anthoring reads, its insert size is 500 ~ 700 \n"
			"Example2: \n"
			"$> anytag -t 32 -X 550 -x 400 -y 700 -c 800 -1 -k 16 -l 32 -m 6 -s 0.95 -2 -s 0.85 -g 200 test2 test2.ins*.s?.fq.gz \n"
			"\n"
			"Note: options end with * are suggested to be customized\n"
			"\n"
	);
	return 1;
}

int usage_aln(){
	printf(
			"Usage: anytag aln [options] <inputXXX.insYYY.s1.fq/fa[.gz]> <inputXXX.insYYY.s2.fq/fa[.gz]> ...\n"
			"Options:\n"
			" -o <string> Output file [stdout]\n"
			"{global and aln options, please see `anytag all`}\n"
			"\n"
		  );
	return 1;
}

int usage_asm(){
	printf(
			"Usage: anytag asm [options] <alignment> <longreads_prefix>\n"
			"Options:\n"
			"{global and asm options, please see `anytag all`}\n"
			"\n"
			"Example1:\n"
			" # read from stdin, and output fis.fasta, fis.sam.gz fis.msa.gz\n"
			" anytag asm -k 4 -l 8 -x 550 -y 650 -t 8 - fis\n"
			"\n"
		  );
	return 1;
}

int usage_lnk(){
	printf(
			"Usage: anytag pair [options]\n"
			"Options:\n"
			" -i <string> FIS file in fasta format [*required]\n"
			" -1 <string> Read1 file of mate-paired [*required]\n"
			" -2 <string> Read2 file of mate-paired [*required]\n"
			" -o <string> Output file, mate-paired UID of FIS [stdout]\n"
			" -t <int>    Number of threads [1]\n"
			" -r <int>    All mate-paired short reads have the same length <-r>, set to 0 if reads' length is various, [0]\n"
			" -n <int>    Number of adject kmer in seed, 2-6 [3]\n"
			" -k <int>    Kmer size 3 - 16 bp, [15]\n"
			" -s <float>  Min similiarity [0.95]\n"
			" -m <int>    Max mismatch 0-7 [6]\n"
			" -g          Trun on gap alignment\n"
			"\n"
			"Example1:\n"
			" anytag pair -x 5000 -i fis.fasta -o fis.fasta.lib.gz -1 read1.fa -2 read2.fq.gz -t 32\n"
			"\n"
		  );
	return 1;
}

int main_aln(int argc, char **argv);

int main_asm(int argc, char **argv);

int main_all(int argc, char **argv);

int main_lnk(int argc, char **argv);

int main(int argc, char **argv){
	if(argc < 2) return usage();
	if(strcasecmp(argv[1], "all") == 0){
		return main_all(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "aln") == 0){
		return main_aln(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "asm") == 0){
		return main_asm(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "lnk") == 0){
		return main_lnk(argc - 1, argv + 1);
	} else {
		return usage();
	}
}
