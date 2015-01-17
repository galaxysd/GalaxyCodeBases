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

int usage(){
	printf(
			"Program: anytag, efficiently converting paired-end short reads (PE) of stepwise descresing insert size\n"
			"                 into pseudo-Sanger(PS) sequences\n"
			"          --- close the internal gap between paired-end reads\n"
			"Author : Jue Ruan <ruanjue@gmail.com>\n"
			"Version: "VERSION"\n"
			"UID    : "UID"\n"
			"Input  : A serial of paired-end short reads files in fq/fa format (supporting .gz).\n"
			"         Short reads MUST have the same read length\n"
			"Output : <output_prefix>.fasta  PS sequences\n"
			"         <output_prefix>.msa short reads layout of PS sequences\n"
			"         <output_prefix>.aln clustered reads of PS sequences\n"
			"         <output_prefix>.dot overlap graph of PS sequences\n"
			"Usage  : anytag [options] <output_prefix> <XXX.insYYY[.varZZZ].s1.fq/fa[.gz]> <XXX.insYYY[.varZZZ].s2.fq/fa[.gz]> ...\n"
			"Options: \n"
			"#global\n"
			" -g <int>    Genome size (Mbp),  REQUIRED *\n"
			" -t <int>    Number of threads [8]*\n"
			" -X <int>    Min insert-size library to be selected as AR [all excepts smallest]\n"
			" -r <int>    Trim all reads into <-r> size, [0]\n"
			//" -R <int>    Clip first <-R> bases, [0]\n"
			" -A          Using reads from AR library as SR too, useful in joining overlapped PE reads [no]\n"
			" -m <int>    Choose modes: default [111]\n"
			"             1, both ends of PE should be unique;\n"
			"             2, either end of PE should be unique;\n"
			"             3, try assemble regardless of PE's uniqueness;\n"
			"             10, don't output the reads' information in <output_prefix>.fasta;\n"
			"             20, will output the reads' information in <output_prefix>.fasta, much larger than mode<10>;\n"
			"             100, will align the PS sequences against all PEs (takes much time), to evalute the structure of PS sequence\n"
			"             200, don't align the PS sequences against all PEs\n"
			"             1000, try assemble very difficult genome, force -m 211 -1 -n 2 -l 60 -s 1.0 -2 -l 40\n"
			"             111 means: 1 + 10 + 100\n"
			" -f <int>    Minimum number of tri-mers in a sequence, filtering low complexity [7]\n"
			" -S <int>    A given number of PS will output intermediate information (msa, dot, aln) [10000]\n"
			"#clustering\n"
			" -1          Begin parse parameters for clustering (only used in options parsing)\n"
			" -n <int>    Number of adject kmer in seed, 2-6, the smaller the less memory\n"
			"             but the less sensitive, please use the default value [4]\n"
			" -N <int>    Skip the n'th seed, -N >= 1 and -N <= -n, can be used multiple times\n"
			"             if -n == 4 and -N == 3, then only index 1,2,4 seeds [NULL]\n"
			" -k <int>    Kmer size 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap, [40]\n"
			" -s <float>  Min similiarity [0.95]\n"
			"#local assembling\n"
			" -2          Begin parse parameters for assembling\n"
			" -k <int>    Kmer size 3 - 16 bp, [9]\n"
			" -l <int>    Min overlap, [20]\n"
			" -s <float>  Min similiarity [0.95]\n"
			" -c          Make consensus again after realign, [no]\n"
			"#debug\n"
			" -D <int|->  Load seqs+index from <prefix>.dump and start from <-D> read pair,\n"
			"             if <-D> is '-', read read ID from stdin and assemble it\n"
			" -d          Dump seqs+index to <prefix>.dump\n"
			//"#gap close\n"
			//" -G <string> Close gaps in scaffolds's fasta file\n"
			""
			"\n"
	);
	return 1;
}

int usage_all(){ return usage(); }

int main_all(int argc, char **argv);

int main(int argc, char **argv){
	return main_all(argc, argv);
}
