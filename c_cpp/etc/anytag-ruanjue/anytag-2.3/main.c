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
			"Input  : A serial of paired-end short reads files in fq/fa format (supporting .gz). \n"
			"Output : <output_prefix>.fasta  PS sequences\n"
			"         <output_prefix>.msa.gz short reads layout of PS sequences\n"
			"         <output_prefix>.aln.gz clustered reads of PS sequences\n"
			"         <output_prefix>.dot.gz overlap graph of PS sequences\n"
			"Usage  : anytag [options] <output_prefix> <XXX.insYYY[.varZZZ].s1.fq/fa[.gz]> <XXX.insYYY[.varZZZ].s2.fq/fa[.gz]> ...\n"
			"Options: \n"
			"#global\n"
			" -t <int>    Number of threads [8]*\n"
			" -c <int>    Max number of SR pairs in local assembly, [400]\n"
			" -X <int>    Min insert-size library to be selected as AR [largest]\n"
			" -x <int>    Min length of PS sequence, [400]*\n"
			" -y <int>    Max length of PS sequence, [800]*\n"
			"             -x and -y give the range of PS sequnece's length\n"
			" -d <int>    Insert size deviation.  When there is no '.varZZZ' in file name, use this value. [100]\n"
			" -r <int>    Trim all reads into <-r> bp, 0: no trim [0]\n"
			" -R <int>    Clip first <-R> bases, useful bases is <-R> ~ <-r or all> [0]\n"
			" -p <int>    To reduce memory peak, split the supporting reads into <-p> parts, [1]\n"
			" -f <int>    Min number of tri-mers in one sequence, filtering low complexity [10]\n"
			" -S <int>    A given number of PS will output intermediate information (msa, dot, aln) [10000]\n"
			"#clustering\n"
			" -1          Begin parse parameters for clustering (only used in options parsing)\n"
			" -C <int>    Max number of hits per query, a larger number will slow the search,\n"
			"             while increase the accuracy of pseudo-sanger sequences. MUST bigger than 1\n"
			"             default: [500]*\n"
			" -n <int>    Number of adject kmer in seed, 2-6, the smaller the less memory\n"
			"             but the less sensitive, please use the default value [4]\n"
			" -a          query ARs againest all reads (SRs + ARs), default: only SRs [no]*\n"
			" -k <int>    Kmer size 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap, [30]\n"
			" -s <float>  Min similiarity [0.97]\n"
			" -m <int>    Max mismatch [4]\n"
			"#local assembling\n"
			" -2          Begin parse parameters for assembling\n"
			" -C <int>    Max number of hits per query, a larger number will slow the search,\n"
			"             while increase the accuracy of pseudo-sanger sequences. MUST bigger than 1\n"
			"             default: [5000]\n"
			" -n <int>    Number of adject kmer in seed [4]\n"
			" -k <int>    Kmer size 3 - 16 bp, [5]\n"
			" -l <int>    Min overlap, [10]\n"
			" -s <float>  Min similiarity [0.90]\n"
			" -m <int>    Max mismatch [6]\n"
			"#realign\n"
			" -3          Begin parse parameters for realign\n"
			" -C <int>    Max number of SR realign to PS sequence, a larger number will significantly\n"
			"             slow the realign, set to 0 to disable realign [1000]\n"
			" -l <int>    Min match region length, [50]\n"
			" -s <float>  Min similarity [0.97]\n"
			" -m <int>    Max mismatch [4]\n"
			" -G          Disable gap alignment [no]\n"
			"#Only used for debugging\n"
			" -i <string> Load alignments from *.aln.gz and skip aln (will disable multi-thread)\n"
			"\n"
			"Example1: \n"
			"$> anytag -t 32 -r 80 my_test1\\ \n"
			"    test1.ins200.s1.fq.gz test1.ins200.s2.fq.gz test1.ins250.s1.fq.gz test1.ins250.s2.fq.gz\\ \n"
			"    test1.ins300.s1.fq.gz test1.ins300.s2.fq.gz test1.ins350.s1.fq.gz test1.ins350.s2.fq.gz\\ \n"
			"    test1.ins400.s1.fq.gz test1.ins400.s2.fq.gz test1.ins450.s1.fq.gz test1.ins450.s2.fq.gz\\ \n"
			"    test1.ins600.s1.fq.gz test1.ins600.s2.fq.gz \n"
			"# use *ins600.s1* and *ins600.s2* as anthor reads, the length of PS will be 500 ~ 700 \n"
			"Example2: \n"
			"$> anytag -x 150 -y 300 -r 100 -1 -l 50 -a -3 -C 0 joint_reads joint.ins200.var70.s1.fq.gz joint.ins200.var70.s2.fq.gz\n"
			"# join reads, -a is enabled\n"
			"Example3:\n"
			"$> anytag -t 32 -x 500 -y 700  -r 100 -p 4 human_fis human_data.*\n"
			"Example4:\n"
			"$> anytag -t 32 -1 -n 2 -s 1 -m 0 het het.*.fq\n"
			"# disable mismatch in clustering\n"
			"\n"
			"Note: options end with * are suggested to be customized\n"
			"\n"
	);
	return 1;
}

int usage_all(){ return usage(); }

int main_all(int argc, char **argv);

int main(int argc, char **argv){
	return main_all(argc, argv);
}
