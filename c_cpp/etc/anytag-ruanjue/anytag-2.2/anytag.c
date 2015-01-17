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

static const char *version = "2.2";

int usage(){
	printf(
			"Program: anytag, efficiently convert paired short reads of descresing insert-sizes into pseudo-Sanger sequences\n"
			"          --- close the internal gap between paired-end reads\n"
			"Author : Jue Ruan <ruanjue@gmail.com>\n"
			"Version: %s\n"
			"Input  : A serial of paired-end short reads files in fq/fa format (supporting .gz). \n"
			"Output : <output_prefix>.fasta  pseudo-Sanger sequences\n"
			"         <output_prefix>.msa    multiple alignments of long sequences\n"
			"Usage  : anytag [options] <output_prefix> <XXX.insYYY[.varZZZ].s1.fq/fa[.gz]> <XXX.insYYY[.varZZZ].s2.fq/fa[.gz]> ...\n"
			"Options: \n"
			"#global\n"
			" -t <int>    Number of threads [1]*\n"
			" -X <int>    Min insert-size library to be selected as AR [largest]*\n"
			" -x <int>    Min length of long sequence (min insert size), [500]*\n"
			" -y <int>    Max length of long sequence (max insert size), [700]*\n"
			"             -x and -y give the range of long sequnece's length\n"
			" -d <int>    Insert size deviation.  When there is no '.varZZZ' in file name, use this value. [50]*\n"
			" -r <int>    Trim all reads into <-r> bp, 0: no trim [0]\n"
			" -p <int>    To reduce memory peak, split the supporting reads into <-p> parts, [1]\n"
			" -c <int>    Max number of reads in one local assembly [500]\n"
			" -f <int>    Min number of tri-mers in one sequence, filtering low complexity [10]\n"
			"#clustering\n"
			" -1          Begin parse parameters for clustering (only used in options parsing)\n"
			" -n <int>    Number of adject kmer in seed, 2-6, the smaller the less memory\n"
			"             but the less sensitive, please use the default value [4]\n"
			" -a          query ARs againest all reads (SRs + ARs), default: only SRs\n"
			" -k <int>    Kmer size 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap, [30]\n"
			" -s <float>  Min similiarity [0.97]\n"
			" -m <int>    Max mismatch [4]\n"
			"#local assembling\n"
			" -2          Begin parse parameters for assembling\n"
			" -n <int>    Number of adject kmer in seed [4]\n"
			" -k <int>    Kmer size 3 - 16 bp, [7]\n"
			" -l <int>    Min overlap, [16]\n"
			" -s <float>  Min similiarity [0.90]\n"
			" -m <int>    Max mismatch [6]\n"
			"#Only used for debugging\n"
			" -o          Output alignments into '<prefix.aln.gz>'\n"
			"             Output dot graphs into '<prefix.dot.gz>'\n"
			" -i <string> Load alignments and skip aln (will disable multi-thread)\n"
			"\n"
			"Example1: \n"
			"$> anytag -t 32 -x 550 -y 650 -r 80 my_test1 test1.ins150.var50.s1.fq.gz test1.ins150.var50.s2.fq.gz\\ \n"
			"    test1.ins200.var50.s1.fq.gz test1.ins200.var50.s2.fq.gz test1.ins250.var50.s1.fq.gz test1.ins250.var50.s2.fq.gz\\ \n"
			"    test1.ins300.var70.s1.fq.gz test1.ins300.var70.s2.fq.gz test1.ins350.var70.s1.fq.gz test1.ins350.var70.s2.fq.gz\\ \n"
			"    test1.ins400.var80.s1.fq.gz test1.ins400.var80.s2.fq.gz test1.ins450.var80.s1.fq.gz test1.ins450.var80.s2.fq.gz\\ \n"
			"    test1.ins500.var100.s1.fq.gz test1.ins500.var100.s2.fq.gz test1.ins550.var100.s1.fq.gz test1.ins550.var100.s2.fq.gz\\ \n"
			"    test1.ins600.var100.s1.fq.gz test1.ins600.var100.s2.fq.gz \n"
			"# use *ins600.s1* and *ins600.s2* as anthoring reads, its insert size is 550 ~ 650 \n"
			"Example2: \n"
			"$> anytag -t 32 -X 550 -x 500 -y 600 -d 60 -c 800 -1 -k 16 -l 50 -m 6 -s 0.95 -2 -s 0.85 test2 test2.* \n"
			"Example3:\n"
			"$> anytag -t 32 -x 550 -y 650 -d 60 -r 75 -p 8 human_fis human_data.*\n"
			"\n"
			"Note: options end with * are suggested to be customized\n"
			"\n"
			, version
	);
	return 1;
}

int usage_all(){ return usage(); }

/*
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
*/

//int main_aln(int argc, char **argv);

//int main_asm(int argc, char **argv);

int main_all(int argc, char **argv);

//int main_lnk(int argc, char **argv);

int main(int argc, char **argv){
	return main_all(argc, argv);
	/*
	if(argc < 2) return usage();
	if(strcasecmp(argv[1], "all") == 0){
		return main_all(argc - 1, argv + 1);
	//} else if(strcasecmp(argv[1], "aln") == 0){
		//return main_aln(argc - 1, argv + 1);
	//} else if(strcasecmp(argv[1], "asm") == 0){
		//return main_asm(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "lnk") == 0){
		return main_lnk(argc - 1, argv + 1);
	} else {
		return usage();
	}
	*/
}
