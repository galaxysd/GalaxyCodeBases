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
 
#define _GNU_SOURCE
#include <unistd.h>
#include "ps_assembly.h"
#include "file_reader.h"

int usage(){
	printf(
	"Fill in gap between the first two reads using fellowing reads\n"
	"Usage: fillin [Options] <fasta_file>\n"
	"Options:\n"
	" -x <int>    Minimum length, [400]\n"
	" -y <int>    Maximum length, [800]\n"
	" -k <int>    Kmer size, [9]\n"
	" -l <int>    Minimum overlap, [20]\n"
	" -s <float>  Minimum similarity, [0.95]\n"
	);
	return 1;
}

int main(int argc, char **argv){
	Graph *g;
	FileReader *inp;
	Sequence *seq;
	float sm;
	uint32_t rid;
	int x, y, ksize, ol, c;
	x = 400;
	y = 800;
	ksize = 9;
	ol = 20;
	sm = 0.95;
	while((c = getopt(argc, argv, "hx:y:k:l:s:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'x': x = atoi(optarg); break;
			case 'y': y = atoi(optarg); break;
			case 'k': ksize = atoi(optarg); break;
			case 'l': ol = atoi(optarg); break;
			case 's': sm = atof(optarg); break;
		}
	}
	if(optind == argc) return usage();
	inp = fopen_filereader(argv[optind]);
	seq = NULL;
	g = init_graph(ksize, ol, sm);
	rid = 0;
	while(fread_fasta(&seq, inp)){
		push_graph(g, rid++, seq->seq.string, seq->seq.size, 0, 0);
	}
	fclose_filereader(inp);
	ready_graph(g);
	index_graph(g);
	align_graph(g);
	set_insert_graph(g, x, y);
	simplify_graph(g);
	shave_graph(g);
	allpaths_graph(g);
	validate_paths_graph(g);
	consensus_graph(g);
	output_contigs_graph(g, NULL, stdout);
	print_dot_graph(g,stdout);
	free_graph(g);
	return 0;
}
