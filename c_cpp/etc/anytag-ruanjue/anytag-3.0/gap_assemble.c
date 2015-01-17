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
 
#include "sr_aln.h"
#include "file_reader.h"

thread_beg_def(mall);
thread_end_def(mall);


SR_SeqDB* build_sdb_from_files(char **files, int n_file, uint32_t n_cpu, uint32_t ksize, uint32_t n_seed, uint32_t low_cpx, uint32_t rd_clip, uint32_t rd_trim){
	SR_SeqDB *sdb;
	FileReader *s1, *s2;
	Sequence *seq1, *seq2;
	int is_fq1, is_fq2, i;
	if(n_file == 0 || (n_file % 2)){
		fprintf(stderr, " -- number of paired-end files is not even in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
	for(i=0;i<n_file;i+=2){
		if((s1 = fopen_filereader(files[i + 0])) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", files[i + 0], __FUNCTION__, __FILE__, __LINE__); abort();
		}
		if((s2 = fopen_filereader(files[i + 1])) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", files[i + 1], __FUNCTION__, __FILE__, __LINE__); abort();
		}
		is_fq1 = (guess_seq_file_type(s1) == 2);
		is_fq2 = (guess_seq_file_type(s2) == 2);
		while(1){
			seq1 = NULL;
			seq2 = NULL;
			if(!((is_fq1)? fread_fastq_adv(&seq1, s1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq1, s1, FASTA_FLAG_NO_NAME))) break;
			if(!((is_fq2)? fread_fastq_adv(&seq2, s2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq2, s2, FASTA_FLAG_NO_NAME))) break;
			if(seq1->seq.size > (int)(rd_clip + rd_trim) && seq2->seq.size > (int)(rd_clip + rd_trim)){
				sr_push_sdb(sdb, seq1->seq.string + rd_clip, seq1->seq.size - rd_clip - rd_trim);
				sr_push_sdb(sdb, seq2->seq.string + rd_clip, seq2->seq.size - rd_clip - rd_trim);
			}
			if(seq1) free_sequence(seq1);
			if(seq2) free_sequence(seq2);
		}
		fclose_filereader(s1);
		fclose_filereader(s2);
	}
	sr_ready_sdb(sdb);
	return sdb;
}

int usage (){
	printf(
"Gap_assemble --  Replace Ns in gap by local assembly\n"
"Usage: gao_assemble [Options] <scaffolds_fasta> <lib1.read1.fa/fq> <lib1.read2.fa/fq> [...]\n"
"Options:\n"
" -k <string> Kmer size, 3 - 16, [15]\n"
			);
	return 1;
}

int main(int argc, char **argv){
}

