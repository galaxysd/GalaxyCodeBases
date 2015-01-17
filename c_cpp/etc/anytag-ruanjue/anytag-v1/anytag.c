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
#include "all_path.h"
#include "thread.h"
#include <regex.h>
#include <unistd.h>

typedef struct {
	char *name;
	char *prefix;
	int ins_size;
	int pair_idx;
	int is_fq;
} SeqFile;

define_list(sfv, SeqFile);

void parse_seqfile(SeqFile *sf, char *name){
	regex_t pat;
	regmatch_t mats[5];
	int ret;
	regcomp(&pat, "(.*?)\\.ins([0-9]+)\\.s([12])\\.f([aq])", REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	ret = regexec(&pat, name, 5, mats, 0);
	if(ret == REG_NOMATCH){
		fprintf(stdout, " -- Wrong file name \"%s\" in %s -- %s:%d --\n", name, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
	}
	sf->name = catstr(1, name);
	strcpy(sf->name, name);
	sf->prefix = malloc(mats[1].rm_eo - mats[1].rm_so + 1);
	memcpy(sf->prefix, sf->name + mats[1].rm_so, mats[1].rm_eo - mats[1].rm_so);
	sf->prefix[mats[1].rm_eo - mats[1].rm_so] = 0;
	sf->ins_size = strtol(sf->name + mats[2].rm_so, NULL, 10);
	sf->pair_idx = strtol(sf->name + mats[3].rm_so, NULL, 10);
	sf->is_fq    = (sf->name[mats[4].rm_so] == 'q');
	regfree(&pat);
}

int cmp_seqfile(const void *e1, const void *e2){
	SeqFile *s1, *s2;
	int cmp;
	s1 = (SeqFile*)e1;
	s2 = (SeqFile*)e2;
	if(s1->ins_size == s2->ins_size){
		cmp = strcmp(s1->prefix, s2->prefix);
		if(cmp == 0){
			if(s1->pair_idx == s2->pair_idx) return 0;
			else if(s1->pair_idx < s2->pair_idx) return -1;
			else return 1;
		} else return cmp;
	} else if(s1->ins_size < s2->ins_size) return 1;
	else return -1;
}

int usage(){
	printf(
		"Usage: anytag <aln|asm|realn>\n"
		" aln    align a set of paired short reads against itself,\n"
		"        find pairs that one end overlap KEY pairs at least\n"
		"        min_overlap basepairs, then output KEY pairs and all\n"
		"        of related pairs. KEY pairs are specified by insert size\n"
		" asm    build string graph for KEY pairs one by one, simpify,\n"
		"        and search a path start from read1 to read2, then call\n"
		"        consensus sequence\n"
		" realn  refine alignment from 'aln' to increase min_overlap,\n"
	);
	return 1;
}

int usage_aln(){
	printf(
			"Usage: anytag aln [options] <prefix##.rankNNN.s1.fq/fa> <prefix##.rankNNN.s2.fq/fa> ...\n"
			"Options:\n"
			" -k <int>    Seed size, 3 - 16 bp, [15]\n"
			" -l <int>    Min overlap [30]\n"
			" -s <float>  Min similiarity [0.97]\n"
			" -m <int>    Max mismatch [6]\n"
			" -M <int>    Max potential hits per query [1024]\n"
			" -f <int>    Min number of tri-mers in one sequence, filtering low complexity [10]\n"
			" -g          Gap alignment (Smith-Waterman), [no]\n"
			" -t <int>    Number of thread [1]\n"
			" -o <string> Prefix of output files [null]\n"
			" -x <int>    Min insert size, used to assemble longreads, [auto]\n"
			" -y <int>    Max insert size, used to assemble longreads, [auto]\n"
			"\n"
			"Example1:\n"
			" anytag -k 15 -l 31 -t 8 -x 500 -y 520 -o drops.aln drops.ins120.s1.fq.gz drops.ins120.s2.fq.gz ... drops.ins520.s2.fq.gz\n"
			"\n"
		  );
	return 1;
}

int usage_asm(){
	printf(
			"Usage: anytag asm [options] <alignment> <longreads_prefix>\n"
			"Options:\n"
			" -x <int>    Min insert size [400]\n"
			" -y <int>    Max insert size [800]\n"
			" -k <int>    Seed size, 3 - 16 bp, [3]\n"
			" -l <int>    Min overlap [6]\n"
			" -s <float>  Min similarity [0.90]\n"
			" -M <int>    Minimum number of reads to discrad local assemble [400]\n"
			" -t <int>    Number of thread [1]\n"
		  );
	return 1;
}

static volatile uint32_t __anytag_n_vis = 0;

int at_cmp_hit_rid(const void *e1, const void *e2){
	uint32_t rid1, rid2;
	rid1 = (*((uint64_t*)e1)) >> 32;
	rid2 = (*((uint64_t*)e2)) >> 32;
	if(rid1 == rid2) return 0;
	else if(rid1 < rid2) return -1;
	else return 1;
}

int simp_load_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *in){
	uint64_t id;
	if(fread(&id, sizeof(uint64_t), 1, in) == 1){
		hit->rid1 = id >> 32;
		hit->rid2 = id & 0xFFFFFFFFU;
		hit->dir1 = 0;
		hit->dir2 = 0;
		return 1;
	} else return 0;
	sdb = sdb;
}

int simp_dump_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	uint64_t id;
	id = (((uint64_t)hit->rid1) << 32) | hit->rid2;
	return fwrite(&id, sizeof(uint64_t), 1, out);
	sdb = sdb;
}

void output_hit_paired_seqs(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	int n;
	uint64_t id;
	if((hit->rid1 >> 1) < __anytag_n_vis){
		id = (((uint64_t)hit->rid1) << 32) | hit->rid2;
		n = fwrite(&id, sizeof(uint64_t), 1, out);
	}
	if((hit->rid2 >> 1) < __anytag_n_vis){
		id = (((uint64_t)hit->rid2) << 32) | hit->rid1;
		n = fwrite(&id, sizeof(uint64_t), 1, out);
	}
	sdb = sdb;
}

void sort_aln_hit_file(FILE *in, uint32_t max_hits_per_sort, uint32_t max_hits_per_rid, FILE *tmp, FILE *out){
	u64list **alns, *cache;
	u64list *begs, *ends, *offs;
	u32list *idxs;
	uint64_t n_hit, m_hit, l_hit, d_hit;
	uint32_t n_part, c_hit, i, min, midx, rid;
	fseek(in, 0, SEEK_END);
	n_hit = ftell(in) / 8;
	n_part = (n_hit + max_hits_per_sort - 1) / max_hits_per_sort;
	if(n_part == 0) n_part = 1;
	max_hits_per_sort = (n_hit + n_part - 1) / n_part;
	alns = malloc(sizeof(u64list*) * n_part);
	alns[0] = init_u64list(max_hits_per_sort);
	fseek(in, 0, SEEK_SET);
	fseek(tmp, 0, SEEK_SET);
	begs = init_u64list(n_part);
	ends = init_u64list(n_part);
	l_hit = 0;
	for(i=0;i<n_part;i++){
		m_hit = fread(ref_u64list(alns[0], 0), sizeof(uint64_t), max_hits_per_sort, in);
		qsort(ref_u64list(alns[0], 0), m_hit, sizeof(uint64_t), at_cmp_hit_rid);
		m_hit = fwrite(ref_u64list(alns[0], 0), sizeof(uint64_t), m_hit, tmp);
		push_u64list(begs, l_hit);
		push_u64list(ends, l_hit + m_hit);
		l_hit += m_hit;
	}
	free_u64list(alns[0]);
	c_hit  = 1024;
	for(i=0;i<n_part;i++) alns[i] = init_u64list(c_hit);
	idxs = init_u32list(n_part);
	offs = init_u64list(n_part);
	for(i=0;i<n_part;i++) push_u32list(idxs, 1);
	for(i=0;i<n_part;i++) push_u64list(offs, 0);
	cache = init_u64list(max_hits_per_rid);
	rid = 0xFFFFFFFFU;
	while(1){
		min = 0xFFFFFFFFU;
		midx = 0;
		for(i=0;i<n_part;i++){
			if(get_u32list(idxs, i) >= count_u64list(alns[i])){
				if(get_u32list(idxs, i) == 0) continue;
				set_u32list(idxs, i, 0);
				d_hit = get_u64list(ends, i) - get_u64list(begs, i) - get_u64list(offs, i);
				if(d_hit == 0){ set_u64list_size(alns[i], 0); continue; }
				if(d_hit > c_hit) d_hit = c_hit;
				fseek(tmp, (get_u64list(begs, i) + get_u64list(offs, i)) * 8LLU, SEEK_SET);
				d_hit = fread(ref_u64list(alns[i], 0), sizeof(uint64_t), d_hit, tmp);
				set_u64list(offs, i, get_u64list(offs, i) + d_hit);
				set_u64list_size(alns[i], d_hit);
			}
			if((get_u64list(alns[i], get_u32list(idxs, i)) >> 32) < min){
				min = get_u64list(alns[i], get_u32list(idxs, i)) >> 32;
				midx = i;
			}
		}
		if(min == 0xFFFFFFFFU) break;
		if(min != rid){
			if(rid != 0xFFFFFFFFU){
				if(min < rid){
					fprintf(stderr, " -- unexpected error %u < %u in %s -- %s:%d --\n", min, rid, __FUNCTION__, __FILE__, __LINE__);
					fflush(stderr);
					abort();
				} else if(count_u64list(cache) <= max_hits_per_rid){
					dump_u64list(cache, out);
				}
			}
			rid = min;
			clear_u64list(cache);
		}
		push_u64list(cache, get_u64list(alns[midx], get_u32list(idxs, midx)));
		set_u32list(idxs, midx, get_u32list(idxs, midx) + 1);
	}
	fflush(out);
	free_u64list(cache);
	for(i=0;i<n_part;i++) free_u64list(alns[i]);
	free_u64list(begs);
	free_u64list(ends);
	free_u32list(idxs);
	free_u64list(offs);
	free(alns);
}

int main_aln(int argc, char **argv){
	SR_SeqDB *sdb;
	sfv *sfs;
	SeqFile *sf1, *sf2;
	FileReader *s1, *s2;
	Sequence *seq1, *seq2;
	FILE *log, *gzout, *tmp1, *tmp2;
	char *prefix, PREFIX[256], TMP1[256], TMP2[256], str[256];
	uint64_t v, off;
	uint32_t i, rid, rrid, len, dir;
	int ii, c, kmer_size, rd_len, min_ol, max_mm, max_hit, s, gap, n_cpu, min_ins, max_ins, x, y, vis, min_cpx;
	int fd;
	float min_sm;
	kmer_size = 15;
	rd_len = 0;
	min_ol = 31;
	min_sm = 0.9;
	max_mm = 5;
	s = 1;
	max_hit = 1024;
	gap = 0;
	n_cpu = 1;
	prefix = NULL;
	min_ins = -1;
	max_ins = -1;
	min_cpx = 10;
	x = 1;
	y = 1;
	while((c = getopt(argc, argv, "hgk:r:l:s:m:M:f:t:o:x:y:")) != -1){
		switch(c){
			case 'k': kmer_size = atoi(optarg); break;
			case 'r': rd_len = atoi(optarg); break;
			case 'l': min_ol = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'm': max_mm = atoi(optarg); break;
			case 'M': max_hit = atoi(optarg); break;
			case 'g': gap = 1; break;
			case 't': n_cpu = atoi(optarg); break;
			case 'o': prefix = optarg; break;
			case 'x': x = 0; min_ins = atoi(optarg); break;
			case 'y': y = 0; max_ins = atoi(optarg); break;
			case 'f': min_cpx = atoi(optarg); break;
			default : return usage_aln();
		}
	}
	if(optind == argc || (argc - optind) % 2 == 1) return usage_aln();
	if(prefix == NULL) return usage_aln();
	sfs = init_sfv(4);
	for(ii=optind;ii<argc;ii++){
		sf1 = next_ref_sfv(sfs);
		parse_seqfile(sf1, argv[ii]);
	}
	qsort(ref_sfv(sfs, 0), count_sfv(sfs), sizeof(SeqFile), cmp_seqfile);
	for(i=0;2*i<count_sfv(sfs);i++){
		if(strcasecmp(ref_sfv(sfs, 2 * i)->prefix, ref_sfv(sfs, 2 * i + 1)->prefix)
			|| ref_sfv(sfs, 2 * i)->ins_size != ref_sfv(sfs, 2 * i + 1)->ins_size
			|| ref_sfv(sfs, 2 * i)->pair_idx >= ref_sfv(sfs, 2 * i + 1)->pair_idx){
			fprintf(stdout, " -- Wrong pair \"%s\" and \"%s\" in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		if((s1 = fopen_filereader(ref_sfv(sfs, 2 * i)->name)) == NULL){
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i)->name, __FUNCTION__, __FILE__, __LINE__);
			return 1;
		} else fclose_filereader(s1);
		if((s2 = fopen_filereader(ref_sfv(sfs, 2 * i + 1)->name)) == NULL){
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", ref_sfv(sfs, 2 * i + 1)->name, __FUNCTION__, __FILE__, __LINE__);
			return 1;
		} else fclose_filereader(s2);
		if(y && ref_sfv(sfs, 2 * i)->ins_size > max_ins){
			max_ins = ref_sfv(sfs, 2 * i)->ins_size;
		}
	}
	if(x) min_ins = max_ins;
	sdb = sr_init_sdb(prefix, n_cpu, kmer_size, rd_len);
	sr_set_align_parameters(sdb, s, min_ol, min_sm, max_mm, gap);
	sr_set_filter_parameters(sdb, min_cpx, max_hit);
	sprintf(PREFIX, "%s.info", prefix);
	if((log = fopen(PREFIX, "w")) == NULL) log = stdout;
	sprintf(PREFIX, "gzip -1 -c >%s.gz", prefix);
#ifdef DEBUG
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
#endif
	if((gzout = popen(PREFIX, "w")) == NULL){
		fprintf(stdout, " -- Cannot invoke '%s' in %s -- %s:%d --\n", PREFIX, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
	}
	sprintf(TMP1, "%s.tmp.XXXXXX", prefix);
	fd = mkstemp(TMP1);
	close(fd);
	sprintf(TMP2, "%s.tmp.XXXXXX", prefix);
	fd = mkstemp(TMP2);
	close(fd);
	tmp1 = fopen(TMP1, "w+");
	tmp2 = fopen(TMP2, "w+");
	sr_set_output_func(sdb, output_hit_paired_seqs);
	sr_set_load_dump_func(sdb, simp_load_hit, simp_dump_hit);
	rid = 0;
	for(i=0;2*i<count_sfv(sfs);i++){
		seq1 = seq2 = NULL;
		sf1 = ref_sfv(sfs, 2 * i);
		sf2 = ref_sfv(sfs, 2 * i + 1);
		vis = ((sf1->ins_size >= min_ins) && (sf1->ins_size <= max_ins));
		if((s1 = fopen_filereader(sf1->name)) == NULL){
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", sf1->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		if((s2 = fopen_filereader(sf2->name)) == NULL){
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", sf2->name, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		rrid = rid;
		while(1){
			if(!((sf1->is_fq)? fread_fastq_adv(&seq1, s1, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq1, s1, FASTA_FLAG_NO_NAME))) break;
			if(!((sf2->is_fq)? fread_fastq_adv(&seq2, s2, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL) : fread_fasta_adv(&seq2, s2, FASTA_FLAG_NO_NAME))) break;
			sr_push_sdb(sdb, seq1->seq.string, seq1->seq.size, vis);
			sr_push_sdb(sdb, seq2->seq.string, seq2->seq.size, vis);
			rid ++;
		}
		if(seq1) free_sequence(seq1);
		if(seq2) free_sequence(seq2);
		fclose_filereader(s1);
		fclose_filereader(s2);
		fprintf(log, "%s\t%d\t%d\t%d\t%s\n", sf1->prefix, sf1->ins_size, rrid, rid, vis? "yes" : "no");
		fflush(log);
		if(vis) __anytag_n_vis = rid;
	}
	sr_ready_sdb(sdb);
	sr_aln_sdb(sdb);
	fflush(sdb->out);
	sort_aln_hit_file(sdb->out, 4LLU * 1024 * 1024 * 1024 / 8, max_hit * 2, tmp1, tmp2);
	fclose(tmp1);
	unlink(TMP1);
	fseek(tmp2, 0, SEEK_SET);
	rid = 0xFFFFFFFFU;
	while((c = fread(&v, sizeof(uint64_t), 1, tmp2)) == 1){
		if((v >> 33) != rid){
			if(rid != 0xFFFFFFFFU && (v >> 33) < rid){
				fprintf(stderr, " -- unexpected error %u < %u in %s -- %s:%d --\n", (unsigned)(v>>33), rid, __FUNCTION__, __FILE__, __LINE__);
				fflush(stderr);
				abort();
			}
			rid = v >> 33;
			off = sr_rdseq_offset(sdb, 2 * rid);
			len = sr_rdseq_length(sdb, 2 * rid);
			bits2seq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "T\t%u\t+\t%s", rid, str);
			off = sr_rdseq_offset(sdb, 2 * rid + 1);
			len = sr_rdseq_length(sdb, 2 * rid + 1);
			bits2revseq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "\t%s\n", str);
		}
		dir = ((v >> 32) ^ (v & 0xFFFFFFFFU)) & 0x01;
		rrid = (v & 0xFFFFFFFFU) >> 1;
		fprintf(gzout, "S\t%u\t%c", rrid, "+-"[dir]);
		if(dir){
			off = sr_rdseq_offset(sdb, 2 * rrid + 1);
			len = sr_rdseq_length(sdb, 2 * rrid + 1);
			bits2seq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "\t%s", str);
			off = sr_rdseq_offset(sdb, 2 * rrid);
			len = sr_rdseq_length(sdb, 2 * rrid);
			bits2revseq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "\t%s\n", str);
		} else {
			off = sr_rdseq_offset(sdb, 2 * rrid);
			len = sr_rdseq_length(sdb, 2 * rrid);
			bits2seq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "\t%s", str);
			off = sr_rdseq_offset(sdb, 2 * rrid + 1);
			len = sr_rdseq_length(sdb, 2 * rrid + 1);
			bits2revseq(str, ref_u64list(sdb->rd_seqs, 0), off, len);
			fprintf(gzout, "\t%s\n", str);
		}
	}
	fclose(tmp2);
	unlink(TMP2);
	sr_free_sdb(sdb);
	unlink(prefix);
	for(i=0;i<count_sfv(sfs);i++){ free(ref_sfv(sfs, i)->name); free(ref_sfv(sfs, i)->prefix); }
	free_sfv(sfs);
	if(log != stdout) fclose(log);
	pclose(gzout);
	return 0;
}

thread_begin_def(mpasm);
MyPath *mp;
String *txt;
int asm_ret;
thread_end_def(mpasm);

thread_begin_func(mpasm);
thread_begin_loop(mpasm);
aln_mypath(mpasm->mp);
mpasm->asm_ret = solve_mypath(mpasm->mp);
if(mpasm->asm_ret){
	do_consensus_mypath(mpasm->mp);
}
thread_end_loop(mpasm);
thread_end_func(mpasm);

int main_asm(int argc, char **argv){
	FileReader *fr;
	FILE *fails, *lr_seqs, *lr_stcs;
	int kmer_size, rd_len, min_ol, max_mm, gap, n_cpu, min_ins, max_ins, lmt, c, n, dir, f;
	int output_aln, output_dot, output_cns, simp;
	float min_sm;
	char flag, *seq1, *seq2, *prefix, cmd[256];
	uint32_t pid, len1, len2;
	thread_preprocess(mpasm);
	kmer_size = 3;
	rd_len = 0;
	min_ol = 6;
	min_sm = 0.90;
	max_mm = 10;
	gap    = 1;
	n_cpu  = 1;
	min_ins = 400;
	max_ins = 800;
	lmt = 400;
	f = 1;
	simp = 0;
	output_aln = output_dot = output_cns = 0;
	while((c = getopt(argc, argv, "hfdeg:k:r:l:s:M:t:x:y:")) != -1){
		switch(c){
			case 'k': kmer_size = atoi(optarg); break;
			case 'r': rd_len = atoi(optarg); break;
			case 'l': min_ol = atoi(optarg); break;
			case 'g': gap    = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'M': lmt    = atoi(optarg); break;
			case 't': n_cpu  = atoi(optarg); break;
			case 'x': min_ins = atoi(optarg); break;
			case 'y': max_ins = atoi(optarg); break;
			case 'f': f = 0; break;
			case 'd': output_aln = output_dot = output_cns = 1; break;
			case 'e': simp = 1; break;
			default : return usage_asm();
		}
	}
	if(optind + 1 >= argc) return usage_asm();
	if((fr = fopen_filereader(argv[optind])) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", argv[optind], __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	prefix = argv[optind+1];
	sprintf(cmd, "gzip -1 -c >%s.failed.gz", prefix);
	if((fails = popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot invoke %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "%s.fasta", prefix);
	if((lr_seqs = fopen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	sprintf(cmd, "gzip -1 -c >%s.danguo.gz", prefix);
	if((lr_stcs= popen(cmd, "w")) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
		fflush(stderr);
		return 1;
	}
	thread_begin_init(mpasm, n_cpu);
	mpasm->mp = init_mypath(1, kmer_size, rd_len, min_ol, min_sm, max_mm, (gap > 0), min_ins, max_ins, lmt);
	mpasm->mp->output_aln = output_aln;
	mpasm->mp->output_dot = output_dot;
	mpasm->mp->output_cns = output_cns;
	mpasm->mp->exec_simp  = simp;
	mpasm->txt = init_string(1024);
	mpasm->asm_ret = -1;
	thread_end_init(mpasm);
	n = 0;
	thread_begin_operate(mpasm, 0);
	while((c = fread_table(fr)) != -1){
		if(c < 5) continue;
		flag = get_col_str(fr, 0)[0];
		pid  = atol(get_col_str(fr, 1));
		dir  = (get_col_str(fr, 2)[0] == '-');
		seq1 = get_col_str(fr, 3);
		len1 = get_col_len(fr, 3);
		seq2 = get_col_str(fr, 4);
		len2 = get_col_len(fr, 4);
		if(flag == 'T'){
			thread_wake(mpasm);
			thread_waitfor_one_idle(mpasm);
			n = thread_index(mpasm);
			if(mpasm->asm_ret == 1) output_consensus_mypath(mpasm->mp, lr_seqs, lr_stcs);
			else if(mpasm->asm_ret == 0 && f)  fprintf(fails, "%s", mpasm->txt->string);
			reset_mypath(mpasm->mp);
			if(f){
				clear_string(mpasm->txt);
				append_string(mpasm->txt, fr->line->string, fr->line->size);
				add_char_string(mpasm->txt, '\n');
			}
			if(dir){
				push_mypath(mpasm->mp, 2 * pid + 1, seq1, len1, 0);
				push_mypath(mpasm->mp, 2 * pid + 0, seq2,  len2, 1);
			} else {
				push_mypath(mpasm->mp, 2 * pid + 0, seq1, len1, 0);
				push_mypath(mpasm->mp, 2 * pid + 1, seq2, len2, 1);
			}
		} else if(flag == 'S'){
			if(dir){
				push_mypath(mpasm->mp, 2 * pid + 1, seq1, len1, 0);
				push_mypath(mpasm->mp, 2 * pid + 0, seq2, len2, 1);
			} else {
				push_mypath(mpasm->mp, 2 * pid + 0, seq1, len1, 0);
				push_mypath(mpasm->mp, 2 * pid + 1, seq2, len2, 1);
			}
			if(f){
				append_string(mpasm->txt, fr->line->string, fr->line->size);
				add_char_string(mpasm->txt, '\n');
			}
		}
	}
	fclose_filereader(fr);
	thread_wake(mpasm);
	thread_waitfor_all_idle(mpasm);
	thread_begin_close(mpasm);
	if(mpasm->asm_ret == 1) output_consensus_mypath(mpasm->mp, lr_seqs, lr_stcs);
	else if(mpasm->asm_ret == 0 && f) fprintf(fails, "%s", mpasm->txt->string);
	free_mypath(mpasm->mp);
	free_string(mpasm->txt);
	thread_end_close(mpasm);
	pclose(fails);
	fclose(lr_seqs);
	pclose(lr_stcs);
	return 0;
}

int main(int argc, char **argv){
	if(argc < 2) return usage();
	if(strcasecmp(argv[1], "aln") == 0){
		return main_aln(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "asm") == 0){
		return main_asm(argc - 1, argv + 1);
	} else {
		return usage();
	}
}
