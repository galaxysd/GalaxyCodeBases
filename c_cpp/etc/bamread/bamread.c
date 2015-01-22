#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

int sam_format2(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
	int i;
	uint8_t *s;
	const bam1_core_t *c = &b->core;

	str->l = 0;
	kputsn(bam_get_qname(b), c->l_qname-1, str); kputc('\t', str); // query name
	kputw(c->flag, str); kputc('\t', str); // flag
	if (c->tid >= 0) { // chr
		kputs(h->target_name[c->tid] , str);
		kputc('\t', str);
	} else kputsn("*\t", 2, str);
	kputw(c->pos + 1, str); kputc('\t', str); // pos
	kputw(c->qual, str); kputc('\t', str); // qual
	if (c->n_cigar) { // cigar
		uint32_t *cigar = bam_get_cigar(b);
		for (i = 0; i < c->n_cigar; ++i) {
			kputw(bam_cigar_oplen(cigar[i]), str);
			kputc(bam_cigar_opchr(cigar[i]), str);
		}
	} else kputc('*', str);
	kputc('\t', str);
	if (c->mtid < 0) kputsn("*\t", 2, str); // mate chr
	else if (c->mtid == c->tid) kputsn("=\t", 2, str);
	else {
		kputs(h->target_name[c->mtid], str);
		kputc('\t', str);
	}
	kputw(c->mpos + 1, str); kputc('\t', str); // mate pos
	kputw(c->isize, str); kputc('\t', str); // template len
	if (c->l_qseq) { // seq and qual
		uint8_t *s = bam_get_seq(b);
		for (i = 0; i < c->l_qseq; ++i) kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], str);
		kputc('\t', str);
		s = bam_get_qual(b);
		if (s[0] == 0xff) kputc('*', str);
		else for (i = 0; i < c->l_qseq; ++i) kputc(s[i] + 33, str);
	} else kputsn("*\t*", 3, str);
	s = bam_get_aux(b); // aux
	while (s+4 <= b->data + b->l_data) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s++;
		kputc('\t', str); kputsn((char*)key, 2, str); kputc(':', str);
		if (type == 'A') {
			kputsn("A:", 2, str);
			kputc(*s, str);
			++s;
		} else if (type == 'C') {
			kputsn("i:", 2, str);
			kputw(*s, str);
			++s;
		} else if (type == 'c') {
			kputsn("i:", 2, str);
			kputw(*(int8_t*)s, str);
			++s;
		} else if (type == 'S') {
			if (s+2 <= b->data + b->l_data) {
				kputsn("i:", 2, str);
				kputw(*(uint16_t*)s, str);
				s += 2;
			} else return -1;
		} else if (type == 's') {
			if (s+2 <= b->data + b->l_data) {
				kputsn("i:", 2, str);
				kputw(*(int16_t*)s, str);
				s += 2;
			} else return -1;
		} else if (type == 'I') {
			if (s+4 <= b->data + b->l_data) {
				kputsn("i:", 2, str);
				kputuw(*(uint32_t*)s, str);
				s += 4;
			} else return -1;
		} else if (type == 'i') {
			if (s+4 <= b->data + b->l_data) {
				kputsn("i:", 2, str);
				kputw(*(int32_t*)s, str);
				s += 4;
			} else return -1;
		} else if (type == 'f') {
			if (s+4 <= b->data + b->l_data) {
				ksprintf(str, "f:%g", *(float*)s);
				s += 4;
			} else return -1;

		} else if (type == 'd') {
			if (s+8 <= b->data + b->l_data) {
				ksprintf(str, "d:%g", *(double*)s);
				s += 8;
			} else return -1;
		} else if (type == 'Z' || type == 'H') {
			kputc(type, str); kputc(':', str);
			while (s < b->data + b->l_data && *s) kputc(*s++, str);
			if (s >= b->data + b->l_data)
				return -1;
			++s;
		} else if (type == 'B') {
			uint8_t sub_type = *(s++);
			int32_t n;
			memcpy(&n, s, 4);
			s += 4; // no point to the start of the array
			if (s + n >= b->data + b->l_data)
				return -1;
			kputsn("B:", 2, str); kputc(sub_type, str); // write the typing
			for (i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
				kputc(',', str);
				if ('c' == sub_type)	  { kputw(*(int8_t*)s, str); ++s; }
				else if ('C' == sub_type) { kputw(*(uint8_t*)s, str); ++s; }
				else if ('s' == sub_type) { kputw(*(int16_t*)s, str); s += 2; }
				else if ('S' == sub_type) { kputw(*(uint16_t*)s, str); s += 2; }
				else if ('i' == sub_type) { kputw(*(int32_t*)s, str); s += 4; }
				else if ('I' == sub_type) { kputuw(*(uint32_t*)s, str); s += 4; }
				else if ('f' == sub_type) { ksprintf(str, "%g", *(float*)s); s += 4; }
			}
		}
	}
	return str->l;
}

#define PRINT_OPAQUE_STRUCT(p)  print_mem((p), sizeof(*(p)))
// From http://stackoverflow.com/questions/5349896/print-a-struct-in-c
void print_mem(void const *vp, size_t n) {
	unsigned char const *p = vp;
	for (size_t i=0; i<n; i++) {
		printf("%02x", p[i]);
		if ( (i&7) == 7 ) {
			putchar('\n');
		} else {
			putchar(' ');
		}
	}
	putchar('\n');
};

int main(int argc, char *argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <in.bam>\n",argv[0]);
		return EXIT_FAILURE;
	}
	samFile* fp = sam_open(argv[1], "r");	// hts_open
	if (fp == NULL) { fprintf(stderr, "[%s] fail to open BAM.\n", __func__); return 1; }
	bam_hdr_t* header = sam_hdr_read(fp);
	bam1_t *b = bam_init1();
	kstring_t ks = { 0, 0, NULL };
	int i;
	for (i = 0; i < header->n_targets; ++i) {
		// Print out contig name and length
		printf("%s\t%d\n", header->target_name[i], header->target_len[i]);
	}
	PRINT_OPAQUE_STRUCT(header);
	int r;
	while ((r = sam_read1(fp, header, b)) >= 0) { // read one alignment from `in'
		sam_format2(header, b, &ks);	// The same as sam_format1 in sam.c of htslib.
		fprintf(stdout, "[%s]\n", ks.s);
	}
	if (r < -1) {
		fprintf(stderr, "[!] truncated file.\n");
	}
	PRINT_OPAQUE_STRUCT(b);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	return 0;
}

// ./bamread /share/users/xuxiao/catwork/samtool/4079_n.uni.sort.mge.rmd.bam
