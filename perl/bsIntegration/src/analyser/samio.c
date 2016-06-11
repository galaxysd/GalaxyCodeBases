#include "functions.h"
#include <htslib/sam.h>

#ifdef DEBUGa
#include "getch.h"
#endif

int do_grep() {
#ifdef DEBUGa
	printf("[!]do_grep\n");
#endif
	BamInfo_t *pbam;
	kh_cstr_t BamID;
	khiter_t ki;

	samFile *in, *in2;
	bam_hdr_t *h;
	hts_idx_t *idx;
	bam1_t *b, *d;
	htsFile *out;
	//hts_opt *in_opts = NULL, *out_opts = NULL;
	int r = 0, exit_code = 0;

	for (ki = kh_begin(bamNFOp); ki != kh_end(bamNFOp); ++ki) {
		if (kh_exist(bamNFOp, ki)) {
			BamID = kh_key(bamNFOp, ki);
			pbam = &kh_value(bamNFOp, ki);
#ifdef DEBUGa
			fprintf(stderr, "%u [%s]=%s\t%u %u\n",ki,BamID,pbam->fileName,pbam->insertSize,pbam->SD);
#endif
			in = sam_open(pbam->fileName, "r");
			in2 = sam_open(pbam->fileName, "r");
			if (in == NULL) {
				fprintf(stderr, "[x]Error opening \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			}
			h = sam_hdr_read(in);
			if (h == NULL) {
				fprintf(stderr, "[x]Couldn't read header for \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			}
			h->ignore_sam_err = 0;
			b = bam_init1();
			d = bam_init1();
			out = hts_open("-", "w");
			if (out == NULL) {
				fprintf(stderr, "[x]Error opening standard output\n");
				return EXIT_FAILURE;
			}
			if (sam_hdr_write(out, h) < 0) {
				fprintf(stderr, "[!]Error writing output header.\n");
				exit_code = 1;
			}
			if ((idx = sam_index_load(in, pbam->fileName)) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
				return 1;
			}
			while ((r = sam_read1(in, h, b)) >= 0) {
				bool flag = false;
				const bam1_core_t *c = &b->core;
				char *qname = bam_get_qname(b);
				if (c->n_cigar) {
					uint32_t *cigar = bam_get_cigar(b);
					for (int i = 0; i < c->n_cigar; ++i) {
						if (bam_cigar_opchr(cigar[i])=='S') {	// soft clipping
							if ( bam_cigar_oplen(cigar[i]) >= myConfig.minGrepSlen ) {
								flag = true;
							}
						}
					}
				}
				if (flag) {
					kstring_t ks = { 0, 0, NULL };
					if (sam_format1(h, b, &ks) < 0) {
						fprintf(stderr, "Error writing output.\n");
						exit_code = 1;
						break;
					} else {
						printf(">[%s]\n",ks_str(&ks));
					}
					free(ks.s);
					if (c->mtid < 0) {
						printf("-[*]\n");
					} else if (c->mtid == c->tid) {	// Only grep those mapped on same ChrID. <---须加上一方在病毒的情况
						//char *nrname = h->target_name[c->mtid];
						//int32_t mpos = c->mpos;	// from 0
						hts_itr_t *iter;
						if ((iter = sam_itr_queryi(idx, c->mtid, c->mpos, c->mpos+1)) == 0) {
							fprintf(stderr, "[E::%s] fail to parse region '%s(%d):%d'\n", __func__, h->target_name[c->mtid], c->mtid, c->mpos);
							continue;
						}
						while ((r = sam_itr_next(in2, iter, d)) >= 0) {
							// 应该只取第一个左端点符合／或qname一致且flag的read12互补的。
							char *qname2 = bam_get_qname(d);
							if (strcmp(qname,qname2) != 0) continue;
							//if (sam_write1(out, h, b) < 0)
							kstring_t ks = { 0, 0, NULL };
							if (sam_format1(h, d, &ks) < 0) {
								fprintf(stderr, "Error writing output.\n");
								exit_code = 1;
								break;
							} else {
								printf("-[%s]\n",ks_str(&ks));
							}
							free(ks.s);
						}
						hts_itr_destroy(iter);
					}
					printf("<--\n");
				}
				/*char *qname = bam_get_qname(b);
				if (sam_write1(out, h, b) < 0) {
					fprintf(stderr, "[x]Error writing output.\n");
					exit_code = 1;
					break;
				}*/
			}
			//r = sam_close(out);	// stdout can only be closed once
			if (r < 0) {
				fprintf(stderr, "Error closing output.\n");
				exit_code = 1;
			}
#ifdef DEBUGa
			fflush(NULL);
			pressAnyKey();
#endif
			bam_destroy1(b);
			bam_destroy1(d);
			bam_hdr_destroy(h);
			r = sam_close(in);
			r = sam_close(in2);
			if (r < 0) {
				fprintf(stderr, "Error closing input.\n");
				exit_code = 1;
			}
		}
	}
	return exit_code;
}
