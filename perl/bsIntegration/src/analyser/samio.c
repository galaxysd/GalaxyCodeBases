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
	khiter_t ki, bami;
	kstring_t ks1 = { 0, 0, NULL };
	kstring_t ks2 = { 0, 0, NULL };

	samFile *in, *in2;
	bam_hdr_t *h;
	hts_idx_t *idx;
	bam1_t *b, *d;
	htsFile *out;
	//hts_opt *in_opts = NULL, *out_opts = NULL;
	int r = 0, exit_code = 0;

	for (bami = kh_begin(bamNFOp); bami != kh_end(bamNFOp); ++bami) {
		if (kh_exist(bamNFOp, bami)) {
			BamID = kh_key(bamNFOp, bami);
			pbam = &kh_value(bamNFOp, bami);
#ifdef DEBUGa
			fprintf(stderr, "%u [%s]=%s\t%u %u\n",bami,BamID,pbam->fileName,pbam->insertSize,pbam->SD);
#endif
			in = sam_open(pbam->fileName, "r");
			in2 = sam_open(pbam->fileName, "r");
			if (in == NULL) {
				fprintf(stderr, "[x]Error opening \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			}
			h = sam_hdr_read(in);

			out = hts_open("-", "w");
			if (out == NULL) {
				fprintf(stderr, "[x]Error opening standard output\n");
				return EXIT_FAILURE;
			}
			if (sam_hdr_write(out, h) < 0) {
				fprintf(stderr, "[!]Error writing output header.\n");
				exit_code = 1;
			}

			int8_t *ChrIsHum = malloc(h->n_targets * sizeof(int8_t));
			if (h == NULL) {
				fprintf(stderr, "[x]Couldn't read header for \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			} else {
				for (int32_t i=0; i < h->n_targets; ++i) {
					//ChrIsHum[i] = -1;
					ki = kh_get(chrNFO, chrNFOp, h->target_name[i]);
					if (ki == kh_end(chrNFOp)) {
						errx(4,"[x]Cannot find ChrID for [%s] !",h->target_name[i]);
					} else {
						ChrInfo_t * tmp = &kh_value(chrNFOp, ki);
						ChrIsHum[i] = tmp->isHum;
						//printf(">>> %d Chr:%s %d\n",i,h->target_name[i],ChrIsHum[i]);
					}
				}
			}
			h->ignore_sam_err = 0;
			b = bam_init1();
			d = bam_init1();
			if ((idx = sam_index_load(in, pbam->fileName)) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
				return 1;
			}
			while ((r = sam_read1(in, h, b)) >= 0) {
				int8_t flag = false;
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
					flag = 0;	// recycle
					//kstring_t ks = { 0, 0, NULL };
					if (sam_format1(h, b, &ks1) < 0) {
						fprintf(stderr, "Error writing output.\n");
						exit_code = 1;
						break;
					} else {
						//printf(">[%s]\n",ks_str(&ks1));
						flag |= 1;
					}
					//free(ks.s);
					if (c->mtid < 0) {
						//printf("-[*]\n");
						//printf("~~~ %d ",flag);
						flag &= ~1;
						//printf("%d\n",flag);
					} else if ((c->mtid == c->tid && ChrIsHum[c->tid]) || (ChrIsHum[c->tid] ^ ChrIsHum[c->mtid])) {	// Only grep those mapped on same Human ChrID, or diff species/一方在病毒的情况.
						//char *nrname = h->target_name[c->mtid];
						//int32_t mpos = c->mpos;	// from 0
						uint16_t flag1 = c->flag & (BAM_FREAD1 | BAM_FREAD2);
						hts_itr_t *iter;
						if ((iter = sam_itr_queryi(idx, c->mtid, c->mpos, c->mpos+1)) == 0) {
							fprintf(stderr, "[E::%s] fail to parse region '%s(%d):%d'\n", __func__, h->target_name[c->mtid], c->mtid, c->mpos);
							continue;
						}
						while ((r = sam_itr_next(in2, iter, d)) >= 0) {	// 存在 第一个左端点符合的
							uint16_t flag2 = (d->core).flag & (BAM_FREAD1 | BAM_FREAD2);
							if ( (flag1 + flag2) != (flag1 | flag2) ) continue;	// Read1 with Read2 only. Read12互补
							if (c->mpos != (d->core).pos) continue;
							char *qname2 = bam_get_qname(d);
							if (strcmp(qname,qname2) != 0) continue;	// qname一致
							//if (sam_write1(out, h, b) < 0)
							//kstring_t ks = { 0, 0, NULL };
							if (sam_format1(h, d, &ks2) < 0) {
								fprintf(stderr, "Error writing output.\n");
								exit_code = 1;
								break;
							} else {
								//printf("-[%s]\n",ks_str(&ks2));
								if (flag & 2) fprintf(stdout, "Multiple Read2 found for [%s] !\n",qname);	// can find itself.
								flag |= 2;
							}
							//free(ks.s);
						}
						hts_itr_destroy(iter);
					}
					if (flag == 3) {
						printf(">[%s]\n",ks_str(&ks1));
						printf("-[%s]\n",ks_str(&ks2));
						printf("<--\n");
					}
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
			hts_idx_destroy(idx);
			bam_destroy1(b);
			bam_destroy1(d);
			bam_hdr_destroy(h);
			r = sam_close(in);
			r = sam_close(in2);
			if (r < 0) {
				fprintf(stderr, "Error closing input.\n");
				exit_code = 1;
			}
			free(ChrIsHum);
#ifdef DEBUGa
			fflush(NULL);
			pressAnyKey();
#endif
		}
	}
	free(ks1.s);
	free(ks2.s);
	return exit_code;
}
