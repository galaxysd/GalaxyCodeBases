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
	kstring_t ks3 = { 0, 0, NULL };

	samFile *in;
	bam_hdr_t *h;
	hts_idx_t *idx;
	bam1_t *b, *d, *d2;
	//htsFile *out;
	//hts_opt *in_opts = NULL, *out_opts = NULL;
	int r = 0, exit_code = 0;

	for (bami = kh_begin(bamNFOp); bami != kh_end(bamNFOp); ++bami) {
		//printf(">[%d]:\n",bami);
		if (kh_exist(bamNFOp, bami)) {
			//printf("-[%d]:\n",bami);
			BamID = kh_key(bamNFOp, bami);
			pbam = &kh_value(bamNFOp, bami);
			fprintf(stderr, "%u [%s]=%s\t%u %u\n",bami,BamID,pbam->fileName,pbam->insertSize,pbam->SD);

			in = sam_open(pbam->fileName, "r");
			if (in == NULL) {
				fprintf(stderr, "[x]Error opening \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			}
			h = sam_hdr_read(in);
/*			out = hts_open("-", "w");
			if (out == NULL) {
				fprintf(stderr, "[x]Error opening standard output\n");
				return EXIT_FAILURE;
			}
			if (sam_hdr_write(out, h) < 0) {
				fprintf(stderr, "[!]Error writing output header.\n");
				exit_code = 1;
			}
*/
			int8_t *ChrIsHum;
			if (h == NULL) {
				fprintf(stderr, "[x]Couldn't read header for \"%s\"\n", pbam->fileName);
				return EXIT_FAILURE;
			} else {
				ChrIsHum = malloc(h->n_targets * sizeof(int8_t));
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
			d2 = bam_init1();
			if ((idx = sam_index_load(in, pbam->fileName)) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
				return 1;
			}
			while ((r = sam_read1(in, h, b)) >= 0) {
				int8_t flag = false;
				const bam1_core_t *c = &b->core;
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
					if (getPairedSam(in, idx, b, d) != 0) {
						flag &= ~1;
						continue;
					} else {
						flag |= 2;
						if (c->flag & BAM_FSECONDARY) {
							if (getPairedSam(in, idx, d, d2) == 0) {
								sam_format1(h, d2, &ks3);
								flag |= 4;
							}
						}
					}
/*
对于 BAM_FSECONDARY(256) 的 Read，跳两次 与 读 SA 项，效果一样。
>[sf95_Ref_48245009_48245108_48245208_Vir_-_2000_2044_R_100_90	353	chr2	13996555	0	50S40M	chr18	48245109	0ACACAACAATGTTCCGGAGACTCTAAGGCCTCCCGATACAGAGCAGAGGCCACACACACACACACCATGGAATACTATTCAGCCAAAAAA	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	NM:i:0	MD:Z:40	AS:i:40	XS:i:40	RG:Z:Fsimout_mB	SA:Z:rgi|59585|emb|X04615.1|,2000,-,40S46M4S,60,0;	YC:Z:CT	YD:Z:f]
-[sf95_Ref_48245009_48245108_48245208_Vir_-_2000_2044_R_100_90	177	chr18	48245109	9	40S50M	gi|59585|emb|X04615.1|2000	0	GTTCCGGAGACTCTAAGGCCTCCCGATACAGAGCAGAGGCCACACACACACACACCATGGAATACTATTCAGCCAAAAAAAGGAATTCAA	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	NM:i:0	MD:Z:50	AS:i:50	XS:i:46	RG:Z:Fsimout_mB	SA:Z:rgi|59585|emb|X04615.1|,2000,+,50S40M,9,0;	YC:Z:GA	YD:Z:f]
+[sf95_Ref_48245009_48245108_48245208_Vir_-_2000_2044_R_100_90	113	gi|59585|emb|X04615.1|	2000	60	40S46M4S	chr18	48245109	0	TTTTTTGGCTGAATAGTATTCCATGGTGTGTGTGTGTGTGGCCTCTGCTCTGTATCGGGAGGCCTTAGAGTCTCCGGAACATTGTTGTGT	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	NM:i:0	MD:Z:46	AS:i:46	XS:i:27	RG:Z:Fsimout_mB	SA:Z:fchr2,13996555,+,50S40M,0,0;	YC:Z:CT	YD:Z:r]
*/
					if (sam_format1(h, d, &ks2) < 0) {
						fprintf(stderr, "Error writing output.\n");
						exit_code = 1;
						break;
					}
					if ((flag & 3) == 3) {
						printf(">[%s]\n",ks_str(&ks1));
						printf("-[%s]\n",ks_str(&ks2));
						if (flag & 4) {
							printf("+[%s]\n",ks_str(&ks3));
						}
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
/*			r = sam_close(out);   // stdout can only be closed once
			if (r < 0) {
				fprintf(stderr, "Error closing output.\n");
				exit_code = 1;
			}
*/
			hts_idx_destroy(idx);
			bam_destroy1(b);
			bam_destroy1(d);
			bam_destroy1(d2);
			bam_hdr_destroy(h);
			r = sam_close(in);
			free(ChrIsHum);
#ifdef DEBUGa
			fflush(NULL);
			pressAnyKey();
#endif
		}
		//printf("<[%d]:\n",bami);
	}
	getPairedSam(NULL, NULL, NULL, NULL);	// sam_close(fp2);
	//printf("---[%d]---\n",exit_code);
	ks_release(&ks1);
	ks_release(&ks2);
	ks_release(&ks3);
	return exit_code;
}
