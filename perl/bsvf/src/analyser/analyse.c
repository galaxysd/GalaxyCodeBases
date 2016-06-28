#include "functions.h"

int do_analyse() {
#ifdef DEBUGa
	printf("[!]do_analyse\n");
#endif
	kstring_t kstr = { 0, 0, NULL };
	khiter_t ki, bami;
	BamInfo_t *pbam;
	kh_cstr_t BamID;
	for (bami = kh_begin(bamNFOp); bami != kh_end(bamNFOp); ++bami) {
		if (kh_exist(bamNFOp, bami)) {
			BamID = kh_key(bamNFOp, bami);
			pbam = &kh_value(bamNFOp, bami);
			ksprintf(&kstr, "%s.grep", pbam->fileName);
			fprintf(stderr, "%u [%s]=%s\t%u %u\n",bami,BamID,kstr.s,pbam->insertSize,pbam->SD);
			// 一作说用旧的每次exec的方案就行，嘛，既然大家都不关心效率，咱还坚持啥。反正CPU够快、Chobits也不会喊累。
		}
	}
	return 0;
}
