ifndef STAGE0_MK
STAGE0_MK=stage0.mk

include init.mk

stage0: tmp/stage0
tmp/stage0: tmp/_stage0-Samplelst \
	tmp/_stage0-cutadapter

tmp/_stage0-Samplelst:
	mkdir -p tmp
	find ${FQ_PATH} -name '*.gz'|sort > tmp/rawfq.lst
	./utils/rawfq2list.pl tmp/rawfq.lst Sample.lst
	touch $@

tmp/_stage0-cutadapter:
	touch $@


final: stage0

endif
