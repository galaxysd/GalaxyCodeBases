MDAPATH := mda
MLBACPATH := mlbac
OUTPUTPREFIX := out

ALNARG:=aln -l 17 -q 10
SAMPEARG:=sampe -a 800
SAMSORTARG:=sort -l 9
REF:=ref/hg19p10XYM

INF := $(patsubst %_1.fastq.gz,%,$(wildcard $(MDAPATH)/*_1.fastq.gz $(MLBACPATH)/*_1.fastq.gz))

OUTF := $(addprefix $(OUTPUTPREFIX)/,$(INF))
SAI1 := $(addsuffix _1.sai,$(OUTF))
SAI2 := $(addsuffix _2.sai,$(OUTF))
SAMS := $(addsuffix .sam.gz,$(OUTF))
BAMS := $(addsuffix .bam,$(OUTF))
BAMSORT := $(addsuffix .sort.bam,$(OUTF))
BAMRMDUP := $(addsuffix .rmdup.bam,$(OUTF))


PATHS := $(addprefix $(OUTPUTPREFIX)/,$(MDAPATH) $(MLBACPATH))

NUMPROC := $(shell grep -c ^processor /proc/cpuinfo)
SAICOUNT := $(words $(SAI1))
SAITHREADS := $(shell echo $(NUMPROC)/$(SAICOUNT)/2 |bc)
ifneq ($(SAITHREADS),)
	ifneq ($(SAITHREADS),0)
		ifneq ($(SAITHREADS),1)
			ALNARG += -t $(SAITHREADS)
		 endif
	endif
endif
FREEMEM := $(shell free -m|grep -e '-/+ buffers/cache'|awk '{print $$NF}')
SAMMEM := $(shell echo $(FREEMEM)/$(SAICOUNT) |bc)
ifeq ($(SAMMEM),)
	SAMMEM := 768
endif

all: $(OUTF)
	@echo "[$(OUTF)]" "[$(INF)]" "[$@]" "[$(SAMMEM) $(FREEMEM)]"
	date > $(OUTPUTPREFIX)/_alldone.log

$(PATHS):
	mkdir -p $@

$(OUTF): $(PATHS) $(BAMRMDUP)
	$(eval IN := $(patsubst $(OUTPUTPREFIX)/%,%,$@))
	echo "$@" "$<" "$(IN)"
	touch $@

$(SAI1) $(SAI2): $(PATHS)
	$(eval IN := $(patsubst $(OUTPUTPREFIX)/%,%,$(@:%.sai=%)))
	bwa $(ALNARG) $(REF) $(IN).fastq.gz >$@ 2>$@.log
	@date >>$@.log
	@echo done. >>$@.log

$(SAMS): $(SAI1) $(SAI2)
	$(eval IN := $(@:%.sam.gz=%))
	$(eval FQ := $(patsubst $(OUTPUTPREFIX)/%,%,$(IN)))
	bwa $(SAMPEARG) $(REF) $(IN)_1.sai $(IN)_2.sai $(FQ)_1.fastq.gz $(FQ)_2.fastq.gz·2>$(IN).sam.log |gzip -9c >$@
	@date >>$(IN).sam.log
	@echo done. >>$(IN).sam.log

$(BAMS): $(SAMS)
	samtools view -bS $< >$@ 2>$@.log
	@date >>$@.log
	@echo "done (sam.gz to bam)">>$@.log

$(BAMSORT): $(BAMS)
	if [ "$(SAMMEM)" -gt "768" ]; then \
		SAMSORTMEM="-m $(SAMMEM)";\
	else \
		SAMSORTMEM=;\
	fi;\
	samtools $(SAMSORTARG) $$SAMSORTMEM $< $(@:%.bam=%) >>$<.log
	@date >>$<.log
	@echo "done (bam sort)">>$<.log

$(BAMRMDUP): $(BAMSORT)
	$(eval LOG := $(@:%.rmdup.bam=%.bam.log))
	samtools rmdup $< $@ >>$(LOG)
	@date >>$(LOG)
	@echo "done (bam rmdup)">>$(LOG)

help:
	@echo "Usage: make -j |tee make.log"

clean:
	rm -fr out/*

.PHONY: help clean
