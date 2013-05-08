INPUTPATH := mda mlbac
OUTPUTPREFIX := out

ALNARG:=aln -l 17 -q 10
SAMPEARG:=sampe -a 800
SAMSORTARG:=sort -l 9
REF:=ref/hg19p10XYM

INF := $(patsubst %_1.fastq.gz,%,$(foreach ONEINPUTPATH,$(INPUTPATH),$(wildcard $(ONEINPUTPATH)/*_1.fastq.gz)))

OUTF := $(addprefix $(OUTPUTPREFIX)/,$(INF))
SAI1 := $(addsuffix _1.sai,$(OUTF))
SAI2 := $(addsuffix _2.sai,$(OUTF))
SAMS := $(addsuffix .sam.gz,$(OUTF))
BAMS := $(addsuffix .bam,$(OUTF))
BAMSORT := $(addsuffix .sort.bam,$(OUTF))
BAMRMDUP := $(addsuffix .rmdup.bam,$(OUTF))

PATHS := $(addprefix $(OUTPUTPREFIX)/,$(MDAPATH) $(MLBACPATH))

NUMPROC := $(shell grep -c ^processor /proc/cpuinfo)
SAMCOUNT := $(words $(SAI1))
SAITHREADS := $(shell echo 1.5*$(NUMPROC)/$(SAMCOUNT) |bc)
ifneq ($(SAITHREADS),)
	ifneq ($(SAITHREADS),0)
		ifneq ($(SAITHREADS),1)
			ALNARG += -t $(SAITHREADS)
		 endif
	endif
endif
FREEMEM := $(shell free -m|grep -e '-/+ buffers/cache'|awk '{print $$NF}')
SAMMEM := $(shell echo $(FREEMEM)/$(SAMCOUNT)/$(SAITHREADS) |bc)

NEEDED_COMMANDS := bc bwa samtools grep mkdir free gzip

all: $(OUTF)
	@echo "[$(OUTF)]" "[$(INF)]" "[$@]" "[$(SAMMEM) $(FREEMEM)]"
	date > $(OUTPUTPREFIX)/_alldone.log

check:
	@for thecmd in $(NEEDED_COMMANDS); do \
	if ! command -v "$${thecmd%% *}" >/dev/null 2>&1; then \
			checkok="0"; \
			echo "[x]'$${thecmd%% *}' not found."; \
		fi; \
	done; \
	if [ "$${checkok}" == "0" ]; then \
		echo "[!]Please install missing cmd(s) above."; \
		exit -1; \
	fi;

$(PATHS): check
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
		SAMSORTMEM="-@ $(SAITHREADS) -m $(SAMMEM)M";\
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

help: check
	@if [ "$(NUMPROC)" -gt "$(SAMCOUNT)"  ]; then \
		JOBCNT="$(SAMCOUNT)";\
	fi;\
	echo -e "Usage: make -j $(SAMCOUNT) |tee make.log\n\n$(NUMPROC) core(s) & $(FREEMEM) mb free memory found now.\nBWA will run sai with [$(ALNARG)] for your $(SAMCOUNT) bam file(s)."
	@if [ "$(SAMMEM)" -gt "768" ]; then \
		echo -e "SAMTOOLS will sort with [-@ $(SAITHREADS) -m $(SAMMEM)M], if free memory remains same when you type make."; \
	fi

clean:
	rm -fr out/*

.PHONY: help clean check
