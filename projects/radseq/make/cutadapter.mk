RADSEQFQPATH := work/radseq
WGSFQPATH := work/parents

ADAPTER13:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER15:=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
ADAPTER23:=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
ADAPTER25:=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

CUTADAPTARG:=-e 0.08 -n 2 -m 0 -O 5
CUTADRADARG:=-e 0.08 -n 1 -m 0 -O 5

FQEXTS:=.fq.gz
CUTADAPTCMD:=cutadapt

RADSEQFQS1 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
RADSEQFQS2 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
WGSFQS1 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
WGSFQS2 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
COMMONFQS := $(RADSEQFQS1) $(RADSEQFQS2) $(WGSFQS1) $(WGSFQS2)

.PHONY: all clean

#all: $(RADSEQFQS1:.1.fq.gz=.rad) $(COMMONFQS:.fq.gz=)
all: $(COMMONFQS:.fq.gz=)
	@echo all [$@] [$<] [${MAKEOPTS}]

%.1: %.1.fq.gz
#	@echo "1 $(CUTADAPTCMD) -a $(ADAPTER13) -g $(ADAPTER15) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz $< > $@.cut.log"
	$(CUTADAPTCMD) -a $(ADAPTER13) -g $(ADAPTER15) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz $< > $@.cut.log
	touch $@

%.2: %.2.fq.gz
#	@echo "2 $(CUTADAPTCMD) -a $(ADAPTER23) -g $(ADAPTER25) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz $< > $@.cut.log"
	$(CUTADAPTCMD) -a $(ADAPTER23) -g $(ADAPTER25) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz $< > $@.cut.log
	touch $@

#%.rad: %.1.fq.gz
#	@echo "RAD1 $(CUTADAPTCMD) -a $(ADAPTER13) $(CUTADRADARG) -o $(@:.rad=).cut1.gz -r $(@:.rad=).rest1.gz $< > $(@:.rad=).cut1.log"
#	@echo perl $(@:.rad=).cut1.gz $(@:.rad=).cutM.gz
#	@echo "RAD2 $(CUTADAPTCMD) -g $(ADAPTER15) $(CUTADRADARG) -o $(@:.rad=).cut2.gz -r $(@:.rad=).rest2.gz $(@:.rad=).cutM.gz > $(@:.rad=).cut2.log"
#	$(CUTADAPTCMD) -a $(ADAPTER13) $(CUTADRADARG) -o $(@:.rad=).cut1.gz -r $(@:.rad=).rest1.gz $< > $(@:.rad=).cut1.log

clean:
	-rm $(COMMONFQS:.fq.gz=) $(COMMONFQS:.fq.gz=.rest.gz) $(COMMONFQS:.fq.gz=.cut.log)
	-rm $(COMMONFQS:.fq.gz=.cut.gz)
