RADSEQFQPATH := fqrad
WGSFQPATH := fqwgs

ADAPTER13:=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC
ADAPTER15:=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
ADAPTER23:=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
ADAPTER25:=GATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT

CUTADAPTARG:=-e 0.08 -n 2 -m 0 -O 5
CUTADRADARG:=-e 0.08 -n 1 -m 0 -O 5

FQEXTS:=.fq.gz
CUTADAPTCMD:=cutadapt

RADSEQFQS1 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
RADSEQFQS2 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
WGSFQS1 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
WGSFQS2 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
COMMONFQS := $(RADSEQFQS2) $(WGSFQS1) $(WGSFQS2)
ALLITEMS := $(RADSEQFQS1:.1.fq.gz=.rad) $(COMMONFQS:.fq.gz=)

.PHONY: all clean

all: $(ALLITEMS)
	@echo all [$@] [$<] [${MAKEOPTS}]

%.1: %.1.fq.gz
	$(CUTADAPTCMD) -a $(ADAPTER13) -g $(ADAPTER15) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz --info-file=$@.nfo.gz $< > $@.cut.log
	touch $@

%.2: %.2.fq.gz
	$(CUTADAPTCMD) -a $(ADAPTER23) -g $(ADAPTER25) $(CUTADAPTARG) -o $@.cut.gz -r $@.rest.gz --info-file=$@.nfo.gz $< > $@.cut.log
	touch $@

%.rad: %.1.fq.gz
	$(CUTADAPTCMD) -a $(ADAPTER13) -g $(ADAPTER15) $(CUTADAPTARG) -o $(@:.rad=.1.cut.gz) -r $(@:.rad=.1.rest.gz) --info-file=$(@:.rad=.1.nfo.gz) $< > $(@:.rad=.1.cut.log)
	touch $@

clean:
	-rm $(ALLITEMS) $(COMMONFQS:.fq.gz=.rest.gz) $(COMMONFQS:.fq.gz=.cut.log)
	-rm $(COMMONFQS:.fq.gz=.cut.gz)
