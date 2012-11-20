FQPATH := out

ALNARG:=aln -l 17 -q 10
SAMPEARG:=sampe -a 800
REF:=cat62

FQEXTS:=.fq.gz
BWACMD:=bwa

FQS1 = $(foreach d,$(FQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
FQS2 = $(foreach d,$(FQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
ALLITEMS := $(FQS1:.1.fq.gz=.sam)

.PHONY: all clean

all: $(ALLITEMS)
	@echo all [$@] [$<] [${MAKEOPTS}]

%.sai: %.fq.gz
	@echo "1 $(BWACMD) $(ALNARG) $(REF) $< > $@ 2>$@.log"
	$(BWACMD) $(ALNARG) $(REF) $< > $@ 2>$@.log

%.sam: %.1.sai %.2.sai
	@echo "2 $(BWACMD) $(SAMPEARG) $(REF) $(@:.sam=).1.sai $(@:.sam=).2.sai $(@:.sam=).1.fq.gz $(@:.sam=).2.fq.gz 2>$@.log |gzip -9c >$@.gz"
	$(BWACMD) $(SAMPEARG) $(REF) $(@:.sam=).1.sai $(@:.sam=).2.sai $(@:.sam=).1.fq.gz $(@:.sam=).2.fq.gz 2>$@.log |gzip -9c >$@.gz
	touch $@


