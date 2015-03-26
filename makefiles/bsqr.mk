BAMS = $(sort $(wildcard *.bam))
BAMF = $(patsubst %.bam,%,$(BAMS))
OUT1 = $(addsuffix .recal.tbl,$(BAMF))
OUT2 = $(addsuffix .post.tbl,$(BAMF))
OUT3 = $(addsuffix .recal.tbl,$(BAMF))
OUT5 = $(addsuffix .cmp.pdf,$(BAMF))
OUT4 = $(addsuffix .recal.bam,$(BAMF))

REF := /bak/seqdata/genomes/Felis_catus_80_masked/Felis_catus80_chr.fa

.PHONY: all

all: $(OUT1) $(OUT2) $(OUT3) $(OUT4) $(OUT5)
	echo $(BAMS)

%.recal.tbl: %.bam
	gatk -T BaseRecalibrator -nct 8 -mte -R $(REF) -I $(*).bam -knownSites filtered2.vcf -o $@ 2>$@.log

%.post.tbl: %.recal.tbl %.bam
	gatk -T BaseRecalibrator -nct 8 -mte -R $(REF) -I $(*).bam -knownSites filtered2.vcf -BQSR $(@:.post.tbl=.recal.tbl) -o $@ 2>$@.log

%.cmp.pdf: %.recal.tbl %.post.tbl
	gatk -T AnalyzeCovariates -R $(REF) -before $(@:.cmp.pdf=.recal.tbl) -after $(@:.cmp.pdf=.post.tbl) -plots $@ 2>$@.log

%.recal.bam: %.recal.tbl %.bam
	gatk -T PrintReads -R $(REF) -I $(*).bam -BQSR $(@:.recal.bam=.recal.tbl) -o $@ 2>$@.log
