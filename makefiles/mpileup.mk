OUTMAIN := mpileup
OUT = $(addsuffix $(shell date +%Y%m%d%Z%H%M%S).bcf,$(OUTMAIN)_)

CMD := samtools mpileup -g -d 1000 -t DP,DPR,DV,DP4,SP -f /bak/seqdata/genomes/Felis_catus_80_masked/Felis_catus80_chr.fa

#BAMS := pti096_clean_aln_pe_rmdup.bam pti183_clean_aln_pe_rmdup.bam pti301_clean_aln_pe_rmdup.bam pti332_clean_aln_pe_rmdup.bam
BAMS := $(sort $(wildcard *.bam))
#CHRS := scaffold1000 scaffold979 scaffold982
RAWCHRS := $(sort $(shell samtools view -H $(firstword $(BAMS)) | sed -n 's/^@SQ\tSN:\([^\t]*\).*/\1/p'))

p+ = $(subst |,+,$1)
+p = $(subst +,|,$1)

CHRS := $(call p+,$(RAWCHRS))

BYCHR = $(addsuffix .bcf,$(addprefix bychr/_,$(CHRS)))
BAIS = $(addsuffix .bai,$(BAMS))
#GATKBAIS = $(BAMS:.bam=.bai)

define newline

endef

NOOP :=
space = $(NOOP) $(NOOP)

.PHONY: all clean cleanall list

all: list $(OUT) bcfbychr.lst

list:
	@echo -e '$(subst $(space),\n,$(BAMS))' >_bams.lst
	@echo -e '$(subst $(space),\n,$(CHRS))' >_chrs.lst
	@echo -e '\n$(subst $(space),\n,$(RAWCHRS))' >>_chrs.lst
	@echo -e '$(subst $(space),\n,$(BYCHR))' >_bcfbychr.lst

bai: $(BAIS)

$(OUT): $(BYCHR) bcfbychr.lst
	bcftools concat -O b -f bcfbychr.lst -o $(OUT)

bychr/:
	mkdir bychr

bcfbychr.lst: $(BYCHR)
	@echo dbg2 [$@] [$*] [$^] [$|]
	echo -e '$(subst $(space),\n,$(BYCHR))' >bcfbychr.lst

%.bai: $*
	samtools index $*

bychr/_%.bcf: $(BAMS) | bychr/
	$(CMD) -r "$(call +p,$(*))" $(BAMS) >$@

clean:
	-rm -fr bychr

cleanall:
	-rm -fr bychr
	-rm bcfbychr.lst
	-rm $(OUT)
