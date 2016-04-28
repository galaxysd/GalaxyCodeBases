OUTMAIN := mpileup
DATESTR := $(shell date +%Y%m%d%Z%H%M%S)
OUTP := $(OUTMAIN)_$(DATESTR)
OUT := $(OUTP).bcf
VCF := $(OUTP).vcf
SAMTOOLS := samtools
BCFTOOLS := bcftools

REF := /bak/seqdata/genomes/Felis_catus_80_masked/Felis_catus80_chr.fa
CMD := $(SAMTOOLS) mpileup -g -d 1000 -t AD,ADF,DP,DPR,DV,DP4,SP -f $(REF)

#BAMS := pti096_clean_aln_pe_rmdup.bam pti183_clean_aln_pe_rmdup.bam pti301_clean_aln_pe_rmdup.bam pti332_clean_aln_pe_rmdup.bam
BAMS := $(sort $(wildcard *.bam))
#CHRS := scaffold1000 scaffold979 scaffold982
RAWCHRS := $(sort $(shell $(SAMTOOLS) view -H $(firstword $(BAMS)) | sed -n 's/^@SQ\tSN:\([^\t]*\).*/\1/p'))

p+ = $(subst |,+,$1)
+p = $(subst +,|,$1)

CHRS := $(call p+,$(RAWCHRS))

BYCHR = $(addsuffix .bcf,$(addprefix bychr/_,$(CHRS)))
#BYCHRINX = $(addsuffix .csi,$(BYCHR))
BAIS = $(addsuffix .bai,$(BAMS))
#GATKBAIS = $(BAMS:.bam=.bai)

define newline

endef

NOOP :=
space = $(NOOP) $(NOOP)

.PHONY: all clean cleanall list

all: list $(OUT) bcfbychr.lst $(VCF)
	$(BCFTOOLS) index $(VCF).gz &
	$(BCFTOOLS) stats -F $(REF) -s - -d 0,1500,1 $(VCF).gz >$(VCF).stats
	$(BCFTOOLS) norm -Df $(REF) -c e -m+both $(VCF).gz | $(BCFTOOLS) filter -sLowQual -e'%QUAL<10' | $(BCFTOOLS) filter -m+ -sDepthHigh -e'DP>650' | $(BCFTOOLS) filter -m+ -sDepthLow -e'DP<2' | $(BCFTOOLS) filter -m+ -sBadSites -e'%QUAL<10 && RPB<0.1' | tee >($(BCFTOOLS) view -Oz -o $(VCF).filtering.gz) | $(BCFTOOLS) view -f .,PASS -Oz -o $(VCF).filtered.gz
	$(BCFTOOLS) index $(VCF).filtered.gz &
	$(BCFTOOLS) index $(VCF).filtering.gz

list:
	@echo -e '$(subst $(space),\n,$(BAMS))' >_bams.lst
	@echo -e '$(subst $(space),\n,$(CHRS))' >_chrs.lst
	@echo -e '\n$(subst $(space),\n,$(RAWCHRS))' >>_chrs.lst
	@echo -e '$(subst $(space),\n,$(BYCHR))' >_bcfbychr.lst

bai: $(BAIS)

$(OUT): $(BYCHR) bcfbychr.lst
	ulimit -n 2048
	$(BCFTOOLS) concat -a -O b -f bcfbychr.lst -o $(OUT)

$(VCF): $(OUT)
	$(BCFTOOLS) call -vmO z -f GQ,GP -o $(@).gz $<

bychr/:
	mkdir bychr

bcfbychr.lst: $(BYCHR)
	@echo dbg2 [$@] [$*] [$^] [$|]
	echo -e '$(subst $(space),\n,$(BYCHR))' >bcfbychr.lst

%.bai: $*
	$(SAMTOOLS) index $*

bychr/_%.bcf: $(BAMS) bai | bychr/
	$(CMD) -r "$(call +p,$(*))" $(BAMS) >$@
	$(BCFTOOLS) index $@

clean:
	-rm -fr bychr

cleanall:
	-rm -fr bychr
	-rm bcfbychr.lst
	-rm $(OUT)
