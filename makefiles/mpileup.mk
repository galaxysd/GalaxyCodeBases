OUT := test.bcf
CMD := samtools mpileup -g -d 1000 -t DP,DPR,DV,DP4,SP -f /bak/seqdata/genomes/TigerRefGenome/P_tigris.scaffold.fa
BAMS := pti096_clean_aln_pe_rmdup.bam pti183_clean_aln_pe_rmdup.bam pti301_clean_aln_pe_rmdup.bam pti332_clean_aln_pe_rmdup.bam
CHRS := scaffold1000 scaffold979 scaffold982

BYCHR = $(addsuffix .bcf,$(addprefix bychr/_,$(CHRS)))
BAIS = $(addsuffix .bai,$(BAMS))

define newline

endef

NOOP :=
space = $(NOOP) $(NOOP)

.PHONY: all clean cleanall

all: $(OUT) $(BAIS) bcfbychr.lst

$(OUT): $(BYCHR) bcfbychr.lst
	bcftools concat -f bcfbychr.lst -o $(OUT)

bychr/:
	mkdir bychr

bcfbychr.lst: $(BYCHR)
	@echo dbg2 [$^] [$|]
	echo -e '$(subst $(space),\n,$(BYCHR))' >bcfbychr.lst

%.bai: $*
	samtools index $*

bychr/_%.bcf: $(BAMS) bychr/ $(BAIS)
	@echo dbg1 [$@] [$*]
	$(CMD) -r $* $(BAMS) >$@

clean:
	-rm -fr bychr

cleanall:
	-rm -fr bychr
	-rm bcfbychr.lst
	-rm $(OUT)
