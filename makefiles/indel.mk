BAMS = $(sort $(wildcard *.recal.bam))
BAMF = $(patsubst %.recal.bam,%,$(BAMS))
OUT1 = $(addsuffix .realignbam,$(BAMF))

REF := /bak/seqdata/genomes/Felis_catus_80_masked/Felis_catus80_chr.fa

.PHONY: all

all: $(OUT1)
	echo [$(BAMS)] -> [$(OUT1)]

%.realignbam: %.recal.bam
	gatk -T IndelRealigner -LOD 1 -model USE_READS -R $(REF) -I $(*).recal.bam -targetIntervals indelzones.intervals -o  $@ 2>$@.log

# gatk -T RealignerTargetCreator -R ../ref/Felis_catus80_chr.fa -I bam0.list -o indelzones.intervals -nt 24 -et STDOUT -K /opt/jar/huxs_pku.edu.cn.key -log indelzones.log &
