SAMPLEID := mbrain mkidney
LANEID_Illumina := L001
LANEID_Salus := L001 L002 L003 L004
PLATFORMS := Illumina Salus

SRBIN := /share/result/spatial/test_huxs/10x/space-ranger/spaceranger-2.1.1/bin/spaceranger

define platform_rule
TAGS_$1 := $$(foreach s, $$(SAMPLEID), $1/$$(s)_S1_L001)
endef
$(foreach platform,$(PLATFORMS),$(eval $(call platform_rule,$(platform))))
ALLTAGS := $(TAGS_Illumina) $(TAGS_Salus)

RAW_Salus := $(foreach a, $(SAMPLEID), $(foreach b, $(LANEID_Salus), $(foreach r,R1 R2,raw4_Salus/$(a)_S1_$(b)_$(r)_001.fastq.gz )))
MERGE_Salus := $(foreach tag,$(TAGS_Salus),raw_$(tag) )
MERGE_Illumina := $(foreach a, $(SAMPLEID), $(foreach r,R1 R2,raw_Illumina/$(a)_S1_L001 ))
FASTP_TAGS := $(foreach tag,$(ALLTAGS),fastp_$(tag) )
SUBFQ_TAGS := $(foreach tag,$(ALLTAGS),sub_$(tag) )
FQ_DIRS := $(foreach p,raw fastp sub,$(foreach pf,$(PLATFORMS),$(p)_$(pf)))
SRLOGS := $(foreach smp, $(SAMPLEID), $(foreach pf,$(PLATFORMS),sr_$(pf)/sub_$(pf)_$(smp).log ))

#all:	$(MERGE_Salus) $(FASTP_TAGS) $(foreach d,$(FQ_DIRS),$(d).stats) $(SRLOGS)
all:	$(foreach d,$(FQ_DIRS),$(d).stats) $(SRLOGS)
	@echo $^

define smp_rule
raw_Salus/$1_S1_L001: $$(filter raw4_Salus/$1_S1_%,$$(RAW_Salus))
	mkdir -p $$(dir $$@)
	@cat $$(filter %_R1_001.fastq.gz,$$^) > $$@_R1_001.fastq.gz
	@cat $$(filter %_R2_001.fastq.gz,$$^) > $$@_R2_001.fastq.gz
	touch $$@
raw_Illumina/$1_S1_L001:
	touch $$@
sub_Illumina/$1_S1_L001: $$(filter fastp_Illumina/$1_S1_L001,$$(FASTP_TAGS))
	mkdir -p $$(dir $$@)
	cp --reflink=auto -l $$^_R1_001.cut.gz $$@_R1_001.fastq.gz
	cp --reflink=auto -l $$^_R2_001.cut.gz $$@_R2_001.fastq.gz
	touch $$@
endef
$(foreach sample, $(SAMPLEID), $(eval $(call smp_rule,$(sample))))

define pf_smp_rule
fastp_$1/$2_S1_L001: $$(filter raw_$1/$2_S1_L001,$$(MERGE_Salus) $$(MERGE_Illumina))
	mkdir -p $$(dir $$@)
	micromamba run -n salus fastp --thread 4 -z 4 -A --max_len1 28 --max_len2 0 --dont_eval_duplication -q 20 -u 82 -n 4 --average_qual 20 --length_required 28 -y -Y 30 -g -x \
    -i $$<_R1_001.fastq.gz -I $$<_R2_001.fastq.gz -o $$@_R1_001.fastp.gz -O $$@_R2_001.fastp.gz -j $$@.fastp.json -h $$@.fastp.html 2>$$@.fastp.log
	micromamba run -n salus cutadapt -Z -j8 --pair-filter=any --minimum-length 28:50 --length 100 --times 0 --report full --json=$$@.cut.json \
	  -o $$@_R1_001.cut.gz -p $$@_R2_001.cut.gz $$@_R1_001.fastp.gz $$@_R2_001.fastp.gz >$$@.cut.log
	touch $$@
endef
$(foreach pf, $(PLATFORMS), $(foreach sample, $(SAMPLEID), $(eval $(call pf_smp_rule,$(pf),$(sample)))) )

define pf_rule
sr_$(1)/sub_$(1)_mbrain.log: $$(filter sub_$1/mbrain_S1_L001,$$(SUBFQ_TAGS))
	mkdir -p $$(dir $$@)
	micromamba run -n salus $$(SRBIN) count --id=srp_$(1)_mbrain --transcriptome=refdata-gex-mm10-2020-A --fastqs=sub_$1 --sample=mbrain \
    --slide=V12A25-031 --slidefile=V12A25-031.gpr --localcores=36 --localmem=64 --vdrmode=disable \
    --area=A1 --image=mbrain.tif >$$@
sr_$(1)/sub_$(1)_mkidney.log: $$(filter sub_$1/mkidney_S1_L001,$$(SUBFQ_TAGS))
	mkdir -p $$(dir $$@)
	micromamba run -n salus $$(SRBIN) count --id=srp_$(1)_mkidney --transcriptome=refdata-gex-mm10-2020-A --fastqs=sub_$1 --sample=mkidney \
    --slide=V12A25-031 --slidefile=V12A25-031.gpr --localcores=36 --localmem=64 --vdrmode=disable \
    --area=A1 --image=mkidney.tif >$$@
endef
$(foreach pf, $(PLATFORMS), $(eval $(call pf_rule,$(pf))))

sub_Salus/mbrain_S1_L001: $(filter fastp_Salus/mbrain_S1_L001,$(FASTP_TAGS))
	mkdir -p $(dir $@)
	micromamba run -n salus seqkit --compress-level 4 -o $@_R1_001.fastq.gz -w 200 -j 1 sample -p 0.72 -s 123456 $^_R1_001.cut.gz
	micromamba run -n salus seqkit --compress-level 4 -o $@_R2_001.fastq.gz -w 200 -j 1 sample -p 0.72 -s 123456 $^_R2_001.cut.gz
	touch $@

sub_Salus/mkidney_S1_L001: $(filter fastp_Salus/mkidney_S1_L001,$(FASTP_TAGS))
	mkdir -p $(dir $@)
	micromamba run -n salus seqkit --compress-level 4 -o $@_R1_001.fastq.gz -w 200 -j 1 sample -p 0.76 -s 123456 $^_R1_001.cut.gz
	micromamba run -n salus seqkit --compress-level 4 -o $@_R2_001.fastq.gz -w 200 -j 1 sample -p 0.76 -s 123456 $^_R2_001.cut.gz
	touch $@

%.stats: $(foreach sample, $(SAMPLEID), %/$(sample)_S1_L001)
	micromamba run -n salus seqkit stats -a -o $@ -j8 $(dir $<)*.gz

.PHONY scp:
	scp makefile s0:/share/result/spatial/data/BoAo_sp/sub/
