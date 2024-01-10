SAMPLEID := mbrain mkidney
LANEID_Illumina := L001
LANEID_Salus := L001 L002 L003 L004
PLATFORMS := Illumina Salus

define platform_rule
TAGS_$1 := $$(foreach s, $$(SAMPLEID), $1/$$(s)_S1_L001)
endef
$(foreach platform,$(PLATFORMS),$(eval $(call platform_rule,$(platform))))
ALLTAGS := $(TAGS_Illumina) $(TAGS_Salus)

RAW_Salus := $(foreach a, $(SAMPLEID), $(foreach b, $(LANEID_Salus), $(foreach r,R1 R2,raw4_Salus/$(a)_S1_$(b)_$(r)_001.fastq.gz )))
MERGE_Salus := $(foreach tag,$(TAGS_Salus),raw_$(tag) )
MERGE_Illumina := $(foreach a, $(SAMPLEID), $(foreach r,R1 R2,raw_Illumina/$(a)_S1_L001 ))
FASTP_TAGS := $(foreach tag,$(ALLTAGS),fastp_$(tag) )
FQ_DIRS := $(foreach p,raw fastp,$(foreach pf,$(PLATFORMS),$(p)_$(pf)))

all:	$(MERGE_Salus) $(FASTP_TAGS) $(foreach d,$(FQ_DIRS),$(d).stats)
	@echo $^

define merge_rule
raw_Salus/$1_S1_L001: $$(filter raw4_Salus/$1_S1_%,$$(RAW_Salus))
	mkdir -p $$(dir $$@)
	@cat $$(filter %_R1_001.fastq.gz,$$^) > $$@_R1_001.fastq.gz
	@cat $$(filter %_R2_001.fastq.gz,$$^) > $$@_R2_001.fastq.gz
	touch $$@
raw_Illumina/$1_S1_L001:
	touch $$@
endef
$(foreach sample, $(SAMPLEID), $(eval $(call merge_rule,$(sample))))

define pf_smp_rule
fastp_$1/$2_S1_L001: $$(filter raw_$1/$2_S1_L001,$$(MERGE_Salus) $$(MERGE_Illumina))
	mkdir -p $$(dir $$@)
	micromamba run -n salus fastp --thread 4 -z 4 -A --max_len1 28 --max_len2 0 --dont_eval_duplication -q 20 -u 82 -n 4 --average_qual 20 --length_required 28 -y -Y 30 -g -x \
	        -i $$<_R1_001.fastq.gz -I $$<_R2_001.fastq.gz \
	        -o $$@_R1_001.fastq.gz -O $$@_R2_001.fastq.gz \
	        -j $$@.json -h $$@.html 2>$$@.log
endef
$(foreach pf, $(PLATFORMS), $(foreach sample, $(SAMPLEID), $(eval $(call pf_smp_rule,$(pf),$(sample)))) )

define pf_rule
endef
$(foreach pf, $(PLATFORMS), $(eval $(call pf_rule,$(pf))))

%.stats: %
	seqkit stats -a -o $@ -j8 $^/*.fastq.gz

.PHONY scp:
	scp makefile s0:/share/result/spatial/data/BoAo_sp/sub/
