RADSEQFQPATH := work/radseq
WGSFQPATH := work/parents

ADAPTER13:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
ADAPTER15:=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
ADAPTER23:=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
ADAPTER25:=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

CUTADAPTARG:=-e 0.08 -n 2 -m 0 -O 5

FQEXTS:=.fq.gz
CUTADAPTCMD:=cutadapt

RADSEQFQS1 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
RADSEQFQS2 = $(foreach d,$(RADSEQFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
WGSFQS1 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.1$(FQEXTS))))
WGSFQS2 = $(foreach d,$(WGSFQPATH),$(wildcard $(addprefix $(d)/*,.2$(FQEXTS))))
COMMONFQS := $(RADSEQFQS2) $(WGSFQS1) $(WGSFQS2)

.PHONY: all clean

all: $(RADSEQFQS1:.1.fq.gz=.rad) $(COMMONFQS:.fq.gz=)
	@echo all [$@] [$<]

%.1: %.1.fq.gz
	@echo 1 $(CUTADAPTCMD) [$@] [$<]

%.2: %.2.fq.gz
	@echo 2 $(CUTADAPTCMD) [$@] [$<]

%.rad: %.1.fq.gz
	@echo RAD1 $(CUTADAPTCMD) [$@] [$<]
