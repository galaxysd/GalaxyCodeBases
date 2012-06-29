ifndef INIT_MK
INIT_MK=init.mk

include helpers.mk
include cfg.mk

N_PROCESSORS:=$(shell grep '^processor' /proc/cpuinfo | wc -l)
#MAKEOPTS:=-j$(shell echo ${N_PROCESSORS}+1 | bc) -l${N_PROCESSORS}
MAKEOPTS:=-j${N_PROCESSORS}

#vpath %.fq.gz FQ_PATH
#vpath %.fastq.gz FQ_PATH
FQ_EXTS:=.fq.gz .fastq.gz
FQ_LIST = $(foreach d,$(FQ_PATH),$(wildcard $(addprefix $(d)/*,$(FQ_EXTS))))

endif

