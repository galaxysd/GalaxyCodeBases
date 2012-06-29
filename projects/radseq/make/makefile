#N_PROCESSORS:=$(shell grep '^processor' /proc/cpuinfo | wc -l)
#MAKEOPTS:=-j$(shell echo ${N_PROCESSORS}+1 | bc) -l${N_PROCESSORS}

OS := $(shell uname)
ifeq ($(OS),Linux)
  N_PROCESSORS := $(shell grep -c ^processor /proc/cpuinfo)
else ifeq ($(OS),Darwin)
  N_PROCESSORS := $(shell system_profiler | awk '/Number of CPUs/ {print $$4}{next;}')
endif # $(OS)

MAXLOAD:=$(shell echo ${N_PROCESSORS}-0.5 | bc)
MAKEOPTS := -j$(N_PROCESSORS) -l$(MAXLOAD)

.PHONY: all clean prepare
all: prepare
	$(MAKE) $(MAKEOPTS) -f make/cutadapter.mk > mk.cutadapter.log

clean:
	rm -vfr work

prepare:
	-mkdir lst work
	find /share/users/huxs/work/tiger/fq.p /share/users/huxs/work/tiger/fq.f -name '*.gz'|sort > ./lst/rawfq0.lst
	-mkdir work/radseq work/parents
	./src/rawfq2list.pl ./lst/rawfq0.lst ./lst/rawfq.lst
