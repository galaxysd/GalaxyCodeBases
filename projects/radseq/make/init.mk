ifndef INIT_MK
INIT_MK=init.mk

N_PROCESSORS:=$(shell grep '^processor' /proc/cpuinfo | wc -l)
#MAKEOPTS:=-j$(shell echo ${N_PROCESSORS}+1 | bc) -l${N_PROCESSORS}
MAXLOAD:=$(shell echo ${N_PROCESSORS}-0.5 | bc)
MAKEOPTS:=-j${N_PROCESSORS} -l${MAXLOAD}

endif
