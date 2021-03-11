WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

.PHONY: all

all: $(OUT1)
	@echo "[$(WDLS)] -> [$(OUT1)]"

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@
