WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

.PHONY: all

all: $(OUT1)
	@echo "[$(WDLS)] -> [$(OUT1)]"
	@echo Run with:\"cromwell run -o options.json [-i inputs.json] workflow.wdl\".

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@
