WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

.PHONY: all

all: $(OUT1) fm2.samples.json
	@echo "[$(WDLS)] -> [$(OUT1)]"
	@echo Run with:{JAVA_OPTS=\"-Dbackend.providers.Local.config.concurrent-job-limit=16\" cromwell run -o cromwell_options_no_user.json [-i inputs.json] workflow.wdl}.

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@

fm2.samples.json: fm2.csv
	biowdl-input-converter --skip-file-check fm2.csv |python3 -mjson.tool >fm2.samples.json
