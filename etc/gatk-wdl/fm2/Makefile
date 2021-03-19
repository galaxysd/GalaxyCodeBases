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

.PHONY release: fm2.samples.json
	tar -cLzvf release.tgz fm2.csv fm2.samples.json fm2.json fm2.wdl cromwell_options_no_user.json run.sh tasks bin ref/ForensicM.fa.* dockerImages.yml
	@echo [!]Use \"/share/FGI2017B/prefix/jdk/jdk1.8.0_281/bin/java\" in \"/share/FGI2017B/brew/bin/picard\".
	@echo [!]Remember to execute \"ln -s /share/FGI2017B/users/huxs/2020Q4_fm2/wdl/fq\" for FASTQ files.
	@ls -lh release.tgz
