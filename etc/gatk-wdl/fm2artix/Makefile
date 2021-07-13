WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

all: $(OUT1) fm2test.samples.json
	@echo "[$(WDLS)] -> [$(OUT1)]"
	@echo Run with:{JAVA_OPTS=\"-Dbackend.providers.Local.config.concurrent-job-limit=16\" cromwell run -o cromwell_options_no_user.json [-i inputs.json] workflow.wdl}.

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@

fm2test.samples.json: fm2test.csv
	biowdl-input-converter --skip-file-check $< |python3 -mjson.tool > $@

release: pkg fm2test.samples.json
	tar -cLzvf release.tgz fm2.csv fm2test.samples.json dockerImages.yml bin/ fmtest.json fm2.wdl cromwell_options_no_user.json fm2wdl.zip
	@echo [!]Use \"/share/FGI2017B/prefix/jdk/jdk1.8.0_281/bin/java\" in \"/share/FGI2017B/brew/bin/picard\".
	@echo [!]Remember to execute \"ln -s /share/FGI2017B/users/huxs/2020Q4_fm2/wdl/fq\" for FASTQ files.
	@ls -lh release.tgz

.PHONY pkg:
	-rm fm2wdl.zip
	zip -9r fm2wdl.zip fm2.wdl  tasks/
