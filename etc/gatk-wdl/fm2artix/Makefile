WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

all: $(OUT1) fm2test.samples.json
	@echo "[$(WDLS)] -> [$(OUT1)]"
	@echo Run with:{JAVA_OPTS=\"-Dbackend.providers.Local.config.concurrent-job-limit=16\" cromwell run -o cromwell_options_no_user.json [-i inputs.json] workflow.wdl}.

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@

example/fm2test.samples.json: example/fm2test.csv
	biowdl-input-converter --skip-file-check $< |python3 -mjson.tool > $@

.PHONY pkg:
	-rm fm2wdl.zip fm2pkg.tgz
	zip -9r fm2wdl.zip tasks/
	tar -czvf fm2pkg.tgz fm2.wdl fm2wdl.zip bin report db example

p0FILE := "report/style.html"
p1FILE := "db/annots.bed.gz"
.PHONY etc:
	sed 's/.markdown-body //g' report/github-markdown.css >.build/github-markdownA.css
	css-html-js-minify .build/github-markdownA.css
	echo "<!DOCTYPE html><html><head><style>" >$(p0FILE)
	cat .build/github-markdownA.min.css >>$(p0FILE)
	echo "</style></head><body>" >>$(p0FILE)
	perl bin/MkBed.pl db/snplst.txt db/annots.bed
