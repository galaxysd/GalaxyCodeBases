INCLIST := FullSomaticPipeline

OUTF := $(addsuffix .imports.zip,$(INCLIST))

all: $(OUTF)

$(OUTF): $(OUTF:.imports.zip=)
	@cd $(OUTF:.imports.zip=);echo ------[$@]------; zip -9uvoX ../$@ *.wdl

clean:
	rm $(OUTF)

.PHONY: clean
