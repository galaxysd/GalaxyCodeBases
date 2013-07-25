INPUTPATH := fqrad
OUTPUTPREFIX := out

INF := $(patsubst %.1.cut.gz,%,$(foreach ONEINPUTPATH,$(INPUTPATH),$(wildcard $(ONEINPUTPATH)/*.1.cut.gz)))

OUTF := $(addprefix $(OUTPUTPREFIX)/,$(INF))
SEP1 := $(addsuffix _000000,$(OUTF))
SEP2 := $(addsuffix _210210,$(OUTF))
SEP3 := $(addsuffix _221210,$(OUTF))

PATHS := $(addprefix $(OUTPUTPREFIX)/,$(INPUTPATH))

NEEDED_COMMANDS := bc grep mkdir free gzip ./src/separatebytag.pl

all: $(SEP1) $(SEP2) $(SEP3)
	@echo "[$(OUTF)]" "[$(INF)]" "[$@]" "[$(SAMMEM) $(FREEMEM)]"
	date > $(OUTPUTPREFIX)/_alldone.log

check:
	@for thecmd in $(NEEDED_COMMANDS); do \
	if ! command -v "$${thecmd%% *}" >/dev/null 2>&1; then \
			checkok="0"; \
			echo "[x]'$${thecmd%% *}' not found."; \
		fi; \
	done; \
	if [ "$${checkok}" == "0" ]; then \
		echo "[!]Please install missing cmd(s) above."; \
		exit -1; \
	fi;

$(PATHS): check
	mkdir -p $@

$(SEP1): $(PATHS)
	$(eval IN := $(patsubst $(OUTPUTPREFIX)/%_000000,%,$@))
	echo "$@" "$<" "$(IN)"
	./src/separatebytag.pl tigrad.lst 000000 $(IN).1.cut.gz $(IN).2.cut.gz $(@)
	touch $@

$(SEP2): $(SEP1)
	$(eval IN := $(patsubst %_210210,%,$@))
	./src/separatebytag.pl tigrad.lst 210210 $(IN)_000000.NA.Unknown.1.fq.gz $(IN)_000000.NA.Unknown.2.fq.gz $(@)
	touch $@

$(SEP3): $(SEP2)
	$(eval IN := $(patsubst %_221210,%,$@))
	./src/separatebytag.pl tigrad.lst 221210 $(IN)_210210.NA.Unknown.1.fq.gz $(IN)_210210.NA.Unknown.2.fq.gz $(@)
	touch $@
