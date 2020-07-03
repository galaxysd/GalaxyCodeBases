fCSV = $(wildcard _TabIn/*.csv)
fTXT = $(wildcard _TabIn/*.txt)

report:
	./doAll.r -i $(fCSV) -m $(fTXT)
	@echo =================
	@cat pdf/_QC.txt
	@echo =================

clean:
	-rm *.log
	-rm -fr out

reset: clean
	-rm _TabIn/*
	-rm -fr pdf
