WDLS = $(sort $(wildcard *.wdl))
WDLF = $(patsubst %.wdl,%,$(WDLS))
OUT1 = $(addsuffix .json.wom,$(WDLF))

.PHONY: all release dbsnp dbsnp/GCF_000001405.38.gz

all: $(OUT1) fm2.samples.json
	@echo "[$(WDLS)] -> [$(OUT1)]"
	@echo Run with:{JAVA_OPTS=\"-Dbackend.providers.Local.config.concurrent-job-limit=16\" cromwell run -o cromwell_options_no_user.json [-i inputs.json] workflow.wdl}.

%.json.wom: %.wdl
	womtool validate $(*).wdl && womtool inputs $(*).wdl >$@

fm2.samples.json: fm2.csv
	biowdl-input-converter --skip-file-check fm2.csv |python3 -mjson.tool >fm2.samples.json

release: fm2.samples.json
	tar -cLzvf release.tgz fm2.csv fm2.samples.json fm2.json fm2.wdl cromwell_options_no_user.json run.sh tasks bin ref/ForensicM.fa.* dockerImages.yml
	@echo [!]Use \"/share/FGI2017B/prefix/jdk/jdk1.8.0_281/bin/java\" in \"/share/FGI2017B/brew/bin/picard\".
	@echo [!]Remember to execute \"ln -s /share/FGI2017B/users/huxs/2020Q4_fm2/wdl/fq\" for FASTQ files.
	@ls -lh release.tgz

dbsnp: ref/dbsnp154.vcf.gz
	@echo
	@ls -lh $<

ref/dbsnp154.vcf.gz: dbsnp/GCF_000001405.38.gz dbsnp/GRCh38.chrids
	@# https://www.biostars.org/p/481746/
	@# curl -OJR https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
	-mkdir ref
	bcftools annotate --rename-chrs dbsnp/GRCh38.chrids -Oz -o $@ dbsnp/GCF_000001405.38.gz

dbsnp/GRCh38.chrids:
	-mkdir dbsnp
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt -o dbsnp/GCF_000001405.39_GRCh38.p13_assembly_report.txt
	awk -F'\t' 'x==1 && !/\tna\t/ {print $$7"\t"$$10} /^# Sequence-Name/ {x=1}' dbsnp/GCF_000001405.39_GRCh38.p13_assembly_report.txt >$@

dbsnp/GCF_000001405.38.gz:
	@echo Download from "https://ftp.ncbi.nih.gov/snp/redesign/archive/b154/VCF/GCF_000001405.38.gz[.md5] and https://ftp.ncbi.nih.gov/snp/redesign/archive/b154/VCF/GCF_000001405.38.gz.tbi[.md5]"
