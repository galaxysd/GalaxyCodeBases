# BSVF
[Bisulfite Sequencing Virus integration Finder](https://github.com/BGI-SZ/BSVF)

## Attention

For directional libraries only. PBAT and indirectional libraries are _NOT_ supported.

## Dependencies

[bwa-meth 0.10](https://github.com/brentp/bwa-meth/tree/0a2f9fc7c3fd3c99c4212941c94be73c9c865bb1) depends on 

 + python 2.7+ (including python3)
   - `toolshed` library. can be installed with: 
      * `easy_install toolshed` or
      * `pip install toolshed`

 + samtools command on the `$PATH` (https://github.com/samtools/samtools)

 + bwa mem from: https://github.com/lh3/bwa

 + EMBOSS from: http://emboss.sourceforge.net/

## Install

### Normal

Run `pip install toolshed`.

Run `src/install.sh`.

In case EMBOSS failed to install, you'll need to download the binary from above sites. And put `water` of EMBOSS in to `./bin`. Or, just link `water` to `./bin`.

Your `bsIntegration/bin/` should be like this:
````bash
-rwxr-xr-x  398860 Feb 20 00:48 bwa
-rwxr-xr-x   21892 Sep  1 08:37 bwameth.py
-rwxr-xr-x   27040 Feb 20 01:14 water
-rwxr-xr-x  971772 Feb 20 00:48 samtools
````

### Homebrew/Linuxbrew

````bash
brew tap Ensembl/homebrew-external
brew install emboss bwa samtools python
pip install toolshed

ln -s `which bwa` ./bin/
ln -s `which samtools` ./bin/
ln -s `which water` ./bin/

brew install gcc
cd ./src/analyser/
make
cp -av analyser/bsanalyser ../../bin/
cd ../../bin/
ls -l
````

## Citation

Gao, S., Hu, X., Xu, F., Gao, C., Xiong, K., Zhao, X., … Pedersen, C. N. S. (2017, December 18). BS-virus-finder: virus integration calling using bisulfite sequencing data. GigaScience. <https://doi.org/10.1093/gigascience/gix123>

Repo URL: <https://github.com/BGI-SZ/BSVF>

## Usage

```
./bsuit <command> <config_file>

./bsuit prepare prj.ini
./bsuit aln prj.ini
./bsuit grep prj.ini
./bsuit analyse prj.ini
```

![a Logo](https://raw.githubusercontent.com/BGI-SZ/BSVF/master/logo/BSVFlogo.png)

## Test Run

```
mkdir sim90 && cd sim90 && ./simVirusInserts.pl GRCh38_no_alt_analysis_set.fna.gz X04615.fa.gz s90 && cd ..
mkdir sim50 && cd sim50 && ./simVirusInserts.pl GRCh38_no_alt_analysis_set.fna.gz X04615.fa.gz s50 50 ../sim90/s90.ini && cd ..

./bsuit prepare sim90/s90.ini

./bsuit aln sim90/s90.ini
./run/s90_aln.sh
./bsuit grep sim90/s90.ini
./bsuit analyse sim90/s90.ini

./bsuit aln sim50/s50.ini
./run/s50_aln.sh
./bsuit grep sim50/s50.ini
./bsuit analyse sim50/s50.ini
```

## Reference Files

 * Human: <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>
 * HBV: [gi|59585|emb|X04615.1| Hepatitis B virus genome, subtype ayr](http://www.ncbi.nlm.nih.gov/nuccore/X04615.1?report=GenBank)

## Simulation

    ./simVirusInserts.pl GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz HBV.X04615.fa sim150 150

### GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

A gzipped file that contains FASTA format sequences for the [following](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt):

1. chromosomes from the GRCh38 Primary Assembly unit.  
   Note: the two PAR regions on chrY have been hard-masked with Ns.   
   The chromosome Y sequence provided therefore has the same 
   coordinates as the GenBank sequence but it is not identical to the
   GenBank sequence. Similarly, duplicate copies of centromeric arrays
   and WGS on chromosomes 5, 14, 19, 21 & 22 have been hard-masked 
   with Ns (locations of the unmasked copies are given below). 
2. mitochondrial genome from the GRCh38 non-nuclear assembly unit.
3. unlocalized scaffolds from the GRCh38 Primary Assembly unit.
4. unplaced scaffolds from the GRCh38 Primary Assembly unit.
5. Epstein-Barr virus (EBV) sequence  
   Note: The EBV sequence is not part of the genome assembly but is 
   included in the analysis set as a sink for alignment of reads that
   are often present in sequencing samples.

## Format of `config_file`

### An example

```ini
[RefFiles]
HostRef=/share/HomoGRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
VirusRef=/share/work/bsvir/HBV.AJ507799.2.fa

[DataFiles]
780_T.1=/share/work/bsvir/F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz
780_T.2=/share/work/bsvir/F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz
s01_P.1=/share/work/bsvir/F12HPCCCSZ0010_Upload/s01_P.bs_1.fq.gz
s01_P.2=/share/work/bsvir/F12HPCCCSZ0010_Upload/s01_P.bs_2.fq.gz
;MultiLibExample.1=/test/Lib1/AAAA.1.fq.gz, /test/Lib2/AAAA.1.fq.gz, /test/Lib3/BBBB.1.fq.gz
;MultiLibExample.2=/test/Lib1/AAAA.2.fq.gz, /test/Lib2/AAAA_2.fq , /test/Lib3/BBBB.2.fq.gz
tSE_X.1=/share/work/bsvir/F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz,/share/work/bsvir/F12HPCCCSZ0010_Upload/s01_P.bs_2.fq.gz

[InsertSizes]
780_T=200
780_T.SD=120
s01_P=200
s01_P.SD=30
;MultiLibExample=210
;MultiLibExample.SD=70
tSE_X=90
tSE_X.SD=1

[Output]
WorkDir=/share/work/bsvir/bsI
ProjectID=SZ2015

[Parameters]
Aligner=bwa-meth
MinVirusLength=20
```

## Build

You'll need `cmake` and `autoconf, automake` and devel-libs, as well as `gcc, g++` to compile all sources.

For Mac OS X, install [Homebrew](http://brew.sh/) first. Then:
```bash
xcode-select --install
brew install autoconf automake cmake python
brew install --without-multilib gcc
````

To Build the binaries:
```bash
cd src
./download.sh
./install.sh

pip install toolshed
```

### Details

 + For comment lines, use `;` as the first character.

 + `RefFiles` Section
   - `HostRef` is **Host genome**.
   - `VirusRef` is **Virus sequence**.

 + `DataFiles` Section
   - Each *Sample* need an **unique ID** as *SampleID*. Use `SampleID.1` and `SampleID.2` to specify pair-end sequencing data.
   - For samples with multiple PE sets, join each file with *comma* and keep their order.

 + `InsertSizes` Section
   - For each `SampleID`, use `SampleID` to specify average insert sizes. And use `SampleID.SD` to specify its standard deviation.

 + `Output` Section
   - `WorkDir` is the output directory.
   - `ProjectID` is an **unique ID** for this analyse defined in the `config_file`.

## Description

**BSuit** is a suit to analyse xxx.

### Formats

#### 病毒整合结果文件

````
Chr	breakpoint	virus-start virus-end virusstrand	how-many-reads-support cluster-name
Chr1	3000	200	300	+/-	20 cluster1
````

#### 中间contig信息文件

````
clustername contig-number chrpoint virus-integration
cluster1	contig1	chr1:3000	virus:+:200-300
cluster1	contig2	chr2:4000	viurs:-:300-400
````

## ToDo

Compare with [ViralFusionSeq [VFS]](https://sourceforge.net/projects/viralfusionseq/) and [VirusFinder 2](https://bioinfo.uth.edu/VirusFinder/) on normal WGS data.

## See also

* [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) [0.18.1](https://github.com/FelixKrueger/Bismark/releases)
* [SVDetect](https://sourceforge.net/projects/svdetect/) [r0.8b](https://sourceforge.net/projects/svdetect/files/SVDetect/0.80/)
* [ViralFusionSeq](https://sourceforge.net/projects/viralfusionseq/) and [Virus-Clip](http://web.hku.hk/~dwhho/Virus-Clip.zip)
* [VirusSeq](http://odin.mdacc.tmc.edu/~xsu1/VirusSeq.html), which uses [MOSAIK](https://github.com/wanpinglee/MOSAIK) aligner.
* [VirusFinder 2 (VERSE)](https://bioinfo.uth.edu/VirusFinder/)
* [Vy-PER](http://www.ikmb.uni-kiel.de/vy-per/)
* [seeksv](https://github.com/qiukunlong/seeksv)
* [BS-Seeker2](https://github.com/BSSeeker/BSseeker2)

* [Sherman](https://www.bioinformatics.babraham.ac.uk/projects/sherman/)
  * 模拟时甲基化率设为800.

| Tool | Sequencing Type | Programme Language | 1st Aligenment * | Assembler | 2nd Aligenment # | Epub Date |
|:-----|:---------------:|:------------------:|:-----:|:-----:|:-----:|:-----:|
|[VirusSeq](https://doi.org/10.1093/bioinformatics/bts665)|RNA-Seq, WGS|Perl|MOSAIK to Human|MOSAIK to Virus|MOSAIK to Hybrid|2012 Nov 08|
|[ViralFusionSeq](https://doi.org/10.1093/bioinformatics/btt011)|RNA-Seq, WGS|Perl|BWA-SW to Human|cap3, SSAKE|Blastall to Virus|2013 Jan 12|
|[VERSE(VirusFinder2)](https://doi.org/10.1186/s13073-015-0126-6)|WGS, RNA-Seq|Perl|Bowtie2 to Human, BLAT to Virus, BLASTN to Virus|Trinity|BWA-SW to Hybrid, SVDetect,CREST|2015 Jan 20|
|[Virus-Clip](https://doi.org/10.18632/oncotarget.4187)|RNA-seq|Perl|BWA-MEM to Virus|Virus-Clip|BLASTN to Human|2015 May 19|
|[Vy-PER](https://doi.org/10.1038/srep11534)|WGS, RNA-Seq|Python2|BWA-SW to Human|Vy-PER|BLAT to Virus|2015 Jul 13|
|[seeksv](https://doi.org/10.1093/bioinformatics/btw591)|WGS|C++|BWA to Hybrid|seeksv|seeksv to Hybrid|2016 Sep 14|
|BSVF|WGBS, WGS|Perl,C,C++|BWA-MEM to Hybrid|BSVF|water(EMBOSS) to Hybrid| N/A |

\* for virus-infected reads  
\# for integration infomation

## One More Things

To extract relevant PE reads within 500 bp range from final result, *BS.analyse* for example.

```bash
perl -lane '$a=$F[2]-501;$b=$F[2]+501;print join("\t",$F[1],$a,$b)' ../W2BS_analyse/BS.analyse >zones.bed
vi zones.bed # To remove the first head line
# sort BS.bam to BS.sort.bam and index it.
samtools view -L zones.bed BS.sort.bam > zones.sam
awk '{print $1}' zones.sam | sort | uniq > zones.ids
#samtools view BS.bam | grep -F -f zones.ids >zones.PE.sam
samtools view BS.sort.bam | grep -F -f zones.ids > zones.PEs.sam # sorted one maybe more useful.
```
