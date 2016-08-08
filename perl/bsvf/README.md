# BSVF
Bisulfite Sequencing Virus integration Finder

## Attention

For directional libraries only. PBAT and indirectional libraries are _NOT_ supported.

## Dependencies

`bwa-meth` depends on 

 + python 2.7+ (including python3)
   - `toolshed` library. can be installed with: 
      * `easy_install toolshed` or
      * `pip install toolshed`

   - You can install `bwa-meth` accroding to https://github.com/brentp/bwa-meth

 + samtools command on the `$PATH` (https://github.com/samtools/samtools)

 + bwa mem from: https://github.com/lh3/bwa

 + EMBOSS from: http://emboss.sourceforge.net/

## Install

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

## Usage

```
./bsuit <command> <config_file>

./bsuit prepare prj.ini
./bsuit aln prj.ini
./bsuit grep prj.ini
./bsuit analyse prj.ini
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
