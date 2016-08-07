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

### 代码说明

* 所有输出文件都在`Output.WorkDir`下。
* 无印的`.bam` 是后续分析目前使用的。
* `.sn.bam` 是 sort by read name，暂时没用。
* `.snPstat.bam` 就是把`.sn.bam`用`-F256`过滤下，只是用来数数的。

#### do_pre()

根据`[RefFiles]`信息，生成bwa-meth的参考序列。并将染色体长度及物种来源信息储存于`Ref/Ref.ini`。

#### do_aln()

生成`${ProjectID}_aln.sh`

1. 调用`bwa-meth`比对。结果是按染色体坐标排序的。
2. `samtools merge`同一id的。所以必须先排序才能merge。
3. `samtools sort -m 2415919104 -n` -> `${id}.sn.bam`

#### do_grep()

`config_file`的设定是没种生物学样品取一个名字，内部文件用逗号隔开。  
即上面例子中有`780_T`, `s01_P`, `tSE_X`三个生物学上不同的样品，三者分开分析。

##### 161@BSuitLib.pm: 在**按染色体坐标排序的**bam文件中，寻找满足下列条件的：

* bam文件第7列为"="即比到同一染色体上,且第9列的Ins_Size与`config_file`中均值的差大于3倍SD。	`$flag |= 1`
* bam文件第3列为病毒染色体或者第7列为病毒染色体。	`$flag |= 4`

将符合条件的bam行写到sam文件中，列表信息在`${ProjectID}_grep.ini`中。

##### 同时记录 ％ReadsIndex ，每个生物学样品记录一个。

173: %ReadsIndex: {$ReadID} = [[0,$RefChr,$flag]]

204: push @{$ReadsIndex{$ReadID}},[$sam文件指针,$r12R1,[Chr,Pos,MapQ,CIGAR],$flag];  

##### 对 ％ReadsIndex 的key，按染色体坐标排序，保证病毒排前面。

得到 @IDsorted of $ReadID。

##### 同一$ReadID代表同一对PE，用mergeIn函数按坐标范围一个个叠加，直到断开。

这里只包括符合前面161中条件的那些Read。  
叠加时，某个ReadID的所有hit，是同时加入叠加的。这条今晚还没确认当时是怎么实现的。  
**确认目前的block是桥墩，不是整个桥**

得到的block写入`${ProjectID}__grep/blocks.ini`中。

#### do_analyse

##### 过滤block

* sam文件中比对到多个地方的
* 覆盖深度太低的

##### 组装

* 提取被覆盖的参考序列
* 将Reads按照 原样、tr /Cc/Tt/、tr /Gg/Aa/ 处理成3份。
* 将参考序列同上处理出3份。
* 用修改IDBA/M-Vicuna的在提供参考序列的情况下组装。
* 通过组装合并3份结果。
* 比对判断断点。

## Add on

### Find bridge 墩 / 候选簇

* 人的PE，不同染色体的hit直接扔掉不管。
* Get all soft-cliped, with 10S is enough. 参数＝(人10S，选出来参与归簇)，(病毒20S，不分析)
* 确定unmap的，只保留map上的。
* 3*sd not needed.
* 按PE归簇

### Analyse

* 组装后前两根contig都得判断。两根得aln到人的两端。 （只看前两根是否分别在两端）
* 组装只考虑按甲基化处理后三种碱基的情况。
* 对结果归类到图中那堆分类。对fq用aln到contig的方法归类。 （年后）
* 分析病毒整合的热点区域。 （年后）

### 甲基化处理

看bam的flag，reverse<=>负链
read1 正链， C/T
read1 minus， G/A
read2 plus， G/A
read2 minus C/T

SE的同 read1

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
