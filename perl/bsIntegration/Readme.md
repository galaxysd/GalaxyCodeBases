# Dependencies

`bwa-meth` depends on 

 + python 2.7+ (including python3)
   - `toolshed` library. can be installed with: 
      * `easy_install toolshed` or
      * `pip install toolshed`

   - You can install `bwa-meth` accroding to https://github.com/brentp/bwa-meth

 + samtools command on the `$PATH` (https://github.com/samtools/samtools)

 + bwa mem from: https://github.com/lh3/bwa

 + IDBA-Hybrid from: https://github.com/loneknightpy/idba

# Install

You'll need `cmake` and `aclocal, autoconf, automake` and devel-libs, as well as `gcc, g++` to compile all sources.

```
cd src
./download.sh
./install.sh
pip install toolshed
```

# Usage

```
./bsuit <command> <config_file>

./bsuit prepare prj.ini
./bsuit aln prj.ini
./bsuit grep prj.ini
./bsuit analyse prj.ini
```

# Format of `config_file`

## An example

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

[InsertSizes]
780_T=200
780_T.SD=120
s01_P=200
s01_P.SD=30
;MultiLibExample=210
;MultiLibExample.SD=70

[Output]
WorkDir=/share/work/bsvir/bsI
ProjectID=SZ2015
```

## Details

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

# Description

**BSuit** is a suit to analyse xxx.

## 代码说明

所有输出文件都在`Output.WorkDir`下。

### do_pre()

根据`[RefFiles]`信息，生成bwa-meth的参考序列。并将染色体长度及物种来源信息储存于`Ref/Ref.ini`。

### do_aln()

