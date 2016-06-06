# BSuit

## References

 * sam/bam 文件的访问 <https://github.com/samtools/samtools>, <https://github.com/samtools/htslib>
   * sam/bam 文件的格式 <http://samtools.github.io/hts-specs/>
   * 测试用的 bam 文件 <https://github.com/gaosjlucky/BS-viral-inte/tree/master/example18/simVir4_aln>
 * 考虑参考序列的组装 IDBA-Hybrid from: <https://github.com/loneknightpy/idba>
 * 动态规划比对(align) <https://github.com/yesimon/M-Vicuna>
   * 修改成甲基化的版本 <https://github.com/gaosjlucky/BS-viral-inte/tree/master/src/methlyAln/src>

## Input

 * bam file, sorted and with index.
 * Reference File of _both_ Human and Virus, with samtools fasta index. 
 * Virus ChrIDs split with comma

## Tasks

### task1

 * 从`bam`文件中提取部分 align 到人上，或者部分 align 到病毒上的reads。按出现过的 Read_ID 提取 alignment。
 * 根据 alignment 的坐标，按照最小 overlap＝2，将相同 Read_ID 的所有 alignment 分类成簇。

以下分析对每个“簇”独立处理。

### task2

#### BS-Seq 的 CT_GA 转换的说明

 * Ref: 原始参考序列
 * Watson/Forward: 将`Ref`中大部分的`C`换成`T`.
 * Crick/Reverse: 将`Ref`中大部分的`G`换成`A`.

alignment 后部有`YD:Z:f`/`YD:Z:r`的 tag，其中`f`表示这条记录是 align 到 Forward 链上的；`r`则对应 Reverse 链。

#### 简化的OLC组装

> 根据 alignment 坐标，把每个“簇”中，人和病毒的参考序列被覆盖的区域挑出，两侧延伸 100 bp。

 * 根据 alignment 坐标，把每个“簇”中，覆盖了人、病毒的 reads 分别 layout（摞） 到一起，统计得到 consensus。记录两条 consensus 与参考序列不同，且不符合 CT_GA 转换的点。
 * 把人和病毒的参考序列被覆盖的区域挑出，两侧延伸 100 bp。按照 CT_GA 转换与上一步人与病毒的两个 consensus 结果作align，得到未知区域的来源。
 * 两条 consensus 做align合并。并确定断点。



---

## Input

所有新产生的文件都是INI格式。<https://en.wikipedia.org/wiki/INI_file>

### Sample of Input

````ini
[RefFiles]
HostRef=/bak/seqdata/genomes/HomoGRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
VirusRef=/share/work/bsvir/X04615.fa.gz

[Output]
WorkDir=/share/work/bsvir/bsI
ProjectID=simVir4

[DataFiles]
Fsimout_m13FG.1=/share/bsim/simout_m13FG.1.fq
Fsimout_m13FG.2=/share/bsim/simout_m13FG.2.fq
Fsimout_m2D.1=/share/bsim/simout_m2D.1.fq
Fsimout_m2D.2=/share/bsim/simout_m2D.2.fq
Fsimout_m458AE.1=/share/bsim/simout_m458AE.1.fq
Fsimout_m458AE.2=/share/bsim/simout_m458AE.2.fq
Fsimout_m6.1=/share/bsim/simout_m6.1.fq
Fsimout_m6.2=/share/bsim/simout_m6.2.fq
Fsimout_m7C.1=/share/bsim/simout_m7C.1.fq
Fsimout_m7C.2=/share/bsim/simout_m7C.2.fq
Fsimout_m9.1=/share/bsim/simout_m9.1.fq
Fsimout_m9.2=/share/bsim/simout_m9.2.fq
Fsimout_mB.1=/share/bsim/simout_mB.1.fq
Fsimout_mB.2=/share/bsim/simout_mB.2.fq

[InsertSizes]
Fsimout_m13FG=200
Fsimout_m13FG.SD=1
Fsimout_m2D=200
Fsimout_m2D.SD=1
Fsimout_m458AE=420
Fsimout_m458AE.SD=1
Fsimout_m6=150
Fsimout_m6.SD=1
Fsimout_m7C=150
Fsimout_m7C.SD=1
Fsimout_m9=200
Fsimout_m9.SD=1
Fsimout_mB=100
Fsimout_mB.SD=1
````

* 所有分析要用的文件都在`WorkDir`内，下文中的相对路径都以`WorkDir`为起点。
* 参考序列在`./Ref/`下，信息记录在`./Ref/Ref.ini`中。
* bam文件在`./${ProjectID}/`下，即`./simVir4_aln/`。Main filename 与`[DataFiles]`下的 main name 一致。这些bam文件是sort并建了索引的（*.bam.bai）。如下：

````
./simVir4_aln/Fsimout_m13FG.bam
./simVir4_aln/Fsimout_m2D.bam
./simVir4_aln/Fsimout_m458AE.bam
./simVir4_aln/Fsimout_m6.bam
./simVir4_aln/Fsimout_m7C.bam
./simVir4_aln/Fsimout_m9.bam
./simVir4_aln/Fsimout_mB.bam
````

### Sample of `./Ref/Ref.ini`

````ini
[s8DOq4JYCenEBMcqnHPOXTN0wB4]
VirusChrIDs=gi|86261677|emb|AJ507799.2|
RefChrIDs=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM
Refilename=/share/work/bsvir/bsI/Ref/s8DOq4JYCenEBMcqnHPOXTN0wB4/GH.fa
chr1=248956422
……省略……
chrY=57227415
gi|86261677|emb|AJ507799.2|=171823

[A5XBKg1ErlCK75YDYxJbT9s9xPI]
VirusChrIDs=gi|59585|emb|X04615.1|
RefChrIDs=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM
Refilename=/share/work/bsvir/bsI/Ref/A5XBKg1ErlCK75YDYxJbT9s9xPI/GX.fa
chr1=248956422
chr10=133797422
……省略……
chr9=138394717
chrM=16569
chrX=156040895
chrY=57227415
gi|59585|emb|X04615.1|=3215
````
* `A5XBKg1ErlCK75YDYxJbT9s9xPI`是先后用人和病毒的参考序列的两个文件名通过`Digest::SHA::sha1_base64`出来的，取值范围就是`base64`范围。如果不便计算，可以作为输入参数提供。
   * 例：`Digest::SHA::sha1_base64("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz,X04615.fa.gz")`的值为`"A5XBKg1ErlCK75YDYxJbT9s9xPI"`。
* `RefChrIDs`是所有人的染色体的id，逗号分隔。
* `VirusChrIDs`是所有病毒的染色体的id，逗号分隔。
* `Refilename`是解压合并后的参考序列，用来跑`bwa-meth`。
* 下面跟的是每条染色体的长度。

### Sample of bam files

````
@HD	VN:1.3	SO:coordinate
@SQ	SN:chr1	LN:248956422
……省略……
@SQ	SN:chrM	LN:16569
@SQ	SN:gi|59585|emb|X04615.1|	LN:3215
@RG	ID:Fsimout_m13FG	SM:Fsimout_m13FG
@PG	ID:bwa-meth	PN:bwa-meth	VN:0.10	CL:"/share/bsIntegration/bin/bwameth.py --reference /share/work/bsvir/bsI/Ref/A5XBKg1ErlCK75YDYxJbT9s9xPI/GX.fa -t 24 --read-group Fsimout_m13FG -p /share/work/bsvir/bsI/simVir4_aln/Fsimout_m13FG /share/bsim/simout_m13FG.1.fq /share/bsim/simout_m13FG.2.fq"
……省略……
sf330_Ref_1806938_1807137_1807337_Vir_-_108_287_R_200_90	353	chr1	98374819	0	50S40M	chr18	1807198	0	CATGTTCGGTGCAGGGTCCCCAGTCCTCGAGAAGATTGACGATATGGGTGTCTTGAATATCTCCTGAAGGACTTGCCTGAGGCTGTTTTA	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	NM:i:0	MD:Z:40	AS:i:40	XS:i:40	RG:Z:Fsimout_m13FG	SA:Z:rgi|59585|emb|X04615.1|,108,-,40S50M,60,0;	YC:Z:CT	YD:Z:f
……省略……
````

## Need to do

涉及的常数都作参数。

### Todo 1. Grep Covered Blocks on Each chromosomes

* 用`samtools depth`的方法，统计每条染色体被覆盖的情况。
* bam文件中，PE中一条在人上，一条在病毒上，保留这一对。分类为A1。
* bam文件中，Read自身有多条比对记录，既有若干比到人上的，也有若干比到病毒上的，且二者CIAGR为`M`或`I`的部分中，重叠的部分乘以2小于不重叠的部分相加，且重叠的部分小于20。若干中组合中发生一次，保留这一对中发生了的那些组合。分类为A2。
   * A1,A2 是为了后面方便，算覆盖时不区分。
* bam文件中，发生soft-clip的reads，保留这一对。分类为B。
* bam文件中，PE都比到人的同一染色体上，且bam文件中的insert size与输入文件中`[InsertSizes]`对应记录相比，超过平均数正负3倍SD的，保留这一对。分类为C1。
* bam文件中，PE只有一条比对上去，另一条unmap，保留这一对。分类为C2。
   * 现在没有考虑C类，留下以后再说。
* 扔掉其余的Reads。
* 在人的染色体上，将A类深度大于0，这样的连续区域定义为一个pre-block。
* 若B>=3的地方与A>0的地方overlap>5，且没有碰上其他的pre-block，则延伸pre-block直到B<3。相邻pre-block延伸到互不接触的最小B值，然后按下面的条件判断是否合并。
* 若延伸pre-block过程中遇上其他pre-block，若合并后整个pre-block长度小于`[InsertSizes]`对应记录中平均数加3倍SD，则直接合并。否则保持前面不合并的状态，但记录被哪些B连接。
   * 人最长的一号染色体不到250m，分染色体处理。可以用`struct {uint16_t A_Depth, uint8_t B_Depth, uint8_t C_Depth} []`来分别记录是否unique的reads覆盖的次数，然后`*`过去用`uint32_t Depth []`来统计A+B深度为0，即`Depth[]`大于256的地方。内存开销小于1G。
      * 这样写得注意在编译时按 大/小-endian 分别判断下 byte order。
* 处理完后的pre-block就是block，同时可以得到block在这条人染色体上的覆盖范围，及范围内的平均深度。对前面选择“不合并”的，覆盖范围会有多个。按起点排序，逗号分隔。
* 将block中reads涉及的病毒的部分取出，用来统计在病毒上的覆盖范围，及范围内的平均深度。若有多个范围，按深度降序排列，逗号分隔坐标，分号分隔染色体（如果病毒有多条染色体）。
* 输出调试用文件`./simVir4_grep/blocks.txt`

#### Sample of `blocks.txt`
````ini
[B1]
HostRange=chr1:7603420-7603493=3,7603520-7603593=2
VirusRange=gi|59585|emb|X04615.1|:2722-2799=1,2822-2899=1;Virus2:222-233=2
HostRange1,VirusRange1-1	A1 sf193_Ref_70796427_70796576_70796726_Vir_-_2722_2841_R_150_90	161	chr1	7603420	0	73M17S	gi|59585|emb|X04615.1|	2722	0	CCTTTACCATTATGTAATGGTCTTCTTTGTCTCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACATCATTACTTCAAAACT	HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH	NM:i:0	MD:Z:73	AS:i:73	XS:i:73	RG:Z:Fsimout_m6	YC:Z:GA	YD:Z:r
HostRange1,VirusRange1-2	A2 sf193_Ref_70796427_70796576_70796726_Vir_-_2722_2841_R_150_90	81	gi|59585|emb|X04615.1|	2722	60	13S77M	chr1	7603420	0	TTTTATCAGAGACATCATTACTTCAAAACTAGGCATTATTTACATACTCTGTGGAAGGCTGGCATTCTATATAAGAGAGAAACTACACGC	HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH	NM:i:0	MD:Z:77	AS:i:77	XS:i:28	RG:Z:Fsimout_m6	YC:Z:CT	YD:Z:r
HostRange2,VirusRange2-1	B  sf217_Ref_70796427_70796576_70796726_Vir_-_2722_2841_R_150_90	593	chr1	7603456	0	37M53S	=	7603396	-97	TGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACATCATTACTTCAAAACTAGGCATTATTTACATACTCTGTGGAAGGCTGGCATT	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	NM:i:0	MD:Z:37	AS:i:37	XS:i:39	RG:Z:Fsimout_m6	YC:Z:CT	YD:Z:r
HostRange2,NA	C1 sf207_Ref_70796427_70796576_70796726_Vir_-_2722_2841_R_150_90	593	chr1	7603466	0	27M63S	=	7603406	-87	TGGTTTAAAGTCTGTTTTATCAGAGACATCATTACTTCAAAACTAGGCATTATTTACATACTCTGTGGAAGGCTGGCATTCTATATAAGA	HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH	NM:i:0	MD:Z:27	AS:i:27	XS:i:29	RG:Z:Fsimout_m6	YC:Z:CT	YD:Z:r

[B2]
……省略……
````

### Todo 2. Use Pair-End Information to Link Blocks to Bridges

* 选取A2类Reads中满足前面“重叠部分”长度要求的记录，以及B类及C1类Reads中`MAPQ`大于40的。用来统计blocks间的连接关系。block为节点，PE关系为边。
* 对只有两个block相连的图，记为一个Bridge（定义为只有两个桥墩的桥，如拱桥）。
* 对前面选择“不合并”的block，根据拱桥是否冲突决定是否合并。

### Todo 3. Analyse Junction sites on Each Block

* 用IDBA-Hybrid的方法考虑参考序列组装，将上面Block的Range对应的参考序列两边延伸100bp（Reads长度）作为ref，组装reads。
   * 根据reads判断是有甲基化的四碱基序列，还是没有甲基化的三碱基序列。四碱基的直接装。三碱基的如果能改比对打分矩阵就改它，否则把ref也变成三碱基的。
      * 正负链与AG／CT的关系高升杰补充下。
   * 四碱基的结果与三碱基的结果作alignment合并，可以参考修改的M-Vicuna。
* 将用人作ref得到的组装结果，与用病毒作ref得到的组装结果，作alignment合并。得到完整片段（VirusFragment）。
   * 结果不唯一就选最长的。为了同时考虑block的覆盖深度，可以把A类深度作权重，即覆盖了最多A区的是对的。权重可以加个系数＝0.8让长度优先。
   * 优势不明显的话，下面能合并得较好的，是正确的结果。

### Todo 4. Check whether VirusFragment can be aligned between Blocks of same Bridge

* 将拱桥两个block的结果片段（VirusFragment），检查`[InsertSizes]`对应记录看平均数正负3倍SD范围内是否有可能重叠。能重叠就作alignment合并。
