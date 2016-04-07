# Meerkat 0.185 测试结果

在 Meerkat 的结果中，能找到插入事件，但无法正确报告插入序列为非人类片段，大部分情况报告为其他染色体的片段。  
对断点的判断不准确，在模拟的断点附近有许多 tandem_dup 结果。

按照最标准的模式，参考序列为hg18，比对软件为bwa：  
有95条原始结果。其中没有信息正确的。
只有77条在模拟断点 10bp 附近，但pattern错误。

## 在推荐参数基础上，测试做了如下调整：

* 设置过高的深度上限，没有滤除任何区域（确保`*.blacklist.gz`为空）；
* FastQ的质量值过滤已关闭（-q 0）。
* 灵敏度设为最高（-o 0）。
* 没有进行任何过滤步骤。(包括但不限于，没有用somatic_sv.pl to filter out germline events and other artifacts.)

### 用bwa比对，参考序列为hg18	95 -> 77

````bash
scripts/pre_process.pl -s 20 -k 1500 -q 0 -t 8 -b merged/simout.bam -I hg18/HomoGRCh38 -A hg18/HomoGRCh38.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 0 -t 8 -b merged/simout.bam -F hg18/ -S samtools-0.1.20/
scripts/mechanism.pl -b merged/simout.bam -R hg18rmsk/rmsk-hg18.txt
````

得到 95 条结果，没有未知片段插入报告，将错误比对的其他染色体的情况报告为转座。  
仅考虑是否覆盖chr18上已模拟的断点左右10bp范围，有 77 条。

### 用bwa比对，参考序列为hg18 + X04615.1	6 -> 0

````bash
scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b mergedvir/simout.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 0 -t 8 -b mergedvir/simout.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b mergedvir/simout.bam -R hg18rmsk/rmsk-hg18.txt
````

得到 6 条结果，没有病毒插入报告，将错误比对的其他染色体的情况报告为转座。  
仅考虑是否覆盖chr18上已模拟的断点左右10bp范围，没有结果。放宽到 25bp 时有一条。

### 用bwa-meth比对，参考序列为hg18 + X04615.1	13 -> 0

````bash
samtools merge -p simvir/simVir4.bam simVir4_aln/Fsimout_m*.bam
scripts/pre_process.pl -s 20 -k 15000 -q 0 -t 8 -b simvir/simVir4.bam -I Ref/GX.fa -A Ref/GX.fa.fai -S samtools-0.1.20/
scripts/meerkat.pl -d 5 -s 20 -p 5 -o 0 -t 8 -b simvir/simVir4.bam -F Ref/ -S samtools-0.1.20/
scripts/mechanism.pl -b simvir/simVir4.bam -R hg18rmsk/rmsk-hg18.txt
````

得到 13 条结果，有一例未知片段插入报告，但离模拟的断点有 81bp 距离，且报告有 34bp 的 deletion。

	del_ins	FoSTeS	1300	9	3	chr18	7692039	7692074	34	-	-	-	43	BP:TE_TE

没有病毒插入报告，将错误比对的其他染色体的情况报告为转座。  
仅考虑是否覆盖chr18上已模拟的断点左右10bp范围，没有结果。放宽到 72bp 时有一条。
