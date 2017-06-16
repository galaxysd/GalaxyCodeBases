---
title: 'SMAP Further Analysis Report'
author: FGI
date: 2017-06-15
lang: zh-CN
toc-title: 目录
papersize: A4
...

## 分析结果

### 1 数据基本处理与质控

所有样品的 WGBS 测序数据，将下机数据进行过滤，包括去污染，去测序接 头和低质量碱基比例过高的 reads，得到 clean data。表 1-1 中列出了所有样本的 数据产出概况。图 1.1 显示的是样品 **Simulate01** 的测序碱基含量分布，图 1.2 显 示的是样品 **Simulate01** 的碱基测序质量分布情况。其余样品的测序碱基含量分布 图与碱基测序质量分布情况图可在路径   `/output/sample_name/sample_type/01.Data_Summary_and_QC` 下查询。

: 表 1-1:数据基本处理与质控

| Sample Name | Sample Type | Total read | Total base | Clean read | Clean base | Clean Rate (%) |
| --- |:---:| ---:| ---:| ---:| ---:| ---:|
| Simulate01 | Normal | 664889900 | 83111056500 | 571476968 | 59657917639 | 71.78% |
| Simulate01 | Tumor | 633165724 | 79145540250 | 538312038 | 55353278301 | 69.94% |
| Simulate02 | Normal | 664889900 | 83111056500 | 571476968 | 59657917639 | 71.78% |
| Simulate02 | Tumor | 633165724 | 79145540250 | 538312038 | 55353278301 | 69.94% |

Clean Rate (%) = Clean Data Size (bp)/Raw Data Size (bp)

![图 1.1: 样品 Simulate01 Clean reads 的碱基含量分布图](report/1_1.fq_qc.png)

横坐标表示 reads 上碱基所 在位置，纵坐标表示碱基比例。如果图中碱基分布不平衡则说明测序过程有异常情况 发生。

![](report/1_2.fq_qc.png)

**图 1.2: 样品 Simulate01 Clean reads 的碱基测序质量分布图**  
横坐标为 reads 上碱基位 置;纵坐标为碱基测序质量。

在得到 clean data 之后，使用比对软件 BSMAP 将 reads 比对到参考基因组 上，比对的统计结果如表 1-2 所示;然后根据需要对各个文库的 reads 进行去 duplication 处理。

: 表 1-2:比对结果统计

| Sample Name | Sample Type | Total Reads | Pair Mapped | Pair Mapped Rate |
| --- |:---:| ---:| ---:| ---:|
|Simulate01|Normal|420154240|394936188|94.00%|
|Simulate01|Tumor|398577754|365101772|91.60%|

表 1-3 是所有样品的测序深度和覆盖度情况。 图 1.3 为样品 Simulate01 的测序 深度分布图，理论上，其最高点对应的测序深度与全基因组平均覆盖深度一致或接 近，这个分布图可以用于反映测序是否均匀。其他样品的测序深度分布图可在路径 /output/sample_name/sample_type/ 01.Data_Summary_and_QC 下查询。表 1-4 为所有样 品在全基因组上的 C 位点覆盖度。

![](report/1_3.depth.png){ width=49% } ![](report/1_4.depth.png){ width=49% }

**图 1.3:样品 Simulate01 测序深度覆盖度分布图**

### 2. 全基因组甲基化水平分析

用于分析的 DNA 样品为多细胞样品，因此 C 碱基的甲基化水平是一个 0%~100%范围 内的数值，等于该 C 碱基上覆盖到的支持 mC 的序列数除以有效覆盖的序列总数，通常 CG 甲基化存在于基因和重复序列中，在基因表达调控过程中起到非常重要的作用。非 CG 类型 的序列(CHG 和 CHH)在基因中十分少见，主要存在于基因间区和富含重复序列的区域， 在沉默转座子过程中起关键作用。图 2-1 和表 2-1 为样品 Simulate01 在全基因组各染色 体上的平均甲基化水平，表 2-2 为样品 Simulate01 在全基因组各类型调控元件范围内的 甲基化水平。其余样品在全基因组各染色体上的平均甲基化水平与在全基因组各类型调控元 件 范 围 内 的 甲 基 化 水 平 信 息 表 可 在 路 径 /output/sample_name/sample_type/ 02.Average_Methylation_Level_of_the_Whole_Genome 下查询。(加上点图 所有样品的三张 合一)

### 3. 甲基化C碱基中CG,CHG 与CHH的分布比例

### 4. 甲基化的 CG，CHG，CHH 附近碱基的序列特征分析

### 5. 染色体水平的甲基化 C 碱基密度分布

### 6. 基因组不同转录元件中的 DNA 平均甲基化水平

### 7. ASM

### 8. DMR 的检测

### 9. DMR 相关基因的 GO 和 Pathway 分析

## 分析方法 {#methods}

### 1. 实验流程

### 2. 信息分析流程

### 3. 数据过滤

### 4. 序列比对

### 5. 甲基化水平

### 6. DMR 检测

### 7. 甲基化水平程度差异

### 8. GO 注释

### 9. KEGG 通路富集

## 参考文献

1. xxx
2. yyy

---

* See [Markdown Guide](http://pandoc.org/MANUAL.html#pandocs-markdown) for Writing help.
* Edit `./styles/kultiad-serif.css` to modify `html` styles.
* Use OpenOffice to edit XML file `./templates/odt.template` to modify `docx` styles.
* PDF output is not ready, use `docx` file to generate pdf.

---

<div class="pagebreak"></div>

(@good)  This is a good example.


Here is a literal backtick `` ` ``.

H~2~O is a liquid.  2^10^ is 1024.

As (@good) illustrates, ...

P~a\ cat~

[Small caps]{style="font-variant:small-caps;"}

$1+2\neq3!$

````
pandoc -f markdown  -s --toc -V theme:moon -V toc-title:"目录" -o report.docx report.md
pandoc -f markdown -t html5 -s --toc -V toc-title:"目录" -o report.html report.md
pandoc -f markdown -t html5 -s --toc -V toc-title:"目录" -o report.pdf report.md
pandoc -f markdown -t html5 -s --toc -V theme:moon,toc-title:"目录" -o report.html report.md
````
