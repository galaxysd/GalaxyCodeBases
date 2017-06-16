---
title: 'SMAP Further Analysis Report'
author: FGI
date: 2017-06-15
lang: zh-CN
toc-title: 目录
...

## 分析结果

### 1 数据基本处理与质控

所有样品的 WGBS 测序数据，将下机数据进行过滤，包括去污染，去测序接 头和低质量碱基比例过高的 reads，得到 clean data。表 1-1 中列出了所有样本的 数据产出概况。图 1.1 显示的是样品 **Simulate01** 的测序碱基含量分布，图 1.2 显 示的是样品 **Simulate01** 的碱基测序质量分布情况。其余样品的测序碱基含量分布 图与碱基测序质量分布情况图可在路径   `/output/sample_name/sample_type/01.Data_Summary_and_QC` 下查询。

: 表 1-1:数据基本处理与质控

| Sample Name | Sample Type | Total read | Total base | Clean read | Clean base | Clean Rate (%) |
| --- | --- | ---:| ---:| ---:| ---:| ---:|
| Simulate01 | Normal | 664889900 | 83111056500 | 571476968 | 59657917639 | 71.78% |
| Simulate01 | Tumor | 633165724 | 79145540250 | 538312038 | 55353278301 | 69.94% |
| Simulate02 | Normal | 664889900 | 83111056500 | 571476968 | 59657917639 | 71.78% |
| Simulate02 | Tumor | 633165724 | 79145540250 | 538312038 | 55353278301 | 69.94% |

Clean Rate (%) = Clean Data Size (bp)/Raw Data Size (bp)

![](report/1_1.fq_qc.png)

**图 1.1: 样品 Simulate01 Clean reads 的碱基含量分布图**  
横坐标表示 reads 上碱基所 在位置，纵坐标表示碱基比例。如果图中碱基分布不平衡则说明测序过程有异常情况 发生。

![](report/1_2.fq_qc.png)

**图 1.2: 样品 Simulate01 Clean reads 的碱基测序质量分布图**  
横坐标为 reads 上碱基位 置;纵坐标为碱基测序质量。


### 2. 全基因组甲基化水平分析

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
