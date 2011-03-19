#!/bin/sh
#$ -S /bin/sh
perl /share/raid002/lidong/SNP_filter_bin/add_cn.pl /share/raid002/lidong/40snp/RawSNP /share/raid002/lidong/40snp/population chr1
perl /share/raid002/lidong/SNP_filter_bin/filter_depth_sort_each_chr.pl /share/raid002/lidong/40snp/population /share/raid002/lidong/40snp/population chr1
perl /share/raid002/lidong/SNP_filter_bin/checkeachsite.pl /share/raid002/lidong/silkworm/silkworm-new/ /share/raid002/lidong/40snp/population/chr1.add_cn_cp1.5_q15.filterByChr chr1 >/share/raid002/lidong/40snp/population/chr1.checkEachSite
/share/raid002/lidong/SNP_filter_bin/RankSumTest -i /share/raid002/lidong/40snp/population/chr1.checkEachSite -o /share/raid002/lidong/40snp/population/chr1.rst
perl /share/raid002/lidong/SNP_filter_bin/LCcheck.pl /share/raid002/lidong/silkworm/silkworm-new/ /share/raid002/lidong/40snp/population/chr1.add_cn_cp1.5_q15.filterByChr chr1 >/share/raid002/lidong/40snp/population/chr1.LCcheck
perl /share/raid002/lidong/SNP_filter_bin/check_allSNP.pl /share/raid002/lidong/40snp/population/chr1.add_cn_cp1.5_q15.filterByChr /share/raid002/lidong/40snp/population/chr1.LCcheck >/share/raid002/lidong/40snp/population/chr1.ratioCheck
perl /share/raid002/lidong/SNP_filter_bin/filter_rst_ratio.pl /share/raid002/lidong/40snp/population chr1 >/share/raid002/lidong/40snp/population/chr1.population.snp
