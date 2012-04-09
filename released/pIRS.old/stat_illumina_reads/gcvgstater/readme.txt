illumia_reads_parameter_stator
	by Hu Xuesong @ BGI <huxuesong@genomics.org.cn>

GC-Depth Stat:
1. run soap and soap.coverage to get .depthsingle file(s). gzip is OK to over it.
2. run gc_coverage_bias on all depthsingle files. You will get gc-depth stat by 1 GC% and other files.
3. run gc_coverage_bias_plot on the gc-depth stat file. You'll get PNG plot and a .gc file by 5 GC%.
4. Manually check the .gc file for any abnormal levels due to the lower depth on certain GC% windows.

Error Matrix Stat:
1. run soap or bwa to get .{soap,single} or .sam file(s).
2. run error_matrix_calculator on those file(s). You will get *.{count,ratio}.matrix .
3. You can use error_matrix_merger to merge several .{count,ratio}.matrix files.
   However, it is up to you to keep the read length matches.

Insert size & mapping ratio stat:
1. run soap or bwa to get .{soap,single} or .sam file(s).
2. run alignment_stator *.
* alignment_stator cannot stat. mapping ratio for sam files now.

c7618d69b6a8786715431fbfc33dc20e *alignment_stator
c7618d69b6a8786715431fbfc33dc20e *alignment_stator.pl
6f68d7dfc60a523557953f67b3281f97 *error_matrix_analyzer
6f68d7dfc60a523557953f67b3281f97 *error_matrix_analyzer.pl
8b921ff206f40fc608f06abb126842b0 *error_matrix_calculator
8b921ff206f40fc608f06abb126842b0 *error_matrix_calculator.pl
101728bda732f984c05513671785881b *error_matrix_merger
101728bda732f984c05513671785881b *error_matrix_merger.pl
e5a51d46f5f676057e4123136ba4cae2 *gc_coverage_bias
e0ef280b60477f828bff8d0653dc4e9c *gc_coverage_bias_plot
e0ef280b60477f828bff8d0653dc4e9c *gc_coverage_bias_plot.pl
78f578c83fed400dca2fb23c6e5068d8 *readme.txt
