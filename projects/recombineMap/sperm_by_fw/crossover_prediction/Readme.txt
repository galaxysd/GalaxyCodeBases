简便方法：
makedata.py和makedata_adjust.py(低cut off，结果更多)
参数是分染色体SNP文件路径，如./当前目录，文件识别特征请在程序第7行修改中，现在识别文件名中有"FaMo"的文件。
然后用sort.py 进行排序。
例：
python makedata.py ./
python sort.py all.sum.changepos.standard


HMM方法：
HMMCO.py是核心程序，第一个参数同上。
运行出来大量，result.SXX.chrN文件，然后建个目录mkdir data，把文件移动进去mv result.* data，运行sumData.py进行汇总，第一个参数同上。
得到结果叫Result.all.txt，运行sort.py排序。
其他几个程序是对比简便方法结果的。
可以使用
work.sh：
pypy HMMCO.py ../findpos #改成SNP文件目录
mv result.* data
pypy sumData.py ../findpos
pypy sort.py Result.all.txt
pypy compareData.py raw.data Result.all.txt #第一个参数改成以前的结果，第二个改成HMM结果。
