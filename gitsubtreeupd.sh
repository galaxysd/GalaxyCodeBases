#!/bin/sh

git pull
git diff .

git subtree pull --prefix=c_cpp/lib/cstring git://github.com/cloudwu/cstring.git master
git subtree pull --prefix=c_cpp/lib/htslib git://github.com/samtools/htslib.git master
git fsck
git gc --aggressive --prune=now



# git subtree pull --squash --prefix=c_cpp/lib/gnulib git://git.savannah.gnu.org/gnulib.git master 由于服务器“no common commits”，每次更新都得下载全部历史，所以直接gc会变6G的.git目录。
# git subtree pull --squash --prefix=c_cpp/lib/gnulib git://github.com/coreutils/gnulib.git master 仍然“no common commits”，可能本地必须有上个关键commit，吧？
# git subtree pull --prefix=c_cpp/fromOthers/PEAR git://github.com/xflouris/PEAR.git master 测试数据太大。所以tbz2打包保存。
