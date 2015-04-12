#!/bin/sh

git pull
git diff .

git subtree pull --prefix=c_cpp/lib/cstring git://github.com/cloudwu/cstring.git master
git subtree pull --prefix=c_cpp/lib/htslib git://github.com/samtools/htslib.git master
git subtree pull --squash --prefix=c_cpp/lib/gnulib git://git.savannah.gnu.org/gnulib.git master
git gc --aggressive --prune=now
