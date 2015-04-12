#!/bin/sh

git pull
git submodule update
git submodule foreach git pull
git submodule foreach git checkout master

#git submodule foreach git fsck
#git submodule foreach git gc --aggressive --prune=now

#git submodule foreach git checkout master
git submodule
git diff .

git subtree pull --prefix=c_cpp/lib/cstring git://github.com/cloudwu/cstring.git master
git subtree pull --prefix=c_cpp/lib/htslib git://github.com/samtools/htslib.git master
git subtree pull --squash --prefix=c_cpp/lib/gnulib git://git.savannah.gnu.org/gnulib.git master
git gc --aggressive --prune=now
