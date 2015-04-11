#!/bin/sh

git pull
git submodule update
git submodule foreach git pull
git submodule foreach git checkout master

git fsck
git gc
#git submodule foreach git fsck
git submodule foreach git gc

git submodule foreach git checkout master
git submodule

git subtree pull --prefix=c_cpp/lib/cstring git://github.com/cloudwu/cstring.git master
git subtree pull --prefix=c_cpp/lib/htslib git://github.com/samtools/htslib.git master
