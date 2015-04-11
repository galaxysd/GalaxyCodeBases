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
