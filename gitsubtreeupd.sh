#!/bin/sh

git pull
git diff .

# `git subtree add` when it does not exist.
git subtree pull --prefix=c_cpp/lib/cstring git://github.com/cloudwu/cstring.git master
git subtree pull --prefix=c_cpp/lib/htslib git://github.com/samtools/htslib.git master
git subtree pull --prefix=c_cpp/lib/uthash https://github.com/troydhanson/uthash.git master
git subtree pull --prefix=c_cpp/lib/inih https://github.com/benhoyt/inih.git master
git subtree pull --prefix=c_cpp/lib/klib https://github.com/attractivechaos/klib.git master

# 已被和谐 git subtree pull --prefix=etc/gfwlist2pac git://github.com/clowwindy/gfwlist2pac.git master
git subtree pull --prefix=etc/agfw/genpac https://github.com/JinnLynn/genpac.git master
git subtree pull --squash --prefix=etc/agfw/autoproxy-gfwlist git://github.com/calfzhou/autoproxy-gfwlist.git master
git subtree pull --prefix=etc/agfw/Pcap_DNSProxy git://github.com/chengr28/Pcap_DNSProxy.git master

git subtree pull --prefix=tools/bioawk git://github.com/lh3/bioawk.git master
git subtree pull --prefix=tools/lh3misc https://github.com/lh3/misc.git master

# GitCafe 将于 2016年5月31日 停止所有服务 git subtree pull --prefix=pubdoc/Markdown-Syntax-CN https://gitcafe.com/riku/Markdown-Syntax-CN.git master
git subtree pull --prefix=etc/whichpm https://github.com/mklement0/whichpm.git master
git subtree pull --prefix=tools/create-dmg https://github.com/andreyvit/create-dmg.git master

git subtree pull --prefix=tools/crawler/xeHentai https://github.com/fffonion/xeHentai.git master
# https://github.com/fffonion/httplib2-plus.git mirrored as https://github.com/galaxy001/httplib2-plus

git subtree pull --prefix=tools/Bilibili/danmaku2ass https://github.com/m13253/danmaku2ass.git master
git subtree pull --prefix=tools/Bilibili/danmaku2ass_native https://github.com/typcn/danmaku2ass_native.git master
git subtree pull --prefix=tools/Bilibili/Biligrab https://github.com/cnbeining/Biligrab.git master
git subtree pull --prefix=tools/Bilibili/BiliDan https://github.com/m13253/BiliDan.git master
git subtree pull --prefix=tools/Bilibili/You-Get https://github.com/soimort/you-get.git develop

git subtree pull --prefix=tools/Crypters/lvdo https://github.com/m13253/lvdo.git master
git subtree pull --prefix=tools/wtb https://github.com/mariotaku/wtb.git master

git subtree pull --prefix=js/TagSelect.js https://github.com/graphiq/TagSelect.js.git master
git subtree pull --prefix=js/ThunderLixianExporter https://github.com/binux/ThunderLixianExporter.git master

git subtree pull --prefix=etc/Chrome/BaiduExporter https://github.com/acgotaku/BaiduExporter.git master
# git subtree pull --prefix=etc/Chrome/ThunderLixianAssistant https://github.com/ohsc/ThunderLixianAssistant.git master
# etc/Chrome/ThunderLixianAssistant/src/ThunderLixianExporter -> ../../../../js/ThunderLixianExporter
git subtree pull --prefix=etc/Chrome/115 https://github.com/acgotaku/115.git master

git subtree pull --prefix=etc/Mac/OfficeThinner https://github.com/goodbest/OfficeThinner.git master
git subtree pull --prefix=etc/Mac/adobe_case_sensitive_volumes https://github.com/tzvetkoff/adobe_case_sensitive_volumes.git master
git subtree pull --prefix=etc/Windows/vlmcsd https://github.com/vancepym/vlmcsd.git master

# BioInfo
git subtree pull --prefix=BioInfo/BS-Seq/bwa-meth https://github.com/brentp/bwa-meth.git master


git fsck
git gc --aggressive --prune=now



# git subtree pull --squash --prefix=c_cpp/lib/gnulib git://git.savannah.gnu.org/gnulib.git master 由于服务器“no common commits”，每次更新都得下载全部历史，所以直接gc会变6G的.git目录。
# git subtree pull --squash --prefix=c_cpp/lib/gnulib git://github.com/coreutils/gnulib.git master 仍然“no common commits”，可能本地必须有上个关键commit，吧？
# git subtree pull --prefix=c_cpp/fromOthers/PEAR git://github.com/xflouris/PEAR.git master 测试数据太大。所以tbz2打包保存。

# http://stackoverflow.com/questions/21702996/git-subtree-add-change-prefix-preserving-local-commits?lq=1
# One alternative is splitting from your current commit that include the local commits you want to preserve:
#    git subtree split --prefix=dir1 HEAD
# Create a branch with the printed commit(sha1) just to use it later
#    git branch split_dir_1 <split_commit>
# And then delete the subdirectory and re-add the subtree.
#    git rm -r dir1
#    git commit
#    git subtree add --prefix=dir2 . <split_commit>

# In case `Working tree has modifications. Cannot add.`
# git diff-index HEAD
