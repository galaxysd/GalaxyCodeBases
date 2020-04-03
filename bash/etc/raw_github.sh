#!/bin/bash
# data: 2020-03-31
# author: muzi502
# for: Fuck GFW and download some raw file form github without proxy using jsDelivr CDN 
# usage: save the .she to your local such as /usr/bin/rawg, and chmod +x /usr/bin/rawg
# use rawg https://github.com/ohmyzsh/ohmyzsh/blob/master/tools/install.sh to download

# URL: https://gist.github.com/muzi502/c1bd677c5b4c41115dcfff6e724dee8b

set -xue
# GitHub rul: https://github.com/ohmyzsh/ohmyzsh/blob/master/tools/install.sh
# jsDelivr url: https://cdn.jsdelivr.net/gh/ohmyzsh/ohmyzsh/tools/install.sh

wget $(echo $1 | sed 's/raw.githubusercontent.com/cdn.jsdelivr.net\/gh/' \
               | sed 's/github.com/cdn.jsdelivr.net\/gh/' \
               | sed 's/\/master//' | sed 's/\/blob//' )

# curl $(echo $1 | sed 's/raw.githubusercontent.com/cdn.jsdelivr.net\/gh/' \
#                | sed 's/github.com/cdn.jsdelivr.net\/gh/' \
#                | sed 's/\/master//' | sed 's/\/blob//' )