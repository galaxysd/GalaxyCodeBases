#!/usr/bin/env bash

# https://gist.github.com/Shians/a95c49bda68481fb0adc83600d9eb6b8
# called by
# sh extract_bam_tag.sh input.bam BC
# to print read_id and BC tag value

# two arguments, a bam file and the tag to extract
BAM=$1
TAG=$2

# write a tsv with columns read_id and tag value
echo -e "read_id\t$TAG"
samtools view "$BAM" | grep "$TAG:." | perl -pe 's/(^.+?)\t.*'$TAG':.:(.+?)\t.*/$1\t$2/g'
# regular expression substutes pattern (read_id)*(tag_value)* for (read_id)\t(tag_value)
