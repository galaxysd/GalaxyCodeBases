#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SampleList = 'samples.lst';
my $SNPtsv = 'micro.snp.tsv';

open I,'<',$SampleList or die $?;

