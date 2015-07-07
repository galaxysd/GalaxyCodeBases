#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Ref> <Virus> <Outprefix>\n" if @ARGV <3;

