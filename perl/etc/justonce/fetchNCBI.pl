#!/usr/bin/env perl
use strict;
use warnings;
use Mozilla::CA;
use LWP::Simple;
use Data::Dump qw(ddx);

my $content = get('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=DS362220&db=nuccore&report=fasta&fmt_mask=0&maxplex=1&sendto=t&from=121868&to=122268&maxdownloadsize=1000000');
#die "Couldn't get it!" unless defined $content;

print "[$content]\n";
