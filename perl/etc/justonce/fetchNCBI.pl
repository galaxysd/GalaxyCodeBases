#!/usr/bin/env perl
use strict;
use warnings;
use Mozilla::CA;
use LWP::Simple;
use Data::Dump qw(ddx);

my $deltra = 200;

my $infile = '/Users/galaxy/Downloads/snpcandidatforpcr.out';
open IN,'<',$infile or die $!;
open OUT,'>','pcr.fa' or die $!;
while (<IN>) {
	chomp;
	my @dat = split /\t/;
	next if @dat < 3;
	my (undef,$id) = split /\./,$dat[0];
	print "$id, $dat[2], $dat[3]\n";
	my ($left,$right) = ($dat[2]-$deltra,$dat[2]+$deltra);
	my $url = 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id='.$id.'&db=nuccore&report=fasta&fmt_mask=0&maxplex=1&sendto=t&from='.$left.'&to='.$right.'&maxdownloadsize=1000000';
	my $content = get($url);
	my @ret = split /\n/,$content;
	my $head = shift @ret;
	$head =~ s/^>//;
	my $out = ">${id}_$dat[2] $head\n" . join('',@ret);
	print OUT "$out\n\n";
}
close OUT;
close IN;

__END__

my $content = get('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=DS362220&db=nuccore&report=fasta&fmt_mask=0&maxplex=1&sendto=t&from=121868&to=122268&maxdownloadsize=1000000');
#die "Couldn't get it!" unless defined $content;

print "[$content]\n";
