#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
use Data::Dump qw(ddx);
#use Galaxy::IO;
#use Galaxy::SeqTools;

die "Usage: $0 <vcf.gz with tabix indexed> <out> <regions> (eg.: scaffold75:924209-5441687)\n" if @ARGV<2;
my $vcfs = shift;
my $outfs = shift;
my $regions = shift;

warn "From [$vcfs] to [$outfs]\n";
warn "Regions:[$regions]\n" if $regions;

my $vcf = Vcf->new(file=>$vcfs,region=>$regions);
$vcf->parse_header();
my (@samples) = $vcf->get_samples();
ddx \@samples;

while (my $x=$vcf->next_data_hash()) { 
	next if $$x{QUAL} < 20;
	my ($flag,%GTs,$gtREC)=(0);
	for my $gt (keys %GTs) {
		my ($a1,$a2,$a3) = $vcf->split_gt($gt);
		if ($a3 or ($a1 != $a2)) {
			$flag = 1;
		}
	}
	#ddx [\%GTs,$x] if $flag;
	print $vcf->format_line($x);
}
