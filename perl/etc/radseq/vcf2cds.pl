#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
use Data::Dump qw(ddx);
#use Galaxy::IO;
#use Galaxy::SeqTools;

die "Usage: $0 <vcf.gz with tabix indexed> <out> <regions> (eg.: scaffold75:924209-5441687)\n" if @ARGV<2;
my $fafs='/bak/seqdata/2012/tiger/120512_TigerRefGenome/bychr/scaffold75.fa';
my $gtfs = 'scaffold75.gtf';
my $vcfs = shift;
my $outfs = shift;
my $regions = shift;

warn "From [$vcfs][$gtfs][$fafs] to [$outfs]\n";
warn "Regions:[$regions]\n" if $regions;

my $code=<<Ecode;
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M---------------M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
Ecode
my @t=split /\n/,$code;
my (%t,%gen_code,%start_codes);
map {s/\s//g;@_=split /=/;$t{$_[0]}=[split //,$_[1]]} @t;
#ddx \%t;
for (@{$t{AAs}}) {
#	my $aa=$_;
#	my @bases;
#	map {my $t=shift @{$t{$_}};push @bases,$t} qw/Base1 Base2 Base3/;
#	my $base=join('',@bases);
	my $base=(shift @{$t{Base1}}).(shift @{$t{Base2}}).(shift @{$t{Base3}});
#print "$base";
	my $start=shift @{$t{Starts}};
	++$start_codes{$base} if $start eq 'M';
	$gen_code{$base}=$_;
#print " -> [$_]\n";
}
#ddx [\@t,\%gen_code,\%start_codes];
print "Codon Table:\n$code\n";

my %GeneDat:
open GTF,'<',$gtfs or dir $!;
while(my $in_line = <GTF>){
	next if ($in_line =~ /^\s*\#/);
	chomp;
      #remove leading whitespace and get comments
      my $comments = "";	  
      while($in_line =~ /^(.*)\#(.*)$/){
        $in_line = $1;
        $comments = $2.$comments;
      }
      if($in_line =~ /^\s+(.*)$/){
        $in_line = $1;
      }
      my @data = split /\s+/,$in_line;
      #verify line is correct length
      next if ($#data < 8);
}
close GTG;
warn "GTF loaded.\n";
ddx \%GeneDat;

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
	ddx [\%GTs,$x] if $flag;
	#print $vcf->format_line($x);
}

close GTFERR;
close TX;
