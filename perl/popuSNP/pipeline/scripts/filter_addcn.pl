#!/bin/env perl
use strict;
use warnings;
# original from HeWeiMing;

unless(@ARGV) {
	print "Usage: $0 <config file> <InPut add_cn>\n";
	print " this is used config file!\n";
	print 	"config file:\n";
	print 	"depth	min	max\nquality	value\nlc	value\nrst	value\ncopynumber	value\n";
	exit;
}

my ($conf,$addcn)=@ARGV;

my $depth_min = 20;
my $depth_max = 450;
my $quality = 15;
my $lc = 2;
my $rst = 0.005;
my $cpnum = 1.5;

my (%Dat,@t,$v);
open FC,'<',$conf or die "[x]Error opening $conf: $!\n";
while(<FC>) {
	next if /^[#;]/;
	chomp;
	my @t=split /\t/;
	$v=shift @t;
	$Dat{$v}=\@t;
}
($depth_min,$depth_max) = @{$Dat{depth}};
$quality = $Dat{quality}->[0];
$rst     = $Dat{rst}->[0];
$cpnum   = $Dat{copynumber}->[0];

warn "depth: $depth_min $depth_max\n";
warn "quality: $quality\n";
warn "rst: $rst\n";
warn "copynumber: $cpnum\n";

close FC;

my $num=0;
open ADDCN,'<',$addcn or die "[x]Error opening $addcn: $!\n";
open OUT,'>',"${addcn}.filter"|| die "$!\n" ;
open OUTdb,'>',"${addcn}.dbsnp"|| die "$!\n" ;
while(<ADDCN>) {
	chomp ;
	my @inf=split;
	next if($inf[3]<$depth_min || $inf[3]>$depth_max);
	next if($inf[9]<$quality);
	next if($inf[10]>$cpnum);
	next if($inf[12]==0);
	next if($inf[13]==0);
#	next if($inf[-2]>$lc);
	next if ($inf[-1]<$rst);
	print OUT "$_\n";
#	print $_ ; exit;
	my @l;
	my %hash;
	$num++;
	push @l,$inf[0],$inf[1];
	$l[2]=1;
	$l[3]=0;
	$l[4]=0;
	$hash{A}=0;
	$hash{C}=0;
	$hash{T}=0;
	$hash{G}=0;
	$hash{$inf[5]}=0.5;
	$hash{$inf[6]}=0.5;
	my $ID = '0000000'.$num;
	$ID=substr($ID,-8);	# Don't you know sprintf %08d ...
	$ID="$ID";
	my $out=join "\t",@l,$hash{A},$hash{C},$hash{T},$hash{G},$ID;
	print OUTdb $out,"\n";
}

close  OUT;
close  OUTdb,
close  ADDCN;
