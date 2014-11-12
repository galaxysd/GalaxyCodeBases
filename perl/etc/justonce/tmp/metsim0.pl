#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $HomSNPrate = 0.0005;
my $HetSNPrate = 0.0005;
my $TranStoV = 4;	# 转换比颠换，transitions，transversions
my ($minRefLen,$maxRefLen) = (40,220);
my $infilename = 'hg19chr18.bed.frag.trim';

open I,'<',$infilename or die;

my $lines = 0;
while (<I>) {
	chomp;
	my ($chr,$s,$e,$len) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	++$lines;
}
my $each = int($lines/10);

print "Lines of [$minRefLen,$maxRefLen]: $lines\nEach of 10 win: $each\n\n";

seek I,0,0;

my ($len0,$s0) = 0;
while ($len0 < $minRefLen or $len0 > $maxRefLen) {
	$_=<I>;
	(undef,$s0,undef,$len0) = split /\t/;
}
my $i=1;

my @Rate = (
	[.1,.1],[1,0],[0,1],[.5,.1],[.1,.5],
	[.7,.3],[.3,.7],[.1,.5],[.5,.1],[.5,.5]
);
my $k=0;

open O,'>','zone.lst' or die;
print STDERR "[!]Writing zone.lst:\n";
my ($chr,$s,$e,$len,$ee);
while (<I>) {
	chomp;
	($chr,$s,$e,$len) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	++$i;
	unless ($i % $each) {
		#warn "[$k]\n";
		if ($k <= $#Rate) {
			print join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
			print O join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
		}
		do {
			$_=<I>;
			last unless defined $_;
			(undef,$s0,$ee,$len) = split /\t/;
		} while ($len < $minRefLen or $len > $maxRefLen);
		++$i;
		++$k;
	}
}
$k = $#Rate;
$e = $ee if $e < $ee;
print join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
print O join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";

#close I;
close O;
warn "[!]done !\n\n";

#########

sub getSNP($) {
	my $inbase = $_[0];
	my $outBase = $inbase;
	my $SorV = rand(1);
	if ($SorV > $TranStoV/(1+$TranStoV)) {	# transversions
		my $AB = rand(1);
		if ($AB > 0.5) {
			$outBase =~ tr/AGCTagct/CCAAccaa/;
		} else {
			$outBase =~ tr/AGCTagct/TTGGttgg/;
		}
	} else {	# transitions
		$outBase =~ tr/AGCTagct/GATCgatc/;
	}
	return $outBase;
}

my $SNPremained = 0;
my %usedPos;
sub addSNP($$$$$) {
	my ($chr,$s,$seq,$rate,$type) = @_;
	my $len = length $seq;
	my $SNPcount = $len * $rate + $SNPremained;
	$SNPremained = $SNPcount - int($SNPcount);
	my $newseq = $seq;
	for ( 1 .. int($SNPcount) ) {
		my $Pos;
		do {
			$Pos = int(rand($len));
		} until ( ! exists $usedPos{$chr}{$s}{$Pos} );
		$usedPos{$chr}{$s}{$Pos} = $type;
		my $oldbase = substr $newseq,$Pos,1;
		my $newbase = getSNP($oldbase);
		$newseq = $seq;
		substr $newseq,$Pos,1,$newbase;
		print S join("\t",$type,$chr,$s,$s+$len-1,$s+$Pos,$oldbase,$newbase ,$Pos+1,$SNPcount,$seq,$newseq),"\n";
	}
	return $newseq;
}

seek I,0,0;
print STDERR "[!]Writing snp.lst: ...";
open S,'>','snp.lst';
while(<I>) {
	chomp;
	my ($chr,$s,$e,$len,$seq) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	my $seq1 = addSNP($chr,$s,$seq,$HomSNPrate,'Hom');
	my $seq2 = addSNP($chr,$s,$seq1,$HetSNPrate,'Het');
}

close I;
close S;
warn "\b\b\bdone !\n";

__END__
scp sim0.pl gaoshengjie@192.168.9.100:/ifs4/BC_CANCER/PROJECT/qbe130308_gaosj_Methylation/sim
