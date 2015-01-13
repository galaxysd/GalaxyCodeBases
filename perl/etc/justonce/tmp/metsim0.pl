#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $HomSNPrate = 0.0005;
my $HetSNPrate = 0.0005;
my $TranStoV = 4;	# 转换比颠换，transitions，transversions
my ($minRefLen,$maxRefLen) = (40,220);
my $infilename = 'hg19chr17.bed.frag';

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
	[.1,.1,0,0],[1,0,0,0],[0,1,0,0],[.5,.1,1,0],[.1,.5,0,-1],
	[.7,.3,0,0],[.3,.7,0,0],[.1,.5,0,1],[.5,.1,1,0],[.5,.5,1,1]
);
# 用bisofite处理后，含甲基化的C不会变成T，不含甲基化的C变成T。由于他只能吧C变成T，所以，正链就是C->T，而他的反义互补连则是负连的C->T，反映到正链上就是G->A。
# 比如CG，如果他在watson连上并且没有甲基化，那么bisofte处理后就变成TG，如果是在crick链上变成CA。
# 这里的比例是甲基化率Rm，分别是总SNP的C（癌组织）、N（癌旁组织）、杂合SNP与甲基化相关度（-1=>只有Ref甲基化，+1=>只有Mut甲基化，0=>无关）。后两列只对前面0.5的有效。
#   ref一样的全部非甲基化，variacen的read全部甲基化。这样甲基化绿刚好是0.5。

my $k=0;

open O,'>','zone.lst' or die;
print STDERR "[!]Writing zone.lst:\n";
print O join("\t",'#No','FirstBegin','LastEnd','Rm(C)','Rm(N)','Het2MetR(C)','Het2MetR(N)'),"\n### Het2MetR: -1=>只有Ref甲基化，+1=>只有Mut甲基化，0=>无关。此时不考虑纯合SNP。\n";
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
#$k = $#Rate;
#$e = $ee if $e < $ee;
#print join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
#print O join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";

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

#my $SNPremained = 0;	# 这样意味着Hom与Het共用这个变量。但各自独立的话，每条序列都有两种SNP了。>_< !
my %usedPos;
sub addSNP($$$$$$) {
	my ($chr,$s,$seq,$rate,$type,$SNPremained) = @_;
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
		print S join("\t",$type,$chr,$s,$s+$len-1,$s+$Pos,$oldbase,$newbase,$Pos+1,$SNPcount,$seq,$newseq),"\n";	# multiple lines if more than 1 SNP.
	}
	return [$SNPremained,$newseq];
}

seek I,0,0;
print STDERR "[!]Writing snp.lst: ...";
open S,'>','snp.lst';
my ($stat1,$stat2)=(0,0);
while(<I>) {
	chomp;
	my ($chr,$s,$e,$len,$seq) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	my ($seq1,$seq2);
	($stat1,$seq1) = @{addSNP($chr,$s,$seq,$HomSNPrate,'Hom',$stat1)};
	($stat1,$seq2) = @{addSNP($chr,$s,$seq1,$HetSNPrate,'Het',$stat1)};
}

close I;
close S;
warn "\b\b\bdone !\n";

__END__
scp sim0.pl gaoshengjie@192.168.9.100:/ifs4/BC_CANCER/PROJECT/qbe130308_gaosj_Methylation/sim
