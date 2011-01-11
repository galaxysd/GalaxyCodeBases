#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use GalaxyXS::ChromByte 1.02;# ':all';

unless (@ARGV > 0) {
    print "perl $0 <chr.nfo> <fabychr_path> <marker_len> <marker_file> <parent_snp> <out_file>\n";
    exit 0;
}

my ($chrnfo,$fabychr_path,$marker_len,$marker_file,$psnp_file,$out_file)=@ARGV;
$marker_len=int $marker_len;
die "[x]marker_len must > 10 !\n" if $marker_len < 10;
open M,'<',$marker_file or die $!;
<M>;
my $t=tell M;
$_=<M>;
my $ChrID=(split /\t/)[0];
seek M,$t,0;

my ($ChrLen,$ChrMem,$SNPMem)=(0);
open ChrNFO,'<',$chrnfo or die "[x]Error opening $chrnfo: $!\n";
while (<ChrNFO>) {
	#next if /^#/;
	my ($chrid,$len)=split /\t/;
	next if $chrid ne $ChrID;
	$ChrLen=$len;
}
close ChrNFO;
die "[x]Cannot find correct ChrLen for [$ChrID] !\n" if $ChrLen<1;
$ChrMem=initchr($ChrLen);
$SNPMem=initchr($ChrLen);	# OK, leave the 8th bit and alloc a new byte !

my ($seq,@seq);
open F,'<',$fabychr_path."/$ChrID.fa" or die $!;
while(<F>) {
	chomp;
	my ($id)=split / /,$_,2;
	$id =~ s/^>//;
	die "[x]Wrong RefSeq ! ($ChrID != $id)\n" if $ChrID ne $id;
	$/=">";
	$seq=<F>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
}
close F;

@seq=split //,$seq;
$seq='';

my %bIUB = ( A => 1,
	     C => 4,
	     G => 8,
	     T => 2,
	     M => 5,
	     R => 9,
	     W => 3,
	     S => 12,
	     Y => 6,
	     K => 10,
	     V => 13,
	     H => 7,
	     D => 11,
	     B => 14,
	     N => 15
	     );
my %bIUBex = (
	     U => 2,
	     X => 15,
	     );
my (%ubIUB,%rbIUB);
#$rbIUB{$bIUB{$k}}=$_ for keys %bIUB;
@rbIUB{ values %bIUB } = keys %bIUB;

#@hash1{ keys %hash2 } = values %hash2;
@bIUB{ keys %bIUBex } = values %bIUBex;

for my $k (keys %bIUB) {
	$ubIUB{$k}=16*$bIUB{$k};
}

for my $pos (1..$ChrLen) {
	my $v=0;
	my $base=shift @seq;	# to save running memory as $seq[$pos-1]
	if (exists $bIUB{$base}) {
		$v=$bIUB{$base};
#warn $seq[$pos-1]," $v\t$pos\n";
	} else {
		warn "[!]$pos -> ",$base;
	}
	setbase($ChrMem,$pos,$v) if $v;
}
@seq=();

open PSNP,'<',$psnp_file or die "[x]Error opening $psnp_file: $!\n";
while (<PSNP>) {
	next if /^#/;
	my ($chr,$pos,$ref,$ga,undef,undef,$gb)=split /\s+/;
	my @bases=split //,$ga.$gb;
	my $v;
	$v |= $bIUB{$_} for @bases;
#print STDERR "$ref [$ga.$gb] $v\t";
	#$v = orbase($ChrMem,$pos,$v) if $v && ($v != 240);
#warn "$v\n";
	setbase($SNPMem,$pos,$v) if $v && ($v != 15);
}
close PSNP;

my @poses;
open O,'>',$out_file or die $!;
while (<M>) {
	my $pos=(split /\t/)[1];
=pod
	my ($left,$right)=($pos-$marker_len,$pos+$marker_len);
	my $marker_seq=substr($seq,$left-1,$marker_len).'N'.substr($seq,$pos,$marker_len);
	print O ">${ChrID}_${pos}m\n$marker_seq\n";
=cut
#	my $v=getbase($ChrMem,$pos);
#my $base=$v & 15;
#	my $snp=($v & 240) >> 4;	# $v >> 4 is wrong !
my $base=getbase($ChrMem,$pos);
	my $snp=getbase($SNPMem,$pos);
#warn "[$snp][$base]\t$pos\n";
	if ($snp) {
		setbase($ChrMem,$pos,$snp);
		setbase($SNPMem,$pos,15);
	}
	push @poses,$pos;
}
close M;
for my $pos (@poses) {
	my ($left,$right)=($pos-$marker_len,$pos+$marker_len);
	my $marker_seq;
	for my $ipos ($left..$right) {
		#my $v=getbase($ChrMem,$ipos);
		my $base=getbase($ChrMem,$ipos);
		my $snp=getbase($SNPMem,$ipos);
		$base=$rbIUB{$base};
		if ($snp==15) {	#  or $ipos == $pos
			$base = lc $base;
		}
		$marker_seq .= $base;	#."($snp)";
	}
	print O ">${ChrID}_${pos}m${marker_len}\n$marker_seq\n";
}
close O;

__END__
./exmarkerseq.pl chr.nfo ./fabychr/ 45 ./dat20101214/Chr01.marker ./psnp/parent_Chr01.snp eChr01.marker
cat chrorder | while read a; do ./exmarkerseq.pl chr.nfo ./fabychr/ 45 ./dat20101214/${a}.marker ./psnp/parent_$a.snp ex45${a}.marker;done &

#!/bin/sh
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH,BLASTDB
#$ -cwd -r y -l vf=8g
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-12
SEEDFILE=./chrorder
a=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
eval ./exmarkerseq.pl chr.nfo ./fabychr/ 25 ./dat20101214/${a}.marker ./psnp/parent_$a.snp ./out20110106/ex25${a}.marker
