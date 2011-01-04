#!/bin/env perl
use strict;
use warnings;

unless (@ARGV > 0) {
    print "perl $0 <In_file> <Out_File>\n";
    exit 0;
}

# choose Kosambi, which is slower than Haldane, for density markers
sub corRIL($) {
	my ($bigRth)=@_;
	my $smallr;
	$smallr=0.5*$bigRth/(1-$bigRth);
	# applies the distance function
	my $RILdist = 25 * log((1 + 2 * $smallr) / (1 - 2 * $smallr));
	return $RILdist;
}

open I,'<',$ARGV[0] or die $!;
open O,'>',$ARGV[1] or die $!;
my $t=<I>;
$t=~s/\tPos\t/\tPos\tcM\tΔcM\tkB\/cM\t/;
#print O $t;
print O "#ChrID\tPos\tcM\tΔcM\tkB/cM\tR_RIL\n";
my (@Pos,%Dat);
my ($chr,$pos,$dat);
while (<I>) {
	chomp;
	($chr,$pos,$dat)=split /\t/;
	my @a=split / /,$dat;
	$Dat{$pos}=\@a;
	push @Pos,$pos;
}
my $LastcM=0;
my $ratio;
print O "$chr\t$Pos[0]\t$LastcM\t0\tNA\t0\n";	#,join(' ',@{$Dat{$Pos[0]}}),"\n";
for my $i (1..$#Pos) {
	my $ary_a=$Dat{$Pos[$i-1]};
	my $ary_b=$Dat{$Pos[$i]};
	my ($cross,$cnt)=(0,0);
	for my $p (0..$#$ary_a) {
		next if $ary_a->[$p] eq '-' or $ary_b->[$p] eq '-';
		++$cross if $ary_a->[$p] ne $ary_b->[$p];
		++$cnt;
	}
	my $cross0 = $cross/$cnt;
	$cross = &corRIL($cross0);
	$LastcM += $cross;
	$ratio=$cross? int(0.5+100*($Pos[$i]-$Pos[$i-1])/(1000*$cross)) : 'NA';
	print O "$chr\t$Pos[$i]\t",int(0.5+100*$LastcM)/100,"\t",int(0.5+100*$cross)/100,"\t"
		,$ratio,"\t",100*$cross0,"\n";	#,join(' ',@{$Dat{$Pos[$i]}}),"\n";
}

close O;
close I;
__END__
cat chrorder | while read a; do ./clinkage.pl $a.ccq cm$a.ccq;done &

cat chrorder | while read a; do ./linkage/clinkage.pl dat20101214/$a.marker dat20101214/$a.cm;done &
