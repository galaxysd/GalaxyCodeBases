#!/bin/env perl
use strict;
use warnings;

my $minMAF=0.2;
my $minCOV=0.5;

unless (@ARGV > 0) {
    print "perl $0 <In_file> <Out_prefix> <ChrID>\n";
    exit 0;
}
=pod
my %varToGT=(
	0   => 'A',
	1   => 'B',
	'-' => '-',
	2   => 'H',
);
=cut
my ($Count,$population)=(0);
my (%DAT,%deDup);
open I,'<',$ARGV[0] or die $!;
my $t=<I>;
$t=~s/^#//;
while (<I>) {
    chomp;
    my ($Pos,$D)=split /\s+/,$_,2;
	$DAT{$Pos}=$D;
	++$Count;
}
close I;
$population=split(/\s+/,(values %DAT)[0]);
open O,'>',$ARGV[1].'.intercross' or die $!;
open OO,'>',$ARGV[1].'.order' or die $!;
print O "data type f2 intercross\n$population $Count\n";
for my $pos (sort {$a<=>$b} keys %DAT) {
	my $str=$DAT{$pos};
	$str =~ tr/012/ABH/;
	print O $ARGV[2],"M$pos\t",join("\t",split /\s+/,$str),"\n";
	print OO '*',$ARGV[2],"M$pos\n";
}
close OO;
close O;
__END__
./tmapfmt.pl ./p15/Chr01 Chr01.intercross Chr01
cat chrorder | while read a; do ./tmapfmt.pl ./dat2/${a}_num $a $a;done &

./tmap-1.0/pedmerge ril1206.Chr01.intercross tt
cat chrorder | while read a; do ./tmap-1.0/pedmerge ./marker/$a.intercross ./marker/$a.pedigree;done &
