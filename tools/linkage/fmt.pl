#!/bin/env perl
use strict;
use warnings;

my $minMAF=0.2;
my $minCOV=0.5;

unless (@ARGV > 0) {
    print "perl $0 <In_file> <Out_File> <ChrID>\n";
    exit 0;
}

open I,'<',$ARGV[0] or die $!;
open O,'>',$ARGV[1] or die $!;
my $t=<I>;
$t=~s/^#//;
print O "#ChrID\tPos\t$t";
while (<I>) {
    chomp;
    my ($Pos,$D)=split /\s+/,$_,2;
    my @Dat=split / /,$D;
    my ($cvg,%cnt)=(0);
    for (@Dat) {
        ++$cnt{$_} if /\d/;
        ++$cvg if /\d/;
    }
    next if $cvg < @Dat * $minCOV;
    ($cvg)=sort {$a<=>$b} values %cnt;
    next if $cvg < @Dat * $minMAF;
    print O "$ARGV[2]\t$Pos\t$D\n"
}
close O;
close I;
__END__
./fmt.pl ./p15/Chr01 Chr01.15 Chr01
cat chrorder | while read a; do ./fmt.pl ./p15/$a $a.15 $a;done &

