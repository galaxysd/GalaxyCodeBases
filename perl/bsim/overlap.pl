#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <SimLst> <AnalyseRes>\n" if @ARGV <2;

my $total=shift;
my $result=shift;
open TT,$total or die $!;
open RT,$result or die $!;
my %hash;

while(<TT>){
        chop;
        my @a=split;
        #print $a[3]."\n";
        next unless(/\w/);
        for( my $kk=$a[3]-10;$kk<$a[3]+10;$kk++){
                $hash{$kk}=1;
        }

}
close TT;
while(<RT>){
        chomp;
        my @a=split;
        #print $a[2]."\n";
        if(/RefCut/){
        if(exists($hash{$a[2]})){
                print $_."\n";
        }
        }


}
close RT;

__END__
./overlap.pl simed.lst simVir4_analyseAll.txt > simgot.lst

grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

