#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
    print "perl $0 <psnp.lst> <out.maf>\n";
    exit;
}

my (%Files,%Dat,%DatAll,@Chr);
open L,'<',$ARGV[0] or die "[!]Error opening $ARGV[0]: $!\n";
while(<L>) {
	chomp;
	my ($chr,$file)=split /\t/;
	push @Chr,$chr;
	$Files{$chr}=$file;
}
close L;

for my $fc (keys %Files) {
	open I,'<',$Files{$fc} or die "[!]Error opening $Files{$fc}: $!\n";
	while (<I>) {
		my ($chr,$nbase1,$nbase2)=(split /\t/)[0,7,8];
		my $maf=$nbase2/($nbase1+$nbase2);
		++$DatAll{$maf};
		++$Dat{$chr}{$maf};
	}
	close I;
}

open O,'>',$ARGV[1] or die "[!]Error opening $ARGV[1]: $!\n";
print O "MAF\tAllCount\t",join("\t",@Chr),"\n";
for my $v (sort {$a <=> $b} keys %DatAll) {
	print O $v,"\t",$DatAll{$v};
	for my $c (@Chr) {
		if (exists $Dat{$c}{$v}) { print O "\t",$Dat{$c}{$v}; }
		else { print O "\t0"; }
	}
	print O "\n";
}
close O;
