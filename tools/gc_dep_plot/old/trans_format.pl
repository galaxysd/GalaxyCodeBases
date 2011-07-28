#!/usr/bin/perl -w

use strict;

my %gc_dep;
my %dep_num;
while (<>){
	chomp;
	my @line = split /\t/;
	push @{$gc_dep{$line[0]}}, @line[1..$#line];
	my $num = $#line;
	$dep_num{$line[0]} = $num;
}

my $gc_num;
for my $gc (sort {$dep_num{$a} <=> $dep_num{$b}}  keys %dep_num) {
	$gc_num = "$gc\t$dep_num{$gc}";

}	
my ($gc, $num) = split /\t/, $gc_num;

for my $i (20..80) {
	print "$i\t";
}
print "\n";

for my $j (0..($num-1)) {
	for my $k (20..80) {
		if (defined ${$gc_dep{$k}}[$j]) {
			print "${$gc_dep{$k}}[$j]\t";
		} else {
			print "NA\t";
		}
	}
	print "\n";
}
