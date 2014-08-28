#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <out> <in1> [in2] [...]\n" if @ARGV < 2;
my ($outf,@infs)=@ARGV;

my (@FH,%Dat,%Max,@Names);
for my $i (@infs) {
	my $t;
	open $t,'<',$i or die "Error opening $i: $!\n";
	push @FH,$t;
}
ddx \@FH;

for my $i (@FH) {
	my $name = readline $i;
	$name =~ s/^# //;
	chomp $name;
	push @Names,$name;
	my $chr = '][';
	my $zeroPos=0;
	while (<$i>) {
		if (/^\[(\w+)\]$/) {
			$chr = $1;
			$zeroPos=0;
			next;
		}
		my $v = (split /\s+/)[0];
		$Dat{$chr}[$zeroPos]{$name} = $v;
		$Max{$chr} = $v if (defined $Max{$chr}) ? $Max{$chr} < $v : 1;
		++$zeroPos;
	}
}
ddx \%Dat;
ddx \%Max;

close $_ for @FH;

open OUT,'>',$outf or die "Error opening $outf: $!\n";
print OUT join("\t",'Chr','Pos',@Names),"\n";
for my $chr (sort { "$a$b"=~/^\d+$/ ? $a<=>$b : $a cmp $b } keys %Dat) {
	my $ArrayRef = $Dat{$chr};
	for my $zeroPos ( 0 .. $#$ArrayRef ) {
		print OUT join("\t",$chr,1+$zeroPos);
		print OUT "\t",$ArrayRef->[$zeroPos]->{$_} for @Names;
		print OUT "\n";
	}
} 
close OUT;

__END__
perl rsmerge.pl t t.out t2.out /bak/seqdata/sperm/t.out
perl rsmerge.pl bamrsplot.tsv xtubam/*.ss
