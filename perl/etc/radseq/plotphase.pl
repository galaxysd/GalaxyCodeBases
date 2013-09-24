#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
Purpose: Read bcf, get tped for p-link
Notes: rad2marker is deprecated.
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::SeqTools;
use Data::Dumper;

die "Usage: $0 <tped prefix> [out].chrID.inp\n" if @ARGV<1;
my $prefix=shift;
my $outfs=shift;
unless (defined $outfs) {
	warn "Using prefix[$prefix] for both input and output.\n";
	$outfs = $prefix;
}
warn "From [$prefix] to [$outfs]\n";

my (@Positions,%HapDat,@IDs,@IDab);
open I,'<',$prefix or die "Error opening $prefix : $!\n";
while (<I>) {
	if (/^Positions of loci: /) {
		@Positions = split /\s+/,$_;
		splice @Positions,0,3;
	}
	if (/^BEGIN BESTPAIRS1/) {
		while (<I>) {
			last if /^END BESTPAIRS1/;
			chomp;
			my $id = (split /\s/,$_)[1];
			$id =~ s/^#//;
			chomp($_ = <I>);
			my @basesA = split /\s+/,$_;
			chomp($_ = <I>);
			my @basesB = split /\s+/,$_;
			for my $base (@basesA,@basesB) {
				if ($base =~ /^\((\w)\)$/) {
					$base = lc $1;
					#last;
				}
				$base = '.'.lc($1) if ($base =~ /^\[(\w)\]$/);
			}
			push @IDs,$id;
			push @IDab,"${id}_a";
			push @IDab,"${id}_b";
			$HapDat{$id} = [\@basesA,\@basesB];
		}
	}
}
#ddx \%HapDat;
close I;

open O,'>',$prefix.'.tsv' or die "Error opening $prefix.tsv : $!\n";
print O "# From: [$prefix]\n",join("\t",'Pos',@IDab),"\n";
for my $i (0 .. $#Positions) {
	my @tmp;
	push @tmp,$Positions[$i];
	for my $id (@IDs) {
		push @tmp,$HapDat{$id}->[0]->[$i];
		push @tmp,$HapDat{$id}->[1]->[$i];
	}
	print O join("\t",@tmp),"\n";
}
close O;

__END__
