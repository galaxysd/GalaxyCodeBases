#!/bin/env perl
#get the DP score matrix for DNA alignment, based on the average alignment identity
use strict;


print "align_identity\tmatch_score\tunmatch_score\n";
for (my $identity = 0.50; $identity<1; $identity+=0.01) {
	my $Pmatch = $identity / 4;
	my $Punmat = (1-$identity) / 12;

	my $Prandom = 1/(4*4);
	my $lamda = 4;

	my $match_score = int($lamda*log($Pmatch/$Prandom)+0.5);
	my $unmat_score = int($lamda*log($Punmat/$Prandom)+0.5);

	print  "$identity\t$match_score\t$unmat_score\n";

}


__END__

