#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my ($flag,$minaln,$ia,$ib,$E,$line,$qSt,$qEd,$sSt,$sEd,@Q,@S)=(1,200);

sub setQS($$$$$$) {
	my ($ia,$ib,$qSt,$qEd,$sSt,$sEd)=@_;
	for my $i ($qSt .. $qEd) {
		#$Q[$i] = $ia if (!defined($Q[$i])) or $Q[$i] < $ia;
		++$Q[$i];
	}
	for my $i ($sSt .. $sEd) {
		#$S[$i] = $ib if (!defined($S[$i])) or $S[$i] < $ib;
		++$S[$i];
	}
}

open I,'<','tig2cat.blast' or die $!;
=pod
 Score = 2.188e+04 bits (11039), Expect = 0.0
 Identities = 12026/12360 (97%), Gaps = 72/12360 (0%)
 Strand = Plus / Plus
=cut
$line = readline(I);
while ($flag) {
	unless ($line =~/^ Score = /) {
		defined($line = readline(I)) or $flag=0;
		next;
	}
	$line =~ / Expect = ([\d\.+-e]+)$/ or die "[$line]";
	$E=$1;
	$E = "1$E" if $E=~/^e/i;
	last if $E>1e-10;
	$line = readline(I);
	$line =~ /^ Identities = (\d+)\/(\d+) / or die;
	($ia,$ib)=($1,$2);
	next if ($ia<$minaln or $ib<$minaln);
	$line = readline(I);
	unless ($line =~ /^ Strand = Plus \/ Plus/) {
		#warn "[$line";
		#print STDERR '.';
		defined($line = readline(I)) or $flag=0;
		next;
	}
	my (@QSt,@QEd,@SSt,@SEd);
	while ($line = readline(I)) {
		last if $line =~/^ Score = /;
		if ($line =~/^Query: (\d+) [^\d]+(\d+)$/) {
			push @QSt,$1;
			push @QEd,$2;
		}
		if ($line =~/^Sbjct: (\d+) [^\d]+(\d+)$/) {
			push @SSt,$1;
			push @SEd,$2;
		}
	}
	$qSt=$QSt[0]; $qEd=$QEd[-1];
	$sSt=$SSt[0]; $sEd=$SEd[-1];
	setQS($ia,$ib,$qSt,$qEd,$sSt,$sEd);
	#print STDERR '+';
#print "$ia,$ib,$E,$qSt,$qEd,$sSt,$sEd\n";
}
close I;
warn "\nDone.\n";

my ($Ql,$Sl,%Qh,%Sh)=(0,0);
for (@Q) {
	unless (defined $_) {
		++$Qh{0};
		next;
	}
	++$Ql if $_ == 1;
	++$Qh{$_};
}
for (@S) {
	unless (defined $_) {
		++$Sh{0};
		next;
	}
	++$Sl if $_ == 1;
	++$Sh{$_};
}
warn "Q:$Ql\nS:$Sl\n";
ddx \%Qh;
ddx \%Sh;

__END__
repeated==1
Q:10557598
S:10504512

repeated<=2
Q:10637014
S:10621036

no repeat filter (directly ++)
Q:12958103
S:16862471

Len:
cat75n1458	16866601
tiger75n1458	12958419

Query= Tiger75n1458

Q
#   "0"   => 1772885,
#   "1"   => 10557598,
#   "2"   => 79416,
#   "3"   => 42748,
#   "4"   => 24955,
#   "5"   => 11571,
S
#   "0"   => 5396572,
#   "1"   => 10504512,
#   "2"   => 116524,
#   "3"   => 54313,
#   "4"   => 36982,
#   "5"   => 34610,

