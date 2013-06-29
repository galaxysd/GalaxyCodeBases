#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

my $ZONE_LENGTH = 200000;
$ZONE_LENGTH = 20000;

die "Usage: $0 <output>\n" if @ARGV < 1;
my ($outf)=@ARGV;

open CHR,'<','chr.lst' or die;
my %ChrGI2ID;
while (<CHR>) {
	my ($gi,$id) = split /\s+/,$_;
	$ChrGI2ID{$gi}=$id;
}
close CHR;
$ChrGI2ID{'='}='=';
ddx \%ChrGI2ID;

open I,'<','xtubam/depth.sh' or die;
our @IDs;
my @tmpids;
while (<I>) {
	if (/^samtools depth /) {
		(undef,undef,@tmpids) = split /\s+/;
		last;
	}
}
close I;
for (@tmpids) {
	if (/\.bam$/) {
		s/\.bam$//;
		push @IDs,$_;
	}
}
ddx \@IDs;
our @IDsAllOne = split //,'1' x @IDs;

sub doaddup($$) {
	my ($adder,$data) = @_;
	for my $i (0 .. $#IDs) {
		$$adder[$i] += $$data[$i];
	}
}
sub doaddifone($$) {
	my ($adder,$data) = @_;
	for my $i (0 .. $#IDs) {
		++$$adder[$i] if $$data[$i] > 0;
	}
}
sub doaddifzero($$) {
	my ($adder,$data) = @_;
	for my $i (0 .. $#IDs) {
		++$$adder[$i] if $$data[$i] == 0;
	}
}
sub doinc($$) {
	my ($adder,$value) = @_;
	for my $i (0 .. $#IDs) {
		$$adder[$i] += $value;
	}
}
sub dosetvalue($$) {
	my ($adder,$value) = @_;
	for my $i (0 .. $#IDs) {
		$$adder[$i] = $value;
	}
}

my $FileName = 'xtubam/depths.xz';
$FileName = 'tmpdep.xz';
open I,'-|',"xz -dc $FileName" or die;
my $posPoint;	# 1-based coordinate
#$posPoint{$_} = 0 for values %ChrGI2ID; 
my (%Dat,$t,%LastmPos);	# [Covered,Empty] x [result, sum, count] => [cov-avg,emptRatio,covsum,covbp,emptbp,lasts]
while (<I>) {
	chomp;
	my ($chrname,$pos,@depthDat) = split /\t/;
	my $chrid = $ChrGI2ID{$chrname};
	my $mPos = int($pos / $ZONE_LENGTH);
	unless (defined $Dat{$chrid}->[$mPos]) {
		$Dat{$chrid}->[$mPos] = [[],[],[],[],[],[]];
		$t = $Dat{$chrid}->[$mPos];
		dosetvalue($t->[0],0);
		dosetvalue($t->[1],0);
		dosetvalue($t->[2],0);
		dosetvalue($t->[3],0);
		dosetvalue($t->[4],0);
		dosetvalue($t->[5],0);
		$posPoint = $mPos*$ZONE_LENGTH -1;
	}
	#if ($posPoint{$chrid} < $pos) {
		doinc($Dat{$chrid}->[$mPos]->[4], $pos - $posPoint - 1);
		#warn "$mPos $pos, $posPoint ",$pos - $posPoint - 1,"\tDep: @depthDat\n";
	#}
	$posPoint = $pos;
	doaddup($Dat{$chrid}->[$mPos]->[2], \@depthDat);
	doaddifone($Dat{$chrid}->[$mPos]->[3], \@depthDat);
	doaddifzero($Dat{$chrid}->[$mPos]->[4], \@depthDat);
	dosetvalue($t->[5],(1+$mPos)*$ZONE_LENGTH -1 - $pos);
	$LastmPos{$chrid} = $mPos;
}
close I;

for my $chrid (keys %Dat) {
	for my $ZoneID (0 .. $#{$Dat{$chrid}}) {
		$t = $Dat{$chrid}->[$ZoneID];
		for my $i (0 .. $#IDs) {
			$t->[0][$i] = $t->[3][$i] + $t->[4][$i];
			unless ($ZoneID == $LastmPos{$chrid}) {
				$t->[0][$i] +=  $t->[5][$i];
			}
		}
	}
}

ddx \%Dat;

__END__
perl rsstat.pl
