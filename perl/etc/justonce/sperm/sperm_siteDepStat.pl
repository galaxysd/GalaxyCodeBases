#!/bin/env perl
use strict;
use warnings;
#use IO::Unread qw(unread);
use Data::Dump qw(ddx);

#my $ZONE_LENGTH = 200000;
#$ZONE_LENGTH = 20000;

die "Usage: $0 <ZONE_LENGTH> <output>\n" if @ARGV < 2;
my ($ZONE_LENGTH,$outf)=@ARGV;
warn "ZONE_LENGTH: $ZONE_LENGTH\n";

=pod
open CHR,'<','chr.lst' or die;
my %ChrGI2ID;
while (<CHR>) {
	my ($gi,$id) = split /\s+/,$_;
	$ChrGI2ID{$gi}=$id;
}
close CHR;
$ChrGI2ID{'='}='=';
#ddx \%ChrGI2ID;
=cut

open I,'<','xtubam/depth.sh' or die;
# $ cat /bak/seqdata/sperm/xtubam/depth.sh
# samtools depth mdaSperm23xtu.bam mdaSperm24xtu.bam mdaSperm28xtu.bam mlbacDonor.bam mlbacSpermS01.bam mlbacSpermS02.bam mlbacSpermS03.bam | xz -9c > depths.xz
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
		s/\.sort$//;
		push @IDs,$_;
	}
}
ddx \@IDs;
#our @IDsAllOne = split //,'1' x @IDs;

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

my $FileName = 'xtubam/depths.gz';
#$FileName = 'tmpdep.xz';
open I,'-|',"gzip -dc $FileName" or die;
my $posPoint;	# 1-based coordinate
my (%Dat,$t,%LastmPos);	# [cov-avg,emptRatio,covsum,covbp,emptbp,lasts]
while (<I>) {
	chomp;
	my ($chrid,$pos,@depthDat) = split /\t/;
	#my $chrid = $ChrGI2ID{$chrname};
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
	doinc($Dat{$chrid}->[$mPos]->[4], $pos - $posPoint - 1);
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
		unless (defined $t) {
			$Dat{$chrid}->[$ZoneID] = [[],[],[],[],[],[]];
			$t = $Dat{$chrid}->[$ZoneID];
			dosetvalue($t->[0],0);
			dosetvalue($t->[1],1);
			dosetvalue($t->[2],0);
			dosetvalue($t->[3],0);
			dosetvalue($t->[4],$ZONE_LENGTH);
			dosetvalue($t->[5],-1);	# Well, just mark it.
			next;
		}
		for my $i (0 .. $#IDs) {
			my $sum = $t->[3][$i] + $t->[4][$i];
			unless ($ZoneID == $LastmPos{$chrid}) {
				$sum +=  $t->[5][$i];
			}
			$t->[0][$i] = $t->[3][$i] ? ($t->[2][$i] / $t->[3][$i]) : 0;
			$t->[1][$i] = $t->[4][$i] / $sum;
		}
	}
}
#ddx \%Dat;

open OUT,'>',$outf or die "Error opening $outf: $!\n";
print OUT "# ZONE_LENGTH: $ZONE_LENGTH\n# Order: ",join(' ',@IDs),"\n#",
	join("\t",'ChrID','zPos','unCovBases','unCovRatio','CovSum','avgDepth','isNUL'),"\n";
for my $chrid (keys %Dat) {
	for my $ZoneID (0 .. $#{$Dat{$chrid}}) {
		$t = $Dat{$chrid}->[$ZoneID];
		print OUT join( "\t",$chrid,$ZoneID,
			join(' ',@{$t->[4]}),join(' ',map { int(.5+1000*$_)/1000 } @{$t->[1]}),
			join(' ',@{$t->[2]}),join(' ',map { int(.5+1000*$_)/1000 } @{$t->[0]}),
			join('',map {($_==-1)?1:0} @{$t->[5]}) ),"\n";
	}
}


__END__
perl sperm_siteDepStat.pl 200000 rss.tsv
perl rsstat.pl 500000 rss5k.tsv
perl rsstat.pl 1000000 rss1m.tsv
