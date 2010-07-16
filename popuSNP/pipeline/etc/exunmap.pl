#!/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(dump ddx);

unless (@ARGV){
	print "perl $0 <fq path/> <fq1,fq2> <soap path/> <soap,single,se>\n";	# soaps.nfo can store the file size of soap. Too late, useless.
	exit;
}

my ($fqpath,$fqs,$soapath,$soaps) = @ARGV;
my @fqs=split /,/,$fqs;
$_=$fqpath.$_ for @fqs;

my @soaps=split /,/,$soaps;
$_=$soapath.$_ for @soaps;
$soaps=$soaps[0];
$soaps=~s/\.(\w+)$/\.unmap/;

#ddx [\@fqs,\@soaps,$soaps];

sub popsoapid($$) {
	my ($fh,$this,$ziki)=@{$_[0]};
	my $theid=$ziki;
	defined(my $line1=<$fh>) or return [$fh,$theid,-1];
	my $id1=(split /\t/,$line1)[0];
	if ($id1 == $theid) {
		defined(my $line2=<$fh>) or return [$fh,-1,-1];;
		my $id2=(split /\t/,$line2)[0];
		die $line2 if $id2 == $theid;	# no time for writing a while
		return [$fh,$theid,$id2];
	} else { return [$fh,$theid,$id1]; }

}

sub popfq($$) {
	my ($FHar,$Pos)=@{$_[0]};
	my $Dat;
	++$Pos;
	for my $f (@$FHar) {
		my @line4;
		for (1..4) {
			defined(my $t=<$f>) or return [$FHar,-1,''];
			push @line4,$t;
		}
		chomp @line4;
		$Dat .= ">$line4[0] $Pos\n$line4[1]\n";	# for Fasta out
		#$Dat .= join("\n",@line4)."\n";	# for FastQ out
	}
	return [$FHar,$Pos,$Dat];
}

my ($fqsFHPosDat,@fqFH,@soapFHpos);
for my $f (@soaps) {
	my $t;
	open $t,'<',$f or die "[x]Error opening $f: $!\n";
	my $a=&popsoapid([$t,-2,-2]);
	$a=&popsoapid($a);
	push @soapFHpos,$a;
}
@soapFHpos = sort { $$a[1] <=> $$b[1] } @soapFHpos;

for my $f (@fqs) {
	my $t;
	open $t,'<',$f or die "[x]Error opening $f: $!\n";
	push @fqFH,$t;
}
$fqsFHPosDat=&popfq([\@fqFH,-1]);	# the 1st fq in soap is '0'

open O,'>',$soaps or die "[x]Error opening $soaps: $!\n";
while ($$fqsFHPosDat[1] != -1) {
	while ($$fqsFHPosDat[1] < $soapFHpos[0][1]) {
		print O $$fqsFHPosDat[2];
		$fqsFHPosDat=&popfq($fqsFHPosDat);
	}
	$fqsFHPosDat=&popfq($fqsFHPosDat) if $$fqsFHPosDat[1] == $soapFHpos[0][1];
	$soapFHpos[0]=&popsoapid($soapFHpos[0]);
	@soapFHpos = sort { $$a[1] <=> $$b[1] } @soapFHpos;
	$fqsFHPosDat=&popfq($fqsFHPosDat) if $$fqsFHPosDat[1] == $soapFHpos[0][1];
}
close O;
__END__
./exunmap.pl ./1fqfilted/WATxozRADDIAAPE=GS-2/100124_I644_FC615G8AAXX_L4_WATxozRADDIAAPE_ 1.fq,2.fq ./v4/2soap/GS-2/WATxozRADDIAAPE/100124_I644_FC615G8AAXX_L4_WATxozRADDIAAPE_ 1.soap,1.single
