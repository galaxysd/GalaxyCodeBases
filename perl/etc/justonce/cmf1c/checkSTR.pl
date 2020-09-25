#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec::Functions;

use Data::Dump qw(ddx);

my $theDir='dat';
my $outF='Checked.tsv';
opendir (DIR, $theDir) or die "[x]Cannot open directory [$theDir], $!";
open O,'>',$outF or die "[x]Cannot open file [$outF], $!";
while (my $file = readdir DIR) {
	next unless $file =~ /\.dat$/i;
	my $fname = catfile($theDir,$file);
	print STDERR "[$file]";
	my ($Fcur,$Fcnt,$ret,@Fdat)=(0,0);
	open FIN,'<',$fname or die "[x]Cannot open file [$fname], $!";
	@Fdat = <FIN>;
	$Fcnt = scalar @Fdat;
	close FIN;
	for (@Fdat) {
		s/\r[\n]*//gm;
	}
	my $fDate = $Fdat[5];
	print STDERR "@[$fDate]: ";
	my $app = $Fdat[7];
	if ($app ne 'GeneMapper ID-X') {
		die "\n[x]Unsupported AppID:[$app].\n"
	}
	my $Scnt = $Fdat[8];
	if ($Scnt < 1) {
		print STDERR "Empty.\n";
		next;
	} else {
		print STDERR "$Scnt";
	}
	$Fcur=9;
	$ret = readCMF1($file,\@Fdat,$Fcur);
	print STDERR "\n";
}
close O;
closedir DIR;

sub readCMF1 {
	my ($file,$pfdat,$Fcur) = @_;
	die if $$pfdat[$Fcur] ne 'DNA Analysis Result';
	++$Fcur;
	my @aSample;
	for (my $i=$Fcur;$i<=$#$pfdat;$i++) {
		if ($$pfdat[$i] eq 'DNA Analysis Result') {
			parseAsample(\@aSample);
			@aSample=();
			next;
		}
		push @aSample,$$pfdat[$i];
	}
}

sub parseAsample {
	my ($pSdat) = @_;
	#ddx \$pSdat;
	die if $$pSdat[0] ne '1.0';
	die if $$pSdat[1] ne 'PCR';
	my $Sid = $$pSdat[2];
	my $STRcnt = $$pSdat[8];
	my ($STRcur,$STRc)=(9,0);
	#print "> $STRcur\n";
	while ($STRcur <= $#$pSdat) {
		#print ">>>$Sid $STRcnt $STRcur $$pSdat[$STRcur]\n";
		my $strID = $$pSdat[$STRcur];
		die if $$pSdat[$STRcur+1] != 1;
		my $strAcnt = $$pSdat[$STRcur+5];
		die $strAcnt if $strAcnt < 0;
		my @strA=();
		if ($strAcnt > 0) {
			for my $x (($STRcur+6) .. ($STRcur+5+$strAcnt)) {
				push @strA,$$pSdat[$x];
			}
		}
		++$STRc;
		#ddx [$STRc,$strID,@strA];
		$STRcur += $strAcnt+6;
	}
}
