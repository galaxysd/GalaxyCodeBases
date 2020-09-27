#!/usr/bin/env perl

use strict;
use warnings;
use File::Spec::Functions;

#use Data::Dump qw(ddx);

my %toCheck = (
	'DYS643' => 16,
	'DYS557' => 23,
);
my $theDir='dat';
my $outF='Checked.tsv';
if (-d $theDir) {
	my $Ocnt=0;
	# $pathname exists and is a directory
	opendir (DIR, $theDir) or die "[x]Cannot open directory [$theDir], $!";
	open O,'>',$outF or die "[x]Cannot open file [$outF], $!";
	while (my $file = readdir DIR) {
		next unless $file =~ /\.dat$/i;
		my $fname = catfile($theDir,$file);
		print STDERR "[$file]: ";
		my ($Fcur,$Fcnt,$ret,@Fdat)=(0,0);
		open FIN,'<',$fname or die "[x]Cannot open file [$fname], $!";
		@Fdat = <FIN>;
		$Fcnt = scalar @Fdat;
		close FIN;
		for (@Fdat) {
			s/\r*[\n]*//gm;
		}
		my $fDate = $Fdat[5];
		#print STDERR "@[$fDate]: ";
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
		#ddx $ret;
		if (keys(%{$ret})>0) {
			for my $k (keys %{$ret}) {
				my $str;
				for my $kk (keys %{$ret->{$k}}) {
					$str .= ' '.$kk.':'.$ret->{$k}->{$kk};
				}
				print O join("\t",$file,$k,$str),"\n";
				++$Ocnt;
			}
		}
	}
	close O;
	closedir DIR;
	if ($Ocnt > 0) {
		die "[!]See [$outF] for [$Ocnt] results.\n";
	} else {
		unlink $outF;
		die "[!]No results found.\n";
	}
} elsif (-e $theDir) {
	# $pathname exists, but is NOT a directory
	unlink $theDir;
	mkdir $theDir or die '[x]Please create folder [$theDir] beside me !\n';
} else {
	mkdir $theDir or die '[x]Please create folder [$theDir] beside me !\n';
}
die "[!]Put .dat file in [$theDir] and run me again.\n";

sub readCMF1 {
	my ($file,$pfdat,$Fcur) = @_;
	die if $$pfdat[$Fcur] ne 'DNA Analysis Result';
	++$Fcur;
	my @aSample;
	my %retS;
	for (my $i=$Fcur;$i<=$#$pfdat;$i++) {
		if ($$pfdat[$i] eq 'DNA Analysis Result') {
			my ($sid,$retp) = parseAsample(\@aSample);
			if (keys(%{$retp})>0) {
				$retS{$sid}=$retp;
			}
			@aSample=();
			next;
		}
		push @aSample,$$pfdat[$i];
	}
	return \%retS;
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
	my %ret;
	while ($STRcur <= $#$pSdat) {
		#print ">>>$Sid $STRcnt $STRcur $$pSdat[$STRcur]\n";
		my $strID = $$pSdat[$STRcur];
		die if $$pSdat[$STRcur+1] != 1;
		my $strAcnt = $$pSdat[$STRcur+5];
		die $strAcnt if $strAcnt < 0;
		my @strA=();
		my $flag=0;
		my @tmp=();
		if ($strAcnt > 0) {
			for my $x (($STRcur+6) .. ($STRcur+5+$strAcnt)) {
				push @strA,$$pSdat[$x];
				if (exists $toCheck{$strID}) {
					if ($$pSdat[$x] > $toCheck{$strID}) {
						$flag=1;
						#print STDERR " ($strID:$$pSdat[$x])";
						print STDERR '.';
						push @tmp,$$pSdat[$x];
					}
				}
			}
		}
		if ($flag) {
			$ret{$strID} = join(',',@tmp);
		}
		++$STRc;
		#ddx [$STRc,$strID,@strA];
		$STRcur += $strAcnt+6;
	}
	return ($Sid,\%ret);
}
