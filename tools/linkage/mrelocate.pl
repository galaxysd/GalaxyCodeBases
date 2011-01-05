#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use Galaxy::ShowHelp;
use Time::HiRes qw ( gettimeofday tv_interval );
use Statistics::Regression;
use Data::Dump qw(ddx);

$main::VERSION=0.0.1;

our $opts='i:o:n:g:v:r:c:l:h:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_c, $opt_l, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_h);

our $help=<<EOH;
\t-i Blast filted list (./markerblastf.lst)
\t-l linkage map list (./linkagemap.lst)
\t-n N zone file (Pa64.Nzone)
\t-o Output file
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./markerblastf.lst' if ! $opt_i;
$opt_l='./linkagemap.lst' if ! $opt_l;
$opt_n='Pa64.Nzone' if ! $opt_n;
$opt_v=0 if ! $opt_v;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_n][$opt_l]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN

my (%LinkageMap,%cMcluster,%Marker);

open L,'<',$opt_l or die "Error: $!\n";
while (<L>) {
	chomp;
	my ($ChrID,$File)=split /\t/;
	open LM,'<',$File or die "Error: $!\n";
	while (<LM>) {
		next if /^#/;
		my ($Chr,$Pos,$cM)=split /\t/;
		die "[x][$ChrID] not match [$File] !\n" if $Chr ne $ChrID;
		$LinkageMap{$Chr}{$Pos}=$cM;
	}
}
close L;

sub splitMarkerid($) {
	my $MarkerID=$_[0];
	my ($mChr,$mPos,$mSiderLen)=split /[_m]/,$MarkerID,3;
	return [$mChr,$mPos,$mSiderLen];
}
sub getRel($$$$$) {
	my ($Qid,$Qs,$Qe,$Ss,$Se,$BTOP)=@_;
	my ($mChr,$mPos,$mSiderLen)=@{&splitMarkerid($Qid)};
	my $cM=$LinkageMap{$mChr}{$mPos};
	unless (defined $cM) {
		warn "[!]Cannot find Marker @ $mChr,$mPos !\n";
		$cM=0;
	}
	my $LeftBPalongQ=$mSiderLen-$Qs+1;
	my $WalkingOn=$LeftBPalongQ;
	my @btop0=split /(\D+)/,$BTOP;
	# $ perl -le '$a="45YT9-Ac-c-c-11TC4";@b=split /(\D+)/,$a;print "[",join("|",@b),"]"'
	# [45|YT|9|-Ac-c-c-|11|TC|4]
	my @btop=();
	for (@btop0) {
		if (/\d+/) {
			push @btop,$_;
		} else {
			my @bin=split /([\w-]{2})/;
			# $ perl -le '$a="-Ac-c-c-";@b=split /([\w-]{2})/,$a,0;print "[",join("|",@b),"]"'
			# [|-A||c-||c-||c-]
			for (@bin) {
				next unless $_;
				push @btop,$_;
			}
		}
	}
	my $LeftBPalongS=0;
	for (@btop) {
		if ($WalkingOn <= 0) {
			print STDERR "$Qid [$_] " if $opt_v and $WalkingOn == 0;
			last;
		}
		print STDERR "-$Qid [$_]-" if $opt_v>1;
		if (/\d/) {
			$LeftBPalongS += $_;
			$WalkingOn -= $_;
		} else {
			my @op=split //;
			if (/-/) {
				--$WalkingOn if $op[1] eq '-' and $op[0] ne '-';
				++$LeftBPalongS if $op[0] eq '-' and $op[1] ne '-';
			} else {
				--$WalkingOn;
				++$LeftBPalongS;
			}
		}
		print STDERR "-$WalkingOn $LeftBPalongS\n" if $opt_v>1;
		my $NOP;
	}
	my $strand;
	if ($Ss < $Se) {
		$strand=1;
		#$mChr .= "\t+";
	} else {
		$strand=-1;
		#$mChr .= "\t-";
	}
	my $Spos=$Ss + $strand*$LeftBPalongS;
	warn "$mChr,$Spos,$cM\n" if $opt_v;
	return [$mChr,$Spos,$cM,$strand];
}
open L,'<',$opt_i or die "Error: $!\n";
while (<L>) {
	chomp;
	my ($ChrID,$File)=split /\t/;
	open LM,'<',$File or die "Error: $!\n";
	my $lastcM=-1;
	while (<LM>) {
		chomp;
		my @Dat=split /\t/;
		my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP,$Hit)=@Dat;
		next if $Sid !~ /^chr/i;
		my ($mChr,$mPos,$mSiderLen)=@{&splitMarkerid($Qid)};
		next if $mChr ne $Sid;	# Well, since we caltuate cM in such a way ... So, no duplication on same cM
		next unless exists $LinkageMap{$mChr}{$mPos};
		my $cM=$LinkageMap{$mChr}{$mPos};
		if ($lastcM != $cM) {
			$lastcM=$cM;
			my $weight=-log($E)/$Hit;	# in case order problem, we will have to choose by weight
			my ($chr,$pos,$cM,$strand)=@{&getRel($Qid,$Qs,$Qe,$Ss,$Se,$BTOP)};
			#push @{$cMcluster{$Sid}{$cM}},[$pos,$strand,$weight];
			$cMcluster{$Sid}{$cM}=[$pos,$strand,$weight,$Qid];
			#warn "$Qid\t$Sid\t$pos,$strand,$weight\n" if scalar @{$cMcluster{$Sid}{$cM}} > 1;
		}
	}
}
close L;
ddx \%cMcluster;
ddx \%LinkageMap;
=pod
for my $Sid (sort keys %cMcluster) {
	my $reg = Statistics::Regression->new( "cM", [ "const", "BP" ] );
	for my $cM (keys %{$cMcluster{$Sid}}) {
		for (@{$cMcluster{$Sid}{$cM}}) {
			my ($pos,$strand,$weight)=@$_;
			$weight /= 100 if $strand == -1;	# Well, ...
			$reg->include( $cM, [ 1.0, $pos ], $weight );
		}
	}
	my ($b,$k) = $reg->theta;
	warn "$Sid,$k,$b\n";
	$reg->print();
}
=cut
for my $Sid (sort keys %cMcluster) {
	for my $cM (sort {$a<=>$b} keys %{$cMcluster{$Sid}}) {
	}
}












__END__
