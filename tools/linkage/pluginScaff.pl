#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use Galaxy::ShowHelp;
use Time::HiRes qw ( gettimeofday tv_interval );
use Data::Dump qw(ddx);

$main::VERSION=0.0.1;

our $opts='i:o:n:g:v:r:c:l:p:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_c, $opt_l, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_h, $opt_p);

our $help=<<EOH;
\t-i Genome FASTA (./Pa64_genome.fa)
\t-p Anchor Pos file (./f3545ChrScaff.pos.n3)
\t-n N zone file (Pa64.Nzone)
\t-o Output prefix
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./Pa64_genome.fa' if ! $opt_i;
$opt_p='./f3545ChrScaff.pos.n3' if ! $opt_p;
$opt_n='Pa64.Nzone' if ! $opt_n;
$opt_v=0 if ! $opt_v;

print STDERR "From [$opt_i] to [$opt_o].{fa,stat} refer to [$opt_n][$opt_p]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%Nzone,%SacffPos,%Genome);
open I,'<',$opt_n or die $!;
while (<I>) {
	#chomp;
	my ($chr,$s,$e,$len)=split /\t/;
	push @{$Nzone{$chr}},[$s,$e];
}
close I;
open I,'<',$opt_p or die $!;
while (<I>) {
	next if /^#/;
	#chomp;
	my ($ScaffID,$sPosRel,$chr,$strand,$scM,$pos,$err)=split /\t/;
	push @{$SacffPos{$chr}},[$ScaffID,$strand,$sPosRel,$pos,$err];
}
close I;
open I,'<',$opt_i or die $!;
{
	local $/=">";
	$_=<I>;
	die "[x]Not a FASTA file ! [$opt_i]\n" unless /^\s*>/;
}
while (<I>) {
	chomp;
	my ($id,$desc)=split / /,$_,2;
	$desc='' unless $desc;
	$/=">";
	my $seq=<I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	print STDERR ">$id,\t[$desc]: ",length $seq,"\n";
	$Genome{$id}=$seq;
}
close I;

sub searchtoPlug($$$$) {
	my ($chr,$pos,$left,$right)=@_;
	my $theZone;
	for (@{$Nzone{$chr}}) {
		my ($s,$e)=@$_;
		next if $e < $pos;
		$theZone=$_;
		last if $s > $pos;
	}
	unless ($theZone) {
		$theZone=[1,1];
		for (@{$Nzone{$chr}}) {
			my ($s,$e)=@$_;
			next if $e < $left;
			$theZone=$_;
			last if $s > $right;
		}
	}
	my ($thePos,$cutL,$cutR);
	if ($pos>$theZone->[0] and $pos<$theZone->[1]) {	# in zone
		$thePos=$pos;
		($cutL,$cutR)=($thePos-$theZone->[0]-1,$theZone->[1]-$thePos-1);
	} else {
		$thePos=$theZone->[0];
		$cutL=0;
		$cutR=$theZone->[1]-$theZone->[0]-1;
	}
	return [$thePos,$cutL,$cutR];
}
open S,'>',$opt_o.'.stat' or die $!;
open O,'>',$opt_o.'.fa' or die $!;
for my $chr (sort keys %SacffPos) {
	my $lastPos=1;
	print O ">$chr\n";
	print S "[$chr]\n";
	warn ">$chr\n";
	for my $ref (sort {$a->[3] <=> $b->[3]} @{$SacffPos{$chr}}) {
		my ($ScaffID,$strand,$sPosRel,$pos,$err)=@$ref;
		my $ScaffLen=length $Genome{$ScaffID};
		my ($ScaffLeft,$ScaffRight);
		if ($strand eq '+') {
			$ScaffLeft=1-$sPosRel;
			$ScaffRight=$ScaffLen-$sPosRel;
		} else {
			$ScaffLeft=$sPosRel-$ScaffLen;
			$ScaffRight=$sPosRel-1;
		}
		my ($left,$right)=($pos+$ScaffLeft-$err,$pos+$ScaffRight+$err);	# this should be wrong but should be right.
		my ($thePos,$cutL,$cutR)=@{&searchtoPlug($chr,$pos,$left,$right)};
		print S "$ScaffID\t",length($Genome{$ScaffID}),"\t$strand\t$thePos\t$cutL,$cutR\n";
		print O substr($Genome{$chr},$lastPos-1,$thePos-$lastPos+1-$cutL),"\nx",$Genome{$ScaffID},"x\n";
		$lastPos=$thePos+$cutR;
	}
	print O "\n";
	print S "\n";
}
close O;
close S;

__END__
./pluginScaff.pl -o PA64new.fa
