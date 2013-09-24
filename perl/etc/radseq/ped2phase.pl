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

my $Target = 'scaffold97,scaffold1457';
my %Targets;
my @tTarget = map { s/\s//g;$Targets{$_}=1;$_ } split /\,/,$Target;
#ddx \%Targets,\@tTarget;
#die;

my %mid2pos;
open ID,'<',$prefix.'.dict' or die "Error opening $prefix.dict : $!\n";
open IP,'<',$prefix.'.tped' or die "Error opening $prefix.tped : $!\n";
open IF,'<',$prefix.'.tfam' or die "Error opening $prefix.tfam : $!\n";

my %GTdata;
while (<ID>) {
	chomp;
	my ($chr,$pos,$id) = split /\t/;
	if (exists $Targets{$chr}) {
		$mid2pos{$id} = [$chr,$pos];
		#ddx $mid2pos{$id};
	}
}
close ID;

my ($NumberOfIndividuals,@Individuals);
while (<IF>) {
	my @tmp = split /\t/;
	push @Individuals,$tmp[1];
}
$NumberOfIndividuals = @Individuals;
close IF;

while (<IP>) {
	chomp;
	my ($chrNO,$id,undef,$pos,@dat) = split /\t/;
	if (exists $mid2pos{$id}) {
		#print join(',',@{$mid2pos{$id}},$chrNO,$id,$pos,@dat),"\n";
		my ($chrid,$chrpos) = @{$mid2pos{$id}};
		#$NumberOfIndividuals = @dat;
		for (@dat) {
			s/ //;
			s/0/\?/g;
			#$_ = '??' if $_ eq '00';
		}
		$GTdata{$chrid}{$chrpos} = \@dat;
	}
}
close IP;

#ddx \%GTdata;

for my $chrid (keys %GTdata) {
	my @Locus = sort { $a <=> $b } keys %{$GTdata{$chrid}};
	my $NumberOfLoci = @Locus;
	open O,'>',"$outfs.$chrid.inp" or die "Error opening $outfs.$chrid.inp : $!\n";
	print O "$NumberOfIndividuals\n$NumberOfLoci\nP ";
	print O join(' ',@Locus),"\n",'S' x $NumberOfLoci,"\n";
	for my $i (0 .. $#Individuals) {
		print O "#$Individuals[$i]\n";
		my ($m,$n,@t);
		for (@Locus) {
			@t = split //,$GTdata{$chrid}{$_}->[$i]; 
			$m .= $t[0];
			$n .= $t[1];
		}
		print O "$m\n$n\n";
	}

	close O;
}
__END__
perl ped2phase.pl sw000-18
./hp/phase.2.1.1.source/PHASE sw000-18.scaffold1457.inp p18s1457 >p18s1457.log 2>p18s1457.err &
./hp/phase.2.1.1.source/PHASE sw000-18.scaffold97.inp p18s97 >p18s97.log 2>p18s97.err &
