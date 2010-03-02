#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.3;

our $opts='i:o:s:bvl:f:c:';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_d, $opt_f, $opt_l, $opt_c);

our $help=<<EOH;
\t-i Indel path (./Indel/) for [sample].indel.txt.filter
\t-f fabyChr path (./fabychr/) for [chr].fa
\t-c consensus path (./consensus/) for [sample].[chr].txt
\t-s Samples list (./samples.list) [sample\\s+]
\t-l Length of nearby sequence (200)bp
\t-o Output txt file (./indels.pop)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indels.pop' if ! defined $opt_o;
$opt_i='./Indel/' if ! $opt_i;
$opt_s='./samples.list' if ! $opt_s;
$opt_f='./fabychr/' if ! $opt_f;
$opt_c='./consensus/' if ! $opt_c;
$opt_l=200 if ! $opt_l;

$opt_i =~ s/\/$//;
$opt_f =~ s/\/$//;
print STDERR "From [$opt_i]/ to [$opt_o], with [$opt_s][$opt_f][$opt_l]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

my (@Samples,%Indels);

open( SAMP,'<',$opt_s) or die "Error: $!\n";
while (<SAMP>) {
	chomp;
	my ($sample)=split /\s+/;
	push @Samples,$sample;
}
close SAMP;
if ($opt_v) {print "Samples List:\n";print "[$_]\t" for @Samples;print "\n";}

for my $Samp (@Samples) {
	my $name=$opt_i.'/'.$Samp.'.indel.txt.filter';
	open IN,'<',$name or die "Error: [$name] $!\n";
	warn "Indel:[$Samp]\t...\n";
	while (<IN>) {
		my ($chr,$pos,$nid,$seq,undef,$homhet,undef,$pairs)=split /\t/;
		if ($nid =~ /^I(\d+)$/) {$nid=$1;}
		 elsif ($nid =~ /^D(\d+)$/) {$nid=-$1;}
		 else  {warn "File error @ uid !"; next;}
		if ($homhet eq 'homo') {$seq = uc $seq;}
		 elsif ($homhet eq 'hete') {$seq = lc $seq;}
		 else  {warn "File error @ homhet !"; next;}
		print "$Samp\t$chr,$pos,$nid,$seq,$homhet\n" if $opt_v;
		$Indels{$chr}{$pos}{$Samp}=[$nid,$seq,$pairs];	# $pairs just for test
	}
	close IN;
}
warn "Indels all loaded.\n\nLoading Depth:\n";

my (%Depth,$i,$j);
for my $Samp (@Samples) {
	for my $chr (keys %Indels) {
print STDERR " ${Samp}_$chr\t";
		my $name=$opt_c.'/'.$Samp.'.'.$chr.'.txt';
		$i=$j=0;
		open IN,'<',$name or (warn "Error: [$name] $!\n" and next);
		while (<IN>) {
			++$i;
			my ($chr,$pos,undef,undef,undef,undef,$depth)=split /\t/;
			next unless exists $Indels{$chr}{$pos};
			++$j;
			$Depth{$chr}{$pos}{$Samp}=$depth;
		}
		close IN;
warn "$j/$i\n";
	}
}
warn "Depth all loaded.\n\n";

open( OUT,'>',$opt_o) or die "Error: $!\n";
for my $chr (sort keys %Indels) {
	open FA,'<',$opt_f.'/'.$chr.'.fa' or die "Error: $!\n";
	$_=<FA>;
	chomp;
	/^>(\S+)$/;
	die "[x]FASTA file wrong for $chr & $1" if $1 ne $chr;
	$/=">";
	my $seq = <FA>;
	chomp $seq;
	$seq=~s/\s//g;
	$seq = uc($seq); ## all upper case
	$/="\n";
	my $ref_len = length($seq);
	warn "$chr: $ref_len bp\n";
	close FA;
	for my $pos (sort {$a <=> $b} keys %{$Indels{$chr}}) {
		my ($posU,$lenU,$up);
		if ($pos>=$opt_l+1) {$posU=$pos-$opt_l-1;$lenU=$opt_l;}
		 else {$posU=0;$lenU=$pos-1;}
		if ($lenU>0) {$up=substr($seq,$posU,$lenU);}
		 else {$up='.';}
		my $down=substr($seq,$pos,$opt_l);	# First character is at offset 0
		print OUT "$chr\t$pos\t$up\t$down\t";
		warn "$chr,$pos\t$posU,$lenU\t$pos,$opt_l\n" if $opt_v;
		for my $Samp (@Samples) {
			my $item=$Indels{$chr}{$pos};
			if (defined $$item{$Samp}) {
			    if (@{$$item{$Samp}}==3) {	# if Indel, use depth of indel_file
				$item=join '|',@{$$item{$Samp}}
			    } else {
				$item=join '|',(@{$$item{$Samp}},$Depth{$chr}{$pos}{$Samp});
			    }
			} else {$item='.|.|'.$Depth{$chr}{$pos}{$Samp};}
			print OUT $item,'^';
		}
		print OUT "\n";
	}
}
close OUT;

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
