#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use lib 'E:/BGI/pl/PERL5LIB/';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.3;

our $opts='i:o:s:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_d);

our $help=<<EOH;
\t-i Indel path (./Indel/) for [sample].indel.txt.filter
\t-s Samples list (./samples.list) [sample\\s+]
\t-o Output txt file (./indels.pop)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indels.pop' if ! defined $opt_o;
$opt_i='./Indel/' if ! $opt_i;
$opt_s='./samples.list' if ! $opt_s;

$opt_i =~ s/\/$//;
print STDERR "From [$opt_i]/ to [$opt_o], with [$opt_s]\n";
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

open( OUT,'>',$opt_o) or die "Error: $!\n";
for my $chr (sort keys %Indels) {
	for my $pos (sort {$a <=> $b} keys %{$Indels{$chr}}) {
		print OUT "$chr\t$pos\t";
		for my $Samp (@Samples) {
			my $item=$Indels{$chr}{$pos};
			if (defined $$item{$Samp}) {
			    if (@{$$item{$Samp}}==3) {	# if Indel, use depth of indel_file
				$item=join '|',@{$$item{$Samp}}
			    } else {
				$item=join '|',(@{$$item{$Samp}},'NA');
			    }
			} else {$item='.|.|NA'}
			print OUT $item,'^';
		}
		print OUT "\n";
	}
}
close OUT;

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
