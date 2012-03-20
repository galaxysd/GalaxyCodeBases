#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my ($in,$out)=("refGene.filter.gff", "ymy20120320.lst");

my (%GFF,%TransID,%Dat);

sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}
open IN,'<',$in or die "Error: $!\n";
while (<IN>) {
	next if /^\s+$/ or /^[;#]+/;
	chomp;
	s/\r//g;
	my ($seqname, $source, $primary, $start, $end,
	$score, $strand, $frame, $groups) = split /\t/;
	my @groups = split(/\s*;\s*/, $groups);
	my (%groups,$name);
	for my $group (@groups) {
		my ($tag,$value) = split /=/,$group;
		$tag             = unescape($tag);
		my @values       = map {unescape($_)} split /,/,$value;
		#$groups{$tag}=\@values;	# patch for those alter-splices
		$groups{$tag}=$values[0];
		die if @values>1;
	}
	if ($primary eq 'mRNA') {
		push @{$GFF{$groups{'name'}}},$groups{'ID'};
	}
}

for my $k (keys %GFF) {
	if (@{$GFF{$k}} != 2) {
		delete $GFF{$k};
	} else {
		for (@{$GFF{$k}}) {
			++$TransID{$_};
		}
	}
}
ddx \%GFF;
close IN;

open IN,'<',$in or die "Error: $!\n";	# It is a small file after all, ...
while (<IN>) {
	next if /^\s+$/ or /^[;#]+/;
	chomp;
	s/\r//g;
	my ($seqname, $source, $primary, $start, $end,
	$score, $strand, $frame, $groups) = split /\t/;
	my @groups = split(/\s*;\s*/, $groups);
	my (%groups,$name);
	for my $group (@groups) {
		my ($tag,$value) = split /=/,$group;
		$tag             = unescape($tag);
		my @values       = map {unescape($_)} split /,/,$value;
		#$groups{$tag}=\@values;	# patch for those alter-splices
		$groups{$tag}=$values[0];
	}
	if ($primary =~ /CDS|3-UTR/) {
		next unless exists $TransID{$groups{'Parent'}};
		push @{$Dat{$groups{'Parent'}}{$primary}},[$start, $end];
	} elsif ($primary eq 'mRNA') {
		next unless exists $TransID{$groups{'ID'}};
		$Dat{$groups{'ID'}}{'Strand'}=$strand;
	}
}
close IN;

open OUT,'>',$out or die "Error: $!\n";
for my $k (keys %GFF) {
	my ($IDA,$IDB)=@{$GFF{$k}};
	my $RefA=$Dat{$IDA};
	my $RefB=$Dat{$IDB};
	my $Strand = $$RefA{"Strand"};
	die if $Strand ne $$RefB{"Strand"};
	my ($ua,$ub,$uc,$ud);
	unless (exists $$RefA{"3-UTR"} and exists $$RefB{"3-UTR"}) {
		#ddx $RefA; ddx $RefB;
		warn "$k: [$IDA,$IDB]\n";
		if (exists $$RefA{"3-UTR"}) {
			($ua,$ub)=@{$$RefA{"3-UTR"}->[0]};
			($uc,$ud)=(0,0);
		} elsif (exists $$RefB{"3-UTR"}) {
			($ua,$ub)=(0,0);
			($uc,$ud)=@{$$RefB{"3-UTR"}->[0]};
		} else {
			warn "$k: [$IDA,$IDB]\n";
			next;
		}
	} else {
		($ua,$ub)=@{$$RefA{"3-UTR"}->[0]};
		($uc,$ud)=@{$$RefB{"3-UTR"}->[0]};
	}
	next if $ua==$uc and $ub==$ud;
	my ($a,$b,$c,$d);
	if ($Strand eq '+') {
		($a,$b)=@{$$RefA{"CDS"}->[-1]};	# sorted
		($c,$d)=@{$$RefB{"CDS"}->[-1]};
	} elsif ($Strand eq '-') {
		($a,$b)=@{$$RefA{"CDS"}->[0]};
		($c,$d)=@{$$RefB{"CDS"}->[0]};
#ddx $RefA; ddx $RefB;
#print "$a,$b $c,$d\n"; die;
	} else {die;}
	next unless $a==$c and $b==$d;
	print OUT "$k\t$IDA, [$ua,$ub]\t$IDB, [$uc,$ud]\t$Strand\t<$a,$b>\n";
}
close OUT;
