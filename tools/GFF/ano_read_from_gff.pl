#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
use Galaxy::IO;
use Data::Dump qw(ddx);

die "Usage: $0 <gff.gz> <gene_file>\n" if @ARGV<2;
my $infs = shift;
my $outfs = shift;

warn "From [$infs] to [$outfs]\n";

my $IN = openfile($infs);
open O,'>',$outfs or die;
print O "# From [$infs] to [$outfs]\n";

sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}

my ($LN,%GFFdat,$GFFanno)=(0);
while (<$IN>) {
	++$LN;
	next if /^\s+$/ or /^[;#]+/;
	chomp;
	s/\r//g;	# There DO be some file with \r.
	my ($seqname, $source, $primary, $start, $end,
	$score, $strand, $frame, $groups) = split /\t/;	# well, reading file no need to Optimize
	my @groups = split(/\s*;\s*/, $groups);
	my (%groups,$name);
	for my $group (@groups) {
		my ($tag,$value) = split /=/,$group;
		$tag             = unescape($tag);
		my @values       = map {unescape($_)} split /,/,$value;
		$groups{$tag}=\@values;
	}
	my @name_order=qw/Parent ID/;
	@name_order=qw/ID Parent/ if $primary =~ /mRNA/i;
	for (@name_order) {
		if ($groups{$_}) {$name=$groups{$_};last;}
	}
	for (@$name) {
		#@dat=($seqname,$primary,$start,$end,$strand,$frame,$groups,$_);
		push @{$GFFdat{$_}{$strand}},[$seqname, $primary, $start, $end, $frame, $LN, $groups];
		#print "$seqname,$primary,$start,$end,$strand,$frame,$_\n" if $opt_v;
	}
}
close $IN;

for my $geneid (sort keys %GFFdat) {
	my @Strands = sort keys %{$GFFdat{$geneid}};
	die if @Strands != 1;
	for my $strand ( @Strands ) {
		my $AnoArray = $GFFdat{$geneid}{$strand};
		print O "\n[$geneid] $strand\n";
		#ddx $AnoArray;
		for (sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @$AnoArray) {
			print O join("\t",@$_),"\n";
		}
	}
}
close O;

__END__
perl ano_read_from_gff.pl ref_Dasnov3.0_gnomon_scaffolds.gff3.gz Dasnov3.gene1 &

