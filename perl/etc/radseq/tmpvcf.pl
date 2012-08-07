#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
#use Galaxy::IO::VCF;
use Vcf;

die "Usage: $0 <bcf/vcf> <out>\n" if @ARGV<2;
my $bcfs=shift;
my $outfs=shift;

my ($fh,$vcf);
if ($bcfs =~ /\.bcf$/) {
	$fh = openpipe('bcftools view',$bcfs);
	$vcf = Vcf->new(fh=>$bcfs);
} else {
	$fh = openpipe('bcftools view -S',$bcfs);
	$vcf = Vcf->new(file=>$bcfs);
}
$vcf->parse_header();
my (@samples) = $vcf->get_samples();
ddx \@samples;

while (my $x=$vcf->next_data_hash()) { 
#	for my $gt (keys %{$$x{gtypes}}) {
#		my ($al1,$sep,$al2) = $vcf->parse_alleles($x,$gt);
#		print "\t$gt: $al1$sep$al2\n";
#	}
#	print "\n";
	next unless $$x{REF} =~ /^[ATCG]$/;
	next if $$x{QUAL} < 20;
	my ($flag,%GTs)=(0);
	for my $gt (keys %{$$x{gtypes}}) {
		if ($$x{gtypes}{$gt}{DP} > 0) {
			++$GTs{$$x{gtypes}{$gt}{GT}};
		}
	}
	next if (keys %GTs) > 1;	# one GT;
	for my $gt (keys %GTs) {
		$flag = 1 if $gt =~ /0/;
		my ($a1,$a2,$a3) = $vcf->split_gt($gt);
		if ($a3 or ($a1 != $a2)) {
			$flag = 1;
		}
	}
	#ddx [\%GTs,$x] if $flag;
	next if $flag;
	print $vcf->format_line($x);
}
=pod
        BHX011.bam: A/A
        BHX019.bam: A/A
        JHH001.bam: A/A

# tmpvcf.pl:35: {
#   ALT    => ["A"],
#   CHROM  => "scaffold1001",
#   FILTER => ["."],
#   FORMAT => ["GT", "PL", "DP", "SP", "GQ"],
#   gtypes => {
#               "BHX011.bam" => { DP => 8, GQ => 32, GT => "1/1", PL => "131,24,0", SP => 0 },
#               "BHX019.bam" => { DP => 6, GQ => 26, GT => "1/1", PL => "107,18,0", SP => 0 },
#               "JHH001.bam" => { DP => 0, GQ => 9, GT => "1/1", PL => "0,0,0", SP => 0 },
#             },
#   ID     => ".",
#   INFO   => { AC1 => 6, AF1 => 1, DP => 14, DP4 => "0,0,9,5", FQ => -33.3, MQ => 35, VDB => 0.0256 },
#   POS    => 5193,
#   QUAL   => 202,
#   REF    => "C",
# }
=cut

#while (my $x=$vcf->next_data_array()) {
#	print join('|',@$x),"\n"; #"$$x[0]:$$x[1]\n";
#	ddx $x;
#}
# tmpvcf.pl:62: [
#   "scaffold1001",
#   5193,
#   ".",
#   "C",
#   "A",
#   202,
#   ".",
#   "DP=14;VDB=0.0256;AF1=1;AC1=6;DP4=0,0,9,5;MQ=35;FQ=-33.3",
#   "GT:PL:DP:SP:GQ",
#   "1/1:0,0,0:0:0:9",
#   "1/1:131,24,0:8:0:32",
#   "1/1:107,18,0:6:0:26",
# ]

$vcf->close();

__END__
./tmpvcf.pl JB19.bcvgL.vcf.gz c > tmp
perl -lane 'BEGIN {my %a;} next if /^#/; ++$a{$F[0]}; END {print "$_\t$a{$_}" for sort {$a{$b} <=> $a{$a}} keys %a;}' tmp|cat -n > tmp.s
