#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
use Data::Dump qw(ddx);
use Galaxy::IO::FASTAQ;
use Galaxy::SeqTools qw(translate revcom);
use Galaxy::Data;

die "Usage: $0 <vcf.gz> <minGQ> <out>\n" if @ARGV<3;
my $vcfs = shift;
my $qual = shift;
my $outfs = shift;

warn "From [$vcfs] with [minGQ=$qual] to [$outfs]\n";

open O,'>',$outfs."q$qual.snp.lst" or die;
open OI,'>',$outfs."q$qual.indel.lst" or die;

my $vcf = Vcf->new(file=>$vcfs);
$vcf->parse_header();
my (@samples) = $vcf->get_samples();
#ddx \@samples;
for (@samples) {
	s/\.bam$//i;
}
my $t = 'Samples:'.join(',',@samples);
print O join("\t",'#CHROM','POS','GT','QUAL',$t),"\n";
print OI join("\t",'#CHROM','POS','GT','QUAL',$t),"\n";

while (my $x=$vcf->next_data_hash()) {
	next if $$x{QUAL} < 20;
	my ($Dindel,@GTs)=(0);
	my %GTs = %{$$x{gtypes}};
	my @Bases = ($$x{REF},@{$$x{ALT}});
	$Dindel = length($Bases[1]) - length($Bases[0]);
	for my $sample (keys %GTs) {
		my $sampleID = $sample;
		$sampleID =~ s/\.bam$//i;
		next if $GTs{$sample}{DP}<=0;
		next if $GTs{$sample}{GQ}<$qual;	# 10:219090, 20:88956, 15:150149.
		my ($a1,$a2,$a3) = $vcf->split_gt($GTs{$sample}{GT});
		die "[$a3]" if $a3;
		my ($t,%t)=(0);
		for (($Bases[$a1],$Bases[$a2])) {
			++$t{$_};
			$t = length $_ if $t < length $_;
		}
		if ($t==1) {
			$t = join('',sort keys %t);
			$t = $REV_IUB{$t};
		} else {
			$t = '.';
		}
		push @GTs, join('|',"$a1/$a2","$Bases[$a1]/$Bases[$a2]",$GTs{$sample}{DP},$GTs{$sample}{GQ},$t);
	}
	next unless @GTs;
	if (length($Bases[0]) == 1) {
		print O join("\t",$$x{CHROM},$$x{POS},join('/',@Bases,$Dindel),$$x{QUAL}, join(", ",@GTs) ),"\n";
	} else {
		print OI join("\t",$$x{CHROM},$$x{POS},join('/',@Bases,$Dindel),$$x{QUAL}, join(", ",@GTs) ),"\n";
	}

#	for my $gt (keys %GTs) {
#		my ($a1,$a2,$a3) = $vcf->split_gt($gt);
#		if ($a3 or ($a1 != $a2)) {
#			$flag = 1;
#		}
#	}
	#ddx \%GTs;
	#ddx $x;
# vcf2cds.pl:189: {
#   ALT    => ["A"],
#   CHROM  => "scaffold75",
#   FILTER => ["."],
#   FORMAT => ["GT", "PL", "DP", "SP", "GQ"],
#   gtypes => {
#               "BHX011.bam" => { DP => 19, GQ => 61, GT => "0/0", PL => "0,57,255", SP => 0 },
#               "BHX019.bam" => { DP => 26, GQ => 82, GT => "0/0", PL => "0,78,255", SP => 0 },
#               "JHH001.bam" => { DP => 28, GQ => 99, GT => "0/1", PL => "244,0,255", SP => 9 },
#             },
#   ID     => ".",
#   INFO   => {
#               AC1 => 1,
#               AF1 => 0.1667,
#               DP  => 74,
#               DP4 => "34,26,4,9",
#               FQ  => 999,
#               MQ  => 59,
#               PV4 => "0.13,1,0.12,1",
#               VDB => 0.0365,
#             },
#   POS    => 3572768,
#   QUAL   => 999,
#   REF    => "G",
# }
	#print $vcf->format_line($x);
}
#ddx $vcf;
$vcf->close;

close O;
close OI;

sub cmpstr {
	my ($a, $b) = @_;
	my $c = $a ^ $b;
	my @ret;
	while ($c =~ /[^\0]/g) {
		my $p = pos($c);
		push @ret,[$p,substr($a,$p-1,1),substr($b,$p-1,1)];
	}
	@ret;
}

__END__

bcftools view /bak/archive/projects/Tiger/BCF/WGS/parents.bcgv.bcf |bgzip -c > parents.bcgv.vcf.gz &
#tabix -p vcf parents.bcgv.vcf.gz

perl filtervcf.pl bam2.bcgv.vcf.gz 20 bam2 &

