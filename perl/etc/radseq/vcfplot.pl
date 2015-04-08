#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
#use Galaxy::SeqTools;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

die "Usage: $0 <tfam file> <bcgv bcf> <regions> <list,of,samples to exclude(none)> <Dominant/Recessive> <outprefix>\n" if @ARGV<6;
my $tfamfs=shift;
my $bcfs=shift;
my $regions=shift;
my $samples=shift;
my $DomRec=shift;
my $outfs=shift;

$DomRec='D' if $DomRec =~ /^D/i;
$DomRec='R' if $DomRec =~ /^R/i;

my %SkippedSamples;
my $cmd = 'bcftools query -f\'%CHROM\t%POS\t%REF,%ALT\t%QUAL\t%DP[\t%SAMPLE,%GT,%DP,%GQ]\n\'';
if ($regions ne 'all') {
	$cmd .= " -r '$regions'";
	my @t = split /,/,$samples;
	%SkippedSamples = map { $_ => 1 } @t;
}
if ($samples ne 'none') {
	$cmd .= ' -s^'.$samples;
	my @t = split /,/,$samples;
	%SkippedSamples = map { $_ => 1 } @t;
}
warn "$DomRec [$cmd $bcfs]\n";

my (@tfamSamples,%tfamDat,%tfamSamplePheno,%inFamily,@CaseS,@ControlS);
open L,'<',$tfamfs or die;
while (<L>) {
	next if /^(#|$)/;
	#print OF $_;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	next if exists $SkippedSamples{$ind};
	$tfamDat{$ind} = $_."\n";
	if ($pho == 1 or $pho == 2) {
		$tfamSamplePheno{$ind} = $pho;	# disease phenotype (1=unaff/ctl, 2=aff/case, 0=miss)
		if ($pho == 2) {
			push @CaseS,$ind;
		} else {
			push @ControlS,$ind;
		}
	} else { die; }	# $pho can only be 1 or 2
	push @tfamSamples,$ind;
	push @{$inFamily{$family}},$ind;
}
#ddx \%tfamSamplePheno,\%inFamily;
close L;

my @Samples;
my $fh = openpipe('bcftools query -l',$bcfs);
while (<$fh>) {
	chomp;
	push @Samples,$_ unless exists $SkippedSamples{$_};
}
close $fh;

sub biggerKeyExcept($$) {
	my ($in,$except) = @_;
	my @range1 = grep {$_ ne '.'} keys %{$in};
	my @range2 = grep {$_ ne $except} @range1;
	if (scalar @range2 == 0) {
		return undef;
	} elsif (scalar @range2 == 1) {
		return $range2[0];
	} elsif (scalar @range2 > 1) {
		my @order = sort { $in->{$b} <=> $in->{$a} || $a<=>$b } @range2;
		return $order[0];
	} elsif (scalar @range1 == 1) {
		return $range1[0];
	} else {die;}
}
sub biggerKey($) {
	my $in = $_[0];
	my @order = sort { $in->{$b} <=> $in->{$a} || $a<=>$b } grep {$_ ne '.'} keys %{$in};
	return $order[0];
}
sub decodeGT($$) {
	my ($gt,$array) = @_;
	my @sGTs = split /\//,$gt;
	my %GTs = map { $_ => 1 } @sGTs;
	@sGTs = sort {$a cmp $b} keys %GTs;
	my $ret = '';
	for (@sGTs) {
		if (/^\d+$/) {
			$ret .= $array->[$_];
		}
	}
	$ret = '.' if length $ret == 0;
	return $ret;
}
sub TermColorGT($$) {
	my ($thisGT,$ref) = @_;
	my $ret;
	if (length $thisGT > 1) {
		return GREEN . $thisGT;
	} elsif ($thisGT eq $ref) {
		$ret = YELLOW . $thisGT;
	} elsif ($thisGT eq '.') {
		$ret = WHITE . $thisGT;
	} else {
		$ret = BLUE . $thisGT;
	}
	return $ret.' ';
}

open OUTT,'>',$outfs.'.txt' or die $!;
open OUTS,'>',$outfs.'.svg' or die $!;
print OUTT "### $DomRec [$cmd $bcfs]\n";
print OUTT "## Case_Samples: ",join(', ',@CaseS),"\n";
print OUTT "## Control_Samples: ",join(', ',@ControlS),"\n";
print OUTT join("\t",qw/#Chr Pos Case_Samples/)," . Control_Samples\n";
$fh = openpipe($cmd,$bcfs);
while (<$fh>) {
	chomp;
	my ($chr,$pos,$refalt,$qual,$adp,@sampleDat) = split /\t/;
	next if $qual < 20;
	my @GTs = split /,/,$refalt;
	my $t;
	for (@GTs) {
		$t = length $_;
		last if $t>1;
	}
	next if $t>1;	# skip indels
	#print join(',',$chr,$pos,$qual,$adp),"\t@GTs\t@sampleDat\n";
# gi|753572091|ref|NC_018727.2|,152998276,999,136	A C	FCAP055,1/1,5,18 FCAP056,1/1,8,27 FCAP059,1/1,3,12 FCAP066,1/1,2,9 FCAP067,1/1,6,21 FCAP069,1/1,6,21 FCAP072,1/1,16,51 FCAP075,0/1,7,57 FCAP084,0/1,6,50 FCAP085,1/1,3,12 FCAP086,0/0,6,9 FCAP087,1/1,8,27 FCAP088,1/1,8,27 FCAP089,1/1,5,18 FCAP090,1/1,4,15 FCAP110,1/1,7,24 FCAP111,0/1,8,54 FCAP112,0/1,8,12 FCAP113,0/1,4,30
	my ($skipped,%SampleDat,%caseGT,%ctlGT)=(0);
	die if @Samples != @sampleDat;
	for (@sampleDat) {
		my ($id,$gt,$sdp,$squal) = split /,/,$_;
		if ($sdp==0 or $squal<20) {
			$gt = './.';
			++$skipped;
		}
		my @sGTs = split /\//,$gt;
		$SampleDat{$id} = $gt;
		if ($tfamSamplePheno{$id} == 2) {	# 2=aff/case
			++$caseGT{$_} for @sGTs;
		} elsif ($tfamSamplePheno{$id} == 1) {
			++$ctlGT{$_} for @sGTs;
		} else {die;}
	}
	#next if $skipped > 0.5 * scalar(@sampleDat);
	my $theGT;
	if ($DomRec eq 'D') {
		my $t = biggerKey(\%ctlGT);
		next unless defined $t;
		$theGT = biggerKeyExcept(\%caseGT,$t);
	} elsif ($DomRec eq 'R') {
		$theGT = biggerKey(\%caseGT);
	}
	next unless defined $theGT;
	#ddx \%caseGT,\%ctlGT,$theGT;
	#print decodeGT($theGT,\@GTs);
	print OUTT join("\t",$chr,$pos,
		join('  ',map {TermColorGT(decodeGT($_,\@GTs),decodeGT($theGT,\@GTs));} ($theGT,@SampleDat{@CaseS},'.',@SampleDat{@ControlS}))
	),RESET,"\n";
}
close $fh;
close OUTT;
close OUTS;
#ddx \@CaseS,\@ControlS;

__END__
./vcfplot.pl kinkcats.tfam mpileup_20150403.vcf.filtered.gz 'gi|753572091|ref|NC_018727.2|:151386958-153139134' FCAP114 D plotN114
./vcfplot.pl kinkcats.tfam mpileup_20150403.vcf.filtered.gz 'gi|753572091|ref|NC_018727.2|:151386958-153139134' none D plotN114a

gi|753572091|ref|NC_018727.2|:151586958-152939134
-200k