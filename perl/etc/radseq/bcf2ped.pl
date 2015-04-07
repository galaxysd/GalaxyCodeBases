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

die "Usage: $0 <tfam file> <bcgv bcf> <min Sample Count> <Dominant/Recessive> <out> [sample list]\n" if @ARGV<5;
my $tfamfs=shift;
my $bcfs=shift;
my $minSampleCnt=shift;	# 16 for 16 samples in paper
my $DomRec=shift;
my $outfs=shift;
my $sampleList=shift;

$DomRec='D' if $DomRec =~ /^D/i;
$DomRec='R' if $DomRec =~ /^R/i;

my (%Stat,$t);
open OP,'>',$outfs.'.tped' or die "Error opening $outfs.tped : $!\n";
open OPA,'>',$outfs.'.case.tped' or die $!;
open OPO,'>',$outfs.'.control.tped' or die $!;

open OM,'>',$outfs.'.MinorAllele' or die "Error opening $outfs.MinorAllele : $!\n";
open OD,'>',$outfs.'.dict' or die "Error opening $outfs.dict : $!\n";
if ($tfamfs ne $outfs.'.tfam') {
	open OF,'>',$outfs.'.tfam' or die "Error opening $outfs.tfam : $!\n";
	open OFA,'>',$outfs.'.case.tfam' or die $!;
	open OFO,'>',$outfs.'.control.tfam' or die $!;
} else {die;}
open O,'>',$outfs.'.bcf2pedlog' or die "Error opening $outfs.bcf2pedlog : $!\n";
$t = "# In:[$bcfs], Out:[$outfs]. minSampleCnt=$minSampleCnt Dominant/Recessive:[$DomRec]\n";
print O $t;
print $t;

#open OV,'>',$outfs.'.vcf' or die "Error opening $outfs.vcf : $!\n";

my (%Pheno,@tfamSamples,%tfamDat,%tfamSampleFlag,%inFamily,%ISinFamily);
open L,'<',$tfamfs or die;
while (<L>) {
	next if /^(#|$)/;
	#print OF $_;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	$tfamDat{$ind} = $_."\n";
	if ($ind =~ s/^~//) {
		$tfamSampleFlag{$ind} = 0;
	} elsif ($pho == 1 or $pho == 2) {
		$tfamSampleFlag{$ind} = $pho;
	} else { die; }	# $pho can only be 1 or 2
#	$tfamSampleFlag{$ind} = 3 if $ind !~ /^GZXJ/;
	push @tfamSamples,$ind;
	push @{$inFamily{$family}},$ind;
	$Pheno{$ind} = $pho;	# disease phenotype (1=unaff/ctl, 2=aff/case, 0=miss)
}
for my $ind (@tfamSamples) {
	if ($tfamSampleFlag{$ind} == 0) {
		next;
	} elsif ($tfamSampleFlag{$ind}) {
		if ($tfamSampleFlag{$ind} & 1) {
			print OFO $tfamDat{$ind};
		} 
		if ($tfamSampleFlag{$ind} & 2) {
			print OFA $tfamDat{$ind};
		}
	} else { die; }
	print OF $tfamDat{$ind};
}
close L;
close OF;
close OFA;
close OFO;
for my $family (keys %inFamily) {
	if (scalar @{$inFamily{$family}} == 1) {
		delete $inFamily{$family};
	} else {
		$ISinFamily{$_}=$family for @{$inFamily{$family}};
	}
}
ddx \%inFamily; ddx \%ISinFamily;

my $cmd = 'bcftools view --types snps -m2 -M2';
if (defined $sampleList) {
	$cmd .= "-S $sampleList";
}
print O "Open:[$cmd] [$bcfs]\n";
print "Open:[$cmd] [$bcfs]\n";;
my $th = openpipe($cmd,$bcfs);
my (@Samples,@Parents);
while (<$th>) {
	#print OV $_;
	next if /^##/;
	chomp;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {
			if (/^(\w+)_all\.bam$/) {
				$_ = $1;
			} elsif (/\.\/bam0\//) {
				my $t=(split /_/)[-1]; $t=(split /\./,$t)[0]; $_=$t;
			} else {
				#my $t=(split /\//)[-1];$t=~s/_cut//g;$t=~s/-/./g; $_=join('_',(split /\./,$t)[-5,-6]);
				s/_2$//;
				my $t = $_;
			}
		} splice @data,9;
		# ../5.bam_0000210210_merged/d1_4_merged.D4.JHH001.XTU.sort.rmdup.bam
		#@Parents = grep(!/^GZXJ/,@Samples);
		@Parents=();
		last;
	}
}

die "Sample Count in tfam and bcf not match !\n" if @tfamSamples != @Samples;
my @res;
for (my $i = 0; $i < @Samples; $i++) {
	push @res,join(',',$tfamSamples[$i],$tfamSampleFlag{$Samples[$i]});
	die "Samples in tfam and bcf not match ! [$tfamSamples[$i] <> $Samples[$i]]\n" if $tfamSamples[$i] ne $Samples[$i];
}	# http://stackoverflow.com/questions/2591747/how-can-i-compare-arrays-in-perl
print O "# Samples: [",join('],[',@res),"]\t# 1=control, 2=case, 0=drop\n# Parents: [",join('],[',@Parents),"]\n";
warn "Samples: (1=control, 2=case, 0=drop)\n[",join("]\n[",@res),"]\nParents: [",join('],[',@Parents),"]\n";

my (%ChrID2Num);
sub genChrNumber($) {
	my $id = $_[0];
	if (exists $ChrID2Num{$id}) {
		return $ChrID2Num{$id};
	} else {
		my $lastCnt = scalar keys %ChrID2Num;
		$ChrID2Num{$id} = $lastCnt+1;
		return $ChrID2Num{$id};
	}
}

while (<$th>) {
	next if /^#/;
	chomp;
	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	++$Stat{'VCF_In'};	# also as rs#
	my @Bases = split /\s*,\s*/,"$REF, $ALT";
	my @groups = split(/\s*;\s*/, $INFO);
	my (%INFO,$name);
#	if ($groups[0] eq 'INDEL') {
#		$INFO{'Type'} = 'INDEL';
#		shift @groups;
#	} else {
#		$INFO{'Type'} = 'SNP';
#	}
	for my $group (@groups) {
		my ($tag,$value) = split /=/,$group;
#warn "- $group -> $tag,$value\n";
		my @values = split /,/,$value;
		if (@values == 1) {
			$INFO{$tag}=$values[0];
		} else {
			$INFO{$tag}=\@values;
		}
	}
	my (%GT,%GTcnt);
	my @FMT = split /:/,$FORMAT;
	for my $s (@Samples) {
		my $dat = shift @data or die "bcf file error.";
		my @dat = split /:/,$dat;
		for my $i (@FMT) {
			$GT{$s}{$i} = shift @dat;
		}
	}
	my $SPcnt = 0;
	my (%GTitemCnt,$Mut,@plinkGT);
	for (@Samples) {
		if ($GT{$_}{'DP'} > 0 and $GT{$_}{'GQ'} >= 20) {
			my $gt = $GT{$_}{'GT'};
			++$GTcnt{$gt};
			my @GT = split /[\/|]/,$gt;
			++$SPcnt;
			if ($DomRec eq 'R') {
				if ($Pheno{$_} == 2) {
					++$GTitemCnt{$_} for @GT;
				}
			} elsif ($DomRec eq 'D') {
				if (exists $ISinFamily{$_}) {
					++$GTitemCnt{$_} for @GT;
				}
			} else {die;}
			$GT{$_}{'GTp'} = join ' ',map($Bases[$_], @GT);
			#$GT{$_}{'O_K'} = 1;
		} else {
			$GT{$_}{'GTp'} = '0 0';
			#$GT{$_}{'O_K'} = 0;
		}
		push @plinkGT,$GT{$_}{'GTp'};
	}
	if ($DomRec eq 'R') {
		($Mut) = sort { $GTitemCnt{$b} <=> $GTitemCnt{$a} } keys %GTitemCnt; # Desc
	} elsif ($DomRec eq 'D') {
		($Mut) = sort { $GTitemCnt{$a} <=> $GTitemCnt{$b} } keys %GTitemCnt; # Asc
	} else {die;}
	unless (defined $Mut) {
		++$Stat{'VCF_noAffInd_Skipped'};
		next;
		#$Mut = -1;
	}
	unless ($Mut) {
		++$Stat{'VCF_noMUT_Count'};
		#next;
	}
	#$Mut = $Bases[$Mut];
=pod
++$Stat{'GTcnt'}{$INFO{'FQ'} <=> 0}{scalar(keys %GTcnt)};
ddx $Stat{'GTcnt'};
ddx $CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO,\%INFO,\%GT if scalar(keys %GTcnt) > 1 and $INFO{'FQ'} < 0 and $SPcnt>2;
# rad2marker.pl:135: {
#   -1 => { 1 => 2850, 2 => 526, 3 => 37 },
#   1  => { 1 => 8, 2 => 2507, 3 => 1792 },
# }
=cut
#warn "$CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO\n";
#ddx $CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO,\%GTcnt,\%INFO,\%GT,\%GTitemCnt,$Mut;
	if ($QUAL<20 or scalar(keys %GTcnt)<2 or $SPcnt < $minSampleCnt) {	# No 'FQ' in 'INFO' for v1.1
		++$Stat{'VCF_Skipped'};
		next;
	}
	$Mut = $Bases[$Mut];
	my $SNPid = "r".$Stat{'VCF_In'};
	my $ChrID = genChrNumber($CHROM);
	#print OP join("\t",$1,$SNPid,0,$POS,@plinkGT),"\n";
	print OM join("\t",$SNPid,$Mut),"\n";
	print OD join("\t",${CHROM},${POS},$SNPid),"\n";
	my (@GTall,@GTcase,@GTcontrol);
	for my $i (0 .. $#tfamSamples) {
		my $ind = $tfamSamples[$i];
		if ($tfamSampleFlag{$ind} == 0) {
			next;
		} else {
			push @GTall,$plinkGT[$i];
			if ($tfamSampleFlag{$ind}) {
#				if ($tfamSampleFlag{$ind} == 3) {
#					$plinkGT[$i] = '0 0';
#				}
				if ($tfamSampleFlag{$ind} & 1) {
					push @GTcontrol,$plinkGT[$i];
				}
				if ($tfamSampleFlag{$ind} & 2) {
					push @GTcase,$plinkGT[$i];
				}
			} else { die; }
		}
	}
	print OP join("\t",$ChrID,$SNPid,0,$POS,@GTall),"\n";
	print OPA join("\t",$ChrID,$SNPid,0,$POS,@GTcase),"\n";
	print OPO join("\t",$ChrID,$SNPid,0,$POS,@GTcontrol),"\n";
	++$Stat{'Marker_Out'};
	#print OV $_,"\n";
}
close $th;

close OP;
close OPA;
close OPO;
close OM;
close OD;
#close OV;

print O Dumper(\%ChrID2Num);
print O Dumper(\%Stat);
close O;

ddx \%Stat;

print "[Prepare $outfs.phe.] And then:\np-link --tfile $outfs --reference-allele $outfs.MinorAllele --fisher --out ${outfs}P --model --cell 0 --allow-no-sex\n--pheno $outfs.phe --all-pheno [screen log will be saved by p-link itselt]\n";
__END__

bcftools view -s^FCAP114 -m2 mpileup_20150321HKT071334.vcf.gz \
 | bcftools norm -Df ../ref/Felis_catus80_chr.fa -c e -m+both \
 | bcftools filter -sLowQual -e'%QUAL<10' \
 | bcftools filter -m+ -sDepthHigh -e'DP>650' \
 | bcftools filter -m+ -sDepthLow -e'DP<2' \
 | bcftools filter -m+ -sBadSites -e'%QUAL<10 && RPB<0.1' \
 | tee >(bcftools view -Oz -o filter_20150321HKT071334.vcf.gz) \
 | bcftools view -f .,PASS -Oz -o filtered.vcf.gz &

bcftools query -f '%CHROM,%POS\t%REF|%ALT|%QUAL\t%DP[\t%SAMPLE=%GT,%DP]\n' filtered.vcf.gz |les

./bcf2ped.pl ~/work/catbtail/merged/kinkcats.tfam ~/work/catbtail/merged/mpileup_20150321HKT071334.vcf.gz 18 D xxxxx 2>xxxxx



=== Old ===
grep -hv \# radseq.gt > radseq.tfam

./bcf2bed.pl radseq.tfam radseq.bcgv.bcf radseq 2>&1 | tee radseq.pedlog

grep REC radseq.p.snow.model > radseq.p.snow.model.REC
grep DOM radseq.p.snow.model > radseq.p.snow.model.DOM
sort -nk8 radseq.p.snow.model.REC > radseq.p.snow.model.REC.sortnk8 &
sort -nk8 radseq.p.snow.model.DOM > radseq.p.snow.model.DOM.sortnk8 &
bcftools view -I radseq.bcgv.bcf |grep -v \# |cat -n > radseq.bcgv.bcf.rs &

bcftools view tigers.bcgv.bcf|grep -v \## > tigers.bcgv.vcf &

join -1 3 -2 2 <(sort -k3 radall.dict) <(sort -k2 radallP.model) > tmp &
#sed -n '75,380 p' tmp.s.REC.nk3.scaffold75

cat tmp.s.REC|perl -lane '@a=split /\//,"$F[7]/$F[8]";$sum=$a[0]+$a[1]+$a[2]+$a[3];
$theta=($a[1]+$a[2])/($sum*2);
$p=((1-$theta)**(2*$sum-($a[1]+$a[2])))*($theta**($a[1]+$a[2]))/(0.5**(2*$sum));
$LOD=int(0.5+log($p)*1000/log(10))/1000;
print join("\t",@F,$sum,int(0.5+$theta*100000)/100000,int($LOD),$LOD)' > rec.pa

sort -nk13 -k2 -nk3  rec.pa > rec.pas

grep -P "scaffold75\t" rec.pas > scaffold75.pas
sort -nk3 scaffold75.pas > scaffold75.pas.nk3
perl -lane 'print if $F[10]>16' scaffold75.pas.nk3 > scaffold75.pas.nk3.16

grep -P "scaffold1458\t" rec.pas > scaffold1458.pas
sort -nk3 scaffold1458.pas > scaffold1458.pas.nk3
perl -lane 'print if $F[10]>16' scaffold1458.pas.nk3 > scaffold1458.pas.nk3.16

-------------------

./cat/dojoin.pl radall.dict 3 radallP.model.REC 3 tmp.REC

awk '{print $3" "$1" "$2" "$4" "$6" "$7" "$8" "$9" "$10" "$11}' tmp.REC.out > tmp.REC.out.j

cat tmp.REC.out.j|perl -lane '@a=split /\//,"$F[7]/$F[8]";$sum=$a[0]+$a[1]+$a[2]+$a[3];
$theta=($a[1]+$a[2])/($sum*2);
$p=((1-$theta)**(2*$sum-($a[1]+$a[2])))*($theta**($a[1]+$a[2]))/(0.5**(2*$sum));
$LOD=int(0.5+log($p)*1000/log(10))/1000;
print join("\t",@F,$sum,int(0.5+$theta*100000)/100000,int($LOD),$LOD)' |sort -k2,2 -k3,3n > rec.npa

sort -k13,13n -k2,2 -k3,3n  rec.npa > rec.npas

Tue Aug 13 14:44:52 CST 2013

find *.model|while read a;do grep REC $a > $a.REC.0;done
find *.model|while read a;do perl -lane 'if (@F<8) {splice @F,3,0,".";} die if @F<8; print join("\t",@F)' $a.REC.0 > $a.REC;done
find *.model|while read a;do sort -k 8,8n -k 2,2 $a.REC > $a.REC.nk8;done
find *.dict|perl -lane 's/\.dict//;system "./dojoin.pl $_.dict 3 ${_}P.model.REC.nk8 2 ${_}R.REC"'
find *.dict|perl -lane 's/\.dict//;print "awk ",chr(39),"BEGIN {OFS=\"\\t\"} {print \$3,\$1,\$2,\$4,\$6,\$7,\$8,\$9,\$10,\$11}",chr(39)," ${_}R.REC.out > ${_}R.REC.out.j" '

