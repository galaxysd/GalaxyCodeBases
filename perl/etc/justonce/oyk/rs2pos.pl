#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;

use Data::Dump qw(ddx);

use lib '.';

sub getPE(@) {
	my @AFs = @_;
	my ($p1,$p2)=(0,0);
	for my $i (0 .. $#AFs) {
		my $p1i = $AFs[$i]*(1-$AFs[$i])*(1-$AFs[$i]);
		$p1 += $p1i;
		for my $j (($i+1) .. $#AFs) {
			my $p2i = $AFs[$i]*$AFs[$i]*$AFs[$j]*$AFs[$j]*(4-3*$AFs[$i]-3*$AFs[$j]);
			$p2 += $p2i;
		}
	}
	return $p1-$p2;
}
sub getDPE(@) {
	my @AFs = @_;
	my ($p1,$p2)=(0,0);
	for my $i (0 .. $#AFs) {
		my $p1i = $AFs[$i]*$AFs[$i]*(1-$AFs[$i])*(1-$AFs[$i]);
		$p1 += $p1i;
		for my $j (($i+1) .. $#AFs) {
			my $tb = 1-$AFs[$i]-$AFs[$j];
			my $p2i = 2*$AFs[$i]*$AFs[$j]*$tb*$tb;
			$p2 += $p2i;
		}
	}
	return $p1+$p2;
}
sub get0PE(@) {
	my @AFs = @_;
	my ($p1,$p2)=(0,0);
	for my $i (0 .. $#AFs) {
		my $p1i = $AFs[$i]*(1-$AFs[$i])*(1-$AFs[$i]);
		$p1 += $p1i;
		for my $j (0 .. $#AFs) {
			next if $i == $j;
			my $p2i = $AFs[$i]*$AFs[$i]*$AFs[$j]*$AFs[$j]*(4-3*$AFs[$i]-3*$AFs[$j]);
			$p2 += $p2i;
		}
	}
	return $p1-$p2/2;
}

die "Usage: $0 <in.tsv> ..\n" if @ARGV<1;
my (%Markers,%MarkerAF,$sum);
while (<>) {
	chomp;
	my ($rs,$chr,$pos,@d) = split /\t/,$_;
	$MarkerAF{$rs} = {@d};
	my @AFs = values %{$MarkerAF{$rs}};
	#$sum = eval join '+',@AFs;
	#print "$rs: $sum\n" if $sum != 1;
	my $PE = getPE(@AFs);
	$Markers{$rs} = [$chr,$pos,$PE,'unlocalized',0,'unlocalized',0];
}
#ddx \%Markers;

# zgrep rs112459276 /hwfssz1/pub/database/ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_X.bed.gz
# zgrep rs112459276 /hwfssz1/pub/database/ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/bed_chr_X.bed.gz
my @hg19fs = glob '/hwfssz1/pub/database/ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/bed_chr_*.bed.gz';
my @hg38fs = glob '/hwfssz1/pub/database/ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_*.bed.gz';
#ddx \@hg19fs,\@hg38fs;
sub readPos(@) {
	my $offect = shift @_;
	for (@_) {
		open( GZIN,"-|","gzip -dc $_") or die "Error opening $_: $!\n";
		while(<GZIN>) {
			next if /^track /;
			my ($chrom,$chromStart,$chromEnd,$name,$score,$strand)=split /\t/;
			next unless exists $Markers{$name};
			my $len = $chromEnd - $chromStart;
			warn "[!]Len:$len for [$_]\n" if $len > 1;
			$Markers{$name}->[$offect] = $chrom;
			$Markers{$name}->[$offect+1] = $chromEnd;
			#ddx \$Markers{$name};
			if ($offect == 3) {
				if ($Markers{$name}->[0] ne $chrom or $Markers{$name}->[1] ne $chromEnd) {
					warn "[!]Position Not match to hg19.\n";
					ddx \$Markers{$name};
				}
			}
		}
	}
}
readPos(3,@hg19fs);
readPos(5,@hg38fs);
my @ChrIDs = map {"chr$_"} (1..22,'X','Y','MT','M');
my $i = 0;
my %L = map { $_ => $i++ } @ChrIDs;
$i = 5;
#$i = 0;
my @rsids = sort {
	no warnings 'uninitialized';
	exists $L{$Markers{$a}->[$i]} <=> exists $L{$Markers{$b}->[$i]} ||
	$L{$Markers{$a}->[$i]} <=> $L{$Markers{$b}->[$i]} ||
	$Markers{$a}->[$i+1] <=> $Markers{$b}->[$i+1] ||
	$Markers{$a}->[1] <=> $Markers{$b}->[1] ||
	$a cmp $b
} keys %Markers;
print join("\t",'#SNPid','Chr.hg19','Pos.hg19','Chr.GRCh38','Pos.GRCh38','Alleles with Frequency'),"\n";
for my $id (@rsids) {
	my @d = @{$Markers{$id}};
	my @allels = sort { $MarkerAF{$id}{$a} <=> $MarkerAF{$id}{$b} || $a cmp $b } keys %{$MarkerAF{$id}};
	my @af = map {"$_\t$MarkerAF{$id}{$_}"} @allels;
	print join("\t",$id,@d[3,4,5,6],@af),"\n";
}

=pod
rs149540625 was merged into rs112459276 on August 21, 2014 (Build 142)
rs75611253 was merged into rs3852322 on July 1, 2015 (Build 144)
rs117314748 was merged into rs289279 on July 1, 2015 (Build 144)
rs200370602 was merged into rs10154714 on July 1, 2015 (Build 144)
rs201431024 was merged into rs4617205 on July 1, 2015 (Build 144)
rs113184075 was merged into rs77634512 on July 19, 2016 (Build 147)

rs2484385,rs10453900,rs61800290,rs144913592
unlocalized scaffold in GRCh38

./rs2pos.pl >nippt.tsv 2>nippt.tsv.err &
=cut
