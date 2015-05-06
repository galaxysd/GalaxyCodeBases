#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;

die "Usage: $0 <LD dist tick> <plink LD file> <dict file> <out>\n" if @ARGV<3;
my $disttick=shift;
my $ldfs=shift;
my $dictfs=shift;
my $outfs=shift;

my $r2tick = 0.01;

my %Dict;
open D,'<',$dictfs or die $?;
while (<D>) {
	chomp;
	my ($chr,$pos,$rsid) = split /\t/;
	$Dict{$rsid} = [$chr,$pos];
}
close D;

my (%Heat,%Dist,%R2,%Stat,$DotCnt);
my $ldfh = openfile($ldfs);
chomp(my $firstLine = <$ldfh>);
die unless $firstLine =~ /\sR2 $/;
while (<$ldfh>) {
	chomp;
	my @arr = split /\s+/;
	# ddx \@arr;	# plotLD.pl:29: ["", 0, 21, "r9775647", 0, 8017, "r1680049", 0.25]
	my $snpr1 = $Dict{$arr[3]};
	my $snpr2 = $Dict{$arr[6]};
	if ($snpr1->[0] ne $snpr2->[0]) {
		++$Stat{'Skipped'};
		next;
	}
	++$Stat{'Used'};
	my $dist = abs($snpr1->[1] - $snpr2->[1]);
	$dist = int($dist/$disttick)*$disttick;
	my $r2 = int($arr[7]/$r2tick)*$r2tick;
	++$Heat{$dist}{$r2};
	++$Dist{$dist};
	++$R2{$r2};
	++$DotCnt;
	# ddx \%Heat;
}
close $ldfh;
ddx \%Stat;

open O,'>',$outfs or die $!;
my @DistKeys = sort {$a<=>$b} keys %Dist;
my @R2Keys = sort {$a<=>$b} keys %R2;
print O "# DistMax: $DistKeys[-1], R2Max: $R2Keys[-1]\n";
for my $distT (0 .. $DistKeys[-1]/$disttick) {
	my $dist = $distT*$disttick;
	my $line = "[$dist]";
	for my $r2T (0 .. $R2Keys[-1]/$r2tick) {
		my $r2 = $r2T*$r2tick;
		my $t = 0;
		$t = $Heat{$dist}{$r2} if exists $Heat{$dist}{$r2};
		$line .= "\t$t";
	}
	print O $line,"\n";
}
print O "\n# Dist\n";
for my $dist (@DistKeys) {
	print O $dist,"\t",$Dist{$dist},"\t",$Dist{$dist}/$DotCnt,"\n";
}
print O "\n# R2\n";
for my $r2 (@R2Keys) {
	print O $r2,"\t",$R2{$r2},"\t",$R2{$r2}/$DotCnt,"\n";
}
close O;


__END__
perl plotLD.pl 5000 tiger650.filtered.Amur.r2.ld.bz2 tiger650.filtered.Amur.dict tiger650.filtered.Amur.ld.out

cp plotLD.pl /share/users/huxs/work/tigerLYC/Tiger_SNP_Final/

./bcf2ped.pl samples.tfam tiger650.filtered.gz 3 R tiger650.filtered.Amur Amur.lst
p-link --tfile tiger650.filtered.Amur --r2 --out tiger650.filtered.Amur.r2 --allow-no-sex --ld-window-kb 60000 --ld-window-r2 0 --ld-window 10000
-rw-r--r-- 1 huxs users 2006976400720 May  5 01:22 tiger650.filtered.Amur.r2.ld

perl -lane 'last if /^# Dist$/;next if /^# DistMax/;shift @F;print join("\t",@F);' tiger650.filtered.Amur.ld.out > t.tsv

gnuplot

set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set xrange [0:100]
set cbrange [0:20000]
set cblabel "Cnt"
set xlabel "r2"
set ylabel "Dist5k"
plot 't.tsv' matrix with image

set dgrid3d 50
set cbrange [0:5000]
splot 't.tsv' matrix with pm3d
