#!/bin/env perl
use strict;
use warnings;

die "perl $0 <psnp file(.add_ref)> <out prefix> [minCov,minHom]\n" unless (@ARGV);

my ($minCov,$minHom);
if ($ARGV[2]) {
	($minCov,$minHom)=split /,/,$ARGV[2];
	warn "[!]Cov>=$minCov, Hom>=$minHom\n";
} else { ($minCov,$minHom)=(2,2); }
open PSNP,'<', $ARGV[0] or die "$!\n";
chomp(my $title=<PSNP>);
my @Samples;
open OUT,'>', "$ARGV[1].psnp" or die "$!\n";
if ($title =~ /^#ChrID\t/) {
	print OUT $title,"\n";
	@Samples = split / /,(split /\t/,$title)[-1];
} else {
	seek PSNP,0,0;
}
print STDERR '[!]Samples: [',join('|',@Samples),"]\n";

open NFO,'>', "$ARGV[1].hcnfo" or die "$!\n";
print NFO "Hom\tCov\n";
my ($In,$Out,%Cov,%Hom)=(0,0);
while (<PSNP>) {
	chomp;
	++$In;
	my ($chr,$pos,$ref,$tail)=split /\t/;
	my @indSNP=split / /,$tail;	# /[ACGTRYMKSWHBVDNX-]/
	my ($hom,$cov,$all)=(0,0,0);
	for (@indSNP) {
		++$all;
		++$hom if /[ATCG]/i;
		++$cov unless /-/;
	}
	$hom /= $cov;
	$cov /= $all;
	print NFO "$hom\t$cov\n";
	if ($hom >= $minHom and $cov >= $minCov) {
		print OUT join("\t",$chr,$pos,$ref,$tail),"\n";
		++$Out;
	}
}
close NFO;
close OUT;
warn "\n[!] $In -> $Out, ",int(10000*$Out/$In)/100,"%\n";

__END__
./psnphc.pl psnp/Chr02.add_ref Chr02 0.95,0.99999
cat 9311/chrorder |while read a; do echo ./psnphc.pl ./psnp/$a.add_ref $a 0.95,0.999;done > ff.sh

perl -lane 'next if $F[8]<70;next if $F[8]>170;print join("\t",@F[0,1,2,7,8],$F[8]/($F[7]+$F[8]))' psnp/Chr01_final.snp
cat 9311/chrorder |while read a;do perl -lane 'next if $F[8]<70;next if $F[8]>170;print join("\t",@F[0,1,2,7,8],$F[8]/($F[7]+$F[8]))' psnp/${a}_final.snp > ff$a.snp;done

cat 9311/chrorder |while read a;do perl -e "my %Dat;open IN,'<','ff$a.snp';
while (<IN>) {\$x=(split /\t/)[1]; ++\$Dat{\$x};}
close IN;
open T,'<','psnp/$a.add_ref';
while (<T>) {
	@F=split /\t/;
	next unless \$Dat{\$F[1]};
	print \$_;
}
close T;
" > aa$a.snp;
done
