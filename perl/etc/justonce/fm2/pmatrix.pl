#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dump qw(ddx);

my %idF;
while (<DATA>) {
	chomp;
	my ($id,$fname) = split /\s+/;
	next unless defined $fname;
	$idF{$id} = $fname;
}
#ddx \%idF;

my %GT;
for my $id (keys %idF) {
	open I,'<',$idF{$id} or die $?;
	while (<I>) {
		chomp;
		my ($rsid,$depth,$qual,$tgt) = split /\t/;
		my @GTs = split /\//,$tgt;
		my $gt = join('',sort @GTs);
		$GT{$rsid}{$id} = $gt;
	}
	close I;
}
#ddx \%GT;

my @ids = sort keys %idF;
my %Cnt;
for my $rsid (keys %GT) {
	my $oneRS = $GT{$rsid};
	for my $id1 (@ids) {
		for my $id2 (@ids) {
			if ($oneRS->{$id1} eq $oneRS->{$id2}) {
				++$Cnt{$id1}{$id2};
			}
		}
	}
}
ddx \%Cnt;

print join("\t",('ID',@ids)),"\n";
for my $id1 (@ids) {
	print $id1;
	for my $id2 (@ids) {
		print "\t",$Cnt{$id1}{$id2};
	}
	print "\n";
}

# find test-output2 -name snp0.txt|perl -lane '@x=split /\//,$_;print "$x[2]\t$_"' >pmatrix.tsv

__DATA__
19M940F	test-output2/samples/19M940F/SNP/snp0.txt
NA12878	test-output2/samples/NA12878/SNP/snp0.txt
19M1202FN	test-output2/samples/19M1202FN/SNP/snp0.txt
19M1196F	test-output2/samples/19M1196F/SNP/snp0.txt
YH-DNA	test-output2/samples/YH-DNA/SNP/snp0.txt
2800M	test-output2/samples/2800M/SNP/snp0.txt
LD-DNA1	test-output2/samples/LD-DNA1/SNP/snp0.txt
LD-DNA2	test-output2/samples/LD-DNA2/SNP/snp0.txt

