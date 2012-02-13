#!/bin/env perl
use strict;
use warnings;

my %dat=(
'Drone1' => './bam/LWX-D1.bam.d.vcf',
'Scouts1' => './bam/LWX-S1.bam.d.vcf',
'Scouts2' => './bam/LWX-S2.bam.d.vcf',
'Recruit2' => './bam/LWX-R2.bam.d.vcf',
'Recruit3' => './bam/LWX-R3.bam.d.vcf',
);

my %SNP; # chr,pos,{GT with ref}
my @Samples = sort keys %dat;

for my $s (@Samples) {
	open VCF,'<',$dat{$s} or die "Error: $!\n";
	print "Loading $s\n";
	while (<VCF>) {
		next if /^#/;
		my ($chr,$pos,$ID,$ref,$alt,$QUAL,$FILTER,$INFO,$FORMAT,$data)=split /\t/;
		next if $INFO =~ /INDEL;/;
		next if $QUAL < 20;
		die unless $FORMAT eq 'GT:PL:GQ';
		my ($GT,$PL,$GQ)=split /:/,$data;
		my @GTs=split /[\/|]/,$GT;
		if ($s eq 'Drone1') {
			next if $GTs[0] ne $GTs[1];
		}
		my @GeneTypes=split /,/,$alt;
		unshift @GeneTypes,$ref;
		$SNP{$chr}{$pos}{'Ref'}=$ref . $ref;
		$SNP{$chr}{$pos}{$s}=$GeneTypes[$GTs[0]] . $GeneTypes[$GTs[1]];
#print "$s:{$chr}{$pos}{$s} $ref $GeneTypes[$GTs[0]] $GeneTypes[$GTs[1]] $GT <@GTs>\n";
	}
	close VCF;
}

open O,'>','bees.mfa' or die "Error: $!\n";
for my $s ('Ref',@Samples) {
	print O ">$s\n";
	print "Writing $s\n";
	for my $chr (sort keys %SNP) {
		for my $pos (sort {$a<=>$b} keys %{$SNP{$chr}}) {
			if (exists $SNP{$chr}{$pos}{$s}) {
				print O $SNP{$chr}{$pos}{$s};
			} else {
				print O $SNP{$chr}{$pos}{'Ref'};
			}
		}
	}
	print O "\n";
}
close O;
