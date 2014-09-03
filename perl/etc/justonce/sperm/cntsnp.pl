#!/bin/env perl
use strict;
use warnings;
use Data::Dump;

die "Usage: $0 <blood vcf.gz> <sperm vcf.gz>\n" if @ARGV < 2;
my ($inbd,$insp)=@ARGV;

open IBD, "-|", "zcat $inbd" or die "$!";
open ISP, "-|", "zcat $insp" or die "$!";
open O,'>',"$insp.stat";

my %SNP;

sub getGT($) {
	my $str = $_[0];
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$INFO,undef,$data) = split /\t/,$str;
	my @GeneTypes=split /,/,$alt;
	unshift @GeneTypes,$ref;
	my $GT = (split /:/,$data)[0];
	my @GTs=split /[\/|]/,$GT;
	my @ret = ( $chr, $pos, $GeneTypes[$GTs[0]] . $GeneTypes[$GTs[1]] );
	return \@ret
}

while (<ISP>) {
	next if /^#/;
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$INFO,undef,$data) = split /\t/;
	next if $INFO =~ /INDEL;/;
	next unless $chr =~ /^chr\d+/;
	my $GT;
	($chr,$pos,$GT) = @{getGT($_)};
	$SNP{$chr}{$pos}{'SP'} = $GT;
}
while (<IBD>) {
	next if /^#/;
	my ($chr,$pos,undef,$ref,$alt,undef,undef,$INFO,undef,$data) = split /\t/;
	next if $INFO =~ /INDEL;/;
	next unless $chr =~ /^chr\d+/;
	my $GT;
	($chr,$pos,$GT) = @{getGT($_)};
	$SNP{$chr}{$pos}{'BD'} = $GT;
}
#ddx \%SNP;
my (%Results,%subResults);
for my $chr (keys %SNP) {
	for my $pos (keys %{$SNP{$chr}}) {
		my ($GT1,$GT2) = ( $SNP{$chr}{$pos}{'BD'}, $SNP{$chr}{$pos}{'SP'} );
		unless ( defined $GT1 ) {
			++$Results{$chr}{'SpermOnly'};
			++$Results{'_All_'}{'SpermOnly'};
			next;
		}
		unless ( defined $GT2 ) {
			++$Results{$chr}{'BloodOnly'};
			++$Results{'_All_'}{'BloodOnly'};
			next;
		}
		my @gt1 = split //,$GT1;
		my @gt2 = split //,$GT2;
		my ($flag1,$flag2);
		if ($gt1[0] eq $gt1[1]) { $flag1='Hom'; } else { $flag1='Het'; }
		if ($gt2[0] eq $gt2[1]) { $flag2='Hom'; } else { $flag2='Het'; }
#print "$gt1[0] $gt1[1] $gt2[0] $gt2[1] $flag1,$flag2\n";
		my ($str,$subtype);
		if ($GT1 eq $GT2) {
			$str = 'Same';
		} else {
			$str = $flag1 . '-' . $flag2;
			if ( $flag1 eq $flag2 ) {
				if ($flag1 eq 'Hom') {
					$subtype = 'Hom' . $gt1[1] . $gt2[1];
				} else {
					$subtype = 'Het' . $GT1 . $GT2;
				}
			}
		}
		++$Results{$chr}{$str};
		++$Results{'_All_'}{$str};
		if ( defined $subtype ) {
			++$subResults{$chr}{$subtype};
			++$subResults{'_All_'}{$subtype};
		}
	}
}
ddx \%subResults;
ddx \%Results;

for my $chr (sort keys %Results) {
	print O "$chr\t";
	for my $type (sort keys %{$Results{$chr}}) {
		print O "$type: $Results{$chr}{$type}\t";
	}
	if (exists $subResults{$chr}) {
		print O "|\t";
		for my $type (sort keys %{$subResults{$chr}}) {
			print O "$type: $subResults{$chr}{$type}\t";
		}
	}
	print O "\n";
}

close IBD;
close ISP;
close O;

__END__
perl cntsnp.pl blood.gz sperm23.gz

perl cntsnp.pl blood.mda.filter.gz sperm23.vcf.filter.gz



