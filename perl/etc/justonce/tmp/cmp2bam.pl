#!/usr/bin/env perl
use strict;
use warnings;
#use IO::Unread;

die "Usage: $0 <bigger_bam> <sorted_bam> <out>\n" if @ARGV < 3;
my ($bam1,$bam2,$out)=@ARGV;

open IN1,"-|","samtools view $bam1" or die "Error opening $bam1: $!\n";
open IN2,"-|","samtools view $bam2" or die "Error opening $bam2: $!\n";

my %Dat;
while (<IN1>) {
	chomp;
	my ($id,$flag,$chr,$pos) = split /\t/;
	#$Dat{$id} = [$chr,$pos];
	$Dat{$id} = 1;
}
close IN1;

my ($lastCP,@CurrArray)=('');
open OUT,'>',$out or die "Error opening $out: $!\n";
while (<IN2>) {
	#chomp;
	my ($id,$flag,$chr,$pos) = split /\t/;
	my $t = "$chr,$pos";
	if ($lastCP eq $t) {
		push @CurrArray,[$id,$_];
	} else {
		$lastCP = $t;
		my $flag=0;
		for my $t (@CurrArray) {
			unless (exists $Dat{$t->[0]}) {
				$flag=1;
				$t->[1] = '*'.$t->[1];
				#print $t->[0],"\n";
			}
		}
		if ($flag==1) {
			for my $t (@CurrArray) {
				print OUT $t->[1];
			}
			print OUT "\n";
		}
		@CurrArray=([$id,$_]);
	}
}
close IN2;
close OUT;

__END__
./cmp.pl Tiger_aln_bam/pti096_clean_aln_pe.bam Tiger_aln_rmdup/pti096_clean_aln_pe_rmdup.bam pti096.cmp.sam &

# http://www.nature.com/nature/journal/v399/n6737/full/399682a0.html
# Cultures in chimpanzees
