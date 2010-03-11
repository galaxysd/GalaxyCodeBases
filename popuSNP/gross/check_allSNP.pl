#!/usr/bin/perl -w
use strict;
print STDERR $0,"\n";

my %genotyping;

# perl /share/raid002/lidong/SNP_filter_bin/check_allSNP.pl /share/raid002/lidong/40snp/population/chr4.add_cn_cp1.5_q15.filterByChr /share/raid002/lidong/40snp/population/chr4.LCcheck >/share/raid002/lidong/40snp/population/chr4.ratioCheck

open (A,$ARGV[0]) || die $!; ##genotyping result
while(<A>){
	chomp;
	my @words = split;
#	$genotyping{$words[0]."\t".$words[1]} = $words[5]."\t".$words[6];
	$genotyping{$words[0]."\t".$words[1]} = join "\t",@words[0..$#words];
}
close A;

my %check;
	open (A,$ARGV[1]) || die $!; ## LCcheck result
	while(<A>){
		chomp;
		my @words = split(/\t/);
		next if (! exists $genotyping{$words[0]."\t".$words[1]});
		$check{$words[0]."\t".$words[1]} = $words[-1];
	}
	close A;

foreach my $key (sort {my ($a,$b),(split(/\t/,$a))[1]<=>(split(/\t/,$b))[1]}keys %genotyping){
	my @snp = (split(/\t/,$genotyping{$key}))[5,6];
	my %snp;
	foreach (@snp){
		$snp{$_} = 1;
	}
	my %stat;my %het=();
	my @reads = split(/\,/,$check{$key});
#	print $key,"\t";
	print $genotyping{$key},"\t";
	for (my $i=0;$i<@reads;$i++){
		if ($reads[$i] =~ /^~$/){
#			print "- ";
		}else{
			my $seq = (split(/\~/,$reads[$i]))[0];my $quality = (split(/\~/,$reads[$i]))[1];
			my @allele = split(//,$seq);my @q=split(//,$quality);
			for (my $k = 0;$k<@q;$k++){
				$q[$k] = ord($q[$k]) - 64;
			}
			my %hash;
			for (my $j = 0;$j<@allele;$j++){
				next if (! exists $snp{$allele[$j]});
				next if ($q[$j]<5);	# SNPed reads with Q>=5
				$hash{$allele[$j]} ++;
				$stat{$allele[$j]} ++;
			}
			my $l = scalar keys %hash;
			if ($l == 1){
				foreach (keys %hash){
#					print $hash{$_},$_;
				}
			}
			if ($l > 1){
				my $output = "";
				foreach (sort keys %hash){
					$het{$_} += $hash{$_} if ($l > 1);	# += should be the same as =
					$output .= $hash{$_}.$_.":";
#					print $hash{$_},$_,":";
				}
				chop $output;	# remove the extra ':'
#				print $output;
			}
#			print " ";
		}
	}
	my @just;my @stat_reverse;
	foreach (sort keys %stat){
##		print $stat{$_},$_,":";
		push @just,[$stat{$_},$_];	# Count, base
	}
	@just = sort{my ($a,$b),$b->[0]<=>$a->[0] || $a->[1] cmp $b->[1]} @just;	# Sort DESC
	print $just[0][0],$just[0][1],":",$just[1][0],$just[1][1];	# all

	print " ";
	if (scalar keys %het > 1){
		print $het{$just[0][1]},$just[0][1],":",$just[1][0],$just[1][1];
		my $rate = $het{$just[0][1]}/$just[1][0];	# major_het / rare, should be 1 ?
		print "\t$rate";
#		my $hom = $just[1][0] - $het{$just[1][1]};
#		my $rate2 = $hom / $het{$just[1][1]};
#		print "\t$rate2";
	}else{
		print "0:0\t0";
	}
	print "\n";
}
