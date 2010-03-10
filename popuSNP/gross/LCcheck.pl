#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my %strain = ();

die "perl $0 <rawdatadir> <infile> <chr> <sample list> \n" unless @ARGV == 4;
print STDERR $0,"\n";

my $rawdatadir = shift;
my $infile = shift;
my $chr = shift;


my $sample_list = shift;


open (S,$sample_list) or die $!;
while(<S>)
{
        chomp;

        my @l = split /\s+/, $_;

        $strain{$l[0]} = $l[1];
#	print STDERR $l[0]."\t".$l[1]."\n";
}
close S;




my %pos;
my $name;
open (G,$infile) || die $!;  ## /share/raid6/jiangl/silkworm_soap/lane118_all_snp/total.snp
while(<G>){
	chomp;
	my @words = split(/\t/);
	my ($chr,$pos) = ($words[0],$words[1]);
	$pos{$chr."_".$pos} = join ("\t", @words[0..$#words]);
}
close G;

my %hash;
foreach my $lib (keys %strain){
	my $pe = `find $rawdatadir/$lib/ -name "*.sort"`;
	print STDERR $lib."\n";
	#my $pe = "$rawdatadir/$lib/${lib}.$chr";
	#print STDERR $pe."\n";
my @PE=split /\n/,$pe;
warn "[$_]\n" for @PE;
	for $pe (@PE) {
		open (A,$pe) || die $!;
		while(<A>){
			chomp;
			my @words = split(/\t/);
			my ($uniq,$chrF,$start) = ($words[3],$words[7],$words[8]);
			next if $chr ne $chrF;
			my $len = length($words[1]);
			my $end = $start + $len -1;
			next if ($uniq != 1);
			foreach my $i($start .. $end){
				if (exists $pos{$chr."_".$i}){
					my $offset = $i - $start;
					my $allele = substr($words[1],$offset,1);
					my $q = substr($words[2],$offset,1);
					$hash{$chr."_".$i} .= $allele.$q;
				}
			}
		}
		close A;
		foreach (keys %pos){
			$hash{$_} .= "~,";
		}
	}
}
#print Dumper %hash;
#die;
my %allele = ();
foreach my $key (sort {my ($a,$b),(split(/_/,$a))[1]<=>(split(/_/,$b))[1]} keys %hash){
#	my ($chr,$pos) = split(/_/,$key);
	my @t = split (/\,/,$hash{$key}); ##if not covered, it will report error
	for (my $j=0;$j<=$#t;$j++){
		my @a = split(//,$t[$j]);
		my $seq = "";
		my $quality = "";
		for (my $k=0;$k<=$#a;$k+=2){
#			next if ($a[$k] eq "~");
			$seq .= $a[$k];
			if ($#a>0 && $a[$k] ne "~"){
				$quality .= $a[$k+1];
			}
		}
		$allele{$key} .= $seq.$quality.",";
	}
#	my $temp;
#	if (exists $allele{'A'}){
#		$temp = $allele{'A'} . "." . "\t";
#	}else{
#		$temp = "." . "\t";
#	}
#	if (exists $allele{'C'}){
#		$temp .= $allele{'C'} . "." . "\t";
#	}else{
#		$temp .= "." . "\t";
#	}
#	if (exists $allele{'T'}){
#		$temp .= $allele{'T'} . "." . "\t";
#	}else{
#		$temp .= "." . "\t";
#	}
#	if (exists $allele{'G'}){
#		$temp .= $allele{'G'} . ".";
#	}else{
#		$temp .= ".";
#	}
#	foreach (keys %count){
#		$temp .= $count{$_}.$_."(".$allele{$_}.")";
#	}
#	%allele = ();
#	%count = ();
#	$hash{$key}{$key2} = $temp;
	print $pos{$key},"\t",$allele{$key},"\n";
}
