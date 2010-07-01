#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage=<<"USAGE";
Usage: <ASN.1 flatfile> | $0 [options] <stdout>
	-genome	symbol of genome version, default=reference for human
	-unval	use all snps including unvalidated ones, default=validated
	-h	this help

USAGE

my %options;
GetOptions (\%options, "genome=s", "unval", "h");
die $usage if(defined $options{h});
my $genome = (defined $options{genome}? $options{genome} : "reference");
my $unval = (defined $options{unval}? 1 : 0);


my $rs_id;  #rs***
my $validation; #YES|NO
my $alleles;
my $het;
my $se;
my $type;  # 1: insertion, 2: exact, 3: deletion, 4: range-insertion, 5: range-exact, 6: range-deletion
my $orient; # +,-
my $chr;
my $chr_loc;
my $count;   #count number of hits on the genome
my @hits;
my $first;

while (<stdin>) {
	if (/^(rs\d+)/) {
	   $rs_id = $1;
		 $validation=0;
		 $alleles="";
		 $het=0;
		 $se=0;
		 $type=0;
		 $count=0;
		 @hits = ();
		 $first = 1;
	}
	elsif (/^SNP/) {
		if(/alleles=\'(\S+)\'/) {
			$alleles=$1;
		}
		if(/het=([.\d]+)/) {
			$het=$1;
		}
		if(/se\(het\)=([\.\d]+)/) {
			$se=$1;
		}
	}
	elsif(/^VAL/) {
		if(/validated=YES/) {
			$validation=1;
		}
	}	
	elsif(/^CTG.*assembly=(\S+).*chr=(\S+).*chr-pos=(\d+).*loctype=(\d+).*orient=(\S)/) {
		next if($1 ne $genome);
		push(@hits, $_);
	}
	elsif((/^SEQ/) && $first)   #output
	{
		$first = 0;
		next if(@hits != 1);  #multiple hits on the genome
		if($hits[0] =~ /^CTG.*assembly=(\S+).*chr=(\S+).*chr-pos=(\d+).*loctype=(\d+).*orient=(\S)/)
		{
			$chr=$2;
			$chr_loc=$3;
			$type=$4;
			$orient=$5;
			
			#output
			next if((!$unval) && (!$validation));      #not validated
			next if(2 != $type);        #not single nucleotide polymorphism
			next if($alleles =~ /-/);   #with gaps
			next if("" eq $alleles);
			my @a=split(/\//, $alleles);
			my $long=0;
			foreach my $i(@a) {
				if((length $i) > 1) {
					$long=1;
					last;
				}
			}
			next if($long);
			if($orient eq '-') {        
				$alleles =~ tr/ACGT/TGCA/;   #convert to the orientation of chr
			}
			print "$rs_id\t$chr\t$chr_loc\t$alleles\t$het\t$se\n";
		}
	}
}
