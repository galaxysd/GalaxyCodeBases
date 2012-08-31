#!/usr/bin/perl

=head1 Name

snp_dist_filter.pl   --  filter snp by neigboring distance

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  snp_dist_filter.pl <*.cns.snp>
  --distance <int>  set the minimum distance between two neighboring SNPs, default=5
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Dist_cutoff, $Verbose,$Help);
GetOptions(
	"distance:i"=>\$Dist_cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Dist_cutoff ||= 5;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $head;
my %data;

while (<>) {
	if (/^\#/) {
		$head .= $_;
	}

	if (/^\w/) {
		my @t = split /\s+/;
		my $ref = shift @t;
		my $cycle = shift @t;
		$data{$ref}{$cycle} = \@t;
	}
	if(/^</){
		last;	
	}
}

$head =~ s/Cycle_number 200/Cycle_number 100/;

print $head;

#print Dumper \%data;

my $outstr;
foreach my $ref (sort keys %data) {
	my $ref_p = $data{$ref};
	for (my $i=1; $i<=100; $i++) {
		$outstr .= "$ref\t$i";
		my $read1_p = $ref_p->{$i};
		my $read2_p = $ref_p->{$i+100};
		for (my $j=0; $j<@$read1_p; $j++) {
			my $sum = $read1_p->[$j] + $read2_p->[$j];
			$outstr .= "\t$sum";
		}
		$outstr .= "\n";
	}
}

print $outstr;

