#!/usr/bin/perl

=head1 Name

matrix_pirs2soapsnp.pl   --  convert pirs matrix to soapsnp format

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
 
  matrix_pirs2soapsnp.pl <*pirs.matrix>
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
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $Ref_num = 0;
my $Qual_num = 0;
my $Cyc_num = 0;
my $Base_num = 0;
my $head;
my %data;

my %BaseNum = ("A",0,"C",1,"G",2,"T",3);

while (<>) {
	if (/^\#/) {
		$head .= $_;
		$Ref_num = $1 if (/Ref_base_number (\d+)/);
		$Qual_num = $1 if (/Quality_number (\d+)/);
		$Cyc_num = $1 if (/Cycle_number (\d+)/);
		$Base_num = $1 if (/Seq_base_num (\d+)/);
	}

	if (/^\w/) {
		my @t = split /\s+/;
		my $ref = shift @t;
		my $cycle = shift @t;
		pop @t;
		$ref = $BaseNum{$ref};
		$cycle -= 1;
		
		if (@t % $Qual_num != 0) {
			print "error #####################\n";
		}
		for (my $i=0; $i<@t; $i++) {
			my $qual = $i%$Qual_num;
			my $base = ($i - $i%$Qual_num)/$Qual_num;
			$data{$qual}{$cycle}{$ref}{$base} = $t[$i];
		}
		
	}
	if(/^</){
		last;	
	}
}

print "#qual\tcycle\tA->A\tA->C\tA->G\tA->T\tC->A\tC->C\tC->G\tC->T\tG->A\tG->C\tG->G\tG->T\tT->A\tT->C\tT->G\tT->T\n";

#print Dumper \%data;

foreach my $qual (sort {$a<=>$b} keys %data) {
	my $qual_p = $data{$qual};
	foreach my $cyc (sort {$a<=>$b} keys %$qual_p) {
		print "$qual\t$cyc";
		my $cyc_p = $qual_p->{$cyc};
		foreach my $ref (sort {$a<=>$b} keys %$cyc_p) {
			my $ref_p = $cyc_p->{$ref};
			my $total = 0;
			foreach my $base (sort {$a<=>$b} keys %$ref_p) {
				my $freq = $ref_p->{$base};
				$freq = 1 if($freq == 0);
				$total += $freq;
			}
			foreach my $base (sort {$a<=>$b} keys %$ref_p) {
				my $freq = $ref_p->{$base};
				$freq = 1 if($freq == 0);
				printf("\t%e", $freq/$total);
				
			}
		}
		print "\n";
	}
}

